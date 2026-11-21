//! Base-level alignment of a hit's blocks (optional output).
//!
//! The mapper places introns and reports exon blocks; nothing downstream of it needs more for
//! quantification. Tools that read the alignment itself do, so on request each block is
//! aligned base to base and the hit gets a full CIGAR (`M`, `I`, `D` inside blocks, `N`
//! between them, `S` at the ends) and an edit count. Introns are never aligned across: a block
//! is a read stretch against the genome stretch already chosen for it, both ends fixed.
//!
//! The alignment minimises the number of edits (unit costs), by the wavefront method: per edit
//! count, the furthest read position reached on each diagonal. Time and memory grow with the
//! square of the edits, not with the block length, so a clean exon costs little more than
//! reading it.

use crate::index::GenomeIndex;
use crate::map::SplicedHit;

/// CIGAR operation codes (BAM numbering).
const M: u32 = 0;
const I: u32 = 1;
const D: u32 = 2;
const N: u32 = 3;
const S: u32 = 4;

/// Blocks needing more edits than this are written as before: one `M` and the length
/// difference (the wavefronts kept for the traceback take four bytes times the square of this).
const MAX_EDITS: usize = 1500;

/// Blocks longer than twice this many read bases are aligned in pieces of this many.
const PIECE: usize = 512;

const NONE: i32 = i32::MIN / 2;

/// Bases of blocks aligned, and of those the bases aligned again with affine gaps (all threads).
pub static BASES: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static AGAIN: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// How a block is aligned.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Mode {
    /// Fewest edits.
    #[default]
    Edits,
    /// Affine gaps: the fewest-edits alignment first, then every stretch of it with an indel
    /// aligned again with affine gap costs (16 bases kept on each side of an indel; the long
    /// runs without one stay as they are).
    Affine,
}

/// Affine scoring (as minimap2's long-read presets: match 2, mismatch 4, gap 4 + 2 per base).
const MAT: [i8; 25] = {
    let mut m = [-4i8; 25];
    let mut i = 0;
    while i < 4 {
        m[i * 5 + i] = 2;
        i += 1;
    }
    m[24] = 0;
    m
};

/// A run of this many aligned bases without an indel is where the edit alignment is cut into
/// segments (in the middle of the run).
const CUT: u32 = 16;

/// Working memory of the block aligner.
#[derive(Default)]
pub struct Aligner {
    pub mode: Mode,
    /// Bases of blocks, and of those the bases aligned again with affine gaps.
    pub bases: u64,
    pub again: u64,
    tmp: Vec<u32>,
    seg: Vec<u32>,
    ksw: ksw2rs::Aligner,
    /// Wavefronts, one after the other: that of `s` edits holds diagonals `-s..=s` from `s * s`.
    waves: Vec<i32>,
    genome: Vec<u8>,
    codes: Vec<u8>,
}

fn push(ops: &mut Vec<u32>, len: u32, op: u32) {
    if len == 0 {
        return;
    }
    match ops.last_mut() {
        Some(last) if *last & 15 == op => *last += len << 4,
        _ => ops.push(len << 4 | op),
    }
}

/// [`push`] that joins runs only from `start` on (what is before is in the other order).
fn put(ops: &mut Vec<u32>, start: usize, len: u32, op: u32) {
    if ops.len() == start {
        if len > 0 {
            ops.push(len << 4 | op);
        }
    } else {
        push(ops, len, op);
    }
}

/// Align all of `a` (read) to all of `b` (genome) with the fewest edits; or, with `free`, all
/// of `a` to the start of `b`, ending wherever in `b` takes the fewest edits. Appends the CIGAR
/// operations to `ops` and returns the edits and the genome bases used, or `None` (nothing
/// appended) when more than `MAX_EDITS` are needed.
///
/// The sequences are seen through `ext(i, j, max)` alone: how many bases from read position
/// `i` and genome position `j` on are the same, `max` at most.
fn align(
    m: usize,
    n: usize,
    ext: &impl Fn(usize, usize, usize) -> usize,
    free: bool,
    waves: &mut Vec<i32>,
    ops: &mut Vec<u32>,
) -> Option<(u32, usize)> {
    let (m, n) = (m as i32, n as i32);
    let common = |i: i32, j: i32| ext(i as usize, j as usize, (m - i).min(n - j) as usize) as i32;
    if m == 0 || n == 0 {
        let d = if free { 0 } else { n };
        push(ops, m as u32, I);
        push(ops, d as u32, D);
        return Some(((m + d) as u32, d as usize));
    }
    let mut kf = n - m;
    waves.clear();
    waves.push(common(0, 0));
    let mut s = 0usize;
    // (diagonal k = genome position - read position; a wavefront entry is a read position)
    loop {
        if free {
            // the diagonal nearest the main one among those that used up the read
            let w = &waves[s * s..];
            let done = (0..w.len()).filter(|&x| w[x] >= m).min_by_key(|&x| x.abs_diff(s));
            if let Some(x) = done {
                kf = x as i32 - s as i32;
                break;
            }
        } else if kf.unsigned_abs() as usize <= s && waves[s * s + (kf + s as i32) as usize] >= m {
            break;
        }
        s += 1;
        if s > MAX_EDITS {
            return None;
        }
        let base = s * s;
        waves.resize(base + 2 * s + 1, NONE);
        let (prev, cur) = waves.split_at_mut(base);
        let prev = &prev[(s - 1) * (s - 1)..];
        let ps = s as i32 - 1;
        let get = |k: i32| {
            if k < -ps || k > ps {
                NONE
            } else {
                prev[(k + ps) as usize]
            }
        };
        let lo = (-(s as i32)).max(-m);
        let hi = (s as i32).min(n);
        for k in lo..=hi {
            // a mismatch on this diagonal, a genome base skipped, a read base left over
            let mut v = get(k) + 1;
            v = v.max(get(k - 1));
            v = v.max(get(k + 1) + 1);
            if v < 0 || v > m || v + k > n || v + k < 0 {
                // (the best way in overshoots an end: only a shorter one could stay inside,
                // and that happens on the last rows alone; try them in turn)
                v = NONE;
                for c in [get(k) + 1, get(k - 1), get(k + 1) + 1] {
                    if c >= 0 && c <= m && c + k <= n && c + k >= 0 {
                        v = v.max(c);
                    }
                }
                if v == NONE {
                    continue;
                }
            }
            let j = v + k;
            v += common(v, j);
            cur[(k + s as i32) as usize] = v;
        }
    }
    // Walk back from the end, an edit at a time; the matches before each edit are what the
    // diagonal ran on after it was entered.
    let start = ops.len();
    let edits = s as u32;
    let (mut k, mut s) = (kf, s);
    let mut at = m;
    loop {
        if s == 0 {
            put(ops, start, at as u32, M);
            break;
        }
        let ps = s as i32 - 1;
        let prev = &waves[(s - 1) * (s - 1)..s * s];
        let get = |k: i32| {
            if k < -ps || k > ps {
                NONE
            } else {
                prev[(k + ps) as usize]
            }
        };
        let ok = |c: i32, k: i32| c >= 0 && c <= m && c + k <= n && c + k >= 0;
        let (x, d, i) = (get(k) + 1, get(k - 1), get(k + 1) + 1);
        let mut v = NONE;
        for c in [x, d, i] {
            if ok(c, k) {
                v = v.max(c);
            }
        }
        debug_assert!(v != NONE && v <= at);
        put(ops, start, (at - v) as u32, M);
        if ok(x, k) && x == v {
            put(ops, start, 1, M);
            at = v - 1;
        } else if ok(d, k) && d == v {
            put(ops, start, 1, D);
            at = v;
            k -= 1;
        } else {
            put(ops, start, 1, I);
            at = v - 1;
            k += 1;
        }
        s -= 1;
    }
    ops[start..].reverse();
    if start > 0 && ops.len() > start && ops[start - 1] & 15 == ops[start] & 15 {
        ops[start - 1] += ops[start] & !15;
        ops.remove(start);
    }
    Some((edits, (m + kf) as usize))
}

/// Align a read stretch of `m` bases to a genome stretch of `n`, end to end (the sequences
/// seen through `ext`, as in [`align`]): appends the operations, returns the edits. A long block
/// goes a piece of the read at a time, each piece ending where in the genome it fits best, so
/// that the work grows with the edits in a piece.
fn block(
    m: usize,
    n: usize,
    ext: &impl Fn(usize, usize, usize) -> usize,
    waves: &mut Vec<i32>,
    ops: &mut Vec<u32>,
) -> u32 {
    let mut edits = 0;
    let mut whole = true;
    // (where the part still to align starts)
    let (mut a0, mut g0) = (0usize, 0usize);
    while whole && m - a0 > 2 * PIECE {
        let reach = (n - g0).min(PIECE + PIECE / 4);
        let sub = |i: usize, j: usize, max: usize| ext(a0 + i, g0 + j, max);
        match align(PIECE, reach, &sub, true, waves, ops) {
            Some((e, used)) => {
                edits += e;
                a0 += PIECE;
                g0 += used;
            }
            None => whole = false,
        }
    }
    let sub = |i: usize, j: usize, max: usize| ext(a0 + i, g0 + j, max);
    match if whole {
        align(m - a0, n - g0, &sub, false, waves, ops)
    } else {
        None
    } {
        Some((e, _)) => edits += e,
        None => {
            // (what is left of the block, as the block-level output writes it)
            let (lq, lt) = (m - a0, n - g0);
            let mm = lq.min(lt);
            edits += (0..mm).filter(|&x| sub(x, x, 1) == 0).count() as u32 + lq.abs_diff(lt) as u32;
            push(ops, mm as u32, M);
            push(ops, lq.saturating_sub(lt) as u32, I);
            push(ops, lt.saturating_sub(lq) as u32, D);
        }
    }
    edits
}

/// [`align`]'s view of two sequences held as codes, a byte a base (codes above 3, for symbols
/// that are not bases, match nothing).
fn ext_codes<'a>(a: &'a [u8], g: &'a [u8]) -> impl Fn(usize, usize, usize) -> usize + 'a {
    move |i, j, max| {
        let mut l = 0;
        while l < max && a[i + l] == g[j + l] && a[i + l] < 4 {
            l += 1;
        }
        l
    }
}

/// The same for two sequences known to hold bases only, eight at a time.
fn ext_bases<'a>(a: &'a [u8], g: &'a [u8]) -> impl Fn(usize, usize, usize) -> usize + 'a {
    move |i, j, max| {
        let (a, g) = (&a[i..i + max], &g[j..j + max]);
        let mut l = 0;
        while l + 8 <= max {
            let x = u64::from_le_bytes(a[l..l + 8].try_into().unwrap())
                ^ u64::from_le_bytes(g[l..l + 8].try_into().unwrap());
            if x != 0 {
                return l + (x.trailing_zeros() / 8) as usize;
            }
            l += 8;
        }
        while l < max && a[l] == g[l] {
            l += 1;
        }
        l
    }
}

/// Align all of `a` to all of `g` with affine gaps inside a band of `w`: appends the
/// operations. False, with nothing appended, when no end-to-end alignment came back.
fn affine(ksw: &mut ksw2rs::Aligner, a: &[u8], g: &[u8], w: i32, ops: &mut Vec<u32>) -> bool {
    if a.is_empty() || g.is_empty() {
        push(ops, a.len() as u32, I);
        push(ops, g.len() as u32, D);
        return true;
    }
    let ez = ksw.align(&ksw2rs::Extz2Input {
        query: a,
        target: g,
        m: 5,
        mat: &MAT,
        q: 4,
        e: 2,
        w,
        zdrop: -1,
        end_bonus: 0,
        flag: 0,
    });
    let (mut i, mut j) = (0usize, 0usize);
    for &o in &ez.cigar {
        if o & 15 != D {
            i += (o >> 4) as usize;
        }
        if o & 15 != I {
            j += (o >> 4) as usize;
        }
    }
    if (i, j) != (a.len(), g.len()) || ez.cigar.iter().any(|&o| o & 15 > D) {
        return false;
    }
    for &o in &ez.cigar {
        push(ops, o >> 4, o & 15);
    }
    true
}

/// The edits of an alignment of `a` to `g`.
fn count(a: &[u8], g: &[u8], ops: &[u32]) -> u32 {
    let (mut i, mut j, mut e) = (0usize, 0usize, 0u32);
    for &o in ops {
        let len = (o >> 4) as usize;
        match o & 15 {
            M => {
                e += (0..len).filter(|&t| a[i + t] != g[j + t] || a[i + t] > 3).count() as u32;
                i += len;
                j += len;
            }
            I => {
                e += len as u32;
                i += len;
            }
            _ => {
                e += len as u32;
                j += len;
            }
        }
    }
    e
}

impl Aligner {
    /// The full CIGAR of `h` for `read` (ASCII, as sequenced) into `ops` (as `len << 4 | op`),
    /// and the number of edits in its blocks.
    pub fn cigar(&mut self, ix: &GenomeIndex, read: &[u8], h: &SplicedHit, ops: &mut Vec<u32>) -> u32 {
        ops.clear();
        self.codes.clear();
        // (a base that is none of the four gets a code of its own: it matches nothing)
        let code = |b: u8| match b.to_ascii_uppercase() {
            b'A' => 0u8,
            b'C' => 1,
            b'G' => 2,
            b'T' | b'U' => 3,
            _ => 4,
        };
        if h.reverse {
            self.codes.extend(
                read.iter()
                    .rev()
                    .map(|&b| if code(b) < 4 { 3 - code(b) } else { 4 }),
            );
        } else {
            self.codes.extend(read.iter().map(|&b| code(b)));
        }
        // (a read with a symbol that is not a base is compared a base at a time)
        let plain = self.codes.iter().all(|&c| c < 4);
        let runs = ix.n_runs();
        let qlen = read.len() as u32;
        let (first, last) = (h.blocks[0], h.blocks[h.blocks.len() - 1]);
        push(ops, first.qs, S);
        let mut edits = 0;
        for (i, b) in h.blocks.iter().enumerate() {
            if i > 0 {
                ops.push((b.ts - h.blocks[i - 1].te) << 4 | N);
            }
            // The genome's symbols that are not bases (rare) are listed apart; a block that
            // holds one is compared a base at a time, with a code for them that matches nothing.
            let r = runs.partition_point(|&x| x <= b.ts);
            let clean = r == runs.len() || (r % 2 == 0 && runs[r] >= b.te);
            let (m, n) = ((b.qe - b.qs) as usize, (b.te - b.ts) as usize);
            self.genome.clear();
            {
                let mut p = b.ts;
                while p < b.te {
                    let w = ix.word(p);
                    let take = (b.te - p).min(32);
                    self.genome.extend((0..take).map(|k| ((w >> (2 * k)) & 3) as u8));
                    p += take;
                }
                if !clean {
                    for run in runs[r - r % 2..].chunks(2).take_while(|run| run[0] < b.te) {
                        for p in run[0].max(b.ts)..run[1].min(b.te) {
                            self.genome[(p - b.ts) as usize] = 4;
                        }
                    }
                }
            }
            let a = &self.codes[b.qs as usize..b.qe as usize];
            self.bases += m as u64;
            let waves = &mut self.waves;
            let mut first_pass = |ops: &mut Vec<u32>| {
                if plain && clean {
                    block(m, n, &ext_bases(a, &self.genome), waves, ops)
                } else {
                    block(m, n, &ext_codes(a, &self.genome), waves, ops)
                }
            };
            if self.mode != Mode::Edits {
                let start = ops.len();
                first_pass(ops);
                // Cut the edit alignment in the middle of its long runs without an indel,
                // and align again the segments whose indels call for it.
                self.tmp.clear();
                self.tmp.extend_from_slice(&ops[start..]);
                ops.truncate(start);
                self.seg.clear();
                let (mut i0, mut j0, mut i, mut j) = (0usize, 0usize, 0usize, 0usize);
                let n_ops = self.tmp.len();
                for x in 0..=n_ops {
                    // (the part of a long run before its middle closes the segment)
                    let (len, op) = if x < n_ops {
                        (self.tmp[x] >> 4, self.tmp[x] & 15)
                    } else {
                        (0, M)
                    };
                    let cut = x == n_ops || (op == M && len >= CUT);
                    // (of a longer run only the `CUT` bases at each end go to a segment; the
                    // middle stays as it is)
                    let keep = if cut && len > 2 * CUT { len - 2 * CUT } else { 0 };
                    let head = if cut { (len - keep) / 2 } else { len };
                    if head > 0 {
                        push(&mut self.seg, head, op);
                        if op != D {
                            i += head as usize;
                        }
                        if op != I {
                            j += head as usize;
                        }
                    }
                    if cut {
                        let mut runs = 0u32;
                        for &o in &self.seg {
                            if o & 15 != M {
                                runs += 1;
                            }
                        }
                        let wanted = runs > 0;
                        let (sa, sg) = (&a[i0..i], &self.genome[j0..j]);
                        let w = sa.len().abs_diff(sg.len()) as i32 + 16;
                        if wanted && affine(&mut self.ksw, sa, sg, w, ops) {
                            self.again += sa.len() as u64;
                        } else {
                            for &o in &self.seg {
                                push(ops, o >> 4, o & 15);
                            }
                        }
                        self.seg.clear();
                        push(ops, keep, M);
                        i += keep as usize;
                        j += keep as usize;
                        (i0, j0) = (i, j);
                        if len > head + keep {
                            push(&mut self.seg, len - head - keep, op);
                            i += (len - head - keep) as usize;
                            j += (len - head - keep) as usize;
                        }
                    }
                }
                edits += count(a, &self.genome, &ops[start..]);
                continue;
            }
            edits += first_pass(ops);
        }
        if last.qe < qlen {
            ops.push((qlen - last.qe) << 4 | S);
        }
        BASES.fetch_add(
            std::mem::take(&mut self.bases),
            std::sync::atomic::Ordering::Relaxed,
        );
        AGAIN.fetch_add(
            std::mem::take(&mut self.again),
            std::sync::atomic::Ordering::Relaxed,
        );
        edits
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Edit distance by the plain table.
    fn distance(a: &[u8], b: &[u8]) -> u32 {
        let mut row: Vec<u32> = (0..=b.len() as u32).collect();
        for i in 1..=a.len() {
            let mut diag = row[0];
            row[0] = i as u32;
            for j in 1..=b.len() {
                let c = (diag + (a[i - 1] != b[j - 1]) as u32)
                    .min(row[j] + 1)
                    .min(row[j - 1] + 1);
                diag = row[j];
                row[j] = c;
            }
        }
        row[b.len()]
    }

    /// The operations consume both sequences exactly, and their edits are the fewest possible.
    #[test]
    fn alignments_are_optimal() {
        let mut x = 0x9e37_79b9_7f4a_7c15u64;
        let mut rnd = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let mut waves = Vec::new();
        for round in 0..400 {
            // (every tenth round a block long enough to go in pieces)
            let long = round == 9 || round == 19;
            let n = if long { 9000 } else { (rnd() % 300) as usize };
            let b: Vec<u8> = (0..n).map(|_| (rnd() % 4) as u8).collect();
            let rate = [0, 2, 10, 30][round % 4];
            let mut a = Vec::new();
            for &c in &b {
                match rnd() % 100 {
                    r if r < rate / 3 => {}
                    r if r < 2 * rate / 3 => a.extend([c, (rnd() % 4) as u8]),
                    r if r < rate => a.push((c + 1) % 4),
                    _ => a.push(c),
                }
            }
            let mut ops = Vec::new();
            // (both views of the sequences must give the same alignment)
            let mut ops2 = Vec::new();
            let e2 = block(a.len(), b.len(), &ext_bases(&a, &b), &mut waves, &mut ops2);
            let e = block(a.len(), b.len(), &ext_codes(&a, &b), &mut waves, &mut ops);
            assert_eq!((e, &ops), (e2, &ops2));
            if long {
                // (each piece ends where it fits best on its own, which can cost a few edits)
                let d = distance(&a, &b);
                assert!(e <= d + d / 100 + 4);
            } else {
                assert_eq!(e, distance(&a, &b));
            }
            let (mut i, mut j, mut edits) = (0usize, 0usize, 0u32);
            for (k, &op) in ops.iter().enumerate() {
                let len = (op >> 4) as usize;
                assert!(len > 0 && (k == 0 || ops[k - 1] & 15 != op & 15));
                match op & 15 {
                    M => {
                        edits += (0..len).filter(|&t| a[i + t] != b[j + t]).count() as u32;
                        i += len;
                        j += len;
                    }
                    I => {
                        edits += len as u32;
                        i += len;
                    }
                    _ => {
                        edits += len as u32;
                        j += len;
                    }
                }
            }
            assert_eq!((i, j, edits), (a.len(), b.len(), e));
        }
    }
}
