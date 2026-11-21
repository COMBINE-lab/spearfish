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
const PIECE: usize = 2048;

const NONE: i32 = i32::MIN / 2;

/// Working memory of the block aligner.
#[derive(Default)]
pub struct Aligner {
    /// Wavefronts, one after the other: that of `s` edits holds diagonals `-s..=s` from `s * s`.
    waves: Vec<i32>,
    genome: Vec<u8>,
    codes: Vec<u8>,
}

/// Length of the common prefix of `a` and `b`.
#[inline]
fn common(a: &[u8], b: &[u8]) -> usize {
    let n = a.len().min(b.len());
    let mut i = 0;
    while i + 8 <= n {
        let x = u64::from_le_bytes(a[i..i + 8].try_into().unwrap())
            ^ u64::from_le_bytes(b[i..i + 8].try_into().unwrap());
        if x != 0 {
            return i + (x.trailing_zeros() / 8) as usize;
        }
        i += 8;
    }
    while i < n && a[i] == b[i] {
        i += 1;
    }
    i
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
fn align(a: &[u8], b: &[u8], free: bool, waves: &mut Vec<i32>, ops: &mut Vec<u32>) -> Option<(u32, usize)> {
    let (m, n) = (a.len() as i32, b.len() as i32);
    if m == 0 || n == 0 {
        let d = if free { 0 } else { n };
        push(ops, m as u32, I);
        push(ops, d as u32, D);
        return Some(((m + d) as u32, d as usize));
    }
    let mut kf = n - m;
    waves.clear();
    waves.push(common(a, b) as i32);
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
            v += common(&a[v as usize..], &b[j as usize..]) as i32;
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

/// Align read stretch `a` to genome stretch `g`, end to end: appends the operations, returns
/// the edits. A long block goes a piece of the read at a time, each piece ending where in the
/// genome it fits best, so that the work grows with the edits in a piece.
fn block(mut a: &[u8], mut g: &[u8], waves: &mut Vec<i32>, ops: &mut Vec<u32>) -> u32 {
    let mut edits = 0;
    let mut whole = true;
    while whole && a.len() > 2 * PIECE {
        let reach = g.len().min(PIECE + PIECE / 4);
        match align(&a[..PIECE], &g[..reach], true, waves, ops) {
            Some((e, used)) => {
                edits += e;
                a = &a[PIECE..];
                g = &g[used..];
            }
            None => whole = false,
        }
    }
    match if whole {
        align(a, g, false, waves, ops)
    } else {
        None
    } {
        Some((e, _)) => edits += e,
        None => {
            // (what is left of the block, as the block-level output writes it)
            let (lq, lt) = (a.len(), g.len());
            let mm = lq.min(lt);
            edits += (0..mm).filter(|&x| a[x] != g[x]).count() as u32 + lq.abs_diff(lt) as u32;
            push(ops, mm as u32, M);
            push(ops, lq.saturating_sub(lt) as u32, I);
            push(ops, lt.saturating_sub(lq) as u32, D);
        }
    }
    edits
}

impl Aligner {
    /// The full CIGAR of `h` for `read` (ASCII, as sequenced) into `ops` (as `len << 4 | op`),
    /// and the number of edits in its blocks.
    pub fn cigar(&mut self, ix: &GenomeIndex, read: &[u8], h: &SplicedHit, ops: &mut Vec<u32>) -> u32 {
        ops.clear();
        self.codes.clear();
        let code = |b: u8| match b.to_ascii_uppercase() {
            b'C' => 1u8,
            b'G' => 2,
            b'T' | b'U' => 3,
            _ => 0,
        };
        if h.reverse {
            self.codes.extend(read.iter().rev().map(|&b| 3 - code(b)));
        } else {
            self.codes.extend(read.iter().map(|&b| code(b)));
        }
        let qlen = read.len() as u32;
        let (first, last) = (h.blocks[0], h.blocks[h.blocks.len() - 1]);
        push(ops, first.qs, S);
        let mut edits = 0;
        for (i, b) in h.blocks.iter().enumerate() {
            if i > 0 {
                ops.push((b.ts - h.blocks[i - 1].te) << 4 | N);
            }
            // the genome stretch, 32 bases a word
            self.genome.clear();
            let mut p = b.ts;
            while p < b.te {
                let w = ix.word(p);
                let take = (b.te - p).min(32);
                self.genome.extend((0..take).map(|k| ((w >> (2 * k)) & 3) as u8));
                p += take;
            }
            edits += block(
                &self.codes[b.qs as usize..b.qe as usize],
                &self.genome,
                &mut self.waves,
                ops,
            );
        }
        if last.qe < qlen {
            ops.push((qlen - last.qe) << 4 | S);
        }
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
            let e = block(&a, &b, &mut waves, &mut ops);
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
