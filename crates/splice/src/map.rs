//! Spliced mapping of one read: seed, pick loci, chain across introns, place junctions.
//!
//! Nothing here aligns across an intron or through an exon's interior. Level-3 cores give
//! anchors; anchors are grouped into loci and chained with a cost that lets the reference
//! skip ahead cheaply; consecutive chain anchors that the reference separates by more than
//! `min_intron` become a junction, and each junction is placed by a small banded edit
//! distance from both sides of the read bases between the two anchors.

use seeds::BaseSeeder;

use crate::gap::{self, CLIP, Hop, ReadView, Solved, View};
use crate::index::{GenomeIndex, L2_SPEC, L3_SPEC, code};

#[derive(Clone, Debug)]
pub struct MapOptions {
    /// Seeds with more occurrences than this are not used.
    pub cap: usize,
    /// The cap for a second attempt when a read gets no hit.
    pub cap_retry: usize,
    pub max_intron: u32,
    /// A reference skip of at least this many bases between two chained anchors is a junction.
    pub min_intron: u32,
    /// Loci chained per read (best seed weight first).
    pub max_loci: usize,
    /// Only forward-strand transcripts (direct RNA): the read is the transcript's sense strand.
    pub forward_only: bool,
    /// Seed noisy reads with pairs of level-2 cores as well as level-3 cores from the start
    /// (otherwise only reads that level 3 anchors weakly get them, in a second pass).
    pub pairs: bool,
    /// Junction cost per read base (times 64) above which the solver looks for missed exons,
    /// unless the read's own junctions usually cost more.
    pub trigger_x64: u32,
    /// Hits that were aligned and score below this are dropped.
    pub min_score: i32,
    /// Keep each hit's chained seed matches (for the base-level alignment, which then only has
    /// to align between them).
    pub keep_seeds: bool,
}

impl Default for MapOptions {
    fn default() -> Self {
        MapOptions {
            cap: 64,
            cap_retry: 1000,
            max_intron: 200_000,
            min_intron: 40,
            max_loci: 6,
            forward_only: false,
            pairs: true,
            trigger_x64: 16,
            min_score: 100,
            keep_seeds: false,
        }
    }
}

/// One aligned block: read `[qs, qe)` (in the mapped orientation) against genome `[ts, te)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Block {
    pub qs: u32,
    pub qe: u32,
    pub ts: u32,
    pub te: u32,
}

#[derive(Clone, Debug)]
pub struct SplicedHit {
    /// The read maps as its reverse complement.
    pub reverse: bool,
    /// Blocks in genome order, on the concatenated axis; read coordinates are those of the
    /// mapped orientation.
    pub blocks: Vec<Block>,
    /// Transcript strand from the splice sites: `+`, `-`, or `.` when there is no evidence.
    pub tx_strand: u8,
    pub score: i32,
    pub anchors: u32,
    /// Alignment cost summed over the placed junctions and read ends.
    pub edits: u32,
    pub primary: bool,
    /// One letter per junction, in genome order, saying how it was made (for diagnosis).
    pub origins: Vec<u8>,
    /// The chained seed matches as (read position, genome position, length), in chain order;
    /// empty unless [`MapOptions::keep_seeds`]. They are where the chain was built from, not
    /// part of the result: a block need not hold them, nor they match base for base.
    pub seeds: Vec<(u32, u32, u16)>,
}

#[derive(Clone, Copy, Debug)]
struct Anchor {
    t: u32,
    q: u32,
    ql: u16,
    tl: u16,
    /// Read offsets, from `q`, of the end of the core's first homopolymer run and of the start
    /// of its last one. A run at a core's edge can continue into the next exon in the read but
    /// not in the genome, so only the points inside these runs correspond exactly.
    q_in: u16,
    q_out: u16,
    w: f32,
}

/// Per-thread working memory.
#[derive(Default)]
pub struct Scratch {
    codes: [Vec<u8>; 2],
    hits: Vec<Anchor>,
    /// Scratch of the hit sort: keys, and the hits in their new order.
    hit_keys: Vec<u64>,
    /// Pairs of level-2 cores as seeds, per orientation, once they were asked for.
    pair_seeds: [Vec<L2Seed>; 2],
    pairs_done: [bool; 2],
    hits_sorted: Vec<Anchor>,
    f: Vec<f32>,
    back: Vec<u32>,
    chain: Vec<u32>,
    g1: Vec<u8>,
    e1: Vec<u16>,
    e2: Vec<u16>,
    etmp: Vec<u16>,
    gap: crate::gap::Bufs,
    /// Level-2 cores of the read, per orientation, in read order.
    l2: [Vec<L2Seed>; 2],
    /// Level-3 cores of the read, per orientation.
    l3: [Vec<L2Seed>; 2],
    /// Level-3 cores of the read found in the genome, and their number, in the better orientation
    /// (from the last mapping pass).
    found_best: usize,
    l3_best: usize,
    /// Whether the last mapping pass had any seed hit to chain.
    had_hits: bool,
    /// Marks the level-2 label ids wanted by the rescue under way (indexed by id).
    want: Vec<u32>,
    /// Short-word search of a gap: per word, the last read offset with it (plus one); per
    /// read offset, the one before with the same word (plus one); the matches found.
    kfirst: Vec<u16>,
    knext: Vec<u16>,
    khits: Vec<(u32, u16)>,
    /// Counters for reports.
    pub stats: Stats,
}

#[derive(Clone, Copy, Debug)]
struct L2Seed {
    q: u32,
    ql: u16,
    q_in: u16,
    q_out: u16,
    /// Homopolymer runs in the core.
    runs: u16,
    id: u16,
    label: u64,
    /// Where the label's occurrences are in the index (level-3 cores only).
    occ: (u32, u32),
}

#[derive(Default, Clone, Copy, Debug)]
pub struct Stats {
    pub junctions: u64,
    pub rescues_tried: u64,
    pub rescues_kept: u64,
    pub window_bases: u64,
    /// Rescues between anchors: nothing found; found but not cheaper; cheaper but a junction
    /// along the new path is poor; not tried because the junction was too costly to place at all.
    pub empty: u64,
    pub not_cheaper: u64,
    pub poor_junction: u64,
    pub unplaced: u64,
    pub end_tried: u64,
    pub end_kept: u64,
    /// Hits scored by aligning their blocks.
    pub scored: u64,
    pub retried: u64,
    /// Spliced hits ranked above an unspliced rival of slightly higher score.
    pub spliced_preferred: u64,
    /// Hits whose blocks had to be made consistent (a defect upstream when not zero).
    pub tidied: u64,
    /// Splice-site searches: tried, kept, genome bases scanned.
    pub site_tried: u64,
    pub site_kept: u64,
    pub site_bases: u64,
    /// Junctions placed without a table (an exact match from both sides, canonical sites).
    pub fast_placed: u64,
    /// Gap solver work: solves, prefix-table cells, window cells; and cells of other tables.
    pub solves: u64,
    pub prefix_cells: u64,
    pub window_cells: u64,
}

impl Stats {
    pub fn add(&mut self, o: &Stats) {
        self.junctions += o.junctions;
        self.rescues_tried += o.rescues_tried;
        self.rescues_kept += o.rescues_kept;
        self.window_bases += o.window_bases;
        self.empty += o.empty;
        self.not_cheaper += o.not_cheaper;
        self.poor_junction += o.poor_junction;
        self.unplaced += o.unplaced;
        self.end_tried += o.end_tried;
        self.end_kept += o.end_kept;
        self.scored += o.scored;
        self.retried += o.retried;
        self.spliced_preferred += o.spliced_preferred;
        self.tidied += o.tidied;
        self.site_tried += o.site_tried;
        self.site_kept += o.site_kept;
        self.site_bases += o.site_bases;
        self.fast_placed += o.fast_placed;
        self.solves += o.solves;
        self.prefix_cells += o.prefix_cells;
        self.window_cells += o.window_cells;
    }
}

/// A spliced mapper over one index; `Send + Sync` (working memory is per thread, in [`Scratch`]).
pub struct Mapper {
    pub index: std::sync::Arc<GenomeIndex>,
    pub opt: MapOptions,
    seeder: BaseSeeder,
    seeder2: BaseSeeder,
}

pub(crate) const INF: u16 = u16::MAX / 4;

/// Alignment costs, relative to a match: a mismatch, opening a gap, each gap base. With
/// affine gaps one long insertion in a read is cheaper than the scattered matches a junction
/// placed after it would collect, which unit costs get wrong.
pub(crate) const MISMATCH: u16 = 6;
pub(crate) const GAP_OPEN: u16 = 4;
pub(crate) const GAP_EXT: u16 = 2;
/// Each read base in a gap also forfeits the match it would have scored.
pub(crate) const INS_EXT: u16 = GAP_EXT + 2;
/// Splice-site costs per site: GT and AG are free, a GC left end costs `SEMI` (GC..AG), any
/// other `NONCANON`. The two ends of an AT..AC intron cost `ATAC` each: sites are priced one
/// at a time, so a cheap AT or AC would also make GT..AC and AT..AG cheap, which are not
/// introns, and AT..AC itself is rare.
const SEMI: u16 = 11;
const NONCANON: u16 = 26;
const ATAC: u16 = 14;
/// Sort hits by genome position, then read position. The keys (position and index in one
/// word) are sorted, not the hits themselves, which are then moved once.
fn sort_hits(hits: &mut Vec<Anchor>, keys: &mut Vec<u64>, sorted: &mut Vec<Anchor>) {
    if hits.len() < 32 {
        hits.sort_unstable_by_key(|a| (a.t, a.q));
        return;
    }
    keys.clear();
    keys.extend(
        hits.iter()
            .enumerate()
            .map(|(i, a)| ((a.t as u64) << 32) | i as u64),
    );
    if keys.len() < 512 {
        keys.sort_unstable();
    } else {
        // by position in three counting passes of 11, 11 and 10 bits (the index in the low
        // half keeps its order)
        const BITS: [u32; 3] = [11, 11, 10];
        let n = keys.len();
        let mut counts = [[0u32; 2048]; 3];
        for &k in keys.iter() {
            let t = (k >> 32) as usize;
            counts[0][t & 2047] += 1;
            counts[1][(t >> 11) & 2047] += 1;
            counts[2][t >> 22] += 1;
        }
        keys.resize(2 * n, 0);
        let (mut a, mut b) = keys.split_at_mut(n);
        let mut shift = 32;
        for (p, c) in counts.iter_mut().enumerate() {
            let mask = (1u64 << BITS[p]) - 1;
            let mut at = 0u32;
            for x in c.iter_mut() {
                let here = *x;
                *x = at;
                at += here;
            }
            for &k in a.iter() {
                let d = ((k >> shift) & mask) as usize;
                b[c[d] as usize] = k;
                c[d] += 1;
            }
            std::mem::swap(&mut a, &mut b);
            shift += BITS[p];
        }
        // (three passes: the sorted keys are in the second half)
        keys.copy_within(n.., 0);
        keys.truncate(n);
    }
    sorted.clear();
    sorted.extend(keys.iter().map(|&k| hits[(k & 0xffff_ffff) as usize]));
    // hits at one genome position (few): by read position
    let mut s = 0;
    for i in 1..=sorted.len() {
        if i == sorted.len() || sorted[i].t != sorted[s].t {
            if i - s > 1 {
                sorted[s..i].sort_unstable_by_key(|a| a.q);
            }
            s = i;
        }
    }
    std::mem::swap(hits, sorted);
}

/// Banded affine-gap costs of every prefix of `a` against prefixes of `b`.
/// `e[i * w + k]` is the cost of `a[..i]` against `b[..j]` with `j = i + k - band`,
/// `w = 2 * band + 1`. `tmp` holds the gap states. With `drop > 0` the rows stop once the best
/// score of a row (two per base of `a`, minus its cost) falls `drop` below the best seen;
/// returns the number of rows of `a` computed.
#[allow(clippy::too_many_arguments)]
pub(crate) fn band_prefix(
    a: impl Fn(usize) -> u8,
    m: usize,
    b: impl Fn(usize) -> u8,
    n: usize,
    band: usize,
    e: &mut Vec<u16>,
    tmp: &mut Vec<u16>,
    drop: i32,
) -> usize {
    debug_assert!(band >= 1);
    let w = 2 * band + 1;
    e.clear();
    e.resize((m + 1) * w, INF);
    // gap states of the row above (a gap that consumes `a`), per k, and those of the row being
    // filled; each with one spare cell
    tmp.clear();
    tmp.resize(2 * (w + 1), INF);
    let (mut cur, mut new) = tmp.split_at_mut(w + 1);
    // the bases of `b` once each, not once per cell, and only as far as the rows get
    let mut bb: Vec<u8> = Vec::with_capacity(n.min(m + band));
    let open_ins = GAP_OPEN + INS_EXT;
    let mut best = 0i32;
    for i in 0..=m {
        // the row's cells: j = i + k - band within 0..=n
        let k_lo = band.saturating_sub(i);
        if i > n + band {
            // past the end of `b`: nothing reachable in this or any later row
            if drop > 0 {
                return i;
            }
            cur.fill(INF);
            continue;
        }
        let k_hi = (n + band - i).min(w - 1);
        while bb.len() < (i + band).min(n) {
            bb.push(b(bb.len()));
        }
        let (above, here) = e.split_at_mut(i * w);
        let here = &mut here[..w];
        if i == 0 {
            here[band] = 0;
        } else {
            let above = &above[(i - 1) * w..];
            let ai = a(i - 1);
            // Every way into a cell but along the row. The cell above (i - 1, j) is k + 1 in
            // the previous row, the diagonal one the same k. The first cell of the row when
            // it is column 0, and the band's last cell, lack one of them.
            let (mut lo, mut hi) = (k_lo, k_hi);
            if i + k_lo == band {
                let up = (cur[k_lo + 1] + INS_EXT).min(above[k_lo + 1] + open_ins).min(INF);
                here[k_lo] = up;
                new[k_lo] = up;
                lo += 1;
            }
            if k_hi == w - 1 && hi >= lo {
                let j = i + hi - band;
                here[hi] = (above[hi] + if ai != bb[j - 1] { MISMATCH } else { 0 }).min(INF);
                new[hi] = INF;
                hi -= 1;
            }
            if lo <= hi {
                prefix_cells(
                    &above[lo..=hi],
                    &above[lo + 1..=hi + 1],
                    &cur[lo + 1..=hi + 1],
                    &mut here[lo..=hi],
                    &mut new[lo..=hi],
                    &bb[i + lo - band - 1..],
                    ai,
                );
            }
        }
        new[..k_lo].fill(INF);
        new[k_hi + 1..].fill(INF);
        if i == 0 {
            new.fill(INF);
        }
        // along the row: a gap that consumes `b`
        let row = prefix_left(&mut here[k_lo..=k_hi]);
        std::mem::swap(&mut cur, &mut new);
        if drop > 0 {
            let s = 2 * i as i32 - row as i32;
            best = best.max(s);
            if s + drop < best {
                return i;
            }
        }
    }
    m
}

/// The cells of one row of [`band_prefix`] that have both a cell above and a diagonal one:
/// the cheaper of the two ways in, and the gap state for the next row.
#[inline(always)]
fn prefix_cells_body(
    diag: &[u16],
    above: &[u16],
    gap: &[u16],
    here: &mut [u16],
    new: &mut [u16],
    bb: &[u8],
    ai: u8,
) {
    let len = here.len();
    let (diag, above, gap, new, bb) = (
        &diag[..len],
        &above[..len],
        &gap[..len],
        &mut new[..len],
        &bb[..len],
    );
    let open_ins = GAP_OPEN + INS_EXT;
    for x in 0..len {
        let up = (gap[x] + INS_EXT).min(above[x] + open_ins).min(INF);
        let d = diag[x] + ((ai != bb[x]) as u16) * MISMATCH;
        new[x] = up;
        here[x] = d.min(up);
    }
}
avx2_variant! {
    fn prefix_cells = prefix_cells_body, prefix_cells_avx2(
        diag: &[u16], above: &[u16], gap: &[u16], here: &mut [u16], new: &mut [u16], bb: &[u8], ai: u8,
    ) -> ()
}

/// Gaps along one row of [`band_prefix`], over cells that hold every other way in; returns
/// the row's lowest cost. (A gap is never opened from a cell that a gap reached, so the gap
/// state follows from the cells as they were.)
#[inline]
fn prefix_left(here: &mut [u16]) -> u16 {
    #[cfg(target_arch = "x86_64")]
    if here.len() >= 16 && std::is_x86_feature_detected!("avx512bw") {
        // SAFETY: the feature was just detected.
        return unsafe { prefix_left_avx512(here) };
    }
    prefix_left_plain(here)
}

#[inline]
fn prefix_left_plain(here: &mut [u16]) -> u16 {
    let open_del = GAP_OPEN + GAP_EXT;
    let (mut left, mut row) = (INF, INF);
    for c in here.iter_mut() {
        let t = *c;
        let v = t.min(left);
        *c = v;
        left = (left + GAP_EXT).min(t + open_del);
        row = row.min(v);
    }
    row
}

/// `prefix_left_plain`, 32 cells at a time, with the same results. The gap that reaches cell
/// x from cell y < x costs open + (x - 1 - y) * ext, so with a[y] = cell y + open +
/// (len - 1 - y) * ext the cheapest gap into x is the least a[y] before x, less
/// (len - x) * ext: a running minimum, taken inside the 32 by doubling steps and one
/// value carried between them. (Cells are at most INF, so nothing overflows for any row
/// that fits a band.)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw")]
unsafe fn prefix_left_avx512(here: &mut [u16]) -> u16 {
    use std::arch::x86_64::*;
    static IOTA: [u16; 32] = {
        let mut a = [0u16; 32];
        let mut i = 0;
        while i < 32 {
            a[i] = i as u16;
            i += 1;
        }
        a
    };
    let open_del = GAP_OPEN + GAP_EXT;
    let len = here.len();
    debug_assert!(len < 20_000);
    // SAFETY (all of it): loads and stores are masked to the cells that remain in `here`;
    // IOTA is 32 words.
    unsafe {
        let iota = _mm512_loadu_si512(IOTA.as_ptr() as *const _);
        let ramp = _mm512_mullo_epi16(iota, _mm512_set1_epi16(GAP_EXT as i16));
        let none = _mm512_set1_epi16(-1);
        let zero = _mm512_setzero_si512();
        // (bytes 14 and 15 of each 128 bits: its last cell, to all eight of its cells)
        let tops = _mm512_set1_epi16(0x0f0e);
        let last = _mm512_set1_epi16(31);
        let mut carry = _mm512_set1_epi16((INF + len as u16 * GAP_EXT) as i16);
        let mut low = none;
        let mut x = 0;
        while x < len {
            let rem = len - x;
            let k: __mmask32 = if rem >= 32 { !0 } else { (1u32 << rem) - 1 };
            let at = here.as_mut_ptr().add(x) as *mut i16;
            // (cells past the end read as "none"; what they turn into is never stored and
            // never reaches a cell before them)
            let c = _mm512_mask_loadu_epi16(none, k, at);
            let a0 = open_del + (len - 1 - x) as u16 * GAP_EXT;
            let p = _mm512_sub_epi16(_mm512_add_epi16(c, _mm512_set1_epi16(a0 as i16)), ramp);
            // The running minimum of p, as the running maximum of its complement (so that
            // the zeros a shift brings in change nothing): first inside each eight cells,
            // then each eight takes the last cell of the eights before it.
            let mut q = _mm512_xor_si512(p, none);
            q = _mm512_max_epu16(q, _mm512_bslli_epi128::<2>(q));
            q = _mm512_max_epu16(q, _mm512_bslli_epi128::<4>(q));
            q = _mm512_max_epu16(q, _mm512_bslli_epi128::<8>(q));
            let t = _mm512_shuffle_epi8(q, tops);
            let e = _mm512_max_epu16(
                _mm512_max_epu16(
                    _mm512_alignr_epi64::<6>(t, zero),
                    _mm512_alignr_epi64::<4>(t, zero),
                ),
                _mm512_alignr_epi64::<2>(t, zero),
            );
            let run = _mm512_xor_si512(_mm512_max_epu16(q, e), none);
            // (the running minimum includes the cell itself, as a gap of no length that
            // costs GAP_OPEN: never cheaper than the cell)
            let back = _mm512_set1_epi16(((len - x) as u16 * GAP_EXT) as i16);
            let left = _mm512_add_epi16(_mm512_sub_epi16(_mm512_min_epu16(carry, run), back), ramp);
            let v = _mm512_min_epu16(c, left);
            _mm512_mask_storeu_epi16(at, k, v);
            low = _mm512_mask_min_epu16(low, k, low, v);
            carry = _mm512_min_epu16(carry, _mm512_permutexvar_epi16(last, run));
            x += 32;
        }
        let l = _mm256_min_epu16(_mm512_castsi512_si256(low), _mm512_extracti64x4_epi64::<1>(low));
        let l = _mm_min_epu16(_mm256_castsi256_si128(l), _mm256_extracti128_si256::<1>(l));
        (_mm_extract_epi16::<0>(_mm_minpos_epu16(l)) as u16).min(INF)
    }
}

/// Calls `f` with every index of `ids` whose value is one of `want` (at most eight values), in order.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw")]
unsafe fn ids_among_avx512(ids: &[u16], want: &[u16], mut f: impl FnMut(usize)) {
    use std::arch::x86_64::*;
    let mut w = [_mm512_setzero_si512(); 8];
    for (k, &x) in want.iter().enumerate() {
        w[k] = _mm512_set1_epi16(x as i16);
    }
    let w = &w[..want.len()];
    let mut i = 0;
    while i + 32 <= ids.len() {
        // SAFETY: 32 values from `i` on lie in `ids`.
        let v = unsafe { _mm512_loadu_si512(ids.as_ptr().add(i).cast()) };
        let mut m: u32 = 0;
        for x in w {
            m |= _mm512_cmpeq_epi16_mask(v, *x);
        }
        while m != 0 {
            f(i + m.trailing_zeros() as usize);
            m &= m - 1;
        }
        i += 32;
    }
    for (j, id) in ids.iter().enumerate().skip(i) {
        if want.contains(id) {
            f(j);
        }
    }
}

/// The most pieces of a gap's read bases looked for in one approximate search.
const SCAN_PIECES: usize = 16;

/// Approximate search of the genome `[lo, hi)` (`word(p)` gives the 32 bases from `p` on) for pieces of a read (Myers' bit-vector
/// method): `pieces` gives each piece's length and the edits allowed (second and third
/// field), `peq[k][c]` the positions of piece `k` that hold base `c`. For every run of genome
/// positions where a piece ends within its allowed edits, the run's best position is added
/// to `found` as (edits per 64 bases, piece, position after the match).
fn scan(
    word: &impl Fn(u32) -> u64,
    lo: u32,
    hi: u32,
    pieces: &[(usize, usize, u32)],
    peq: &[[u64; 4]; SCAN_PIECES],
    found: &mut Vec<(u32, usize, u32)>,
) {
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx512f") {
        for g in (0..pieces.len()).step_by(8) {
            // SAFETY: the feature was just detected.
            unsafe { scan_avx512(word, lo, hi, pieces, peq, g, found) };
        }
        return;
    }
    scan_plain(word, lo, hi, pieces, peq, found)
}

fn scan_plain(
    word: &impl Fn(u32) -> u64,
    lo: u32,
    hi: u32,
    pieces: &[(usize, usize, u32)],
    peq: &[[u64; 4]; SCAN_PIECES],
    found: &mut Vec<(u32, usize, u32)>,
) {
    let n = pieces.len();
    let mut top = [0u64; SCAN_PIECES];
    let (mut pv, mut mv) = ([0u64; SCAN_PIECES], [0u64; SCAN_PIECES]);
    let mut score = [0u32; SCAN_PIECES];
    // the best end of the run of matching positions each piece is in: (edits, position)
    let mut run = [(u32::MAX, 0u32); SCAN_PIECES];
    for (k, &(_, len, _)) in pieces.iter().enumerate() {
        top[k] = 1 << (len - 1);
        pv[k] = if len == 64 { !0 } else { (1u64 << len) - 1 };
        score[k] = len as u32;
    }
    let mut p = lo;
    while p < hi {
        let w = word(p);
        let cnt = (hi - p).min(32);
        for j in 0..cnt {
            let c = ((w >> (2 * j)) & 3) as usize;
            for k in 0..n {
                let eq = peq[k][c];
                let xv = eq | mv[k];
                let xh = (((eq & pv[k]).wrapping_add(pv[k])) ^ pv[k]) | eq;
                let mut ph = mv[k] | !(xh | pv[k]);
                let mut mh = pv[k] & xh;
                score[k] += (ph & top[k] != 0) as u32;
                score[k] -= (mh & top[k] != 0) as u32;
                ph <<= 1;
                mh <<= 1;
                pv[k] = mh | !(xv | ph);
                mv[k] = ph & xv;
                if score[k] <= pieces[k].2 {
                    if score[k] < run[k].0 {
                        run[k] = (score[k], p + j + 1);
                    }
                } else if run[k].0 != u32::MAX {
                    found.push((run[k].0 * 64 / pieces[k].1 as u32, k, run[k].1));
                    run[k].0 = u32::MAX;
                }
            }
        }
        p += cnt;
    }
    for k in 0..n {
        if run[k].0 != u32::MAX {
            found.push((run[k].0 * 64 / pieces[k].1 as u32, k, run[k].1));
        }
    }
}

/// [`scan_plain`] for the pieces `g..g + 8`, one piece per 64-bit lane, with the same results.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn scan_avx512(
    word: &impl Fn(u32) -> u64,
    lo: u32,
    hi: u32,
    pieces: &[(usize, usize, u32)],
    peq: &[[u64; 4]; SCAN_PIECES],
    g: usize,
    found: &mut Vec<(u32, usize, u32)>,
) {
    use std::arch::x86_64::*;
    let n = (pieces.len() - g).min(8);
    let mut eqs = [[0u64; 8]; 4];
    let (mut tops, mut pvs, mut scores, mut limits) = ([0u64; 8], [0u64; 8], [1u64 << 32; 8], [0u64; 8]);
    for k in 0..n {
        let (_, len, allowed) = pieces[g + k];
        for c in 0..4 {
            eqs[c][k] = peq[g + k][c];
        }
        tops[k] = 1 << (len - 1);
        pvs[k] = if len == 64 { !0 } else { (1u64 << len) - 1 };
        scores[k] = len as u64;
        limits[k] = allowed as u64;
    }
    // SAFETY (all loads and stores below): the arrays hold eight 64-bit values.
    unsafe {
        let eq = [
            _mm512_loadu_si512(eqs[0].as_ptr().cast()),
            _mm512_loadu_si512(eqs[1].as_ptr().cast()),
            _mm512_loadu_si512(eqs[2].as_ptr().cast()),
            _mm512_loadu_si512(eqs[3].as_ptr().cast()),
        ];
        let top = _mm512_loadu_si512(tops.as_ptr().cast());
        let limit = _mm512_loadu_si512(limits.as_ptr().cast());
        let mut pv = _mm512_loadu_si512(pvs.as_ptr().cast());
        let mut mv = _mm512_setzero_si512();
        let mut score = _mm512_loadu_si512(scores.as_ptr().cast());
        let one = _mm512_set1_epi64(1);
        let mut run = [(u32::MAX, 0u32); 8];
        let mut open = 0u8;
        let mut p = lo;
        while p < hi {
            let w = word(p);
            let cnt = (hi - p).min(32);
            for j in 0..cnt {
                let e = eq[((w >> (2 * j)) & 3) as usize];
                let xv = _mm512_or_si512(e, mv);
                let xh = _mm512_or_si512(
                    _mm512_xor_si512(_mm512_add_epi64(_mm512_and_si512(e, pv), pv), pv),
                    e,
                );
                let ph = _mm512_or_si512(
                    mv,
                    _mm512_andnot_si512(_mm512_or_si512(xh, pv), _mm512_set1_epi64(-1)),
                );
                let mh = _mm512_and_si512(pv, xh);
                score = _mm512_mask_add_epi64(score, _mm512_test_epi64_mask(ph, top), score, one);
                score = _mm512_mask_sub_epi64(score, _mm512_test_epi64_mask(mh, top), score, one);
                let ph = _mm512_slli_epi64::<1>(ph);
                let mh = _mm512_slli_epi64::<1>(mh);
                pv = _mm512_or_si512(
                    mh,
                    _mm512_andnot_si512(_mm512_or_si512(xv, ph), _mm512_set1_epi64(-1)),
                );
                mv = _mm512_and_si512(ph, xv);
                let within = _mm512_cmple_epu64_mask(score, limit);
                if within | open != 0 {
                    let mut sv = [0u64; 8];
                    _mm512_storeu_si512(sv.as_mut_ptr().cast(), score);
                    for k in 0..n {
                        if within >> k & 1 == 1 {
                            if (sv[k] as u32) < run[k].0 {
                                run[k] = (sv[k] as u32, p + j + 1);
                            }
                        } else if run[k].0 != u32::MAX {
                            found.push((run[k].0 * 64 / pieces[g + k].1 as u32, g + k, run[k].1));
                            run[k].0 = u32::MAX;
                        }
                    }
                    open = within;
                }
            }
            p += cnt;
        }
        for k in 0..n {
            if run[k].0 != u32::MAX {
                found.push((run[k].0 * 64 / pieces[g + k].1 as u32, g + k, run[k].1));
            }
        }
    }
}

/// The most a read's usual junction cost per base (times 64) may raise the bar for a doubtful
/// junction.
const USUAL_CAP: u32 = 40;

/// Short-word search of a gap: the longest intron searched, and the read bases the matches on
/// one diagonal must cover.
const WORD_SPAN: u32 = 100_000;
const WORD_COVER: u32 = 14;
/// Approximate search of a gap: the edits allowed in a piece of 32 and of 48 read bases.
const SCAN_EDITS: [u32; 2] = [11, 17];
/// The same at a read end, where what remains of the read may not belong to the transcript.
const SCAN_EDITS_END: [u32; 2] = [11, 17];
/// How far from a read end's last aligned base the approximate search looks.
const SCAN_END_SPAN: u32 = 6_000;
/// By how much a read-end solution through an exon found that way must beat clipping.
const SCAN_END_MARGIN: u32 = 20;

/// By how much a read-end solution with an unknown junction must beat clipping: for usual
/// sites and a short intron, usual sites and a long one, other sites and a short intron,
/// other sites and a long one. `END_NEAR` divides short from long.
const END_MARGIN: [u32; 4] = [2, 16, 16, 32];
/// The margin that no read-end solution meets: for an unknown junction with a site that is of
/// none of the kinds of splice site.
const END_REFUSE: u32 = u32::MAX / 4;
const END_NEAR: u32 = 20_000;

/// Cost of a splice site unknown to the annotation, when one was given (on top of its motif):
/// on a clean read, and the most on a noisy one. Errors next to a junction make an unknown
/// site a few bases from a known one fit as well as it, the more often the noisier the read;
/// on a clean read a new site next to a known one is plain to see and must stay cheap.
const UNKNOWN_SITE: u16 = 5;
const UNKNOWN_SITE_NOISY: u16 = 20;
/// A read counts as clean when its chain's anchors cover this share of the read bases the
/// chain spans (in percent), and as fully noisy when they cover no more than the second.
const CLEAN_COVER: u32 = 88;
const NOISY_COVER: u32 = 63;

thread_local! {
    /// The cost of an unknown site on the read being mapped, and whether it was set for it yet
    /// (the first chain built sets it, so that all of a read's hits are priced alike).
    static UNKNOWN: std::cell::Cell<(u16, bool)> = const { std::cell::Cell::new((UNKNOWN_SITE, false)) };
}

/// Motif cost of an intron's left end at `p` (its first base) and of its right end at `p` (the
/// first exon base after it), for a transcript on the forward or reverse strand. On the
/// reverse strand the left end is the reverse complement of the acceptor.
/// Small cost of a canonical site whose neighbouring base is not the usual one: GT is mostly
/// followed by A or G, and AG mostly preceded by C or T (as minimap2's `--splice-flank`).
const FLANK: u16 = 5;

/// Motif cost of a left end whose first three intron bases are `b0`, `b1`, `b2`.
const fn donor_code(minus: bool, b0: u8, b1: u8, b2: u8) -> u16 {
    match (minus, b0, b1) {
        // GT, then A or G; CT (= AG on the other strand), then A or G (= C or T before the AG)
        (false, 2, 3) | (true, 1, 3) => {
            if matches!(b2, 0 | 2) {
                0
            } else {
                FLANK
            }
        }
        (false, 2, 1) => SEMI,                // GC
        (false, 0, 3) | (true, 2, 3) => ATAC, // AT; GT (= AC)
        _ => NONCANON,
    }
}

/// Motif cost of a right end whose last three intron bases are `b3`, `b2`, `b1` (in genome
/// order).
const fn acceptor_code(minus: bool, b3: u8, b2: u8, b1: u8) -> u16 {
    match (minus, b2, b1) {
        // C or T, then AG; C or T, then AC (= GT on the other strand, followed by A or G)
        (false, 0, 2) | (true, 0, 1) => {
            if matches!(b3, 1 | 3) {
                0
            } else {
                FLANK
            }
        }
        (true, 2, 1) => SEMI,                 // GC (= a GC donor)
        (false, 0, 1) | (true, 0, 3) => ATAC, // AC; AT (= an AT donor)
        _ => NONCANON,
    }
}

/// The two functions above for every three bases (the first lowest, two bits each), per
/// strand.
const fn motif_table(donor: bool) -> [[u16; 64]; 2] {
    let mut t = [[0u16; 64]; 2];
    let mut s = 0;
    while s < 2 {
        let mut c = 0;
        while c < 64 {
            let (x, y, z) = ((c & 3) as u8, ((c >> 2) & 3) as u8, ((c >> 4) & 3) as u8);
            t[s][c] = if donor {
                donor_code(s == 1, x, y, z)
            } else {
                acceptor_code(s == 1, x, y, z)
            };
            c += 1;
        }
        s += 1;
    }
    t
}
static DONOR: [[u16; 64]; 2] = motif_table(true);
static ACCEPTOR: [[u16; 64]; 2] = motif_table(false);

#[inline]
pub(crate) fn donor_motif(ix: &GenomeIndex, p: u32, minus: bool) -> u16 {
    if p + 3 > ix.total_len() || p < 1 {
        return NONCANON;
    }
    donor_code(minus, ix.base(p), ix.base(p + 1), ix.base(p + 2))
}

#[inline]
pub(crate) fn acceptor_motif(ix: &GenomeIndex, p: u32, minus: bool) -> u16 {
    if p < 3 {
        return NONCANON;
    }
    acceptor_code(minus, ix.base(p - 3), ix.base(p - 2), ix.base(p - 1))
}

/// [`donor_pen`] and [`acceptor_pen`] at each of `lo..=hi`, from the bases read once:
/// `bases` receives those of `[lo - 3, hi + 3)`. Returns false, with nothing filled, when
/// the range is too near an end of the genome for that (the caller
/// then asks site by site).
pub(crate) fn site_pens(
    ix: &GenomeIndex,
    lo: u32,
    hi: u32,
    minus: bool,
    bases: &mut Vec<u8>,
    donors: &mut Vec<u32>,
    acceptors: &mut Vec<u32>,
) -> bool {
    if lo < 3 || hi as u64 + 3 > ix.total_len() as u64 {
        return false;
    }
    ix.bases(lo - 3, hi + 3, bases);
    let n = (hi - lo) as usize + 1;
    let (dt, at) = (&DONOR[minus as usize], &ACCEPTOR[minus as usize]);
    donors.clear();
    acceptors.clear();
    // (the three bases before position lo + k are bases[k..k + 3], the three from it on
    // bases[k + 3..k + 6])
    let mut c = (bases[0] | (bases[1] << 2) | (bases[2] << 4)) as usize;
    for k in 0..n {
        acceptors.push(at[c] as u32);
        let d = (bases[k + 3] | (bases[k + 4] << 2) | (bases[k + 5] << 4)) as usize;
        donors.push(dt[d] as u32);
        c = (c >> 2) | ((bases[k + 3] as usize) << 4);
    }
    // With known junctions: a known site costs nothing, any other its motif and the price of
    // an unknown site (as `donor_pen` and `acceptor_pen`).
    if let Some(j) = &ix.junctions {
        let unknown = UNKNOWN.get().0 as u32;
        donors.iter_mut().for_each(|c| *c += unknown);
        acceptors.iter_mut().for_each(|c| *c += unknown);
        for &p in j.sites_in(lo, hi, minus, false) {
            donors[(p - lo) as usize] = 0;
        }
        for &p in j.sites_in(lo, hi, minus, true) {
            acceptors[(p - lo) as usize] = 0;
        }
    }
    true
}

/// Cost of an intron's left end at `p`: nothing for a known site; otherwise the motif cost,
/// and a little more when junctions were given and this is not one of them.
#[inline]
pub(crate) fn donor_pen(ix: &GenomeIndex, p: u32, minus: bool) -> u16 {
    match &ix.junctions {
        Some(j) if j.left_site(p, minus) => 0,
        Some(_) => donor_motif(ix, p, minus) + UNKNOWN.get().0,
        None => donor_motif(ix, p, minus),
    }
}

/// Cost of an intron's right end at `p`, as [`donor_pen`].
#[inline]
pub(crate) fn acceptor_pen(ix: &GenomeIndex, p: u32, minus: bool) -> u16 {
    match &ix.junctions {
        Some(j) if j.right_site(p, minus) => 0,
        Some(_) => acceptor_motif(ix, p, minus) + UNKNOWN.get().0,
        None => acceptor_motif(ix, p, minus),
    }
}

/// A junction placed between two fixed points.
#[derive(Clone, Copy, Debug, Default)]
struct Placed {
    /// Read position where the downstream exon starts.
    q: u32,
    /// Last-exon end (first intron base) and next-exon start on the genome.
    donor: u32,
    acceptor: u32,
    cost: u16,
    /// How the junction was made: `d` placed directly, `p` placed directly although the fit
    /// was poor, `w` by the gap solver through an exon that level-2 cores pointed to, `s` through
    /// one found between splice sites, `e` at a read end.
    origin: u8,
}

/// Known exons shorter than this are tried inside any junction that did not fit exactly.
/// Read ends of up to this many unaligned bases are aligned to the genome next to the hit.
const END_CLOSE: u32 = 16;
/// Past a read end, seed hits of one candidate exon lie at most this many read bases apart, and
/// a single hit counts only this close to where the alignment stopped.
const END_HIT_GAP: u32 = 64;
/// The splice-site search past a read end is for fewer unexplained bases than this.
const SITE_END_MAX: u32 = 25;
/// How far back, at most, a further round of the read-end search starts from where the round
/// before it stopped.
const END_BACK: u32 = 12;
/// The approximate search past a read end is for at least this many unexplained bases.
const SCAN_END_MIN: u32 = 40;
/// Price of a known intron at a read end (besides its sites, which cost nothing).
const KNOWN_END: u32 = 8;
/// The most unexplained read bases at an end for which the known introns are tried on their own
/// (more than that is another exon or more, and the seeds find those).
const KNOWN_END_SPAN: u32 = 60;
const MICRO_EXON: u32 = 20;

/// Reference skips shorter than this are written as deletions, not introns.
const MIN_INTRON_OUT: u32 = 30;
/// Skips shorter than this are reported as deletions in the end: among real reads' skips of
/// 30 to 59 bases next to none is a known intron (they are bases the read lacks). Until a hit
/// is chosen they stay junctions, priced as such.
const MIN_INTRON_REPORTED: u32 = 60;

/// Make blocks consistent for output: in genome order without overlap, each starting on the
/// read where the previous one ends. Returns whether anything had to change.
fn tidy(blocks: &mut Vec<Block>, qlen: u32) -> bool {
    let mut out: Vec<Block> = Vec::with_capacity(blocks.len());
    let mut changed = false;
    for &b in blocks.iter() {
        let mut b = b;
        if b.te <= b.ts {
            changed = true;
            if let Some(p) = out.last_mut() {
                p.qe = p.qe.max(b.qe.min(qlen));
            }
            continue;
        }
        if let Some(p) = out.last_mut() {
            if b.ts < p.te + MIN_INTRON_OUT {
                // overlap on the genome, or a skip too short to be an intron (a deletion in
                // the read): one block
                changed |= b.ts < p.te;
                p.te = p.te.max(b.te);
                p.qe = p.qe.max(b.qe.min(qlen));
                continue;
            }
            if b.qs != p.qe {
                changed = true;
                b.qs = p.qe;
            }
        }
        if b.qe < b.qs {
            changed = true;
            b.qe = b.qs;
        }
        if b.qe > qlen {
            changed = true;
            b.qe = qlen;
            b.qs = b.qs.min(qlen);
        }
        out.push(b);
    }
    *blocks = out;
    changed
}

struct EndSolved {
    cost: u32,
    /// Junctions in the direction of travel (away from the chain).
    hops: Vec<Hop>,
    /// Where the alignment stops (read boundary, genome boundary).
    end: (u32, u32),
    /// The transcript strand, when it was chosen at this end.
    minus: Option<bool>,
}

/// A window as the gap solver takes it: position and length in the solver's view, the read
/// offset that goes with the window's start, and the band's half-width (0 for the usual one).
type SolverWindow = (u32, u32, i64, i64);

/// A candidate exon: genome interval and the first and last seed hit it was built from.
#[derive(Clone, Copy, Debug)]
struct Window {
    s: u32,
    e: u32,
    f: Anchor,
    l: Anchor,
    /// Seed hits behind the window (0 for a window that did not come from seeds).
    hits: u32,
}

impl Window {
    /// For the gap solver travelling forwards from read `q`, genome `t`: (position, length,
    /// read offset that goes with the window's start, 0 for the usual band).
    fn forward(&self, q: u32, t: u32) -> SolverWindow {
        (
            self.s - t,
            self.e - self.s,
            (self.f.q as i64 - q as i64) - (self.f.t - self.s) as i64,
            0,
        )
    }
    /// The same travelling backwards from read `q`, genome `t`.
    fn backward(&self, q: u32, t: u32) -> SolverWindow {
        let (lq, lt) = (self.l.q + self.l.ql as u32, self.l.t + self.l.tl as u32);
        (
            t - self.e,
            self.e - self.s,
            (q as i64 - lq as i64) - (self.e - lt) as i64,
            0,
        )
    }
}

/// Candidate exon windows from level-2 hits in genome order. Hits that
/// continue one another on the same diagonal form one group; the groups with the most hits are
/// kept (at most `max_windows`), since a chance hit is usually alone. A window reaches past its
/// hits by the read bases they leave uncovered on that side (the exon's edge can lie anywhere
/// in them), but not past the midpoint to the next group, and stays inside `(t_lo, t_hi)`.
/// Consecutive hits of a group lie at most `q_gap` read bases apart.
fn windows_of(
    hits: &[Anchor],
    q_lo: u32,
    q_hi: u32,
    t_lo: u32,
    t_hi: u32,
    max_windows: usize,
    q_gap: u32,
) -> Vec<Window> {
    // groups as (first hit, last hit, hits)
    let mut groups: Vec<(Anchor, Anchor, u32)> = Vec::new();
    for &h in hits {
        let mut joined = false;
        for g in groups.iter_mut().rev().take(4) {
            let l = g.1;
            if h.q > l.q && h.t > l.t && h.t - l.t < 2000 && h.q - l.q <= q_gap {
                let d = (h.t - l.t) as i64 - (h.q - l.q) as i64;
                if d.abs() < 30 {
                    g.1 = h;
                    g.2 += 1;
                    joined = true;
                    break;
                }
            }
        }
        if !joined {
            groups.push((h, h, 1));
        }
    }
    let end = |a: Anchor| a.t + a.tl as u32;
    // A group's wish: its hits widened by the read bases they leave uncovered on each side.
    let wish = |g: &(Anchor, Anchor, u32)| {
        let lead = g.0.q.saturating_sub(q_lo).min(150) + 12;
        let trail = q_hi.saturating_sub(g.1.q + g.1.ql as u32).min(150) + 12;
        (g.0.t.saturating_sub(lead), end(g.1) + trail)
    };
    // Strongest groups first; a group whose hits overlap a stronger group's hits is dropped.
    // Windows may overlap one another: the solver takes them in genome order and uses a later
    // one only after it has left an earlier one, so a stray hit costs time, not accuracy.
    let mut by: Vec<usize> = (0..groups.len()).collect();
    by.sort_by_key(|&k| std::cmp::Reverse((groups[k].2, end(groups[k].1) - groups[k].0.t)));
    let mut kept: Vec<(Anchor, Anchor, u32)> = Vec::new();
    for k in by {
        let g = groups[k];
        if kept.iter().any(|o| g.0.t < end(o.1) && o.0.t < end(g.1)) {
            continue;
        }
        kept.push(g);
        if kept.len() == max_windows {
            break;
        }
    }
    kept.sort_by_key(|g| g.0.t);
    let mut out: Vec<Window> = Vec::new();
    for g in &kept {
        let (ws, we) = wish(g);
        let (s, e) = (ws.max(t_lo + 20), we.min(t_hi.saturating_sub(20)));
        if e > s + 4 && e - s <= 2000 {
            out.push(Window {
                s,
                e,
                f: g.0,
                l: g.1,
                hits: g.2,
            });
        }
    }
    out
}

impl Mapper {
    pub fn new(index: std::sync::Arc<GenomeIndex>, opt: MapOptions) -> Self {
        let hpc = index.hpc();
        Mapper {
            index,
            opt,
            seeder: BaseSeeder::new(L3_SPEC, hpc).expect("valid seed spec"),
            seeder2: BaseSeeder::new(L2_SPEC, hpc).expect("valid seed spec"),
        }
    }

    /// Hits of `read` (ASCII bases), best first. Empty when the read does not map.
    pub fn map(&self, read: &[u8], sc: &mut Scratch) -> Vec<SplicedHit> {
        UNKNOWN.set((UNKNOWN_SITE, false));
        let hits = self.map_inner(read, sc);
        sc.stats.solves += std::mem::take(&mut sc.gap.calls);
        sc.stats.prefix_cells += std::mem::take(&mut sc.gap.prefix_cells);
        sc.stats.window_cells += std::mem::take(&mut sc.gap.window_cells);
        hits
    }

    fn map_inner(&self, read: &[u8], sc: &mut Scratch) -> Vec<SplicedHit> {
        let covered = |h: &[SplicedHit]| {
            h.first()
                .map_or(0, |h| h.blocks.iter().map(|b| b.qe - b.qs).sum::<u32>())
        };
        if read.len() < 30 {
            return Vec::new();
        }
        self.seed(read, sc);
        let mut hits = self.map_with(read, self.opt.cap, false, sc);
        // A read that level 3 anchors weakly (no hit, three anchors or fewer, or under 60% of
        // the read covered) is seeded again with pairs of level-2 cores as well.
        let weak = hits.first().is_none_or(|h| h.anchors <= 3) || 10 * covered(&hits) < 6 * read.len() as u32;
        // (a long read with next to none of its level-3 cores in the genome is not from this
        // genome: the second pass would be its most expensive step, for nothing)
        let hopeless = sc.found_best < 2 && read.len() >= 600;
        if weak && !hopeless && self.index.has_pairs() {
            sc.stats.retried += 1;
            // More repetitive seeds are allowed here, but only for a read that has a few of its
            // level-3 cores in the genome: one with next to none is not from this genome, and
            // would only cost time.
            let cap = if 25 * sc.found_best >= sc.l3_best.max(1) && sc.found_best >= 3 {
                self.opt.cap.max(256)
            } else {
                self.opt.cap
            };
            let again = self.map_with(read, cap, true, sc);
            if covered(&again) > covered(&hits) {
                hits = again;
            }
        }
        // Still nothing, and nothing even to chain: allow the repetitive seeds too. (A read whose
        // chains all failed does not get better from more repeats.)
        if !hits.is_empty() || sc.had_hits || self.opt.cap >= self.opt.cap_retry {
            return hits;
        }
        self.map_with(read, self.opt.cap_retry, true, sc)
    }

    /// Seed the read once, in both orientations: level-3 and level-2 cores with raw coordinates.
    fn seed(&self, read: &[u8], sc: &mut Scratch) {
        let ix = &*self.index;
        // Both orientations as codes; other symbols read as A.
        sc.codes[0].clear();
        sc.codes[0].extend(read.iter().map(|&b| code(b) & 3));
        sc.codes[1].clear();
        let fw = std::mem::take(&mut sc.codes[0]);
        sc.codes[1].extend(fw.iter().rev().map(|&c| 3 - c));
        sc.codes[0] = fw;

        for o in 0..2 {
            let ascii: Vec<u8> = sc.codes[o].iter().map(|&c| b"ACGT"[c as usize]).collect();
            // Raw coordinates of a core: its span, the length of its first homopolymer run and
            // the offset of its last one, and (with compression) the number of runs.
            let h = self.seeder.parsed(&ascii);
            let parsed: &[u8] = h.as_ref().map_or(&ascii, |h| &h.seq);
            let coords = |s: u32, e: u32| -> (u32, u32, u16, u16, u16) {
                match &h {
                    Some(h) => {
                        let (qs, qe) = h.to_raw(s, e);
                        (
                            qs,
                            qe,
                            h.run_len(s) as u16,
                            (h.to_raw(e - 1, e).0 - qs) as u16,
                            (e - s) as u16,
                        )
                    }
                    None => {
                        let (mut a, mut b) = (s + 1, e - 1);
                        while a < e && ascii[a as usize] == ascii[s as usize] {
                            a += 1;
                        }
                        while b > s && ascii[b as usize - 1] == ascii[e as usize - 1] {
                            b -= 1;
                        }
                        (s, e, (a - s) as u16, (b - s) as u16, 0)
                    }
                }
            };
            sc.l3[o].clear();
            for s in self.seeder.seed_parsed(parsed) {
                let (qs, qe, q_in, q_out, runs) = coords(s.start, s.end);
                sc.l3[o].push(L2Seed {
                    q: qs,
                    ql: (qe - qs).min(u16::MAX as u32) as u16,
                    q_in,
                    q_out,
                    runs,
                    id: 0,
                    label: s.label,
                    occ: (0, 0),
                });
            }
            sc.l2[o].clear();
            for s in self.seeder2.seed_parsed(parsed) {
                if let Some(id) = ix.l2_id(s.label) {
                    let (qs, qe, q_in, q_out, runs) = coords(s.start, s.end);
                    sc.l2[o].push(L2Seed {
                        q: qs,
                        ql: (qe - qs) as u16,
                        q_in,
                        q_out,
                        runs,
                        id,
                        label: s.label,
                        occ: (0, 0),
                    });
                }
            }
        }
        // Every level-3 core is looked up once per read, all of them a step at a time (see
        // `GenomeIndex::prefetch_dir`); the passes over the read share the result.
        for o in 0..2 {
            for s in &sc.l3[o] {
                ix.prefetch_dir(s.label, false);
            }
        }
        for o in 0..2 {
            for s in &sc.l3[o] {
                ix.prefetch_bucket(s.label, false);
            }
        }
        for o in 0..2 {
            for s in sc.l3[o].iter_mut() {
                s.occ = ix.occ_range(s.label, false);
            }
        }
        sc.pair_seeds[0].clear();
        sc.pair_seeds[1].clear();
        sc.pairs_done = [false; 2];
    }

    /// The read's pairs of level-2 cores in orientation `o` that are in the index, looked up
    /// once per read like the level-3 cores.
    fn pair_seeds(&self, o: usize, sc: &mut Scratch) {
        if std::mem::replace(&mut sc.pairs_done[o], true) {
            return;
        }
        let ix = &*self.index;
        for w in sc.l2[o].windows(2) {
            let l = crate::index::pair_label(w[0].label, w[1].label);
            if w[1].q - w[0].q >= 64 || !crate::index::keep_pair(l) {
                continue;
            }
            ix.prefetch_dir(l, true);
            let ql = w[1].q + w[1].ql as u32 - w[0].q;
            sc.pair_seeds[o].push(L2Seed {
                q: w[0].q,
                ql: ql as u16,
                q_in: w[0].q_in,
                q_out: (w[1].q - w[0].q) as u16 + w[1].q_out,
                runs: 0,
                id: 0,
                label: l,
                occ: (0, 0),
            });
        }
        for s in &sc.pair_seeds[o] {
            ix.prefetch_bucket(s.label, true);
        }
        for s in sc.pair_seeds[o].iter_mut() {
            s.occ = ix.occ_range(s.label, true);
        }
    }

    fn map_with(&self, read: &[u8], cap: usize, pairs: bool, sc: &mut Scratch) -> Vec<SplicedHit> {
        let ix = &*self.index;
        sc.gap.reset();
        let n = read.len();
        let mut out: Vec<SplicedHit> = Vec::new();
        let mut cands: Vec<(f32, bool, Vec<Anchor>)> = Vec::new();
        sc.had_hits = false;
        // Level-3 hits of both orientations first: which one the read maps in shows in them.
        let mut l3_hits: [Vec<Anchor>; 2] = [Vec::new(), Vec::new()];
        let mut found = [0usize; 2];
        for o in 0..2 {
            for s in &sc.l3[o] {
                let (pos, len) = ix.occ(s.occ, false);
                found[o] += !pos.is_empty() as usize;
                if pos.is_empty() || pos.len() > cap {
                    continue;
                }
                let w = (s.ql as u32).min(40) as f32 / (1.0 + (pos.len() as f32).ln());
                for (&t, &tl) in pos.iter().zip(len) {
                    l3_hits[o].push(Anchor {
                        t,
                        q: s.q,
                        ql: s.ql,
                        tl: tl as u16,
                        q_in: s.q_in,
                        q_out: s.q_out,
                        w,
                    });
                }
            }
        }
        for o in 0..2 {
            sc.hits.clear();
            sc.hits.append(&mut l3_hits[o]);
            // Pairs of level-2 cores as seeds too, for a noisy read: one with under a quarter of
            // its level-3 cores in the genome in its better orientation (about 7% error or
            // more). They keep the stretches between anchors short. An orientation that the
            // other one clearly beats is not the one the read maps in and is spared the lookups.
            let best = if found[0] >= found[1] { 0 } else { 1 };
            sc.found_best = found[best];
            sc.l3_best = sc.l3[best].len();
            let noisy = 3 * found[best] < sc.l3[best].len()
                && found[best] >= 3
                && 25 * found[best] >= sc.l3[best].len();
            let beaten = found[1 - o] >= 5 && 4 * found[o] <= found[1 - o];
            if (pairs || (self.opt.pairs && noisy && !beaten)) && ix.has_pairs() {
                self.pair_seeds(o, sc);
                for s in &sc.pair_seeds[o] {
                    let (pos, len) = ix.occ(s.occ, true);
                    if pos.is_empty() || pos.len() > cap {
                        continue;
                    }
                    let wt = (s.ql as u32).min(40) as f32 / (1.0 + (pos.len() as f32).ln());
                    for (&t, &tl) in pos.iter().zip(len) {
                        sc.hits.push(Anchor {
                            t,
                            q: s.q,
                            ql: s.ql,
                            tl: tl as u16,
                            q_in: s.q_in,
                            q_out: s.q_out,
                            w: wt,
                        });
                    }
                }
            }
            if sc.hits.is_empty() {
                continue;
            }
            sc.had_hits = true;
            sort_hits(&mut sc.hits, &mut sc.hit_keys, &mut sc.hits_sorted);
            // Loci: runs of hits no further apart than the longest intron.
            let mut loci: Vec<(f32, usize, usize)> = Vec::new();
            let mut s = 0;
            for i in 1..=sc.hits.len() {
                if i == sc.hits.len() || sc.hits[i].t - sc.hits[i - 1].t > self.opt.max_intron {
                    let w: f32 = sc.hits[s..i].iter().map(|a| a.w).sum();
                    loci.push((w, s, i));
                    s = i;
                }
            }
            let best = loci.iter().map(|l| l.0).fold(0.0, f32::max);
            // (a read whose best locus is one or two seeds deep cannot rank its loci by
            // weight: all of them are looked at, and screened below)
            // (only for a read with a few of its level-3 cores in the genome; one with next to
            // none is not from this genome)
            let of_genome = sc.found_best >= 3 && 25 * sc.found_best >= sc.l3_best.max(1);
            let limit = if best >= 80.0 {
                self.opt.max_loci
            } else if of_genome {
                400
            } else {
                48
            };
            // The heaviest `limit` loci among those within a quarter of the best, heaviest
            // first (a read of chance hits has thousands of loci of one hit each: they are
            // not all sorted).
            loci.retain(|l| l.0 >= 0.25 * best);
            if loci.len() > limit {
                loci.select_nth_unstable_by(limit, |a, b| b.0.total_cmp(&a.0));
                loci.truncate(limit);
            }
            loci.sort_unstable_by(|a, b| b.0.total_cmp(&a.0));
            for &(_, s, e) in loci.iter() {
                if let Some((score, chain)) = self.chain_locus(s, e, sc) {
                    cands.push((score, o == 1, chain));
                }
            }
        }
        // Only chains close to the best are worth building.
        cands.sort_by(|a, b| b.0.total_cmp(&a.0));
        // When even the best chain is one or two anchors, chance hits score like the true one.
        // Each candidate is then screened by whether the read goes on matching next to its
        // anchors (two word comparisons a side), and the best few of those are built.
        let mut screened_ok = false;
        if cands.first().is_some_and(|c| c.2.len() <= 2) && cands.len() > 1 {
            screened_ok = true;
            let mut screened: Vec<(u32, f32, bool, Vec<Anchor>)> = Vec::new();
            for (score, reverse, chain) in cands.drain(..) {
                let read = &sc.codes[reverse as usize];
                let (f, l) = (chain[0], chain[chain.len() - 1]);
                let (lq, lt) = (l.q + l.ql as u32, l.t + l.tl as u32);
                let ml = self.diagonal_matches(read, f.q, f.t, true);
                let mr = self.diagonal_matches(read, lq, lt, false);
                let m = ml + mr;
                // (chance gives 6-10 a side; a true anchor of a read with 10% errors about 14-20;
                // an anchor at a read end has one side only)
                if m >= 24 || ml.max(mr) >= 14 {
                    screened.push((m, score, reverse, chain));
                }
            }
            screened.sort_by(|a, b| {
                (b.0, b.1)
                    .partial_cmp(&(a.0, a.1))
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            screened.truncate(16);
            // Among those, the ones the read really continues from: each is extended a little
            // way to both sides (the extension stops as soon as the score drops), and the
            // candidates are ranked by how far that got. A chance hit gets a few bases.
            let mut reached: Vec<(u32, bool, Vec<Anchor>)> = Vec::new();
            for (_, _, reverse, chain) in screened {
                let read = std::mem::take(&mut sc.codes[reverse as usize]);
                let (f, l) = (chain[0], chain[chain.len() - 1]);
                let far = self.extend(&read, f.q, f.t, true, sc).0
                    + self
                        .extend(&read, l.q + l.ql as u32, l.t + l.tl as u32, false, sc)
                        .0;
                sc.codes[reverse as usize] = read;
                if far >= 30 {
                    reached.push((far, reverse, chain));
                }
            }
            reached.sort_by_key(|r| std::cmp::Reverse(r.0));
            cands = reached
                .into_iter()
                .take(4)
                .map(|(far, r, c)| (far as f32, r, c))
                .collect();
        }
        let top = cands.first().map_or(0.0, |c| c.0);
        // The best eight chains; and, past those, up to two that cross an intron: a gene with
        // many processed pseudogenes has more unspliced look-alikes than that, and each scores
        // above the gene itself, whose chain pays for its introns.
        let crosses = |chain: &[Anchor]| {
            chain
                .windows(2)
                .any(|w| (w[1].t - w[0].t) as i64 - (w[1].q - w[0].q) as i64 >= self.opt.min_intron as i64)
        };
        let mut more = 2;
        for (k, (score, reverse, chain)) in cands.into_iter().enumerate() {
            if score < 0.7 * top || (k >= 8 && more == 0) {
                break;
            }
            if k >= 8 {
                if !crosses(&chain) {
                    continue;
                }
                more -= 1;
            }
            let read = std::mem::take(&mut sc.codes[reverse as usize]);
            sc.gap.reset(); // another chain, possibly the other orientation of the read
            out.push(self.build(reverse, &chain, &read, n as u32, score, screened_ok, sc));
            sc.codes[reverse as usize] = read;
        }
        // Rival hits, and a lone hit resting on few anchors or with more edits found between its
        // anchors than they are worth, are scored by aligning their blocks; a lone well-anchored
        // hit keeps its chain score and is never aligned inside its exons.
        if out.len() > 1 || out.first().is_some_and(|h| h.anchors <= 3 || h.score <= 0) {
            for h in out.iter_mut() {
                let read = std::mem::take(&mut sc.codes[h.reverse as usize]);
                h.score = self.align_score(&read, h, sc);
                sc.codes[h.reverse as usize] = read;
                if std::env::var_os("GMAP_DEBUG").is_some() {
                    eprintln!("scored {} blocks={:?}", h.score, h.blocks);
                }
                sc.stats.scored += 1;
            }
            out.retain(|h| h.score >= self.opt.min_score);
            // A spliced read also aligns, unspliced, to a processed pseudogene of its gene, and
            // there it can cover a few more bases: no read ends lost at junctions, and the
            // read's poly(A) tail matches the one the pseudogene was copied with. When a
            // spliced hit is within 10% of the best unspliced one and aligns at least as well
            // (per base; or, where the unspliced hit covers all the spliced one does, in what
            // its alignment costs in all), it is the gene the read came from.
            let cover = |h: &SplicedHit| h.blocks.iter().map(|b| b.qe - b.qs).sum::<u32>().max(1) as i64;
            // (the read interval of a hit, on the read as given)
            let span = |h: &SplicedHit| {
                let (a, b) = (h.blocks[0].qs, h.blocks[h.blocks.len() - 1].qe);
                if h.reverse {
                    (n as u32 - b, n as u32 - a)
                } else {
                    (a, b)
                }
            };
            // (blocks a short gap apart are one exon with a deletion)
            let introns = |h: &SplicedHit| {
                h.blocks
                    .windows(2)
                    .filter(|w| w[1].ts - w[0].te >= MIN_INTRON_REPORTED)
                    .count()
            };
            if let Some(u) = out.iter().filter(|h| introns(h) == 0).max_by_key(|h| h.score) {
                let (us, uc, (ua, ub)) = (u.score as i64, cover(u), span(u));
                for h in out.iter_mut().filter(|h| introns(h) > 0) {
                    let (hs, hc, (ha, hb)) = (h.score as i64, cover(h), span(h));
                    let within = ua <= ha + 3 && hb <= ub + 3;
                    let as_well = hs * uc >= us * hc || (within && 2 * hc - hs <= 2 * uc - us);
                    if hs <= us && 10 * hs >= 9 * us && as_well {
                        h.score = us as i32 + 1;
                        sc.stats.spliced_preferred += 1;
                    }
                }
            }
        }
        out.sort_by_key(|h| std::cmp::Reverse(h.score));
        if let Some(h) = out.first_mut() {
            h.primary = true;
        }
        for h in out.iter_mut() {
            let mut k = 1;
            while k < h.blocks.len() {
                if h.blocks[k].ts < h.blocks[k - 1].te + MIN_INTRON_REPORTED {
                    let b = h.blocks.remove(k);
                    h.blocks[k - 1].te = b.te;
                    h.blocks[k - 1].qe = b.qe;
                    if k - 1 < h.origins.len() {
                        h.origins.remove(k - 1);
                    }
                } else {
                    k += 1;
                }
            }
        }
        out
    }

    /// Chain the hits `sc.hits[s..e]` (sorted by genome position): the best chain and its score.
    fn chain_locus(&self, s: usize, e: usize, sc: &mut Scratch) -> Option<(f32, Vec<Anchor>)> {
        const H: usize = 40;
        let hits = &sc.hits[s..e];
        let m = hits.len();
        sc.f.clear();
        sc.back.clear();
        let (mut best, mut best_i) = (f32::MIN, 0usize);
        for j in 0..m {
            let a = hits[j];
            let mut fj = a.ql as f32;
            let mut bj = u32::MAX;
            for i in (j.saturating_sub(H)..j).rev() {
                let p = hits[i];
                if p.q >= a.q || p.t >= a.t {
                    continue;
                }
                let v = sc.f[i] + self.step(p, a);
                if v > fj {
                    fj = v;
                    bj = i as u32;
                }
            }
            sc.f.push(fj);
            sc.back.push(bj);
            if fj > best {
                best = fj;
                best_i = j;
            }
        }
        sc.chain.clear();
        let mut i = best_i as u32;
        while i != u32::MAX {
            sc.chain.push(i);
            i = sc.back[i as usize];
        }
        sc.chain.reverse();
        Some((best, sc.chain.iter().map(|&i| hits[i as usize]).collect()))
    }

    /// Anchors for long stretches of the read that a chain leaves without one (several exons
    /// in a burst of errors): the level-2 hits between the two chain neighbours are chained,
    /// and those that another hit of that chain continues on the same diagonal are added to
    /// the chain. A chance hit is alone on its diagonal; two that agree are not chance. The
    /// stretch is then several shorter gaps, each within the solver's reach.
    fn refine(&self, reverse: bool, chain: &[Anchor], qlen: u32, sc: &mut Scratch) -> Vec<Anchor> {
        const LONG: u32 = 200;
        const NEAR: i64 = 16;
        const RUN: usize = 3;
        let min_intron = self.opt.min_intron as i64;
        let mut out: Vec<Anchor> = Vec::with_capacity(chain.len());
        let (mut f, mut back, mut path): (Vec<f32>, Vec<u32>, Vec<usize>) = Default::default();
        for (k, &a) in chain.iter().enumerate() {
            out.push(a);
            let Some(&b) = chain.get(k + 1) else { break };
            let diag = |x: Anchor, y: Anchor| (y.t as i64 - x.t as i64) - (y.q as i64 - x.q as i64);
            if b.q < a.q + a.ql as u32 + LONG || diag(a, b) < min_intron {
                continue;
            }
            let mut hits = self.rescue(reverse, Some(a), Some(b), 0, qlen, sc);
            // only hits that can follow `a` and precede `b`
            hits.retain(|&h| diag(a, h) >= -NEAR && diag(h, b) >= -NEAR);
            if hits.len() < 2 {
                continue;
            }
            f.clear();
            back.clear();
            let (mut best, mut best_j) = (f32::MIN, 0usize);
            for j in 0..hits.len() {
                let h = hits[j];
                let (mut fj, mut bj) = (h.ql as f32, u32::MAX);
                for i in (j.saturating_sub(64)..j).rev() {
                    let p = hits[i];
                    if p.q + p.ql as u32 > h.q || p.t + p.tl as u32 > h.t {
                        continue;
                    }
                    let d = diag(p, h);
                    let cost = if d >= min_intron {
                        8.0
                    } else if d.abs() <= NEAR + (h.q - p.q) as i64 / 8 {
                        0.5 * d.abs() as f32
                    } else {
                        continue;
                    };
                    if f[i] + h.ql as f32 - cost > fj {
                        fj = f[i] + h.ql as f32 - cost;
                        bj = i as u32;
                    }
                }
                f.push(fj);
                back.push(bj);
                if fj > best {
                    (best, best_j) = (fj, j);
                }
            }
            path.clear();
            let mut i = best_j as u32;
            while i != u32::MAX {
                path.push(i as usize);
                i = back[i as usize];
            }
            path.reverse();
            let agree =
                |x: usize, y: usize| diag(hits[x], hits[y]).abs() <= NEAR && hits[y].q - hits[x].q <= 200;
            // runs of hits that continue one another; a run of fewer than `RUN` is left out
            let mut s0 = 0;
            for n in 1..=path.len() {
                if n == path.len() || !agree(path[n - 1], path[n]) {
                    if n - s0 >= RUN {
                        out.extend(path[s0..n].iter().map(|&x| hits[x]));
                    }
                    s0 = n;
                }
            }
        }
        out
    }

    /// Blocks from a chain: junctions placed, exons without an anchor recovered, read ends solved.
    #[allow(clippy::too_many_arguments)]
    fn build(
        &self,
        reverse: bool,
        chain: &[Anchor],
        read: &[u8],
        qlen: u32,
        score: f32,
        screened: bool,
        sc: &mut Scratch,
    ) -> SplicedHit {
        // A chain is solid enough to look for more exons from when it has three anchors or more
        // and they cover a fortieth of the read (50 bases at least): a few chance hits strung
        // together across huge skips are not worth the solver.
        let anchored: u32 = chain.iter().map(|a| a.ql as u32).sum();
        // (a thin chain that passed the screen of thin candidates is trusted the same way: the
        // read was seen to go on matching next to its anchors)
        let solid = screened || (chain.len() >= 3 && anchored >= (qlen / 40).max(50));
        let refined;
        let chain = if solid && chain.windows(2).any(|w| w[1].q >= w[0].q + 200) {
            refined = self.refine(reverse, chain, qlen, sc);
            &refined[..]
        } else {
            chain
        };
        if self.index.junctions.is_some() && !UNKNOWN.get().1 {
            let (mut cover, mut end) = (0u32, 0u32);
            for a in chain {
                let e = a.q + a.ql as u32;
                cover += e.saturating_sub(a.q.max(end));
                end = end.max(e);
            }
            let span = end.saturating_sub(chain[0].q).max(1);
            let pct = (100 * cover / span).clamp(NOISY_COVER, CLEAN_COVER);
            let price = UNKNOWN_SITE
                + ((CLEAN_COVER - pct) * (UNKNOWN_SITE_NOISY - UNKNOWN_SITE) as u32
                    / (CLEAN_COVER - NOISY_COVER)) as u16;
            UNKNOWN.set((price, true));
        }
        if std::env::var_os("GMAP_DEBUG").is_some() {
            let v: Vec<_> = chain.iter().map(|a| (a.q, a.ql, a.t, a.w)).collect();
            eprintln!("chain reverse={reverse} score={score} solid={solid} {v:?}");
        }
        // Each pair of chain neighbours that the genome separates by an intron, placed directly.
        let mut pairs: Vec<(Anchor, Anchor, [Placed; 2], bool)> = Vec::new();
        // the prefix tables of the junctions that may go to the solver, to not compute them twice
        let mut tables: Vec<Option<gap::Tables>> = Vec::new();
        for w in chain.windows(2) {
            let (a, b) = (w[0], w[1]);
            let d = (b.t - a.t) as i64 - (b.q - a.q) as i64;
            if d < self.opt.min_intron as i64 {
                continue;
            }
            sc.stats.junctions += 1;
            let (q0, t0, q1, t1) = self.fixed(a, b);
            let p = self.place(read, q0, t0, q1, t1, sc);
            let c = p[0].cost.min(p[1].cost);
            let doubtful = c < INF && c as u32 >= 18 + q1.saturating_sub(q0) * self.opt.trigger_x64 / 64;
            tables.push((solid && doubtful).then(|| sc.gap.stash()));
            pairs.push((a, b, p, false));
        }
        // Two levels of doubt about a junction. One that costs more than a clean read's would
        // gets the candidates backed by two or more seed hits. Only one that costs clearly more
        // than this read's junctions usually do (the median cost per base follows the read's
        // error rate) also gets the lone hits and the seedless search, so that a noisy read
        // does not send every junction through all of it.
        let per_base = |x: &(Anchor, Anchor, [Placed; 2], bool)| {
            let (q0, _, q1, _) = self.fixed(x.0, x.1);
            let c = x.2[0].cost.min(x.2[1].cost);
            (c < INF && q1 > q0).then(|| 64 * c as u32 / (q1 - q0))
        };
        let mut rates: Vec<u32> = pairs.iter().filter_map(per_base).collect();
        rates.sort_unstable();
        // (capped: on a read whose junctions are mostly bad the median says nothing about a
        // good one, and with no cap nothing on such a read would count as bad)
        let usual = if rates.len() >= 3 {
            rates[rates.len() / 2].min(USUAL_CAP)
        } else {
            0
        };
        let rate = (self.opt.trigger_x64).max(usual * 3 / 2); // cost per base, times 64
        let mut vote = [0u32; 2];
        for x in pairs.iter_mut() {
            let (q0, _, q1, _) = self.fixed(x.0, x.1);
            let (c, m) = (x.2[0].cost.min(x.2[1].cost), q1.saturating_sub(q0));
            x.3 = solid && (c >= INF || c as u32 >= 18 + m * self.opt.trigger_x64 / 64);
            // every placed junction votes on the strand (a noisy read's junctions are mostly
            // doubtful, and leaving them out left the strand undecided and the work doubled)
            if x.2[0].cost < INF && x.2[1].cost < INF {
                vote[0] += x.2[0].cost as u32;
                vote[1] += x.2[1].cost as u32;
            }
        }
        // The transcript strand, when the well-placed junctions agree on one.
        let strands: Vec<bool> = if self.opt.forward_only {
            vec![reverse]
        } else if vote[0] + 8 <= vote[1] {
            vec![false]
        } else if vote[1] + 8 <= vote[0] {
            vec![true]
        } else {
            vec![false, true]
        };
        // Where a junction leaves read bases unexplained, look for level-2 cores of those bases
        // inside that intron and solve the gap through the exons they point to.
        let mut cand: Vec<Option<Vec<Window>>> = vec![None; pairs.len()];
        let mut seedless: Vec<Option<Vec<SolverWindow>>> = vec![None; pairs.len()];
        let mut best: Option<(u32, bool, Vec<Placed>)> = None;
        for &minus in &strands {
            let mut total = 0u32;
            let mut juncs: Vec<Placed> = Vec::new();
            for (k, &(a, b, p, poor)) in pairs.iter().enumerate() {
                let direct = p[minus as usize];
                let (q0, t0, q1, t1) = self.fixed(a, b);
                if !poor {
                    // A junction that fits well enough may still hide a known exon of a few
                    // bases (written as a short insertion, it costs too little to raise doubt):
                    // with known junctions, try the short known exons inside the intron.
                    if direct.cost > 0 && direct.cost < INF && q1 > q0 && t1 > t0 {
                        let m = (q1 - q0) as usize;
                        let mut vw = self.known_exons(t0, t1, m);
                        vw.retain(|w| w.1 < MICRO_EXON);
                        if !vw.is_empty() {
                            vw.sort_by_key(|w| w.0);
                            sc.gap.known_only = true;
                            let t = gap::solve(
                                &self.index,
                                read,
                                ReadView { q: q0, rev: false },
                                m,
                                View {
                                    origin: t0,
                                    rev: false,
                                },
                                Some(t1),
                                &vw,
                                minus,
                                self.opt.min_intron / 2,
                                &mut sc.gap,
                            );
                            sc.gap.known_only = false;
                            if std::env::var_os("GMAP_DEBUG").is_some() {
                                eprintln!(
                                    "  micro q0={q0} t0={t0} m={m} windows={vw:?} cost={} hops={:?} direct={direct:?}",
                                    t.cost, t.hops
                                );
                            }
                            if t.hops.len() > 1
                                && t.cost < direct.cost as u32
                                && self.hop_margin(&t, minus) <= 4
                            {
                                total += t.cost;
                                juncs.extend(t.hops.iter().map(|h| Placed {
                                    q: h.q,
                                    donor: h.leave,
                                    acceptor: h.enter,
                                    cost: 0,
                                    origin: b'k',
                                }));
                                continue;
                            }
                        }
                    }
                    total += direct.cost as u32;
                    juncs.push(direct);
                    continue;
                }
                if let Some(t) = &tables[k] {
                    sc.gap.restore(t);
                }
                let hard = {
                    let c = p[0].cost.min(p[1].cost);
                    c >= INF || c as u32 >= 18 + q1.saturating_sub(q0) * rate / 64
                };
                if cand[k].is_none() {
                    sc.stats.rescues_tried += 1;
                    let path = self.rescue(reverse, Some(a), Some(b), 0, qlen, sc);
                    sc.stats.empty += path.is_empty() as u64;
                    cand[k] = Some(windows_of(&path, q0, q1, t0, t1, 8, u32::MAX));
                }
                let wins = cand[k].as_ref().unwrap();
                let m = q1.saturating_sub(q0) as usize;
                let mut done = false;
                // the solution through candidate exons so far, and how its exons were found
                let mut sol: Option<(Solved, u8)> = None;
                // what a junction of this length usually costs on this read: an answer above it
                // still leaves something unexplained
                let fine = 18 + m as u32 * self.opt.trigger_x64.max(usual) / 64;
                if (!wins.is_empty() || self.index.junctions.is_some()) && q1 > q0 && t1 > t0 && m <= 600 {
                    let rv = ReadView { q: q0, rev: false };
                    let gv = View {
                        origin: t0,
                        rev: false,
                    };
                    let half = self.opt.min_intron / 2;
                    let known = self.known_exons(t0, t1, m);
                    // Candidates backed by two or more hits first; the lone hits (mostly
                    // strays) only if those explain nothing.
                    let strong = wins.iter().filter(|w| w.hits >= 2).count();
                    let mut s = Solved {
                        cost: u32::MAX,
                        ..Default::default()
                    };
                    // Then, if junctions were given, the known exons inside this intron on their
                    // own, believed only when every junction of the solution is a known one.
                    let mut found = false;
                    for (min_hits, from_known) in [(2u32, false), (0, false), (0, true)] {
                        // the lone hits are tried when the junction is clearly bad, and also when
                        // the well-supported candidates gave an answer that is still a poor fit
                        // (a second exon may be missing in the same gap)
                        let lone_too = hard || (found && s.cost >= fine);
                        if (from_known && (known.is_empty() || found))
                            || (!from_known && wins.is_empty())
                            || (!from_known && min_hits == 2 && strong == 0)
                            || (!from_known && min_hits == 0 && (!lone_too || strong == wins.len()))
                        {
                            continue;
                        }
                        let mut vw: Vec<SolverWindow> = if from_known {
                            known.clone()
                        } else {
                            wins.iter()
                                .filter(|w| w.hits >= min_hits)
                                .map(|w| w.forward(q0, t0))
                                .collect()
                        };
                        // (known exons inside the intron are candidates in every round: where
                        // the seed hits point to one of two missing exons, the other is often
                        // a known one)
                        if !from_known {
                            // (but not one that a seed window already covers)
                            let seeded = vw.len();
                            for k in &known {
                                if !vw[..seeded].iter().any(|w| w.0 <= k.0 && k.0 + k.1 <= w.0 + w.1) {
                                    vw.push(*k);
                                }
                            }
                        }
                        vw.sort_by_key(|w| w.0);
                        sc.gap.known_only = from_known;
                        let t = gap::solve(
                            &self.index,
                            read,
                            rv,
                            m,
                            gv,
                            Some(t1),
                            &vw,
                            minus,
                            half,
                            &mut sc.gap,
                        );
                        sc.gap.known_only = false;
                        let limit = if found {
                            s.cost
                        } else {
                            (direct.cost as u32).min(1 << 20)
                        };
                        if std::env::var_os("GMAP_DEBUG").is_some() {
                            eprintln!(
                                "  try min_hits={min_hits} known={from_known} windows={vw:?} cost={} hops={:?} limit={limit} fine={fine} hard={hard}",
                                t.cost, t.hops
                            );
                        }
                        if t.hops.len() > 1
                            && t.cost < limit
                            && (!from_known || self.hop_margin(&t, minus) <= 4)
                        {
                            found = true;
                            s = t;
                            if s.cost < fine {
                                break;
                            }
                        }
                    }
                    if found {
                        sc.stats.rescues_kept += (minus == strands[0]) as u64;
                        sol = Some((s, b'w'));
                    }
                }
                // No level-2 hit led to an explanation, or only to one that is still a poor fit
                // (a lone hit can point next to an exon, and some stretch near it will fit the
                // read's bases better than nothing): look for an exon between a pair of splice
                // sites that matches the unexplained bases, and for one that shares short words
                // with them.
                let poor_fit = sol.as_ref().is_none_or(|x| x.0.cost >= 2 * fine);
                if poor_fit && hard && q1 > q0 && t1 > t0 && m <= 600 {
                    let (sites, mid) = self.site_windows(read, q0, t0, q1, t1, minus, sc);
                    // (the two searches below do not depend on the strand)
                    if seedless[k].is_none() {
                        let mut w = self.word_windows(read, q0, t0, q1, t1, sc);
                        if let Some(mid) = mid {
                            w.extend(self.scan_windows(read, q0, t0, q1, t1, mid));
                        }
                        seedless[k] = Some(w);
                    }
                    let words = seedless[k].clone().unwrap();
                    if !sites.is_empty() || !words.is_empty() {
                        sc.stats.site_tried += (minus == strands[0]) as u64;
                        let mut vw: Vec<SolverWindow> = sites.iter().map(|w| w.forward(q0, t0)).collect();
                        let by_words = !words.is_empty();
                        vw.extend(words);
                        if sol.is_some() {
                            vw.extend(wins.iter().map(|w| w.forward(q0, t0)));
                        }
                        vw.sort_by_key(|w| w.0);
                        let rv = ReadView { q: q0, rev: false };
                        let gv = View {
                            origin: t0,
                            rev: false,
                        };
                        let half = self.opt.min_intron / 2;
                        let s = gap::solve(
                            &self.index,
                            read,
                            rv,
                            m,
                            gv,
                            Some(t1),
                            &vw,
                            minus,
                            half,
                            &mut sc.gap,
                        );
                        let limit = sol
                            .as_ref()
                            .map_or((direct.cost as u32).min(1 << 20), |x| x.0.cost);
                        if std::env::var_os("GMAP_DEBUG").is_some() {
                            eprintln!(
                                "  seedless windows={vw:?} cost={} hops={:?} limit={limit}",
                                s.cost, s.hops
                            );
                        }
                        if s.hops.len() > 1 && s.cost < limit {
                            sc.stats.site_kept += (minus == strands[0]) as u64;
                            sol = Some((s, if by_words { b'k' } else { b's' }));
                        }
                    }
                }
                if let Some((s, origin)) = sol {
                    total += s.cost;
                    juncs.extend(s.hops.iter().map(|h| Placed {
                        q: h.q,
                        donor: h.leave,
                        acceptor: h.enter,
                        cost: 0,
                        origin,
                    }));
                    done = true;
                }
                if std::env::var_os("GMAP_DEBUG").is_some() {
                    eprintln!(
                        "gap q0={q0} q1={q1} t0={t0} t1={t1} minus={minus} direct={} windows={wins:?} done={done}",
                        direct.cost
                    );
                }
                if !done {
                    sc.stats.not_cheaper += (minus == strands[0]) as u64;
                    total += direct.cost.min(2000) as u32;
                    juncs.push(Placed {
                        origin: b'p',
                        ..direct
                    });
                }
            }
            if best.as_ref().is_none_or(|b| total < b.0) {
                best = Some((total, minus, juncs));
            }
        }
        let (mut edits, minus, juncs) = best.unwrap();
        let undecided = !self.opt.forward_only && juncs.is_empty();
        let tx_strand =
            if undecided || (strands.len() == 2 && pairs.iter().all(|p| p.2[0].cost == p.2[1].cost)) {
                b'.'
            } else if minus {
                b'-'
            } else {
                b'+'
            };
        let mut blocks: Vec<Block> = Vec::new();
        let first = chain[0];
        let (mut bq, mut bt) = (first.q, first.t);
        let mut origins: Vec<u8> = juncs.iter().map(|p| p.origin).collect();
        for p in &juncs {
            blocks.push(Block {
                qs: bq,
                qe: p.q,
                ts: bt,
                te: p.donor,
            });
            bq = p.q;
            bt = p.acceptor;
        }
        let last = chain[chain.len() - 1];
        blocks.push(Block {
            qs: bq,
            qe: last.q + last.ql as u32,
            ts: bt,
            te: last.t + last.tl as u32,
        });

        // Read ends, from a point inside the outermost anchors.
        // (A seed can run a few bases past an exon's edge when the intron begins as the next
        // exon does, so each end is solved from up to `BACK` bases further in, as far as the
        // read matches the genome base for base and the outermost block reaches.)
        const BACK: u32 = 8;
        let (_, _, mut q1, mut t1) = self.fixed(last, last);
        let rescue = solid;
        if rescue {
            let bl = blocks[blocks.len() - 1];
            let (rv, gv) = (
                ReadView { q: q1, rev: true },
                View {
                    origin: t1,
                    rev: true,
                },
            );
            let room = BACK
                .min(q1.saturating_sub(bl.qs + 1))
                .min(t1.saturating_sub(bl.ts + 1));
            let mut k = 0;
            while k < room && rv.base(read, k as usize) == gv.base(&self.index, k) {
                k += 1;
            }
            (q1, t1) = (q1 - k, t1 - k);
        }
        let right = self.end(read, reverse, q1, t1, false, minus, undecided, rescue, qlen, sc);
        let (mut q0, mut t0, _, _) = self.fixed(first, first);
        if rescue {
            let b0 = blocks[0];
            let (rv, gv) = (
                ReadView { q: q0, rev: false },
                View {
                    origin: t0,
                    rev: false,
                },
            );
            let room = BACK
                .min(b0.qe.saturating_sub(q0 + 1))
                .min(b0.te.saturating_sub(t0 + 1));
            let mut k = 0;
            while k < room && rv.base(read, k as usize) == gv.base(&self.index, k) {
                k += 1;
            }
            (q0, t0) = (q0 + k, t0 + k);
        }
        let left = self.end(read, reverse, q0, t0, true, minus, undecided, rescue, qlen, sc);
        edits += left.cost + right.cost;
        // left end: hops come in travel order (away from the anchor); in genome order the
        // furthest exon is first
        let b0 = blocks[0];
        let mut pre: Vec<Block> = Vec::new();
        let (mut eq, mut et) = left.end;
        for h in left.hops.iter().rev() {
            // travelling backwards, `enter` is the end of the exon further left
            pre.push(Block {
                qs: eq,
                qe: h.q,
                ts: et,
                te: h.enter,
            });
            eq = h.q;
            et = h.leave;
        }
        blocks[0] = Block {
            qs: eq,
            qe: b0.qe,
            ts: et,
            te: b0.te,
        };
        origins.splice(0..0, std::iter::repeat_n(b'e', pre.len()));
        origins.extend(std::iter::repeat_n(b'e', right.hops.len()));
        blocks.splice(0..0, pre);
        let nb = blocks.len() - 1;
        let bl = blocks[nb];
        let (mut sq, mut st) = (bl.qs, bl.ts);
        blocks.pop();
        for h in &right.hops {
            blocks.push(Block {
                qs: sq,
                qe: h.q,
                ts: st,
                te: h.leave,
            });
            sq = h.q;
            st = h.enter;
        }
        blocks.push(Block {
            qs: sq,
            qe: right.end.0,
            ts: st,
            te: right.end.1,
        });
        let fixed = tidy(&mut blocks, qlen);
        if fixed || origins.len() + 1 != blocks.len() {
            origins = vec![b'?'; blocks.len().saturating_sub(1)];
        }
        sc.stats.tidied += fixed as u64;
        // A hit lies on one sequence, that of its longest block: blocks that a chain or a read
        // end reached on a neighbouring one (the sequences lie end to end) are dropped.
        if !blocks.is_empty() {
            let k = (0..blocks.len())
                .max_by_key(|&x| blocks[x].te - blocks[x].ts)
                .unwrap_or(0);
            let o = self.index.offsets();
            let rid = self.index.locate(blocks[k].ts).0;
            let (lo, hi) = (o[rid], o[rid + 1]);
            if blocks[0].ts < lo || blocks[blocks.len() - 1].te > hi {
                let a = (0..=k)
                    .rev()
                    .take_while(|&x| blocks[x].ts >= lo)
                    .last()
                    .unwrap_or(k);
                let z = (k..blocks.len())
                    .take_while(|&x| blocks[x].te <= hi)
                    .last()
                    .unwrap_or(k);
                blocks.drain(z + 1..);
                blocks.drain(..a);
                origins = vec![b'?'; blocks.len() - 1];
            }
            // A few read bases left over at an end are aligned too when they fit the genome
            // next to the hit: where a read ends then shows in where its hit ends.
            let nb = blocks.len() - 1;
            if let Some(d) = self.close(read, blocks[0].qs, blocks[0].ts - lo, true, blocks[0].ts, sc) {
                blocks[0].qs = 0;
                blocks[0].ts -= d;
            }
            if let Some(d) = self.close(
                read,
                qlen - blocks[nb].qe,
                hi - blocks[nb].te,
                false,
                blocks[nb].te,
                sc,
            ) {
                blocks[nb].qe = qlen;
                blocks[nb].te += d;
            }
        }
        let tx_strand = if tx_strand == b'.' && blocks.len() > 1 && !self.opt.forward_only {
            if left.minus.or(right.minus).unwrap_or(minus) {
                b'-'
            } else {
                b'+'
            }
        } else {
            tx_strand
        };
        SplicedHit {
            reverse,
            blocks,
            tx_strand,
            score: score as i32 - (edits / 3) as i32,
            anchors: chain.len() as u32,
            edits,
            primary: false,
            origins,
            seeds: if self.opt.keep_seeds {
                chain.iter().map(|a| (a.q, a.t, a.ql.min(a.tl))).collect()
            } else {
                Vec::new()
            },
        }
    }

    /// The `c` read bases left over at a read end (`left`: before the hit), against the genome
    /// next to the hit at `t` with `room` bases to the end of the sequence: the number of genome
    /// bases they align to with the fewest edits (on a tie, the number closest to `c`), if those
    /// edits are few enough to be sequencing errors. A low-complexity end (a tail) stays clipped.
    fn close(&self, read: &[u8], c: u32, room: u32, left: bool, t: u32, sc: &mut Scratch) -> Option<u32> {
        if c == 0 || c > END_CLOSE || room == 0 {
            return None;
        }
        let (c, n) = (c as usize, (c + 3).min(room) as usize);
        let rb = |i: usize| if left { read[c - 1 - i] } else { read[read.len() - c + i] } & 3;
        if c >= 6 {
            let mut count = [0usize; 4];
            (0..c).for_each(|i| count[rb(i) as usize] += 1);
            if 4 * count.iter().max().unwrap() >= 3 * c {
                return None;
            }
        }
        if left {
            self.index.bases(t - n as u32, t, &mut sc.g1);
        } else {
            self.index.bases(t, t + n as u32, &mut sc.g1);
        }
        let g = &sc.g1;
        let gb = |j: usize| if left { g[n - 1 - j] } else { g[j] };
        let mut row = [0u8; END_CLOSE as usize + 4];
        (0..=n).for_each(|j| row[j] = j as u8);
        for i in 0..c {
            let mut diag = row[0];
            row[0] = i as u8 + 1;
            for j in 0..n {
                let v = (diag + (rb(i) != gb(j)) as u8)
                    .min(row[j] + 1)
                    .min(row[j + 1] + 1);
                diag = row[j + 1];
                row[j + 1] = v;
            }
        }
        let (e, _, j) = (0..=n).map(|j| (row[j], j.abs_diff(c), j)).min().unwrap();
        (3 * e as usize <= c + 2).then_some(j as u32)
    }

    /// Solve a read end: the bases before read position `q` (`left`) or from it on, starting
    /// at genome position `t`. Works in rounds until the read is used up or nothing more fits:
    /// run along the genome while it matches, then solve what follows (up to 400 bases) with
    /// the candidate exons found for it. Returns the junctions found, where the alignment
    /// stops, and the strand used when it had to be chosen here.
    #[allow(clippy::too_many_arguments)]
    fn end(
        &self,
        read: &[u8],
        reverse: bool,
        q: u32,
        t: u32,
        left: bool,
        minus: bool,
        undecided: bool,
        rescue: bool,
        qlen: u32,
        sc: &mut Scratch,
    ) -> EndSolved {
        let ix = &*self.index;
        let (mut q, mut t) = (q, t);
        let mut out = EndSolved {
            cost: 0,
            hops: Vec::new(),
            end: (q, t),
            minus: None,
        };
        let half = self.opt.min_intron / 2;
        let mut strands: Vec<bool> = if undecided { vec![false, true] } else { vec![minus] };
        if !rescue {
            // a thin chain: only run along the genome while it matches (stops early on junk)
            let (dq, dt, _) = self.extend(read, q, t, left, sc);
            out.end = if left { (q - dq, t - dt) } else { (q + dq, t + dt) };
            return out;
        }
        for _round in 0..8 {
            let mut m = if left { q } else { qlen - q };
            if m == 0 {
                break;
            }
            // A long stretch first runs along the genome until it stops matching; only what is
            // left (with a margin, since a run can overshoot an exon edge by chance) is solved.
            if m > 100 {
                let (dq, dt, _) = self.extend(read, q, t, left, sc);
                if m - dq < 10 {
                    out.end = if left { (q - dq, t - dt) } else { (q + dq, t + dt) };
                    break;
                }
                let back = dq.min(dt).min(25);
                let (aq, at) = (dq - back, dt - back);
                (q, t) = if left { (q - aq, t - at) } else { (q + aq, t + at) };
                m -= aq;
            }
            let whole = m <= 400;
            let m = m.min(400) as usize;
            let rv = ReadView { q, rev: left };
            let gv = View { origin: t, rev: left };
            let plain = gap::solve(ix, read, rv, m, gv, None, &[], strands[0], half, &mut sc.gap);
            let reached = plain.end.0.abs_diff(q);
            let plain_cost = plain.cost;
            let mut best = (plain, None);
            let known = self.index.junctions.is_some();
            // What is left of the read may be a poly(A) tail (poly(T) on the other strand) or
            // other low-complexity sequence: that is clipped, never searched for.
            let tail = {
                let n = (m - reached as usize).min(24);
                let mut c = [0usize; 4];
                for i in 0..n {
                    c[rv.base(read, reached as usize + i) as usize & 3] += 1;
                }
                n > 0 && 4 * c.iter().max().unwrap() >= 3 * n
            };
            // With known junctions, first the known introns that start where the alignment
            // stops: a read end that goes on in the exon after one of them is believed on a
            // few bases, and then nothing else is searched.
            let mut settled = false;
            if rescue && !tail && known && m as u32 - reached >= 4 && m as u32 - reached <= KNOWN_END_SPAN {
                let p_end = best.0.end.1.abs_diff(t);
                let mut kh = self.known_hops(gv, p_end, reached as usize, m - reached as usize);
                kh.sort_by_key(|w| w.0);
                if !kh.is_empty() {
                    let j = self.index.junctions.as_ref().unwrap();
                    for &mi in &strands {
                        sc.gap.known_only = true;
                        sc.gap.known_extra = KNOWN_END;
                        let s = gap::solve(ix, read, rv, m, gv, None, &kh, mi, half, &mut sc.gap);
                        sc.gap.known_only = false;
                        sc.gap.known_extra = 0;
                        let all_known = !s.hops.is_empty()
                            && s.hops
                                .iter()
                                .all(|h| j.intron(h.leave.min(h.enter), h.leave.max(h.enter)));
                        if all_known && s.cost < best.0.cost && s.cost <= plain_cost {
                            settled = m as u32 - s.end.0.abs_diff(q) < 4;
                            best = (s, Some(mi));
                        }
                    }
                }
            }
            if rescue && !tail && !settled && m as u32 - reached >= 10 {
                sc.stats.end_tried += 1;
                // level-2 cores of the unexplained bases of this round (the bases the genome
                // next to the hit explains have their hits there), within reach of the longest
                // intron
                let slack = reached.saturating_sub(12);
                let here = Anchor {
                    t,
                    q: if left { q - slack } else { q + slack },
                    ql: 0,
                    tl: 0,
                    q_in: 0,
                    q_out: 0,
                    w: 0.0,
                };
                let hits = if left {
                    self.rescue(reverse, None, Some(here), q - m as u32, qlen, sc)
                } else {
                    self.rescue(reverse, Some(here), None, 0, q + m as u32, sc)
                };
                // windows start at the unexplained bases, not at the start of this round
                let windows = |hits: &[Anchor]| {
                    if left {
                        windows_of(
                            hits,
                            q.saturating_sub(m as u32),
                            q - slack,
                            t.saturating_sub(self.opt.max_intron),
                            t,
                            4,
                            END_HIT_GAP,
                        )
                    } else {
                        let hi = t.saturating_add(self.opt.max_intron).min(ix.total_len());
                        windows_of(hits, q + slack, q + m as u32, t, hi, 4, END_HIT_GAP)
                    }
                };
                // Most introns are short, and within the longest intron's reach a noisy read
                // has chance groups of hits with more hits than its next exon has: the best
                // candidates nearby are kept as well as the best overall.
                let near: Vec<Anchor> = hits
                    .iter()
                    .copied()
                    .filter(|h| h.t.abs_diff(t) <= END_NEAR)
                    .collect();
                // A read end that is no part of the transcript (an adapter, a primer) has
                // chance hits too, one at a time. A single hit is a candidate only nearby, and
                // only when it lies in the first unexplained bases: a read that went on in an
                // exon further off, or for longer, would have more.
                let mut wins = windows(&hits);
                wins.retain(|w| w.hits >= 2);
                if near.len() < hits.len() {
                    for w in windows(&near) {
                        if wins.iter().all(|o| w.e <= o.s || o.e <= w.s) {
                            wins.push(w);
                        }
                    }
                    wins.sort_by_key(|w| w.s);
                }
                wins.retain(|w| {
                    w.hits >= 2
                        || (if left {
                            t - w.e <= END_NEAR
                        } else {
                            w.s - t <= END_NEAR
                        } && if left {
                            (q - reached).saturating_sub(w.l.q + w.l.ql as u32) <= END_HIT_GAP
                        } else {
                            w.f.q.saturating_sub(q + reached) <= END_HIT_GAP
                        })
                });
                let vw: Vec<SolverWindow> = if left {
                    wins.iter().rev().map(|w| w.backward(q, t)).collect()
                } else {
                    wins.iter().map(|w| w.forward(q, t)).collect()
                };
                let p_end = best.0.end.1.abs_diff(t);
                let unexplained = m as u32 - reached;
                // (the approximate search needs a long stretch to tell an exon from chance)
                let scanned = if unexplained >= SCAN_END_MIN {
                    self.scan_windows_end(read, rv, m, reached as usize, gv, p_end)
                } else {
                    Vec::new()
                };
                // the candidates from seed hits, plus known exons, exons that start at a
                // splice site and match the unexplained bases, and those of the approximate
                // search
                let mut cands: Vec<Vec<SolverWindow>> = Vec::new();
                for &mi in &strands {
                    let mut all = vw.clone();
                    all.extend(self.known_hops(gv, p_end, reached as usize, m - reached as usize));
                    // (an exon found by its splice site and a few matching bases alone is
                    // believed for a short read end only: those have no seed to be found by)
                    if unexplained < SITE_END_MAX {
                        all.extend(self.site_windows_end(
                            read,
                            rv,
                            m,
                            reached as usize,
                            gv,
                            p_end,
                            mi,
                            true,
                            sc,
                        ));
                    }
                    all.extend(scanned.iter().copied());
                    all.sort_by_key(|w| w.0);
                    cands.push(all);
                }
                // With the strand undecided, first one solve with every site priced as on its
                // cheaper strand and both strands' candidates: no strand can do better, so
                // when that does not beat what is in hand, neither strand is solved.
                let mut hopeless = false;
                if strands.len() == 2 && cands.iter().any(|c| !c.is_empty()) {
                    let mut both = cands[0].clone();
                    both.extend(cands[1].iter().copied());
                    both.sort_by_key(|w| w.0);
                    both.dedup();
                    sc.gap.set_both_strands(true);
                    let s = gap::solve(ix, read, rv, m, gv, None, &both, strands[0], half, &mut sc.gap);
                    sc.gap.set_both_strands(false);
                    hopeless = s.cost >= best.0.cost;
                }
                for (&mi, all) in strands.iter().zip(&cands) {
                    if hopeless {
                        break;
                    }
                    let mut s = Solved {
                        cost: u32::MAX,
                        ..Default::default()
                    };
                    if !all.is_empty() {
                        s = gap::solve(ix, read, rv, m, gv, None, all, mi, half, &mut sc.gap);
                    }
                    if s.cost == u32::MAX {
                        continue;
                    }
                    if std::env::var_os("GMAP_DEBUG").is_some() {
                        eprintln!(
                            "end left={left} q={q} t={t} m={m} minus={mi} reached={reached} hits={} windows={:?} plain={} cost={} hops={:?} end={:?}",
                            hits.len(),
                            all,
                            best.0.cost,
                            s.cost,
                            s.hops,
                            s.end
                        );
                    }
                    // (an exon found by the approximate search alone has no seed behind it)
                    let by_scan = s.hops.iter().any(|h| {
                        let p = h.enter.abs_diff(t);
                        scanned.iter().any(|w| w.0 <= p && p <= w.0 + w.1)
                            && !vw.iter().any(|w| w.0 <= p && p <= w.0 + w.1)
                    });
                    let margin = self
                        .hop_margin(&s, mi)
                        .max(if by_scan { SCAN_END_MARGIN } else { 0 });
                    if !s.hops.is_empty() && s.cost < best.0.cost && s.cost + margin <= plain_cost {
                        best = (s, Some(mi));
                    }
                }
                sc.stats.end_kept += best.1.is_some() as u64;
            }
            let (s, chosen) = best;
            let got = s.end.0.abs_diff(q);
            out.cost += s.cost.saturating_sub(CLIP * (m as u32 - got));
            out.end = s.end;
            let hopped = !s.hops.is_empty();
            out.hops.extend(s.hops);
            if let Some(mi) = chosen {
                out.minus = out.minus.or(Some(mi));
                strands = vec![mi];
            }
            // Go on only if this round used all it was given and more of the read remains.
            // (after a junction, what is still left of the read may lie in a further exon)
            let more = hopped && got > 0 && m as u32 - got >= if known { 4 } else { 10 };
            if !more && (whole || (m as u32 - got >= 10 && !hopped) || got < 30) {
                break;
            }
            // The next round starts a few bases back, on the exact matches this one ended with:
            // an exon's last bases can match the intron beside it by chance, and a junction
            // there is placed by the round that crosses it.
            let seg = out.hops.last().map_or(got, |h| s.end.0.abs_diff(h.q).min(got));
            let back = self
                .exact_run(read, s.end.0, s.end.1, !left)
                .min(END_BACK)
                .min(seg);
            (q, t) = if left {
                (s.end.0 + back, s.end.1 + back)
            } else {
                (s.end.0 - back, s.end.1 - back)
            };
        }
        out
    }

    /// How many read bases next to read position `q` (before it when `left`) match the genome
    /// next to `t` exactly, one after another on the same diagonal; at most 24.
    fn exact_run(&self, read: &[u8], q: u32, t: u32, left: bool) -> u32 {
        const LOW: u64 = 0x5555_5555_5555_5555;
        let ix = &*self.index;
        let avail = if left { q as usize } else { read.len() - q as usize };
        let gv = View { origin: t, rev: left };
        if avail < 24 || gv.room(ix) <= 32 {
            return 0;
        }
        let Some(g) = gv.word(ix, 0) else { return 0 };
        let rw = (0..24).fold(0u64, |w, i| {
            let b = if left {
                read[q as usize - 1 - i]
            } else {
                read[q as usize + i]
            };
            w | (b as u64 & 3) << (2 * i)
        });
        let x = (g ^ rw) & ((1u64 << 48) - 1);
        (((x | (x >> 1)) & LOW).trailing_zeros() / 2).min(24)
    }

    /// Alignment score of a hit: per block, two per read base minus the banded affine cost of
    /// the block against its genome interval; splice sites off the canonical pair cost as in
    /// junction placement.
    fn align_score(&self, read: &[u8], h: &SplicedHit, sc: &mut Scratch) -> i32 {
        let ix = &*self.index;
        let minus = h.tx_strand == b'-';
        let mut score = 0i32;
        for (k, b) in h.blocks.iter().enumerate() {
            let (m, n) = ((b.qe - b.qs) as usize, (b.te - b.ts) as usize);
            if m == 0 {
                continue;
            }
            if m == n {
                // equal lengths: when the block matches with a few mismatches and no shift, that
                // is its alignment, and no table is needed
                let (q, t) = (b.qs as usize, b.ts);
                let allowed = 2 + m / 64;
                let mut mm = 0;
                for i in 0..m {
                    mm += (read[q + i] != ix.base(t + i as u32)) as usize;
                    if mm > allowed {
                        break;
                    }
                }
                if mm <= allowed {
                    score += 2 * m as i32 - (mm as u16 * MISMATCH) as i32;
                    if k > 0 && h.tx_strand != b'.' {
                        score -= donor_pen(ix, h.blocks[k - 1].te, minus) as i32
                            + acceptor_pen(ix, b.ts, minus) as i32;
                    }
                    continue;
                }
            }
            if m.abs_diff(n) > 64.max(m / 4) {
                // lengths too far apart to align in a band: the block scores as unaligned
                score -= 20;
                continue;
            }
            // (a long block of a noisy read drifts: its two ends can agree in length while the
            // middle is tens of bases off the diagonal)
            let band = (8 + m / 16).min(48) + m.abs_diff(n);
            let w = 2 * band + 1;
            let (q, t) = (b.qs as usize, b.ts);
            band_prefix(
                |i| read[q + i],
                m,
                |j| ix.base(t + j as u32),
                n,
                band,
                &mut sc.e1,
                &mut sc.etmp,
                0,
            );
            let cost = sc.e1[m * w + (n + band - m)];
            score += 2 * m as i32 - cost.min(2 * m as u16 + 20) as i32;
            if k > 0 && h.tx_strand != b'.' {
                score -=
                    donor_pen(ix, h.blocks[k - 1].te, minus) as i32 + acceptor_pen(ix, b.ts, minus) as i32;
            }
        }
        score
    }

    /// Candidate exons for the read bases `[q0, q1)` inside the intron `(t0, t1)`, found by
    /// short words: every 8-base word of those read bases is looked for along the intron, and
    /// matches that agree on a diagonal and together cover `WORD_COVER` read bases or more make
    /// a candidate. Read bases with many errors keep too few level-2 cores for the usual search,
    /// but still share short words with their exon; by chance, a diagonal collects one match
    /// at most. For junctions that nothing else explained, and introns up to `WORD_SPAN`.
    fn word_windows(
        &self,
        read: &[u8],
        q0: u32,
        t0: u32,
        q1: u32,
        t1: u32,
        sc: &mut Scratch,
    ) -> Vec<SolverWindow> {
        const K: usize = 8;
        const TOP: u32 = 2 * (K as u32 - 1);
        let ix = &*self.index;
        let m = (q1 - q0) as usize;
        let (lo, hi) = (t0 + 20, t1.saturating_sub(20));
        if !(24..=600).contains(&m) || hi <= lo + 2 * K as u32 || hi - lo > WORD_SPAN {
            return Vec::new();
        }
        if sc.kfirst.is_empty() {
            sc.kfirst.resize(1 << (2 * K), 0);
        }
        sc.knext.clear();
        sc.knext.resize(m, 0);
        // (words of one base or of two alternating are left out)
        let dull = |c: u32| (c ^ (c >> 2)) & 0x3fff == 0 || (c ^ (c >> 4)) & 0x0fff == 0;
        let (mut code, mut run) = (0u32, 0usize);
        for i in 0..m {
            let b = read[q0 as usize + i] as u32;
            run = if b > 3 { 0 } else { run + 1 };
            code = (code >> 2) | ((b & 3) << TOP);
            if run >= K && !dull(code) {
                let at = i + 1 - K;
                sc.knext[at] = sc.kfirst[code as usize];
                sc.kfirst[code as usize] = at as u16 + 1;
            }
        }
        // matches as (genome position of read offset 0, plus m to stay positive; read offset)
        sc.khits.clear();
        let (mut g, mut n, mut p) = (0u32, 0usize, lo);
        'scan: while p < hi {
            let w = ix.word(p);
            let cnt = (hi - p).min(32);
            for j in 0..cnt {
                g = (g >> 2) | ((((w >> (2 * j)) & 3) as u32) << TOP);
                n += 1;
                if n < K {
                    continue;
                }
                let mut x = sc.kfirst[g as usize];
                while x != 0 {
                    let at = x as u32 - 1;
                    sc.khits.push((p + j + 1 - K as u32 + m as u32 - at, at as u16));
                    x = sc.knext[at as usize];
                }
                if sc.khits.len() > 8192 {
                    break 'scan;
                }
            }
            p += cnt;
        }
        // (the table is cleared by walking the read's words again)
        let (mut code, mut run) = (0u32, 0usize);
        for i in 0..m {
            let b = read[q0 as usize + i] as u32;
            run = if b > 3 { 0 } else { run + 1 };
            code = (code >> 2) | ((b & 3) << TOP);
            if run >= K {
                sc.kfirst[code as usize] = 0;
            }
        }
        if sc.khits.len() > 8192 || sc.khits.len() < 2 {
            return Vec::new();
        }
        sc.khits.sort_unstable();
        // groups of matches on diagonals no more than 3 apart from one to the next
        let mut found: Vec<(u32, SolverWindow)> = Vec::new();
        let mut offs: Vec<u16> = Vec::new();
        let hits = &sc.khits;
        let mut a = 0;
        for b in 1..=hits.len() {
            if b < hits.len() && hits[b].0 - hits[b - 1].0 <= 3 {
                continue;
            }
            if b - a >= 2 {
                offs.clear();
                offs.extend(hits[a..b].iter().map(|h| h.1));
                offs.sort_unstable();
                let (mut cover, mut end) = (0u32, 0u32);
                for &o in &offs {
                    let (s, e) = (o as u32, o as u32 + K as u32);
                    cover += e - s.max(end).min(e);
                    end = end.max(e);
                }
                if cover >= WORD_COVER {
                    // where read offset 0 falls on the group's middle diagonal, and the window:
                    // the matched bases, and on either side the read bases still unexplained
                    // (150 at most), as for windows from level-2 hits
                    let zero = hits[(a + b) / 2].0 as i64 - m as i64;
                    let (first, last) = (offs[0] as i64, offs[offs.len() - 1] as i64 + K as i64);
                    let gs = (zero + first - first.min(150) - 8).max(lo as i64);
                    let ge = (zero + last + (m as i64 - last).min(150) + 8).min(hi as i64);
                    if ge > gs + K as i64 {
                        found.push((cover, ((gs - t0 as i64) as u32, (ge - gs) as u32, gs - zero, 0)));
                    }
                }
            }
            a = b;
        }
        // the best covered first; a few, apart from one another
        found.sort_by_key(|f| std::cmp::Reverse(f.0));
        let mut out: Vec<SolverWindow> = Vec::new();
        for (_, w) in found {
            if out.len() < 4 && out.iter().all(|o| w.0 + w.1 <= o.0 || o.0 + o.1 <= w.0) {
                out.push(w);
            }
        }
        out
    }

    /// Candidate exons for the read bases `[q0, q1)` inside the intron `(t0, t1)`, found by
    /// approximate matching: pieces of those read bases (32 and 48 long, overlapping) are
    /// looked for along the intron with Myers' bit-vector method, which gives, at every genome
    /// position, the fewest edits between a piece and any stretch of the genome ending there.
    /// A piece that matches within about a third of its length makes a candidate. This finds
    /// exons whose read bases are too noisy to share even a short word with the genome. For
    /// junctions that nothing else explained, and introns up to `WORD_SPAN`.
    #[allow(clippy::too_many_arguments)]
    fn scan_windows(
        &self,
        read: &[u8],
        q0: u32,
        t0: u32,
        q1: u32,
        t1: u32,
        mid: (usize, usize),
    ) -> Vec<SolverWindow> {
        let ix = &*self.index;
        let m = (q1 - q0) as usize;
        let (lo, hi) = (t0 + 20, t1.saturating_sub(20));
        // the stretch neither side explains, with a few bases more on either side
        let (a, b) = (mid.0.saturating_sub(8), (mid.1 + 8).min(m));
        if !(32..=600).contains(&m) || b < a + 32 || hi <= lo + 32 || hi - lo > WORD_SPAN {
            return Vec::new();
        }
        // pieces: (offset in the gap, length, edits allowed)
        let mut pieces: Vec<(usize, usize, u32)> = Vec::new();
        for (len, allowed) in [(32usize, SCAN_EDITS[0]), (48, SCAN_EDITS[1])] {
            if b - a < len {
                continue;
            }
            let mut o = a;
            while o + len < b {
                pieces.push((o, len, allowed));
                o += 16;
            }
            pieces.push((b - len, len, allowed));
        }
        if pieces.len() > SCAN_PIECES {
            // too many for one pass: every other one
            let mut k = 0;
            pieces.retain(|_| {
                k += 1;
                k % 2 == 1
            });
            pieces.truncate(SCAN_PIECES);
        }
        let mut peq = [[0u64; 4]; SCAN_PIECES];
        for (k, &(o, len, _)) in pieces.iter().enumerate() {
            for i in 0..len {
                let b = read[q0 as usize + o + i];
                if b < 4 {
                    peq[k][b as usize] |= 1 << i;
                }
            }
        }
        // (edits per 64 bases, piece, genome position after the match)
        let mut found: Vec<(u32, usize, u32)> = Vec::new();
        scan(&|p| ix.word(p), lo, hi, &pieces, &peq, &mut found);
        // the best matches first; a few, apart from one another
        found.sort_unstable();
        let mut out: Vec<SolverWindow> = Vec::new();
        for (_, k, end) in found {
            let (o, len, _) = pieces[k];
            // where read offset 0 of the gap falls on the match's diagonal, and the window: the
            // stretch neither side explains, on that diagonal, and a few bases more
            let zero = end as i64 - (o + len) as i64;
            let gs = (zero + a as i64 - 16).max(lo as i64);
            let ge = (zero + b as i64 + 16).min(hi as i64);
            let w: SolverWindow = ((gs - t0 as i64) as u32, (ge - gs).max(0) as u32, gs - zero, 0);
            if ge > gs + 16 && out.len() < 4 && out.iter().all(|o| w.0 + w.1 <= o.0 || o.0 + o.1 <= w.0) {
                out.push(w);
            }
        }
        out
    }

    /// Candidate exons for read bases that neither side of a junction explains, found without
    /// any seed: genome intervals inside the intron that start after an acceptor site, end
    /// before a donor site, are about as long as the unexplained bases, and match them.
    #[allow(clippy::too_many_arguments)]
    fn site_windows(
        &self,
        read: &[u8],
        q0: u32,
        t0: u32,
        q1: u32,
        t1: u32,
        minus: bool,
        sc: &mut Scratch,
    ) -> (Vec<Window>, Option<(usize, usize)>) {
        let ix = &*self.index;
        let m = (q1 - q0) as usize;
        // How far the read follows the genome from either side.
        let band = 4 + m / 10;
        let n = (m + band).min((t1 - t0) as usize);
        let reach = |e: &Vec<u16>, rows: usize| {
            let w = 2 * band + 1;
            let mut best = (0i32, 0usize);
            for i in 1..=rows {
                let row = e[i * w..(i + 1) * w].iter().copied().min().unwrap_or(INF);
                let s = 2 * i as i32 - row as i32;
                if s > best.0 {
                    best = (s, i);
                }
            }
            best.1
        };
        let rows = band_prefix(
            |i| read[q0 as usize + i],
            m,
            |j| ix.base(t0 + j as u32),
            n,
            band,
            &mut sc.e1,
            &mut sc.etmp,
            30,
        );
        let up = reach(&sc.e1, rows);
        let rows = band_prefix(
            |i| read[q1 as usize - 1 - i],
            m,
            |j| ix.base(t1 - 1 - j as u32),
            n,
            band,
            &mut sc.e2,
            &mut sc.etmp,
            30,
        );
        let down = reach(&sc.e2, rows);
        if up + down + 6 > m {
            return (Vec::new(), None);
        }
        // The unexplained bases, with a margin on both sides (the read can follow an intron a
        // few bases by chance), and the exon lengths that fit them.
        const MARGIN: usize = 4;
        let (up, down) = (up.saturating_sub(MARGIN), down.saturating_sub(MARGIN));
        let k = m - up - down;
        let mid = &read[q0 as usize + up..q1 as usize - down];
        let (lo, hi) = (t0 + up as u32 + 20, (t1 - down as u32).saturating_sub(20));
        let mid_at = Some((up, m - down));
        if hi <= lo || k > 200 + 2 * MARGIN {
            return (Vec::new(), mid_at);
        }
        let slack = 2 + k / 10;
        let (len_lo, len_hi) = (k.saturating_sub(2 * MARGIN + slack).max(4), k + slack);
        // an exon under ten bases is only looked for close by: further off, chance matches win
        if len_lo < 10 && hi - lo > 20_000 {
            return (Vec::new(), mid_at);
        }
        let mut found: Vec<(usize, Window)> = Vec::new();
        let mut scanned = 0u64;
        scanned += (hi - lo) as u64;
        // acceptors: AG before the exon on the forward strand, AC on the reverse
        let (a0, a1) = if minus { (0, 1) } else { (0, 2) };
        let mut starts: Vec<u32> = Vec::new();
        ix.pair_scan(lo.saturating_sub(2), hi.saturating_sub(2), a0, a1, |s| {
            starts.push(s + 2)
        });
        // where an exon can end: before a usual donor site
        // (GT on the forward strand, CT on the reverse)
        let mut ends = vec![0u64; ((hi - lo) as usize >> 6) + 1];
        ix.pair_scan(lo, hi + 1, if minus { 1 } else { 2 }, 3, |e| {
            ends[(e - lo) as usize >> 6] |= 1 << ((e - lo) & 63)
        });
        // the 32 read bases from each offset of the unexplained stretch, two bits each, to
        // compare 32 bases with the genome at once
        const LOW: u64 = 0x5555_5555_5555_5555;
        let mut packed = vec![0u64; k + 1];
        for o in (0..k).rev() {
            packed[o] = (packed[o + 1] << 2) | (mid[o] as u64 & 3);
        }
        for p in starts {
            if p < lo {
                continue;
            }
            let gw = ix.word(p);
            for len in len_lo..=len_hi {
                let e = p + len as u32;
                if e > hi || ends[(e - lo) as usize >> 6] >> ((e - lo) & 63) & 1 == 0 {
                    continue;
                }
                // the exon against the unexplained bases at every offset it can have in them
                let allowed = len / 8;
                let mut best = usize::MAX;
                for o in 0..=k.saturating_sub(len) {
                    let n = len.min(k - o);
                    // the first 32 bases at once: most offsets fail here
                    let x = gw ^ packed[o];
                    let x = (x | (x >> 1)) & LOW & if n >= 32 { !0 } else { (1u64 << (2 * n)) - 1 };
                    if x.count_ones() as usize > allowed {
                        continue;
                    }
                    let mut d = 0;
                    for i in 0..n {
                        d += (mid[o + i] != ix.base(p + i as u32)) as usize;
                        if d > allowed {
                            break;
                        }
                    }
                    best = best.min(d);
                }
                if best <= allowed {
                    let at = Anchor {
                        t: p,
                        q: q0 + up as u32,
                        ql: 0,
                        tl: 0,
                        q_in: 0,
                        q_out: 0,
                        w: 0.0,
                    };
                    let w = Window {
                        s: p.saturating_sub(3).max(lo),
                        e: (e + 3).min(hi),
                        f: at,
                        l: at,
                        hits: 0,
                    };
                    found.push((best * 64 / len.max(1), w));
                }
            }
        }
        sc.stats.site_bases += scanned;
        // Non-overlapping, best match first; a few at most.
        found.sort_by_key(|f| f.0);
        let mut out: Vec<Window> = Vec::new();
        for (_, w) in found {
            if out.len() < 4 && out.iter().all(|o| w.e + 20 <= o.s || o.e + 20 <= w.s) {
                out.push(w);
            }
        }
        out.sort_by_key(|w| w.s);
        (out, mid_at)
    }

    /// How much cheaper than clipping a read-end solution with junctions must be to be believed:
    /// little when every junction is a known one, otherwise enough that a short chance match
    /// (a level-2 core is 11-13 bases) cannot pay for an exon.
    fn hop_margin(&self, s: &Solved, minus: bool) -> u32 {
        let known = self.index.junctions.as_ref().is_some_and(|j| {
            s.hops.iter().all(|h| {
                let (l, r) = (h.leave.min(h.enter), h.leave.max(h.enter));
                j.left_site(l, minus) && j.right_site(r, minus)
            })
        });
        if known {
            return 4;
        }
        // An unknown junction at a read end is believed on little when both its sites are the
        // usual ones and the intron is short, and only on much more otherwise: among real
        // reads' end junctions, those are the ones that turn out to be annotated.
        let ix = &*self.index;
        s.hops
            .iter()
            .map(|h| {
                let (l, r) = (h.leave.min(h.enter), h.leave.max(h.enter));
                let (d, a) = (donor_motif(ix, l, minus), acceptor_motif(ix, r, minus));
                let plain = d + a == 0;
                // (a site of none of the known kinds: among simulated reads' unknown end
                // junctions, hardly one in ten of those is right)
                if d.max(a) >= NONCANON {
                    return END_REFUSE;
                }
                match (plain, r - l < END_NEAR) {
                    (true, true) => END_MARGIN[0],
                    (true, false) => END_MARGIN[1],
                    (false, true) => END_MARGIN[2],
                    (false, false) => END_MARGIN[3],
                }
            })
            .max()
            .unwrap_or(END_MARGIN[0])
    }

    /// Known exons inside the intron `(t0, t1)` that the read bases of a gap of `m` bases could
    /// hold, as windows for the gap solver (searched over the whole gap, since nothing says
    /// where in the read they fall).
    fn known_exons(&self, t0: u32, t1: u32, m: usize) -> Vec<SolverWindow> {
        let Some(j) = &self.index.junctions else {
            return Vec::new();
        };
        if t1 < t0 + 60 {
            return Vec::new();
        }
        let mut out: Vec<SolverWindow> = Vec::new();
        for &(s, e) in j.exons_in(t0 + 20, t1 - 20) {
            if e + 20 <= t1 && (e - s) as usize <= m + 20 && out.len() < 12 {
                // (the exon lies somewhere in the m read bases: its first base goes with a read
                // position from 0 to m less its length, give or take a few for indels)
                let (l, m) = ((e - s) as i64, m as i64);
                out.push((s - t0, e - s, (m - l) / 2, (m - l).abs() / 2 + 8 + l / 8));
            }
        }
        out
    }

    /// Windows past a read end that known junctions lead to: introns that start within a few
    /// bases of where the plain extension stopped (`p_end`, in the view's coordinates).
    fn known_hops(&self, gv: View, p_end: u32, reached: usize, k: usize) -> Vec<SolverWindow> {
        let Some(j) = &self.index.junctions else {
            return Vec::new();
        };
        let x = gv.real(p_end);
        let vl = (k as u32 + 12).min(300);
        let mut out: Vec<SolverWindow> = Vec::new();
        if gv.rev {
            // travelling backwards: introns whose right end is near x; the next exon ends at
            // their left end
            for &(right, left) in j.from_right(x.saturating_sub(12), x + 30) {
                let over = right as i64 - x as i64; // > 0: the extension ran past the site
                if left >= vl {
                    out.push((gv.origin - left, vl, reached as i64 - over, 24));
                }
            }
        } else {
            for &(left, right) in j.from_left(x.saturating_sub(30), x + 12) {
                let over = x as i64 - left as i64;
                if right >= gv.origin {
                    out.push((right - gv.origin, vl, reached as i64 - over, 24));
                }
            }
        }
        out.truncate(8);
        out
    }

    /// Candidate exons past a read end, found without any seed: places where an intron can end
    /// (in the direction of travel) and the genome then matches the read bases that the plain
    /// extension left unexplained. Returned as the gap solver's windows.
    #[allow(clippy::too_many_arguments)]
    fn site_windows_end(
        &self,
        read: &[u8],
        rv: ReadView,
        m: usize,
        reached: usize,
        gv: View,
        p_end: u32,
        minus: bool,
        near: bool,
        sc: &mut Scratch,
    ) -> Vec<SolverWindow> {
        const MARGIN: usize = 4;
        let ix = &*self.index;
        let start = reached.saturating_sub(MARGIN);
        let k = m - start;
        if k < 10 + MARGIN {
            return Vec::new();
        }
        // a short tail is only looked for close by: further off, chance matches win
        // Seedless candidates are looked for nearby only: a scan of the longest intron's reach
        // costs more than every other step and found an exon for under 2% of real read ends.
        let _ = near;
        let reach = 20_000u32.min(self.opt.max_intron);
        let hi = (p_end + reach).min(gv.room(ix).saturating_sub(k as u32 + 8));
        let probe: Vec<u8> = (0..k.min(28 + 2 * MARGIN))
            .map(|i| rv.base(read, start + i))
            .collect();
        let mut found: Vec<(usize, u32, SolverWindow)> = Vec::new();
        // Where an intron can end in the direction of travel: travelling forwards, after an
        // acceptor (AG, or AC for a reverse-strand transcript); backwards, before a donor (GT,
        // or CT).
        let from = p_end + 20;
        let mut cands: Vec<u32> = Vec::new();
        if from < hi {
            if gv.rev {
                let (b0, b1) = if minus { (1, 3) } else { (2, 3) };
                let (a, b) = (gv.real(hi), gv.real(from));
                ix.pair_scan(a + 1, b + 1, b0, b1, |s| cands.push(gv.origin - s));
            } else {
                let (b0, b1) = if minus { (0, 1) } else { (0, 2) };
                let (a, b) = (gv.real(from), gv.real(hi));
                ix.pair_scan(a.saturating_sub(2), b.saturating_sub(2), b0, b1, |s| {
                    cands.push(s + 2 - gv.origin)
                });
            }
        }
        // The probe at each offset the exon can start at in it, as packed words, so that a
        // candidate is compared by XOR and a bit count.
        const LOW: u64 = 0x5555_5555_5555_5555;
        let offsets = (2 * MARGIN).min(probe.len().saturating_sub(10));
        let packed: Vec<(u64, u64, u32)> = (0..=offsets)
            .map(|o| {
                let len = (probe.len() - o).min(20);
                let w = (0..len).fold(0u64, |w, i| w | (probe[o + i] as u64) << (2 * i));
                (w, (1u64 << (2 * len)) - 1, (len / 8) as u32)
            })
            .collect();
        for p in cands {
            if p < from || p >= hi || gv.enter_motif(ix, p, minus) > 5 {
                continue;
            }
            let Some(g) = gv.word(ix, p) else { continue };
            // the exon can start a few read bases into the probe
            for (o, &(w, mask, allowed)) in packed.iter().enumerate() {
                let x = (g ^ w) & mask;
                let d = ((x | (x >> 1)) & LOW).count_ones();
                if d <= allowed {
                    let vs = p.saturating_sub(3);
                    let vl = (k as u32 + 12).min(300);
                    found.push((d as usize, p, (vs, vl, (start + o) as i64 - (p - vs) as i64, 0)));
                    break;
                }
            }
        }
        sc.stats.site_bases += (hi.saturating_sub(p_end)) as u64;
        // best match first, then nearest; a few at most
        found.sort_by_key(|f| (f.0, f.1));
        found.truncate(3);
        found.into_iter().map(|f| f.2).collect()
    }

    /// Candidate exons for the bases of a read end that the genome next to it does not explain
    /// (from `reached` bases into the `m` of the view `rv`), found by approximate matching as in
    /// [`Self::scan_windows`]: the first of those bases, in pieces, are looked for within
    /// `END_NEAR` of the view's origin, past `p_end`.
    fn scan_windows_end(
        &self,
        read: &[u8],
        rv: ReadView,
        m: usize,
        reached: usize,
        gv: View,
        p_end: u32,
    ) -> Vec<SolverWindow> {
        let ix = &*self.index;
        let start = reached.saturating_sub(8);
        let from = p_end + 20;
        let hi = (p_end + SCAN_END_SPAN.min(self.opt.max_intron)).min(gv.room(ix).saturating_sub(64));
        if m < start + 32 || hi <= from + 64 {
            return Vec::new();
        }
        // pieces as (offset in the view, length, edits allowed)
        let mut pieces: Vec<(usize, usize, u32)> = Vec::new();
        for (len, allowed, n) in [(32usize, SCAN_EDITS_END[0], 3), (48, SCAN_EDITS_END[1], 2)] {
            for k in 0..n {
                if start + 16 * k + len <= m {
                    pieces.push((start + 16 * k, len, allowed));
                }
            }
        }
        // (a piece is matched in genome order: backwards in the view when it runs backwards)
        let mut peq = [[0u64; 4]; SCAN_PIECES];
        for (k, &(o, len, _)) in pieces.iter().enumerate() {
            for i in 0..len {
                let b = rv.base(read, if rv.rev { o + len - 1 - i } else { o + i });
                if b < 4 {
                    peq[k][b as usize] |= 1 << i;
                }
            }
        }
        let (lo_real, hi_real) = if gv.rev {
            (gv.origin - hi, gv.origin - from)
        } else {
            (gv.origin + from, gv.origin + hi)
        };
        let mut found: Vec<(u32, usize, u32)> = Vec::new();
        scan(&|p| ix.word(p), lo_real, hi_real, &pieces, &peq, &mut found);
        found.sort_unstable();
        let mut out: Vec<SolverWindow> = Vec::new();
        for (_, k, end) in found {
            let (o, len, _) = pieces[k];
            // where the piece starts in the view, and the window: from a few bases before
            // it, for the read bases that remain
            let ps = if gv.rev {
                gv.origin - end
            } else {
                (end - gv.origin).saturating_sub(len as u32)
            };
            let vs = ps.saturating_sub(16).max(from);
            let vl = ((m - o) as u32 + 28)
                .min(300)
                .min(gv.room(ix).saturating_sub(vs + 8));
            if ps >= vs && vl >= 16 && out.len() < 3 && out.iter().all(|w| vs + vl <= w.0 || w.0 + w.1 <= vs)
            {
                out.push((vs, vl, o as i64 - (ps - vs) as i64, 0));
            }
        }
        out
    }

    /// Score change of chaining anchor `a` after `p`: the read bases `a` adds, minus the price of
    /// the gap between them.
    #[inline]
    fn step(&self, p: Anchor, a: Anchor) -> f32 {
        let (dq, dt) = ((a.q - p.q) as i64, (a.t - p.t) as i64);
        let d = dt - dq;
        let cost = if d >= self.opt.min_intron as i64 {
            // an intron: a flat price, slightly more for long ones
            6.0 + 0.5 * (d as f32).log2()
        } else {
            0.5 * d.abs() as f32 + if dq > 3000 { 0.01 * dq as f32 } else { 0.0 }
        };
        (a.ql as i64).min(dq).min(dt.max(0)) as f32 - cost
    }

    /// Level-2 anchors for read bases the chain leaves unexplained, searched only where they can
    /// lie: between chain neighbours `a` and `b`, or, when one of them is missing, within the
    /// longest intron's reach past the chain's first or last anchor, for the read bases from
    /// `q_lo` on or up to `q_hi`. Returns every hit, in genome order.
    fn rescue(
        &self,
        reverse: bool,
        a: Option<Anchor>,
        b: Option<Anchor>,
        q_lo: u32,
        q_hi: u32,
        sc: &mut Scratch,
    ) -> Vec<Anchor> {
        let ix = &*self.index;
        let (qa, ta) = match a {
            Some(a) => (a.q + a.ql as u32, a.t + a.tl as u32),
            None => (q_lo, b.map_or(0, |b| b.t.saturating_sub(self.opt.max_intron))),
        };
        let (qb, tb) = match b {
            Some(b) => (b.q, b.t),
            None => (q_hi, ta.saturating_add(self.opt.max_intron).min(ix.total_len())),
        };
        if qb < qa + 8 || tb <= ta {
            return Vec::new();
        }
        let seeds = &sc.l2[reverse as usize];
        let lo = seeds.partition_point(|s| s.q < qa);
        let hi = lo + seeds[lo..].partition_point(|s| s.q + s.ql as u32 <= qb);
        if lo == hi {
            return Vec::new();
        }
        let (pos, ids) = ix.l2_window(ta, tb);
        sc.stats.window_bases += (tb - ta) as u64;
        let mut hits: Vec<Anchor> = Vec::new();
        let mut at = |k: usize, first: usize| {
            let (t, id) = (pos[k], ids[k]);
            for s in seeds[first..hi].iter().filter(|s| s.id == id) {
                // the genome length of the core: as many runs as in the read
                // (without compression the genome core is as long as the read's)
                let mut e = if s.runs == 0 { t + s.ql as u32 } else { t };
                for _ in 0..s.runs {
                    let c = ix.base(e);
                    e += 1;
                    while e < tb && ix.base(e) == c {
                        e += 1;
                    }
                }
                if e > tb || hits.len() >= 8192 {
                    continue;
                }
                hits.push(Anchor {
                    t,
                    q: s.q,
                    ql: s.ql,
                    tl: (e - t).min(u16::MAX as u32) as u16,
                    q_in: s.q_in,
                    q_out: s.q_out,
                    w: 0.0,
                });
            }
        };
        // A few labels are looked for by comparing many genome labels at a time; more, through
        // a table.
        let mut few = [0u16; 8];
        let mut n = 0;
        for s in &seeds[lo..hi] {
            if !few[..n].contains(&s.id) {
                if n == 8 {
                    n = 9;
                    break;
                }
                few[n] = s.id;
                n += 1;
            }
        }
        #[cfg(target_arch = "x86_64")]
        let wide = n <= 8 && std::is_x86_feature_detected!("avx512bw");
        #[cfg(not(target_arch = "x86_64"))]
        let wide = false;
        if wide {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: the feature was just detected.
            unsafe {
                ids_among_avx512(ids, &few[..n], |k| at(k, lo))
            };
            return hits;
        }
        if sc.want.is_empty() {
            sc.want.resize(1 << 16, u32::MAX);
        }
        // want[id] = index of the first read seed with that id (others are found by scanning on)
        for (k, s) in seeds[lo..hi].iter().enumerate().rev() {
            sc.want[s.id as usize] = (lo + k) as u32;
        }
        for (k, &id) in ids.iter().enumerate() {
            let first = sc.want[id as usize];
            if first != u32::MAX {
                at(k, first as usize);
            }
        }
        for s in &seeds[lo..hi] {
            sc.want[s.id as usize] = u32::MAX;
        }
        hits
    }

    /// The corresponding read/genome points a junction after anchor `a` and before anchor `b`
    /// is placed between: the end of `a`'s first homopolymer run and the start of `b`'s last
    /// one. The anchors themselves lie between the two points on purpose: after homopolymer
    /// compression an anchor can run a few bases past the true exon edge (the intron's first
    /// bases often look like the next exon's), so the edge is searched for inside them too.
    fn fixed(&self, a: Anchor, b: Anchor) -> (u32, u32, u32, u32) {
        let ix = &*self.index;
        let mut t0 = a.t + 1;
        let c = ix.base(a.t);
        while t0 < a.t + a.tl as u32 && ix.base(t0) == c {
            t0 += 1;
        }
        let mut t1 = b.t + b.tl as u32 - 1;
        let c = ix.base(t1);
        while t1 > b.t && ix.base(t1 - 1) == c {
            t1 -= 1;
        }
        (a.q + a.q_in as u32, t0, b.q + b.q_out as u32, t1)
    }

    /// The best split of read `[q0, q1)` into a part that continues the genome from `t0` and a
    /// part that ends at `t1`, for a transcript on the forward strand (`[0]`) and on the reverse
    /// strand (`[1]`): the gap solver without candidate exons. The prefix tables are shared by
    /// the two strands and stay in `sc.gap` for whoever solves this gap next.
    fn place(&self, read: &[u8], q0: u32, t0: u32, q1: u32, t1: u32, sc: &mut Scratch) -> [Placed; 2] {
        let none = Placed {
            q: q0,
            donor: t0,
            acceptor: t1,
            cost: INF,
            origin: b'd',
        };
        if q1 <= q0 || t1 <= t0 || (q1 - q0) > 600 {
            // anchors that overlap on the read, or too many read bases to explain from two
            // sides (an unanchored exon; rescue comes later)
            return [none; 2];
        }
        let m = (q1 - q0) as usize;
        let half = self.opt.min_intron / 2;
        // Most junctions of a clean read need no table: the read bases match the genome exactly
        // from the left up to some point and from the right back to some point, and where the
        // two stretches overlap there is a split at canonical sites. That split costs nothing,
        // so no alignment can beat it.
        {
            let ix = &*self.index;
            let (r, span) = (&read[q0 as usize..q1 as usize], (t1 - t0) as usize);
            let mut l = 0;
            while l < m && l < span && r[l] == ix.base(t0 + l as u32) {
                l += 1;
            }
            let mut rr = 0;
            while rr < m && rr < span && r[m - 1 - rr] == ix.base(t1 - 1 - rr as u32) {
                rr += 1;
            }
            if l + rr >= m && m + (half as usize) <= span {
                let mut fast = [none; 2];
                let mut free = false;
                for (s, o) in fast.iter_mut().enumerate() {
                    let mut best = (u16::MAX, 0usize);
                    for i in (m - rr)..=l {
                        let c = donor_pen(ix, t0 + i as u32, s == 1)
                            + acceptor_pen(ix, t1 - (m - i) as u32, s == 1);
                        if c < best.0 {
                            best = (c, i);
                        }
                    }
                    free |= best.0 == 0;
                    *o = Placed {
                        q: q0 + best.1 as u32,
                        donor: t0 + best.1 as u32,
                        acceptor: t1 - (m - best.1) as u32,
                        cost: best.0,
                        origin: b'd',
                    };
                }
                if free {
                    sc.gap.reset(); // the tables in there are not this junction's
                    sc.stats.fast_placed += 1;
                    return fast;
                }
            }
        }
        let rv = ReadView { q: q0, rev: false };
        let gv = View {
            origin: t0,
            rev: false,
        };
        let mut out = [none; 2];
        for (s, o) in out.iter_mut().enumerate() {
            let solved = gap::solve(
                &self.index,
                read,
                rv,
                m,
                gv,
                Some(t1),
                &[],
                s == 1,
                half,
                &mut sc.gap,
            );
            if let (true, Some(h)) = (solved.cost < INF as u32, solved.hops.first()) {
                *o = Placed {
                    q: h.q,
                    donor: h.leave,
                    acceptor: h.enter,
                    cost: solved.cost as u16,
                    origin: b'd',
                };
            }
        }
        out
    }

    /// Of the 24 read bases next to read position `q` (before it when `left`), how many match
    /// the genome next to `t` on the best of the five nearest diagonals. Chance gives about 6-10.
    fn diagonal_matches(&self, read: &[u8], q: u32, t: u32, left: bool) -> u32 {
        const LOW: u64 = 0x5555_5555_5555_5555;
        let ix = &*self.index;
        let avail = if left { q as usize } else { read.len() - q as usize };
        if avail < 24 {
            return 0;
        }
        let rw = (0..24).fold(0u64, |w, i| {
            let b = if left {
                read[q as usize - 1 - i]
            } else {
                read[q as usize + i]
            };
            w | (b as u64 & 3) << (2 * i)
        });
        let gv = View { origin: t, rev: left };
        let room = gv.room(ix);
        let mut best = 0u32;
        for shift in 0..5u32 {
            // read offset 0 against genome offsets 0..=2, and read offsets 1..=2 against 0
            let (gs, rs) = if shift < 3 { (shift, 0) } else { (0, shift - 2) };
            if gs + 32 >= room {
                continue;
            }
            if let Some(g) = gv.word(ix, gs) {
                let n = 24 - rs;
                let x = (g ^ (rw >> (2 * rs))) & ((1u64 << (2 * n)) - 1);
                best = best.max(n - ((x | (x >> 1)) & LOW).count_ones());
            }
        }
        best
    }

    /// Extend from read position `q` / genome position `t` towards the read's start (`left`) or
    /// end. Returns the read and genome bases added and the edits among them.
    fn extend(&self, read: &[u8], q: u32, t: u32, left: bool, sc: &mut Scratch) -> (u32, u32, u32) {
        let ix = &*self.index;
        let m = if left { q as usize } else { read.len() - q as usize };
        if m == 0 {
            return (0, 0, 0);
        }
        // A quick look before any table: if the next 24 read bases match the genome on none of
        // the five nearest diagonals better than chance allows, there is nothing to extend
        // into (junk past a stray anchor, or the bases after an exon's edge).
        if m >= 24 && self.diagonal_matches(read, q, t, left) < 12 {
            return (0, 0, 0);
        }
        let band = (8 + m / 20).min(32);
        let w = 2 * band + 1;
        let n = if left {
            (m + band).min(t as usize)
        } else {
            (m + band).min((ix.total_len() - t) as usize)
        };
        if left {
            ix.bases(t - n as u32, t, &mut sc.g1);
        } else {
            ix.bases(t, t + n as u32, &mut sc.g1);
        }
        let g = &sc.g1;
        if left {
            band_prefix(
                |i| read[q as usize - 1 - i],
                m,
                |j| g[n - 1 - j],
                n,
                band,
                &mut sc.e1,
                &mut sc.etmp,
                60,
            );
        } else {
            band_prefix(
                |i| read[q as usize + i],
                m,
                |j| g[j],
                n,
                band,
                &mut sc.e1,
                &mut sc.etmp,
                60,
            );
        }
        // The furthest point whose score (one per read base, minus three per edit) is best.
        let mut best = (0i32, 0usize, 0usize, 0u16);
        for i in 1..=m {
            // the row's cells within the genome stretch; on a tie, the genome length closest
            // to the read length
            let k_lo = band.saturating_sub(i);
            let mut row = (INF, 0usize);
            if i <= n + band {
                let k_hi = (n + band - i).min(w - 1);
                let best =
                    crate::gap::row_end(&sc.e1[i * w + k_lo..=i * w + k_hi], k_lo as i32 - band as i32);
                if best != u32::MAX {
                    row = ((best >> 16) as u16, i + k_lo + (best & 0xff) as usize - band);
                }
            }
            let s = 2 * i as i32 - row.0 as i32;
            if s > best.0 {
                best = (s, i, row.1, row.0);
            }
            if s + 60 < best.0 {
                break;
            }
        }
        (best.1 as u32, best.2 as u32, best.3 as u32)
    }
}

/// Mapping quality of a read's best hit, from the scores of the best two: 60 for a lone hit,
/// falling to 0 as the second approaches the first.
pub fn mapq(hits: &[SplicedHit]) -> u8 {
    match (hits.first(), hits.get(1)) {
        (Some(a), Some(b)) => {
            (60.0 * (1.0 - b.score.max(0) as f32 / a.score.max(1) as f32)).clamp(0.0, 60.0) as u8
        }
        _ => 60,
    }
}

/// PAF line of a hit (with a block-level `cg` tag: `M` per block, `I`/`D` for a length
/// difference inside it, `N` for introns). With `base`, the base-level CIGAR and edit count
/// from [`crate::base::Aligner::cigar`], `cg` holds that alignment and an `NM` tag is added.
pub fn paf_line(
    ix: &GenomeIndex,
    name: &str,
    qlen: u32,
    h: &SplicedHit,
    mapq: u8,
    base: Option<(&[u32], u32)>,
    out: &mut Vec<u8>,
) {
    use std::fmt::Write as _;
    use std::io::Write as _;
    let (first, last) = (h.blocks[0], h.blocks[h.blocks.len() - 1]);
    let (rid, ts) = ix.locate(first.ts);
    let off = ix.offsets()[rid];
    let tlen = ix.offsets()[rid + 1] - off;
    let (qs, qe) = if h.reverse {
        (qlen - last.qe, qlen - first.qs)
    } else {
        (first.qs, last.qe)
    };
    let mut matches = 0u32;
    let mut cg = String::new();
    for (i, b) in h.blocks.iter().enumerate() {
        if i > 0 {
            let _ = write!(cg, "{}N", b.ts - h.blocks[i - 1].te);
        }
        let (lq, lt) = (b.qe - b.qs, b.te - b.ts);
        let mm = lq.min(lt);
        matches += mm;
        if mm > 0 {
            let _ = write!(cg, "{mm}M");
        }
        if lq > lt {
            let _ = write!(cg, "{}I", lq - lt);
        } else if lt > lq {
            let _ = write!(cg, "{}D", lt - lq);
        }
    }
    let mut nm = String::new();
    if let Some((ops, edits)) = base {
        cg.clear();
        let mut aligned = 0;
        for &op in ops.iter().filter(|&&op| op & 15 != 4) {
            let _ = write!(cg, "{}{}", op >> 4, b"MIDNS"[(op & 15) as usize] as char);
            aligned += if op & 15 == 0 { op >> 4 } else { 0 };
        }
        matches = aligned.saturating_sub(edits);
        let _ = write!(nm, "\tNM:i:{edits}");
    }
    let _ = writeln!(
        out,
        "{name}\t{qlen}\t{qs}\t{qe}\t{}\t{}\t{tlen}\t{ts}\t{}\t{matches}\t{}\t{mapq}\ttp:A:{}\tAS:i:{}\tcm:i:{}\tts:A:{}\tjo:Z:{}{nm}\tcg:Z:{cg}",
        if h.reverse { '-' } else { '+' },
        ix.names()[rid],
        last.te - off,
        last.te - first.ts,
        if h.primary { 'P' } else { 'S' },
        h.score,
        h.anchors,
        h.tx_strand as char,
        if h.origins.is_empty() {
            "-"
        } else {
            std::str::from_utf8(&h.origins).unwrap_or("-")
        },
    );
}

/// SAM header for the index's sequences. `program` is the version and command line of the
/// program writing the file, for the `@PG` line.
pub fn sam_header(ix: &GenomeIndex, program: Option<(&str, &str)>, out: &mut Vec<u8>) {
    use std::io::Write as _;
    out.extend_from_slice(b"@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
    for (i, n) in ix.names().iter().enumerate() {
        let _ = writeln!(out, "@SQ\tSN:{n}\tLN:{}", ix.seq_len(i));
    }
    out.extend_from_slice(b"@PG\tID:spearfish\tPN:spearfish");
    if let Some((version, command_line)) = program {
        // (a tab or a line break would end the field)
        let clean = |s: &str| s.replace(['\t', '\n', '\r'], " ");
        let _ = write!(out, "\tVN:{}", clean(version));
        if !command_line.is_empty() {
            let _ = write!(out, "\tCL:{}", clean(command_line));
        }
    }
    out.push(b'\n');
}

/// SAM record of a hit, with a block-level CIGAR (`M` per block, `I`/`D` for a length
/// difference inside it, `N` for introns, `S` for unaligned read ends). `read` is the read as
/// sequenced; secondary records carry no sequence.
pub fn sam_line(
    ix: &GenomeIndex,
    name: &str,
    read: &[u8],
    h: &SplicedHit,
    mapq: u8,
    base: Option<(&[u32], u32)>,
    out: &mut Vec<u8>,
) {
    use std::io::Write as _;
    let qlen = read.len() as u32;
    let (first, last) = (h.blocks[0], h.blocks[h.blocks.len() - 1]);
    let (rid, ts) = ix.locate(first.ts);
    let flag = if h.reverse { 16 } else { 0 } | if h.primary { 0 } else { 256 };
    let _ = write!(out, "{name}\t{flag}\t{}\t{}\t{mapq}\t", ix.names()[rid], ts + 1);
    if let Some((ops, _)) = base {
        for &op in ops {
            let _ = write!(out, "{}{}", op >> 4, b"MIDNS"[(op & 15) as usize] as char);
        }
    } else {
        if first.qs > 0 {
            let _ = write!(out, "{}S", first.qs);
        }
        for (i, b) in h.blocks.iter().enumerate() {
            if i > 0 {
                let _ = write!(out, "{}N", b.ts - h.blocks[i - 1].te);
            }
            let (lq, lt) = (b.qe - b.qs, b.te - b.ts);
            let mm = lq.min(lt);
            if mm > 0 {
                let _ = write!(out, "{mm}M");
            }
            if lq > lt {
                let _ = write!(out, "{}I", lq - lt);
            } else if lt > lq {
                let _ = write!(out, "{}D", lt - lq);
            }
        }
        if last.qe < qlen {
            let _ = write!(out, "{}S", qlen - last.qe);
        }
    }
    out.extend_from_slice(b"\t*\t0\t0\t");
    if !h.primary {
        out.push(b'*');
    } else if h.reverse {
        out.extend(read.iter().rev().map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            _ => b'N',
        }));
    } else {
        out.extend(read.iter().map(|&b| b.to_ascii_uppercase()));
    }
    let _ = write!(out, "\t*\tAS:i:{}\tcm:i:{}", h.score, h.anchors);
    if let Some((_, edits)) = base {
        let _ = write!(out, "\tNM:i:{edits}");
    }
    if h.tx_strand != b'.' {
        let _ = write!(out, "\tXS:A:{}", h.tx_strand as char);
    }
    out.push(b'\n');
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The prefix table as first written, a cell at a time: what the row kernels must reproduce.
    #[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
    fn reference(
        a: impl Fn(usize) -> u8,
        m: usize,
        b: impl Fn(usize) -> u8,
        n: usize,
        band: usize,
        e: &mut Vec<u16>,
        tmp: &mut Vec<u16>,
        drop: i32,
    ) -> usize {
        let w = 2 * band + 1;
        e.clear();
        e.resize((m + 1) * w, INF);
        // gap states of the row above (a gap that consumes `a`), per k
        tmp.clear();
        tmp.resize(w, INF);
        // the bases of `b` once each, not once per cell, and only as far as the rows get
        let mut bb: Vec<u8> = Vec::with_capacity(n.min(m + band));
        let (open_ins, open_del) = (GAP_OPEN + INS_EXT, GAP_OPEN + GAP_EXT);
        let mut best = 0i32;
        for i in 0..=m {
            // the row's cells: j = i + k - band within 0..=n
            let k_lo = band.saturating_sub(i);
            if i > n + band {
                // past the end of `b`: nothing reachable in this or any later row
                if drop > 0 {
                    return i;
                }
                tmp.fill(INF);
                continue;
            }
            let k_hi = (n + band - i).min(w - 1);
            while bb.len() < (i + band).min(n) {
                bb.push(b(bb.len()));
            }
            let ai = if i > 0 { a(i - 1) } else { 0 };
            let (above, here) = e.split_at_mut(i * w);
            let above = if i > 0 { &above[(i - 1) * w..] } else { &above[..0] };
            let here = &mut here[..w];
            let mut left = INF; // a gap that consumes `b`, ending at the cell to the left
            let mut row = INF;
            for k in k_lo..=k_hi {
                let j = i + k - band;
                // cell above is (i-1, j): k + 1 in the previous row
                let up = if i > 0 && k + 1 < w {
                    (tmp[k + 1] + INS_EXT).min(above[k + 1] + open_ins)
                } else {
                    INF
                };
                left = if j > 0 && k > k_lo {
                    (left + GAP_EXT).min(here[k - 1] + open_del)
                } else {
                    INF
                };
                let v = if i == 0 && j == 0 {
                    0
                } else {
                    let diag = if i > 0 && j > 0 {
                        above[k] + if ai != bb[j - 1] { MISMATCH } else { 0 }
                    } else {
                        INF
                    };
                    diag.min(up).min(left)
                };
                here[k] = v.min(INF);
                tmp[k] = up.min(INF);
                left = left.min(INF);
                row = row.min(v);
            }
            tmp[..k_lo].fill(INF);
            tmp[k_hi + 1..].fill(INF);
            if drop > 0 {
                let s = 2 * i as i32 - row as i32;
                best = best.max(s);
                if s + drop < best {
                    return i;
                }
            }
        }
        m
    }

    fn rng(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    #[test]
    fn prefix_table_matches_reference() {
        let mut st = 0x9E37_79B9_7F4A_7C15u64;
        for round in 0..400 {
            let m = (rng(&mut st) % 90) as usize;
            let band = 1 + (rng(&mut st) % 30) as usize;
            // related sequences, so that paths with gaps and mismatches are both taken
            let a: Vec<u8> = (0..m).map(|_| (rng(&mut st) % 4) as u8).collect();
            let mut b: Vec<u8> = Vec::new();
            for &x in &a {
                match rng(&mut st) % 12 {
                    0 => {}
                    1 => b.extend([x, (rng(&mut st) % 4) as u8]),
                    2 => b.push((rng(&mut st) % 4) as u8),
                    _ => b.push(x),
                }
            }
            b.truncate((rng(&mut st) % 130) as usize);
            let n = b.len();
            let drop = if round % 3 == 0 { 20 } else { 0 };
            let (mut e1, mut t1, mut e2, mut t2) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
            let r1 = band_prefix(|i| a[i], m, |j| b[j], n, band, &mut e1, &mut t1, drop);
            let r2 = reference(|i| a[i], m, |j| b[j], n, band, &mut e2, &mut t2, drop);
            assert_eq!(r1, r2, "rows, round {round}");
            let w = 2 * band + 1;
            let rows = if drop > 0 { (r1 + 1).min(m + 1) } else { m + 1 };
            assert_eq!(
                e1[..rows * w],
                e2[..rows * w],
                "cells, round {round} m {m} n {n} band {band}"
            );
        }
    }

    /// The gap pass of a row, 32 cells at a time, against the cell-by-cell one.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn ids_among_matches_plain() {
        #[cfg(target_arch = "x86_64")]
        {
            if !std::is_x86_feature_detected!("avx512bw") {
                return;
            }
            let mut state = 7u64;
            for len in [0usize, 1, 31, 32, 33, 64, 1000] {
                let ids: Vec<u16> = (0..len).map(|_| (rng(&mut state) % 40) as u16).collect();
                for n in 1..=8usize {
                    let want: Vec<u16> = (0..n).map(|_| (rng(&mut state) % 40) as u16).collect();
                    let plain: Vec<usize> = (0..len).filter(|&k| want.contains(&ids[k])).collect();
                    let mut got = Vec::new();
                    // SAFETY: the feature was just detected.
                    unsafe { ids_among_avx512(&ids, &want, |k| got.push(k)) };
                    assert_eq!(got, plain, "len {len}, {n} values");
                }
            }
        }
    }

    #[test]
    fn scans_match() {
        // a random genome, and pieces that are noisy copies of stretches of it
        let mut st = 7u64;
        let g: Vec<u8> = (0..6100).map(|_| (rng(&mut st) % 4) as u8).collect();
        let word = |p: u32| (0..32).fold(0u64, |w, i| w | (g[p as usize + i] as u64) << (2 * i));
        for round in 0..20 {
            let at = 100 + (rng(&mut st) % 5000) as usize;
            let mut pieces = Vec::new();
            let mut peq = [[0u64; 4]; SCAN_PIECES];
            for k in 0..(1 + round % SCAN_PIECES) {
                let len = [32usize, 48, 64][k % 3];
                for i in 0..len {
                    let c = g[at + 7 * k + i];
                    let c = if rng(&mut st).is_multiple_of(8) {
                        (c + 1) & 3
                    } else {
                        c
                    };
                    peq[k][c as usize] |= 1 << i;
                }
                pieces.push((0, len, len as u32 / 4));
            }
            let (mut a, mut b) = (Vec::new(), Vec::new());
            scan_plain(&word, 20, 5990, &pieces, &peq, &mut a);
            scan(&word, 20, 5990, &pieces, &peq, &mut b);
            a.sort_unstable();
            b.sort_unstable();
            assert!(!a.is_empty());
            assert_eq!(a, b);
        }
    }

    #[test]
    fn row_gaps_match() {
        if !std::is_x86_feature_detected!("avx512bw") {
            return;
        }
        let mut st = 0x2545_F491_4F6C_DD1Du64;
        for round in 0..2000 {
            let len = 1 + (rng(&mut st) % 150) as usize;
            // mostly close costs with a few unreachable cells, or anything up to INF
            let row: Vec<u16> = (0..len)
                .map(|_| match (round % 2, rng(&mut st) % 10) {
                    (0, 0) => INF,
                    (0, _) => (rng(&mut st) % 40) as u16,
                    _ => (rng(&mut st) % (INF as u64 + 1)) as u16,
                })
                .collect();
            let (mut x, mut y) = (row.clone(), row.clone());
            let lx = prefix_left_plain(&mut x);
            // SAFETY: the feature was just detected.
            let ly = unsafe { prefix_left_avx512(&mut y) };
            assert_eq!(x, y, "cells, round {round} {row:?}");
            assert_eq!(lx, ly, "lowest, round {round} {row:?}");
        }
    }
}
