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

pub struct Mapper<'a> {
    pub index: &'a GenomeIndex,
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
/// Splice-site costs per site: GT..AG is free, GC..AG and AT..AC cost `SEMI`, others `NONCANON`.
const SEMI: u16 = 8;
const NONCANON: u16 = 18;
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

/// Cost of a splice site unknown to the annotation, when one was given (on top of its motif).
const UNKNOWN_SITE: u16 = 5;

/// Motif cost of an intron's left end at `p` (its first base) and of its right end at `p` (the
/// first exon base after it), for a transcript on the forward or reverse strand. On the
/// reverse strand the left end is the reverse complement of the acceptor.
/// Small cost of a canonical site whose neighbouring base is not the usual one: GT is mostly
/// followed by A or G, and AG mostly preceded by C or T (as minimap2's `--splice-flank`).
const FLANK: u16 = 5;

#[inline]
pub(crate) fn donor_motif(ix: &GenomeIndex, p: u32, minus: bool) -> u16 {
    if p + 3 > ix.total_len() || p < 1 {
        return NONCANON;
    }
    match (minus, ix.base(p), ix.base(p + 1)) {
        // GT, then A or G
        (false, 2, 3) => {
            if matches!(ix.base(p + 2), 0 | 2) {
                0
            } else {
                FLANK
            }
        }
        // CT (= AG on the other strand), then A or G (= C or T before the AG)
        (true, 1, 3) => {
            if matches!(ix.base(p + 2), 0 | 2) {
                0
            } else {
                FLANK
            }
        }
        (false, 2, 1) | (false, 0, 3) | (true, 2, 3) => SEMI, // GC, AT; GT (= AC)
        _ => NONCANON,
    }
}

#[inline]
pub(crate) fn acceptor_motif(ix: &GenomeIndex, p: u32, minus: bool) -> u16 {
    if p < 3 {
        return NONCANON;
    }
    match (minus, ix.base(p - 2), ix.base(p - 1)) {
        // C or T, then AG
        (false, 0, 2) => {
            if matches!(ix.base(p - 3), 1 | 3) {
                0
            } else {
                FLANK
            }
        }
        // C or T, then AC (= GT on the other strand, followed there by A or G)
        (true, 0, 1) => {
            if matches!(ix.base(p - 3), 1 | 3) {
                0
            } else {
                FLANK
            }
        }
        (false, 0, 1) | (true, 2, 1) | (true, 0, 3) => SEMI, // AC; GC, AT (= GC, AT donors)
        _ => NONCANON,
    }
}

/// Cost of an intron's left end at `p`: nothing for a known site; otherwise the motif cost,
/// and a little more when junctions were given and this is not one of them.
#[inline]
pub(crate) fn donor_pen(ix: &GenomeIndex, p: u32, minus: bool) -> u16 {
    match &ix.junctions {
        Some(j) if j.left_site(p, minus) => 0,
        Some(_) => donor_motif(ix, p, minus) + UNKNOWN_SITE,
        None => donor_motif(ix, p, minus),
    }
}

/// Cost of an intron's right end at `p`, as [`donor_pen`].
#[inline]
pub(crate) fn acceptor_pen(ix: &GenomeIndex, p: u32, minus: bool) -> u16 {
    match &ix.junctions {
        Some(j) if j.right_site(p, minus) => 0,
        Some(_) => acceptor_motif(ix, p, minus) + UNKNOWN_SITE,
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

/// Reference skips shorter than this are written as deletions, not introns.
const MIN_INTRON_OUT: u32 = 30;

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
fn windows_of(
    hits: &[Anchor],
    q_lo: u32,
    q_hi: u32,
    t_lo: u32,
    t_hi: u32,
    max_windows: usize,
) -> Vec<Window> {
    // groups as (first hit, last hit, hits)
    let mut groups: Vec<(Anchor, Anchor, u32)> = Vec::new();
    for &h in hits {
        let mut joined = false;
        for g in groups.iter_mut().rev().take(4) {
            let l = g.1;
            if h.q > l.q && h.t > l.t && h.t - l.t < 2000 {
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

impl<'a> Mapper<'a> {
    pub fn new(index: &'a GenomeIndex, opt: MapOptions) -> Self {
        Mapper {
            index,
            opt,
            seeder: BaseSeeder::new(L3_SPEC, index.hpc()).expect("valid seed spec"),
            seeder2: BaseSeeder::new(L2_SPEC, index.hpc()).expect("valid seed spec"),
        }
    }

    /// Hits of `read` (ASCII bases), best first. Empty when the read does not map.
    pub fn map(&self, read: &[u8], sc: &mut Scratch) -> Vec<SplicedHit> {
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
        let ix = self.index;
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
        let ix = self.index;
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
        let ix = self.index;
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
        for (score, reverse, chain) in cands.into_iter().take(4) {
            if score < 0.7 * top {
                break;
            }
            let read = std::mem::take(&mut sc.codes[reverse as usize]);
            sc.gap.reset(); // another chain, possibly the other orientation of the read
            out.push(self.build(reverse, &chain, &read, n as u32, score, screened_ok, sc));
            sc.codes[reverse as usize] = read;
        }
        // Rival hits, and a lone hit resting on few anchors, are scored by aligning their blocks;
        // a lone well-anchored hit keeps its chain score and is never aligned inside its exons.
        if out.len() > 1 || out.first().is_some_and(|h| h.anchors <= 3) {
            for h in out.iter_mut() {
                let read = std::mem::take(&mut sc.codes[h.reverse as usize]);
                h.score = self.align_score(&read, h, sc);
                sc.codes[h.reverse as usize] = read;
                sc.stats.scored += 1;
            }
            out.retain(|h| h.score >= self.opt.min_score);
            // A spliced read also aligns, unspliced, to a processed pseudogene of its gene, and
            // there it can cover a few more bases (no read ends lost at junctions). When a
            // spliced hit is within 10% of the best unspliced one and aligns at least as well
            // per base, it is the gene the read came from.
            let cover = |h: &SplicedHit| h.blocks.iter().map(|b| b.qe - b.qs).sum::<u32>().max(1) as i64;
            if let Some(u) = out.iter().filter(|h| h.blocks.len() == 1).max_by_key(|h| h.score) {
                let (us, uc) = (u.score as i64, cover(u));
                for h in out.iter_mut().filter(|h| h.blocks.len() > 1) {
                    let (hs, hc) = (h.score as i64, cover(h));
                    if hs < us && 10 * hs >= 9 * us && hs * uc >= us * hc {
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
        let usual = if rates.len() >= 3 {
            rates[rates.len() / 2]
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
        let mut best: Option<(u32, bool, Vec<Placed>)> = None;
        for &minus in &strands {
            let mut total = 0u32;
            let mut juncs: Vec<Placed> = Vec::new();
            for (k, &(a, b, p, poor)) in pairs.iter().enumerate() {
                let direct = p[minus as usize];
                if !poor {
                    total += direct.cost as u32;
                    juncs.push(direct);
                    continue;
                }
                let (q0, t0, q1, t1) = self.fixed(a, b);
                if let Some(t) = &tables[k] {
                    sc.gap.restore(t);
                }
                let hard = {
                    let c = p[0].cost.min(p[1].cost);
                    c >= INF || c as u32 >= 18 + q1.saturating_sub(q0) * rate / 64
                };
                if cand[k].is_none() {
                    sc.stats.rescues_tried += 1;
                    let path = self.rescue(reverse, Some(a), Some(b), qlen, sc);
                    sc.stats.empty += path.is_empty() as u64;
                    cand[k] = Some(windows_of(&path, q0, q1, t0, t1, 8));
                }
                let wins = cand[k].as_ref().unwrap();
                let m = q1.saturating_sub(q0) as usize;
                let mut done = false;
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
                    // what a junction of this length usually costs on this read: an answer above it
                    // still leaves something unexplained
                    let fine = 18 + m as u32 * self.opt.trigger_x64.max(usual) / 64;
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
                        vw.sort_by_key(|w| w.0);
                        let t = gap::solve(
                            self.index,
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
                        let limit = if found {
                            s.cost
                        } else {
                            (direct.cost as u32).min(1 << 20)
                        };
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
                        total += s.cost;
                        juncs.extend(s.hops.iter().map(|h| Placed {
                            q: h.q,
                            donor: h.leave,
                            acceptor: h.enter,
                            cost: 0,
                            origin: b'w',
                        }));
                        done = true;
                    }
                }
                // No level-2 hit led to a cheaper explanation: look for an exon between a pair of
                // splice sites that matches the unexplained bases.
                if !done && hard && q1 > q0 && t1 > t0 && m <= 600 {
                    let sites = self.site_windows(read, q0, t0, q1, t1, minus, sc);
                    if !sites.is_empty() {
                        sc.stats.site_tried += (minus == strands[0]) as u64;
                        let vw: Vec<SolverWindow> = sites.iter().map(|w| w.forward(q0, t0)).collect();
                        let rv = ReadView { q: q0, rev: false };
                        let gv = View {
                            origin: t0,
                            rev: false,
                        };
                        let half = self.opt.min_intron / 2;
                        let s = gap::solve(
                            self.index,
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
                        if s.hops.len() > 1 && s.cost < (direct.cost as u32).min(1 << 20) {
                            sc.stats.site_kept += (minus == strands[0]) as u64;
                            total += s.cost;
                            juncs.extend(s.hops.iter().map(|h| Placed {
                                q: h.q,
                                donor: h.leave,
                                acceptor: h.enter,
                                cost: 0,
                                origin: b's',
                            }));
                            done = true;
                        }
                    }
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
        let (_, _, q1, t1) = self.fixed(last, last);
        let rescue = solid;
        let right = self.end(read, reverse, q1, t1, false, minus, undecided, rescue, qlen, sc);
        let (q0, t0, _, _) = self.fixed(first, first);
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
        }
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
        let ix = self.index;
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
            if rescue && !tail && m as u32 - reached >= if known { 4 } else { 10 } {
                sc.stats.end_tried += 1;
                // level-2 cores of the unexplained bases, within reach of the longest intron
                let here = Anchor {
                    t,
                    q,
                    ql: 0,
                    tl: 0,
                    q_in: 0,
                    q_out: 0,
                    w: 0.0,
                };
                let hits = if left {
                    self.rescue(reverse, None, Some(here), qlen, sc)
                } else {
                    self.rescue(reverse, Some(here), None, qlen, sc)
                };
                // windows start at the unexplained bases, not at the start of this round
                let slack = reached.saturating_sub(12);
                let wins = if left {
                    windows_of(
                        &hits,
                        q.saturating_sub(m as u32),
                        q - slack,
                        t.saturating_sub(self.opt.max_intron),
                        t,
                        4,
                    )
                } else {
                    let hi = t.saturating_add(self.opt.max_intron).min(ix.total_len());
                    windows_of(&hits, q + slack, q + m as u32, t, hi, 4)
                };
                let vw: Vec<SolverWindow> = if left {
                    wins.iter().rev().map(|w| w.backward(q, t)).collect()
                } else {
                    wins.iter().map(|w| w.forward(q, t)).collect()
                };
                let p_end = best.0.end.1.abs_diff(t);
                for &mi in &strands {
                    // plus exons that start at a splice site and match the unexplained bases
                    // Most introns are short: look nearby first, further only if that gives
                    // nothing better than the plain extension.
                    let mut s = Solved {
                        cost: u32::MAX,
                        ..Default::default()
                    };
                    let mut all = vw.clone();
                    all.extend(self.known_hops(gv, p_end, reached as usize, m - reached as usize));
                    let base = all.len();
                    all.extend(self.site_windows_end(read, rv, m, reached as usize, gv, p_end, mi, true, sc));
                    let _ = base;
                    if !all.is_empty() {
                        all.sort_by_key(|w| w.0);
                        s = gap::solve(ix, read, rv, m, gv, None, &all, mi, half, &mut sc.gap);
                    }
                    if s.cost == u32::MAX {
                        continue;
                    }
                    if std::env::var_os("GMAP_DEBUG").is_some() {
                        eprintln!(
                            "end left={left} q={q} t={t} m={m} minus={mi} windows={} plain={} cost={} hops={:?} end={:?}",
                            all.len(),
                            best.0.cost,
                            s.cost,
                            s.hops,
                            s.end
                        );
                    }
                    if !s.hops.is_empty()
                        && s.cost < best.0.cost
                        && s.cost + self.hop_margin(&s, mi) <= plain_cost
                    {
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
            if whole || (m as u32 - got >= 10 && !hopped) || got < 30 {
                break;
            }
            (q, t) = s.end;
        }
        out
    }

    /// Alignment score of a hit: per block, two per read base minus the banded affine cost of
    /// the block against its genome interval; splice sites off the canonical pair cost as in
    /// junction placement.
    fn align_score(&self, read: &[u8], h: &SplicedHit, sc: &mut Scratch) -> i32 {
        let ix = self.index;
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
            if m.abs_diff(n) > 64.max(m / 8) {
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
    ) -> Vec<Window> {
        let ix = self.index;
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
            return Vec::new();
        }
        // The unexplained bases, with a margin on both sides (the read can follow an intron a
        // few bases by chance), and the exon lengths that fit them.
        const MARGIN: usize = 4;
        let (up, down) = (up.saturating_sub(MARGIN), down.saturating_sub(MARGIN));
        let k = m - up - down;
        let mid = &read[q0 as usize + up..q1 as usize - down];
        let (lo, hi) = (t0 + up as u32 + 20, (t1 - down as u32).saturating_sub(20));
        if hi <= lo || k > 200 + 2 * MARGIN {
            return Vec::new();
        }
        let slack = 2 + k / 10;
        let (len_lo, len_hi) = (k.saturating_sub(2 * MARGIN + slack).max(4), k + slack);
        // an exon under ten bases is only looked for close by: further off, chance matches win
        if len_lo < 10 && hi - lo > 20_000 {
            return Vec::new();
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
        for p in starts {
            if p < lo {
                continue;
            }
            for len in len_lo..=len_hi {
                let e = p + len as u32;
                if e > hi || donor_motif(ix, e, minus) > FLANK {
                    continue;
                }
                // the exon against the unexplained bases at every offset it can have in them
                let allowed = len / 8;
                let mut best = usize::MAX;
                for o in 0..=k.saturating_sub(len) {
                    let mut d = 0;
                    for i in 0..len.min(k - o) {
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
        out
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
        if known { 4 } else { 12 }
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
                out.push((s - t0, e - s, 0, m as i64 + (e - s) as i64));
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
        let ix = self.index;
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
    /// longest intron's reach past the chain's first or last anchor. Returns every hit, in genome
    /// order.
    fn rescue(
        &self,
        reverse: bool,
        a: Option<Anchor>,
        b: Option<Anchor>,
        qlen: u32,
        sc: &mut Scratch,
    ) -> Vec<Anchor> {
        let ix = self.index;
        let (qa, ta) = match a {
            Some(a) => (a.q + a.ql as u32, a.t + a.tl as u32),
            None => (0, b.map_or(0, |b| b.t.saturating_sub(self.opt.max_intron))),
        };
        let (qb, tb) = match b {
            Some(b) => (b.q, b.t),
            None => (qlen, ta.saturating_add(self.opt.max_intron).min(ix.total_len())),
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
        if sc.want.is_empty() {
            sc.want.resize(1 << 16, u32::MAX);
        }
        // want[id] = index of the first read seed with that id (others are found by scanning on)
        for (k, s) in seeds[lo..hi].iter().enumerate().rev() {
            sc.want[s.id as usize] = (lo + k) as u32;
        }
        let (pos, ids) = ix.l2_window(ta, tb);
        sc.stats.window_bases += (tb - ta) as u64;
        let mut hits: Vec<Anchor> = Vec::new();
        for (&t, &id) in pos.iter().zip(ids) {
            let first = sc.want[id as usize];
            if first == u32::MAX {
                continue;
            }
            for s in seeds[first as usize..hi].iter().filter(|s| s.id == id) {
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
        let ix = self.index;
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
            let ix = self.index;
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
                self.index,
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
        let ix = self.index;
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
        let ix = self.index;
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

/// PAF line of a hit (with a block-level `cg` tag: `M` per block, `I`/`D` for a length
/// difference inside it, `N` for introns).
pub fn paf_line(ix: &GenomeIndex, name: &str, qlen: u32, h: &SplicedHit, mapq: u8, out: &mut String) {
    use std::fmt::Write;
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
    let _ = writeln!(
        out,
        "{name}\t{qlen}\t{qs}\t{qe}\t{}\t{}\t{tlen}\t{ts}\t{}\t{matches}\t{}\t{mapq}\ttp:A:{}\tAS:i:{}\tcm:i:{}\tts:A:{}\tjo:Z:{}\tcg:Z:{cg}",
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

/// SAM header for the index's sequences.
pub fn sam_header(ix: &GenomeIndex, out: &mut String) {
    use std::fmt::Write;
    out.push_str("@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
    for (i, n) in ix.names().iter().enumerate() {
        let _ = writeln!(out, "@SQ\tSN:{n}\tLN:{}", ix.offsets()[i + 1] - ix.offsets()[i]);
    }
    out.push_str("@PG\tID:spearfish\tPN:spearfish\n");
}

/// SAM record of a hit, with a block-level CIGAR (`M` per block, `I`/`D` for a length
/// difference inside it, `N` for introns, `S` for unaligned read ends). `read` is the read as
/// sequenced; secondary records carry no sequence.
pub fn sam_line(ix: &GenomeIndex, name: &str, read: &[u8], h: &SplicedHit, mapq: u8, out: &mut String) {
    use std::fmt::Write;
    let qlen = read.len() as u32;
    let (first, last) = (h.blocks[0], h.blocks[h.blocks.len() - 1]);
    let (rid, ts) = ix.locate(first.ts);
    let flag = if h.reverse { 16 } else { 0 } | if h.primary { 0 } else { 256 };
    let _ = write!(out, "{name}\t{flag}\t{}\t{}\t{mapq}\t", ix.names()[rid], ts + 1);
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
    out.push_str("\t*\t0\t0\t");
    if !h.primary {
        out.push('*');
    } else if h.reverse {
        out.extend(read.iter().rev().map(|&b| match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' | b'U' => 'A',
            _ => 'N',
        }));
    } else {
        out.extend(read.iter().map(|&b| b.to_ascii_uppercase() as char));
    }
    let _ = write!(out, "\t*\tAS:i:{}\tcm:i:{}", h.score, h.anchors);
    if h.tx_strand != b'.' {
        let _ = write!(out, "\tXS:A:{}", h.tx_strand as char);
    }
    out.push('\n');
}
