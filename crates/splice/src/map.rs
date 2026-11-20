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
    f: Vec<f32>,
    back: Vec<u32>,
    chain: Vec<u32>,
    g1: Vec<u8>,
    g2: Vec<u8>,
    e1: Vec<u16>,
    e2: Vec<u16>,
    etmp: Vec<u16>,
    gap: crate::gap::Bufs,
    /// Level-2 cores of the read, per orientation, in read order.
    l2: [Vec<L2Seed>; 2],
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
    let w = 2 * band + 1;
    e.clear();
    e.resize((m + 1) * w, INF);
    // gap states of the row above (a gap that consumes `a`), per k
    tmp.clear();
    tmp.resize(w, INF);
    let mut best = 0i32;
    for i in 0..=m {
        let mut left = INF; // a gap that consumes `b`, ending at the cell to the left
        let mut row = INF;
        for k in 0..w {
            let j = (i + k) as isize - band as isize;
            if j < 0 || j as usize > n {
                left = INF;
                continue;
            }
            let j = j as usize;
            // cell above is (i-1, j): k + 1 in the previous row
            let up = if i > 0 && k + 1 < w {
                (tmp[k + 1] + INS_EXT).min(e[(i - 1) * w + k + 1] + GAP_OPEN + INS_EXT)
            } else {
                INF
            };
            left = if j > 0 && k > 0 {
                (left + GAP_EXT).min(e[i * w + k - 1] + GAP_OPEN + GAP_EXT)
            } else {
                INF
            };
            let v = if i == 0 && j == 0 {
                0
            } else {
                let diag = if i > 0 && j > 0 {
                    e[(i - 1) * w + k] + if a(i - 1) != b(j - 1) { MISMATCH } else { 0 }
                } else {
                    INF
                };
                diag.min(up).min(left)
            };
            e[i * w + k] = v.min(INF);
            tmp[k] = up.min(INF);
            row = row.min(v);
        }
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

/// Cost of a splice site unknown to the annotation, when one was given (on top of its motif).
const UNKNOWN_SITE: u16 = 5;

/// Motif cost of an intron's left end at `p` (its first base) and of its right end at `p` (the
/// first exon base after it), for a transcript on the forward or reverse strand. On the
/// reverse strand the left end is the reverse complement of the acceptor.
#[inline]
pub(crate) fn donor_motif(ix: &GenomeIndex, p: u32, minus: bool) -> u16 {
    if p + 2 > ix.total_len() {
        return NONCANON;
    }
    match (minus, ix.base(p), ix.base(p + 1)) {
        (false, 2, 3) | (true, 1, 3) => 0,                    // GT; CT (= AG)
        (false, 2, 1) | (false, 0, 3) | (true, 2, 3) => SEMI, // GC, AT; GT (= AC)
        _ => NONCANON,
    }
}

#[inline]
pub(crate) fn acceptor_motif(ix: &GenomeIndex, p: u32, minus: bool) -> u16 {
    if p < 2 {
        return NONCANON;
    }
    match (minus, ix.base(p - 2), ix.base(p - 1)) {
        (false, 0, 2) | (true, 0, 1) => 0,                   // AG; AC (= GT)
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
        let mut hits = self.map_with(read, self.opt.cap, false, sc);
        // A read that level 3 anchors weakly (no hit, three anchors or fewer, or under 60% of
        // the read covered) is seeded again with pairs of level-2 cores as well.
        let weak = hits.first().is_none_or(|h| h.anchors <= 3) || 10 * covered(&hits) < 6 * read.len() as u32;
        if weak && self.index.has_pairs() {
            sc.stats.retried += 1;
            let again = self.map_with(read, self.opt.cap, true, sc);
            if covered(&again) > covered(&hits) {
                hits = again;
            }
        }
        if !hits.is_empty() || self.opt.cap >= self.opt.cap_retry {
            return hits;
        }
        // Still nothing: allow the repetitive seeds too.
        self.map_with(read, self.opt.cap_retry, true, sc)
    }

    fn map_with(&self, read: &[u8], cap: usize, pairs: bool, sc: &mut Scratch) -> Vec<SplicedHit> {
        let ix = self.index;
        sc.gap.reset();
        let n = read.len();
        if n < 30 {
            return Vec::new();
        }
        // Both orientations as codes; other symbols read as A.
        sc.codes[0].clear();
        sc.codes[0].extend(read.iter().map(|&b| code(b) & 3));
        sc.codes[1].clear();
        let fw = std::mem::take(&mut sc.codes[0]);
        sc.codes[1].extend(fw.iter().rev().map(|&c| 3 - c));
        sc.codes[0] = fw;

        let mut out: Vec<SplicedHit> = Vec::new();
        let mut cands: Vec<(f32, bool, Vec<Anchor>)> = Vec::new();
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
            let seeds = self.seeder.seed_parsed(parsed);
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
                    });
                }
            }
            sc.hits.clear();
            for s in &seeds {
                let (pos, len) = ix.occurrences(s.label);
                if pos.is_empty() || pos.len() > cap {
                    continue;
                }
                let (qs, qe, q_in, q_out, _) = coords(s.start, s.end);
                let w = (qe - qs).min(40) as f32 / (1.0 + (pos.len() as f32).ln());
                for (&t, &tl) in pos.iter().zip(len) {
                    sc.hits.push(Anchor {
                        t,
                        q: qs,
                        ql: (qe - qs).min(u16::MAX as u32) as u16,
                        tl: tl as u16,
                        q_in,
                        q_out,
                        w,
                    });
                }
            }
            if pairs && ix.has_pairs() {
                for w in sc.l2[o].windows(2) {
                    let l = crate::index::pair_label(w[0].label, w[1].label);
                    if w[1].q - w[0].q >= 64 || !crate::index::keep_pair(l) {
                        continue;
                    }
                    let (pos, len) = ix.pair_occurrences(l);
                    if pos.is_empty() || pos.len() > cap {
                        continue;
                    }
                    let ql = w[1].q + w[1].ql as u32 - w[0].q;
                    let wt = ql.min(40) as f32 / (1.0 + (pos.len() as f32).ln());
                    for (&t, &tl) in pos.iter().zip(len) {
                        sc.hits.push(Anchor {
                            t,
                            q: w[0].q,
                            ql: ql as u16,
                            tl: tl as u16,
                            q_in: w[0].q_in,
                            q_out: (w[1].q - w[0].q) as u16 + w[1].q_out,
                            w: wt,
                        });
                    }
                }
            }
            if sc.hits.is_empty() {
                continue;
            }
            sc.hits.sort_unstable_by_key(|a| (a.t, a.q));
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
            loci.sort_unstable_by(|a, b| b.0.total_cmp(&a.0));
            let best = loci[0].0;
            for &(w, s, e) in loci.iter().take(self.opt.max_loci) {
                if w < 0.25 * best {
                    break;
                }
                if let Some((score, chain)) = self.chain_locus(s, e, sc) {
                    cands.push((score, o == 1, chain));
                }
            }
        }
        // Only chains close to the best are worth building.
        cands.sort_by(|a, b| b.0.total_cmp(&a.0));
        let top = cands.first().map_or(0.0, |c| c.0);
        for (score, reverse, chain) in cands.into_iter().take(4) {
            if score < 0.5 * top {
                break;
            }
            let read = std::mem::take(&mut sc.codes[reverse as usize]);
            sc.gap.reset(); // another chain, possibly the other orientation of the read
            out.push(self.build(reverse, &chain, &read, n as u32, score, sc));
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
    fn build(
        &self,
        reverse: bool,
        chain: &[Anchor],
        read: &[u8],
        qlen: u32,
        score: f32,
        sc: &mut Scratch,
    ) -> SplicedHit {
        let trigger = |c: u16, m: u32| c >= INF || c as u32 >= 18 + m / 4;
        // Each pair of chain neighbours that the genome separates by an intron, placed directly.
        let mut pairs: Vec<(Anchor, Anchor, [Placed; 2], bool)> = Vec::new();
        let mut vote = [0u32; 2];
        for w in chain.windows(2) {
            let (a, b) = (w[0], w[1]);
            let d = (b.t - a.t) as i64 - (b.q - a.q) as i64;
            if d < self.opt.min_intron as i64 {
                continue;
            }
            sc.stats.junctions += 1;
            let (q0, t0, q1, t1) = self.fixed(a, b);
            let p = self.place(read, q0, t0, q1, t1, sc);
            let poor = trigger(p[0].cost.min(p[1].cost), q1.saturating_sub(q0));
            if !poor {
                vote[0] += p[0].cost as u32;
                vote[1] += p[1].cost as u32;
            }
            pairs.push((a, b, p, poor));
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
                    for min_hits in [2u32, 0] {
                        if (min_hits == 2 && strong == 0 && known.is_empty())
                            || (min_hits == 0 && strong == wins.len() && s.cost != u32::MAX)
                        {
                            continue;
                        }
                        let mut vw: Vec<SolverWindow> = wins
                            .iter()
                            .filter(|w| w.hits >= min_hits)
                            .map(|w| w.forward(q0, t0))
                            .collect();
                        vw.extend(known.iter().copied());
                        vw.sort_by_key(|w| w.0);
                        s = gap::solve(
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
                            break;
                        }
                    }
                    if s.hops.len() > 1 && s.cost < (direct.cost as u32).min(1 << 20) {
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
                if !done && q1 > q0 && t1 > t0 && m <= 600 {
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
        // A chain of one or two anchors is too thin to go looking for more exons from.
        let rescue = chain.len() >= 3;
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
                    let mut tried: Vec<SolverWindow> = Vec::new();
                    for near in [true, false] {
                        let sites =
                            self.site_windows_end(read, rv, m, reached as usize, gv, p_end, mi, near, sc);
                        if !near && sites == tried {
                            break; // looking further found nothing new
                        }
                        all.truncate(base);
                        all.extend(sites.iter().copied());
                        tried = sites;
                        if all.is_empty() {
                            continue;
                        }
                        all.sort_by_key(|w| w.0);
                        s = gap::solve(ix, read, rv, m, gv, None, &all, mi, half, &mut sc.gap);
                        if !s.hops.is_empty() && s.cost + self.hop_margin(&s, mi) <= plain_cost {
                            break;
                        }
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
            let band = 8 + m / 20 + m.abs_diff(n);
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
                if e > hi || donor_motif(ix, e, minus) != 0 {
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
        for p in cands {
            if p < from || p >= hi || gv.enter_motif(ix, p, minus) != 0 {
                continue;
            }
            // the exon can start a few read bases into the probe
            for o in 0..=(2 * MARGIN).min(probe.len().saturating_sub(10)) {
                let len = (probe.len() - o).min(20);
                let allowed = len / 8;
                let mut d = 0;
                for i in 0..len {
                    d += (probe[o + i] != gv.base(ix, p + i as u32)) as usize;
                    if d > allowed {
                        break;
                    }
                }
                if d <= allowed {
                    let vs = p.saturating_sub(3);
                    let vl = (k as u32 + 12).min(300);
                    found.push((d, p, (vs, vl, (start + o) as i64 - (p - vs) as i64, 0)));
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
    /// strand (`[1]`): the edit tables are shared, only the splice-site penalties differ.
    fn place(&self, read: &[u8], q0: u32, t0: u32, q1: u32, t1: u32, sc: &mut Scratch) -> [Placed; 2] {
        let ix = self.index;
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
        let band = 4 + m / 10;
        let w = 2 * band + 1;
        let span = (t1 - t0) as usize;
        let n = (m + band).min(span);
        ix.bases(t0, t0 + n as u32, &mut sc.g1);
        ix.bases(t1 - n as u32, t1, &mut sc.g2);
        let r = &read[q0 as usize..q1 as usize];
        let (g1, g2) = (&sc.g1, &sc.g2);
        band_prefix(|i| r[i], m, |j| g1[j], n, band, &mut sc.e1, &mut sc.etmp, 0);
        band_prefix(
            |i| r[m - 1 - i],
            m,
            |j| g2[n - 1 - j],
            n,
            band,
            &mut sc.e2,
            &mut sc.etmp,
            0,
        );
        let mut out = [none; 2];
        for (s, o) in out.iter_mut().enumerate() {
            let minus = s == 1;
            let dpen: Vec<u16> = (0..=n).map(|j| donor_pen(ix, t0 + j as u32, minus)).collect();
            let apen: Vec<u16> = (0..=n).map(|j| acceptor_pen(ix, t1 - j as u32, minus)).collect();
            let mut best = (INF, 0usize, 0usize, 0usize);
            for i in 0..=m {
                // upstream part: read[..i] against genome[t0 .. t0 + j1]
                let mut l = (INF, 0usize);
                for k in 0..w {
                    let j = (i + k) as isize - band as isize;
                    if j < 0 || j as usize > n {
                        continue;
                    }
                    let c = sc.e1[i * w + k] + dpen[j as usize];
                    if c < l.0 {
                        l = (c, j as usize);
                    }
                }
                let mut rr = (INF, 0usize);
                for k in 0..w {
                    let j = (m - i + k) as isize - band as isize;
                    if j < 0 || j as usize > n {
                        continue;
                    }
                    let c = sc.e2[(m - i) * w + k] + apen[j as usize];
                    if c < rr.0 {
                        rr = (c, j as usize);
                    }
                }
                let c = l.0 + rr.0;
                if c < best.0 && l.1 + rr.1 + (self.opt.min_intron as usize / 2) <= span {
                    best = (c, i, l.1, rr.1);
                }
            }
            if std::env::var_os("GMAP_DEBUG").is_some() {
                let t = |v: &[u8]| v.iter().map(|&c| b"ACGT"[c as usize] as char).collect::<String>();
                eprintln!(
                    "place minus={minus} q0={q0} q1={q1} t0={t0} t1={t1} m={m} n={n} band={band} best={best:?}\n  read={}\n  g1  ={}\n  g2  ={}",
                    t(r),
                    t(&sc.g1),
                    t(&sc.g2)
                );
            }
            if best.0 < INF {
                *o = Placed {
                    q: q0 + best.1 as u32,
                    donor: t0 + best.2 as u32,
                    acceptor: t1 - best.3 as u32,
                    cost: best.0,
                    origin: b'd',
                };
            }
        }
        out
    }

    /// Extend from read position `q` / genome position `t` towards the read's start (`left`) or
    /// end. Returns the read and genome bases added and the edits among them.
    fn extend(&self, read: &[u8], q: u32, t: u32, left: bool, sc: &mut Scratch) -> (u32, u32, u32) {
        let ix = self.index;
        let m = if left { q as usize } else { read.len() - q as usize };
        if m == 0 {
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
            let mut row = (INF, 0usize);
            for k in 0..w {
                let j = (i + k) as isize - band as isize;
                if j < 0 || j as usize > n {
                    continue;
                }
                // on a tie, the genome length closest to the read length
                let c = sc.e1[i * w + k];
                if c < row.0 || (c == row.0 && j.unsigned_abs().abs_diff(i) < row.1.abs_diff(i)) {
                    row = (c, j as usize);
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
