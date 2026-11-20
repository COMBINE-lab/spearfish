//! Spliced mapping of one read: seed, pick loci, chain across introns, place junctions.
//!
//! Nothing here aligns across an intron or through an exon's interior. Level-3 cores give
//! anchors; anchors are grouped into loci and chained with a cost that lets the reference
//! skip ahead cheaply; consecutive chain anchors that the reference separates by more than
//! `min_intron` become a junction, and each junction is placed by a small banded edit
//! distance from both sides of the read bases between the two anchors.

use seeds::BaseSeeder;

use crate::index::{GenomeIndex, L2_SPEC, L3_SPEC, code};

#[derive(Clone, Debug)]
pub struct MapOptions {
    /// Seeds with more occurrences than this are not used.
    pub cap: usize,
    pub max_intron: u32,
    /// A reference skip of at least this many bases between two chained anchors is a junction.
    pub min_intron: u32,
    /// Loci chained per read (best seed weight first).
    pub max_loci: usize,
    /// Only forward-strand transcripts (direct RNA): the read is the transcript's sense strand.
    pub forward_only: bool,
}

impl Default for MapOptions {
    fn default() -> Self {
        MapOptions {
            cap: 64,
            max_intron: 200_000,
            min_intron: 40,
            max_loci: 6,
            forward_only: false,
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
}

#[derive(Default, Clone, Copy, Debug)]
pub struct Stats {
    pub junctions: u64,
    pub rescues_tried: u64,
    pub rescues_kept: u64,
    pub window_bases: u64,
}

pub struct Mapper<'a> {
    pub index: &'a GenomeIndex,
    pub opt: MapOptions,
    seeder: BaseSeeder,
    seeder2: BaseSeeder,
}

const INF: u16 = u16::MAX / 4;

/// Alignment costs, relative to a match: a mismatch, opening a gap, each gap base. With
/// affine gaps one long insertion in a read is cheaper than the scattered matches a junction
/// placed after it would collect, which unit costs get wrong.
const MISMATCH: u16 = 6;
const GAP_OPEN: u16 = 6;
const GAP_EXT: u16 = 2;
/// Splice-site costs per site: GT..AG is free, GC..AG and AT..AC cost `SEMI`, others `NONCANON`.
const SEMI: u16 = 4;
const NONCANON: u16 = 9;
/// Expected cost per read base of a true alignment at long-read error rates (x 8).
const COST_PER_BASE_X8: u32 = 3;

/// Banded affine-gap costs of every prefix of `a` against prefixes of `b`.
/// `e[i * w + k]` is the cost of `a[..i]` against `b[..j]` with `j = i + k - band`,
/// `w = 2 * band + 1`. `tmp` holds the gap states.
fn band_prefix(
    a: impl Fn(usize) -> u8,
    m: usize,
    b: impl Fn(usize) -> u8,
    n: usize,
    band: usize,
    e: &mut Vec<u16>,
    tmp: &mut Vec<u16>,
) {
    let w = 2 * band + 1;
    e.clear();
    e.resize((m + 1) * w, INF);
    // gap states of the row above (a gap that consumes `a`), per k
    tmp.clear();
    tmp.resize(w, INF);
    for i in 0..=m {
        let mut left = INF; // a gap that consumes `b`, ending at the cell to the left
        for k in 0..w {
            let j = (i + k) as isize - band as isize;
            if j < 0 || j as usize > n {
                left = INF;
                continue;
            }
            let j = j as usize;
            // cell above is (i-1, j): k + 1 in the previous row
            let up = if i > 0 && k + 1 < w {
                (tmp[k + 1] + GAP_EXT).min(e[(i - 1) * w + k + 1] + GAP_OPEN + GAP_EXT)
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
        }
    }
}

/// Cost of a donor at `p` (first intron base) and of an acceptor ending at `p` (first exon
/// base after the intron), for a transcript on the forward or reverse strand. On the reverse
/// strand the intron's left end is the reverse complement of the acceptor.
#[inline]
fn donor_pen(ix: &GenomeIndex, p: u32, minus: bool) -> u16 {
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
fn acceptor_pen(ix: &GenomeIndex, p: u32, minus: bool) -> u16 {
    if p < 2 {
        return NONCANON;
    }
    match (minus, ix.base(p - 2), ix.base(p - 1)) {
        (false, 0, 2) | (true, 0, 1) => 0,                   // AG; AC (= GT)
        (false, 0, 1) | (true, 2, 1) | (true, 0, 3) => SEMI, // AC; GC, AT (= GC, AT donors)
        _ => NONCANON,
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
}

impl<'a> Mapper<'a> {
    pub fn new(index: &'a GenomeIndex, opt: MapOptions) -> Self {
        Mapper {
            index,
            opt,
            seeder: BaseSeeder::new(L3_SPEC, true).expect("valid seed spec"),
            seeder2: BaseSeeder::new(L2_SPEC, true).expect("valid seed spec"),
        }
    }

    /// Hits of `read` (ASCII bases), best first. Empty when the read does not map.
    pub fn map(&self, read: &[u8], sc: &mut Scratch) -> Vec<SplicedHit> {
        let ix = self.index;
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
        for o in 0..2 {
            let ascii: Vec<u8> = sc.codes[o].iter().map(|&c| b"ACGT"[c as usize]).collect();
            let h = self.seeder.parsed(&ascii).expect("HPC is on");
            let seeds = self.seeder.seed_parsed(&h.seq);
            sc.l2[o].clear();
            for s in self.seeder2.seed_parsed(&h.seq) {
                if let Some(id) = ix.l2_id(s.label) {
                    let (qs, qe) = h.to_raw(s.start, s.end);
                    sc.l2[o].push(L2Seed {
                        q: qs,
                        ql: (qe - qs) as u16,
                        q_in: h.run_len(s.start) as u16,
                        q_out: (h.to_raw(s.end - 1, s.end).0 - qs) as u16,
                        runs: (s.end - s.start) as u16,
                        id,
                    });
                }
            }
            sc.hits.clear();
            for s in &seeds {
                let (pos, len) = ix.occurrences(s.label);
                if pos.is_empty() || pos.len() > self.opt.cap {
                    continue;
                }
                let (qs, qe) = h.to_raw(s.start, s.end);
                let w = (qe - qs).min(40) as f32 / (1.0 + (pos.len() as f32).ln());
                for (&t, &tl) in pos.iter().zip(len) {
                    sc.hits.push(Anchor {
                        t,
                        q: qs,
                        ql: (qe - qs).min(u16::MAX as u32) as u16,
                        tl: tl as u16,
                        q_in: h.run_len(s.start) as u16,
                        q_out: (h.to_raw(s.end - 1, s.end).0 - qs) as u16,
                        w,
                    });
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
                if let Some(hit) = self.chain_locus(o == 1, s, e, n as u32, sc) {
                    out.push(hit);
                }
            }
        }
        out.sort_by_key(|h| std::cmp::Reverse(h.score));
        if let Some(h) = out.first_mut() {
            h.primary = true;
        }
        out
    }

    /// Chain the hits `sc.hits[s..e]` (sorted by genome position) and build the hit.
    fn chain_locus(
        &self,
        reverse: bool,
        s: usize,
        e: usize,
        qlen: u32,
        sc: &mut Scratch,
    ) -> Option<SplicedHit> {
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
        if best < 40.0 && sc.chain.len() < 2 {
            return None;
        }
        let mut chain: Vec<Anchor> = sc.chain.iter().map(|&i| hits[i as usize]).collect();
        let read = std::mem::take(&mut sc.codes[reverse as usize]);
        let hit = self.build(reverse, &mut chain, &read, qlen, best, sc);
        sc.codes[reverse as usize] = read;
        Some(hit)
    }

    /// Blocks from a chain: junctions placed, read ends extended.
    fn build(
        &self,
        reverse: bool,
        chain: &mut Vec<Anchor>,
        read: &[u8],
        qlen: u32,
        score: f32,
        sc: &mut Scratch,
    ) -> SplicedHit {
        let ix = self.index;
        // Exons without a level-3 anchor: where a junction leaves read bases unexplained, look
        // for level-2 cores of those bases inside that intron and splice the found anchors in.
        let mut i = 0;
        while i + 1 < chain.len() {
            let (a, b) = (chain[i], chain[i + 1]);
            i += 1;
            let d = (b.t - a.t) as i64 - (b.q - a.q) as i64;
            if d < self.opt.min_intron as i64 {
                continue;
            }
            sc.stats.junctions += 1;
            let (q0, t0, q1, t1) = self.fixed(a, b);
            let old = self.place(read, q0, t0, q1, t1, sc);
            let old = old[0].cost.min(old[1].cost);
            let m = q1.saturating_sub(q0);
            if old < INF && (old as u32) < 24 + m * COST_PER_BASE_X8 / 8 {
                continue;
            }
            sc.stats.rescues_tried += 1;
            let sub = self.rescue(reverse, Some(a), Some(b), qlen, sc);
            if sub.is_empty() {
                continue;
            }
            // Keep the found anchors only if the junctions they imply cost less.
            let mut new = 0u32;
            let mut prev = a;
            for &x in sub.iter().chain(std::iter::once(&b)) {
                let d = (x.t - prev.t) as i64 - (x.q - prev.q) as i64;
                if d >= self.opt.min_intron as i64 {
                    let (q0, t0, q1, t1) = self.fixed(prev, x);
                    let p = self.place(read, q0, t0, q1, t1, sc);
                    new += p[0].cost.min(p[1].cost) as u32;
                }
                prev = x;
            }
            let mut path = vec![a];
            path.extend_from_slice(&sub);
            path.push(b);
            if new + 12 < old as u32 && self.junctions_ok(read, &path, sc) {
                sc.stats.rescues_kept += 1;
                let n = sub.len();
                chain.splice(i..i, sub);
                i += n;
            }
        }
        // Read ends: when the bases past the last (before the first) anchor do not continue on
        // the genome, look for further exons the same way, one known side only.
        for _ in 0..4 {
            let last = chain[chain.len() - 1];
            let (q, t) = (last.q + last.ql as u32, last.t + last.tl as u32);
            if qlen - q < 20 || qlen - q - self.extend(read, q, t, false, sc).0 < 12 {
                break;
            }
            sc.stats.rescues_tried += 1;
            let mut path = vec![last];
            path.extend(self.rescue(reverse, Some(last), None, qlen, sc));
            if path.len() < 2 || !self.junctions_ok(read, &path, sc) {
                break;
            }
            sc.stats.rescues_kept += 1;
            chain.extend_from_slice(&path[1..]);
        }
        for _ in 0..4 {
            let first = chain[0];
            if first.q < 20 || first.q - self.extend(read, first.q, first.t, true, sc).0 < 12 {
                break;
            }
            sc.stats.rescues_tried += 1;
            let mut path = self.rescue(reverse, None, Some(first), qlen, sc);
            path.push(first);
            if path.len() < 2 || !self.junctions_ok(read, &path, sc) {
                break;
            }
            sc.stats.rescues_kept += 1;
            path.pop();
            chain.splice(0..0, path);
        }
        let chain: &[Anchor] = chain;
        // Junctions between consecutive anchors, placed for both transcript strands.
        let mut junc: Vec<[Placed; 2]> = Vec::new();
        let mut cost = [0u32; 2];
        for w in chain.windows(2) {
            let (a, b) = (w[0], w[1]);
            let d = (b.t - a.t) as i64 - (b.q - a.q) as i64;
            if d < self.opt.min_intron as i64 {
                continue;
            }
            let (q0, t0, q1, t1) = self.fixed(a, b);
            let p = self.place(read, q0, t0, q1, t1, sc);
            cost[0] += p[0].cost as u32;
            cost[1] += p[1].cost as u32;
            junc.push(p);
        }
        let minus = if self.opt.forward_only {
            reverse
        } else {
            cost[1] < cost[0]
        };
        let tx_strand = if junc.is_empty() || (!self.opt.forward_only && cost[0] == cost[1]) {
            b'.'
        } else if minus {
            b'-'
        } else {
            b'+'
        };
        let mut edits = 0u32;
        let mut blocks: Vec<Block> = Vec::new();
        let first = chain[0];
        let (mut bq, mut bt) = (first.q, first.t);
        for p in &junc {
            let p = p[(minus && !self.opt.forward_only) as usize];
            edits += p.cost.min(1000) as u32;
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

        // Read ends: extend while the bases keep matching.
        let (dq, dt, e0) = self.extend(read, blocks[0].qs, blocks[0].ts, true, sc);
        blocks[0].qs -= dq;
        blocks[0].ts -= dt;
        let nb = blocks.len() - 1;
        let (dq, dt, e1) = self.extend(read, blocks[nb].qe, blocks[nb].te, false, sc);
        blocks[nb].qe += dq;
        blocks[nb].te += dt;
        edits += e0 + e1;
        debug_assert!(blocks[nb].qe <= qlen);
        blocks.retain(|b| b.te > b.ts || b.qe > b.qs);
        let _ = ix;
        SplicedHit {
            reverse,
            blocks,
            tx_strand,
            score: score as i32 - (edits / 3) as i32,
            anchors: chain.len() as u32,
            edits,
            primary: false,
        }
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
    /// longest intron's reach past the chain's first or last anchor. Returns the best colinear
    /// run, or nothing.
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
                let mut e = t;
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
        if hits.is_empty() {
            return hits;
        }
        // Best path a -> hits... -> b (hits are in genome order).
        let n = hits.len();
        let mut f = vec![0f32; n];
        let mut back = vec![u32::MAX; n];
        let (mut best, mut best_i) = (
            match (a, b) {
                (Some(a), Some(b)) => self.step(a, b),
                _ => 0.0,
            },
            u32::MAX,
        );
        for j in 0..n {
            let h = hits[j];
            f[j] = a.map_or(h.ql as f32, |a| self.step(a, h));
            for i in j.saturating_sub(40)..j {
                let p = hits[i];
                if p.q + p.ql as u32 > h.q || p.t + p.tl as u32 > h.t {
                    continue;
                }
                let v = f[i] + self.step(p, h);
                if v > f[j] {
                    f[j] = v;
                    back[j] = i as u32;
                }
            }
            let v = f[j] + b.map_or(0.0, |b| self.step(h, b));
            if v > best {
                best = v;
                best_i = j as u32;
            }
        }
        let mut path = Vec::new();
        let mut i = best_i;
        while i != u32::MAX {
            path.push(hits[i as usize]);
            i = back[i as usize];
        }
        path.reverse();
        path
    }

    /// Whether every junction along `path` (anchors in order) is placed cheaply.
    fn junctions_ok(&self, read: &[u8], path: &[Anchor], sc: &mut Scratch) -> bool {
        for w in path.windows(2) {
            let d = (w[1].t - w[0].t) as i64 - (w[1].q - w[0].q) as i64;
            if d >= self.opt.min_intron as i64 {
                let (q0, t0, q1, t1) = self.fixed(w[0], w[1]);
                let p = self.place(read, q0, t0, q1, t1, sc);
                let c = p[0].cost.min(p[1].cost);
                if c >= INF || c as u32 >= 18 + q1.saturating_sub(q0) * COST_PER_BASE_X8 / 8 {
                    return false;
                }
            }
        }
        true
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
        band_prefix(|i| r[i], m, |j| g1[j], n, band, &mut sc.e1, &mut sc.etmp);
        band_prefix(
            |i| r[m - 1 - i],
            m,
            |j| g2[n - 1 - j],
            n,
            band,
            &mut sc.e2,
            &mut sc.etmp,
        );
        let mut out = [none; 2];
        for (s, o) in out.iter_mut().enumerate() {
            let minus = s == 1;
            if minus && self.opt.forward_only {
                break;
            }
            let mut best = (INF, 0usize, 0usize, 0usize);
            for i in 0..=m {
                // upstream part: read[..i] against genome[t0 .. t0 + j1]
                let mut l = (INF, 0usize);
                for k in 0..w {
                    let j = (i + k) as isize - band as isize;
                    if j < 0 || j as usize > n {
                        continue;
                    }
                    let c = sc.e1[i * w + k] + donor_pen(ix, t0 + j as u32, minus);
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
                    let c = sc.e2[(m - i) * w + k] + acceptor_pen(ix, t1 - j as u32, minus);
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
        let band = 8 + m / 20;
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
        "{name}\t{qlen}\t{qs}\t{qe}\t{}\t{}\t{tlen}\t{ts}\t{}\t{matches}\t{}\t{mapq}\ttp:A:{}\tAS:i:{}\tcm:i:{}\tts:A:{}\tcg:Z:{cg}",
        if h.reverse { '-' } else { '+' },
        ix.names()[rid],
        last.te - off,
        last.te - first.ts,
        if h.primary { 'P' } else { 'S' },
        h.score,
        h.anchors,
        h.tx_strand as char,
    );
}
