//! Read seeding, anchor lookup and candidate-gene selection (design doc, "Read processing",
//! stages 1–4, single-level seeding per D12b).

use index::{FLAG_MASKED, Index, NO_FAMILY};
use seeds::Seed;

/// A read seed matched to one index node.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Anchor {
    /// Query span `[qs, qe)` in parsed read coordinates.
    pub qs: u32,
    pub qe: u32,
    pub node: u32,
    pub weight: f32,
}

#[derive(Clone, Debug)]
pub struct AnchorParams {
    /// Labels found in more genes than this are ignored for the read.
    pub h_max: usize,
    /// Span cap in the anchor weight.
    pub s_max: u32,
    /// Gene count at which an anchor's weight reaches zero.
    pub g_max: f32,
    /// Keep genes scoring at least this fraction of the best.
    pub alpha_fam: f32,
    /// Keep this many genes, plus any scoring within `f_tie` (relative) of the last kept.
    pub f_keep: usize,
    pub f_tie: f32,
    /// Hard cap on kept genes (ties in large multi-copy families).
    pub f_max: usize,
}

impl Default for AnchorParams {
    fn default() -> Self {
        AnchorParams {
            h_max: 1000,
            s_max: 100,
            // Effectively off: at genome scale, down-weighting multi-gene seeds hid the true
            // locus of reads from multi-copy and repeat-rich genes (reports/genome-v0.md).
            g_max: 1e9,
            alpha_fam: 0.5,
            f_keep: 8,
            f_tie: 0.01,
            f_max: 64,
        }
    }
}

/// Anchors of one candidate gene (family), sorted by (qs, qe, node).
#[derive(Clone, Debug)]
pub struct Candidate {
    pub family: u32,
    pub score: f32,
    pub anchors: Vec<Anchor>,
}

/// Look up read seeds and return all anchors (unsorted).
pub fn anchors(idx: &Index, seeds: &[Seed], p: &AnchorParams) -> Vec<Anchor> {
    let flags = idx.node_flags();
    let mut out = Vec::with_capacity(seeds.len());
    // Lookups miss in cache twice (directory entry, then labels): request the directory entry
    // 2 * AHEAD seeds ahead and the labels AHEAD seeds ahead.
    const AHEAD: usize = 8;
    for (i, s) in seeds.iter().enumerate() {
        if let Some(x) = seeds.get(i + 2 * AHEAD) {
            idx.prefetch_dir(x.label);
        }
        if let Some(x) = seeds.get(i + AHEAD) {
            idx.prefetch_labels(x.label);
        }
        let r = idx.lookup(s.label);
        let genes = r.len();
        if genes == 0 || genes > p.h_max {
            continue;
        }
        let specificity = 1.0 - (genes as f32).ln() / p.g_max.ln();
        let weight = (s.end - s.start).min(p.s_max) as f32 * specificity.max(0.0);
        for node in r {
            if flags[node as usize] & FLAG_MASKED == 0 {
                idx.prefetch_node(node);
                out.push(Anchor {
                    qs: s.start,
                    qe: s.end,
                    node,
                    weight,
                });
            }
        }
    }
    out
}

/// Sum of weights over a greedy, query-increasing, non-overlapping subset of anchors (sorted by
/// query start): the cheap family score of the design doc.
fn greedy_score(anchors: &[Anchor]) -> f32 {
    let mut end = 0;
    let mut score = 0.0;
    for a in anchors {
        if a.qs >= end {
            score += a.weight;
            end = a.qe;
        }
    }
    score
}

/// Group anchors by gene family and keep the best-scoring families.
pub fn candidates(idx: &Index, anchors: Vec<Anchor>, p: &AnchorParams) -> Vec<Candidate> {
    let fam = idx.unitig_family();
    let unitig = idx.node_unitig();
    // Group by family with a sort of packed (family, input index) keys: each anchor's family
    // is looked up once (two dependent random reads). Seeds arrive sorted by (start, end) and
    // a seed's nodes ascend, so within a family the input order is already (qs, qe, node);
    // check it, and sort a family's anchors only if not.
    let mut keys: Vec<u64> = Vec::with_capacity(anchors.len());
    const AHEAD: usize = 8;
    for (i, a) in anchors.iter().enumerate() {
        if let Some(x) = anchors.get(i + AHEAD) {
            idx.prefetch_family(x.node);
        }
        let f = fam[unitig[a.node as usize] as usize];
        if f != NO_FAMILY {
            keys.push(((f as u64) << 32) | i as u64);
        }
    }
    keys.sort_unstable();
    let mut ordered: Vec<Anchor> = keys.iter().map(|&k| anchors[k as u32 as usize]).collect();
    // (family, range in `ordered`, score) per family.
    let mut groups: Vec<(u32, std::ops::Range<usize>, f32)> = Vec::new();
    let mut lo = 0;
    while lo < keys.len() {
        let f = (keys[lo] >> 32) as u32;
        let mut hi = lo + 1;
        while hi < keys.len() && (keys[hi] >> 32) as u32 == f {
            hi += 1;
        }
        let g = &mut ordered[lo..hi];
        if !g.is_sorted_by_key(|a| (a.qs, a.qe, a.node)) {
            g.sort_unstable_by_key(|a| (a.qs, a.qe, a.node));
        }
        groups.push((f, lo..hi, greedy_score(g)));
        lo = hi;
    }
    drop(anchors);
    // Rank and cut families on their scores alone; anchors are copied out only for the kept.
    groups.sort_unstable_by(|a, b| b.2.total_cmp(&a.2).then(a.0.cmp(&b.0)));
    if let Some(best) = groups.first().map(|g| g.2) {
        groups.retain(|g| g.2 >= p.alpha_fam * best);
        // Don't cut through a tie: multi-copy families give many genes the same score.
        if let Some(last) = groups.get(p.f_keep.max(1) - 1).map(|g| g.2) {
            let floor = last - p.f_tie * last.abs();
            let n = groups.iter().take_while(|g| g.2 >= floor).count();
            groups.truncate(n.min(p.f_max.max(p.f_keep)));
        }
    }
    groups
        .into_iter()
        .map(|(family, r, score)| Candidate {
            family,
            score,
            anchors: ordered[r].to_vec(),
        })
        .collect()
}

/// Merge anchors of `c` that lie on one diagonal of the same unitig: anchor `b` extends a run
/// ending at `a` if both map to the same unitig and advance equally in the read and the
/// unitig (`off(b) - off(a) == qs(b) - qs(a) > 0`). In every transcript containing the unitig
/// such anchors are collinear with zero gap discrepancy, so a merged anchor's weight is
/// exactly the score a chain gains by stepping through the run. Chains that would enter or
/// leave a run part-way through are no longer possible, and the chainer sees far fewer
/// anchors. The merged anchor keeps its first node (for the unitig offset) and spans
/// `[qs, qe)` of the whole run, which equals its transcript span.
pub fn merge_collinear(idx: &Index, c: &mut Candidate) {
    if c.anchors.len() < 2 {
        return;
    }
    let unitig = idx.node_unitig();
    let offset = idx.node_offset();
    // Unitig -> (index of the open run in `out`, last anchor's qs, qe, unitig offset).
    let mut open: rustc_hash::FxHashMap<u32, (usize, u32, u32, u32)> = rustc_hash::FxHashMap::default();
    let mut out: Vec<Anchor> = Vec::with_capacity(c.anchors.len());
    for &a in &c.anchors {
        let (u, off) = (unitig[a.node as usize], offset[a.node as usize]);
        if let Some(&(i, qs, qe, o)) = open.get(&u)
            && a.qs > qs
            && off > o
            && off - o == a.qs - qs
        {
            // The chain gain of stepping from the run's last anchor to `a` (no discrepancy).
            let qlen = (a.qe - a.qs) as f32;
            let overlap = qe.saturating_sub(a.qs).min(a.qe - a.qs) as f32;
            out[i].weight += a.weight * (1.0 - overlap / qlen);
            out[i].qe = out[i].qe.max(a.qe);
            open.insert(u, (i, a.qs, a.qe, off));
            continue;
        }
        open.insert(u, (out.len(), a.qs, a.qe, off));
        out.push(a);
    }
    if !out.is_sorted_by_key(|a| (a.qs, a.qe, a.node)) {
        out.sort_unstable_by_key(|a| (a.qs, a.qe, a.node));
    }
    c.anchors = out;
}
