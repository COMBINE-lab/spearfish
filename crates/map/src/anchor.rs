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
            // locus of reads from multi-copy and repeat-rich genes.
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
    crate::stats::add(&crate::stats::ANCHORS, anchors.len() as u64);
    crate::stats::add(&crate::stats::FAMILIES, groups.len() as u64);
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
    crate::stats::add(&crate::stats::KEPT_FAMILIES, groups.len() as u64);
    crate::stats::add(
        &crate::stats::KEPT_ANCHORS,
        groups.iter().map(|g| g.1.len() as u64).sum(),
    );
    groups
        .into_iter()
        .map(|(family, r, score)| Candidate {
            family,
            score,
            anchors: ordered[r].to_vec(),
        })
        .collect()
}

/// [`select`]'s family table entry of a masked node.
pub const MASKED_NODE: u32 = u32::MAX - 1;

/// Each node's family (`NO_FAMILY` if unassigned, [`MASKED_NODE`] if masked), for [`select`]:
/// a label's nodes are contiguous, so this turns two dependent random reads per anchor
/// (unitig, then family) into a sequential scan.
pub fn node_families(idx: &Index) -> Vec<u32> {
    let (flags, unitig, fam) = (idx.node_flags(), idx.node_unitig(), idx.unitig_family());
    unitig
        .iter()
        .zip(flags)
        .map(|(&u, &f)| {
            if f & FLAG_MASKED != 0 {
                MASKED_NODE
            } else {
                fam[u as usize]
            }
        })
        .collect()
}

/// Per-family state of [`select`]: greedy score and end of the last counted anchor (`end == 0`:
/// untouched, as the first counted anchor sets it past 0). After ranking, a kept family's `end`
/// holds `KEPT - slot`.
#[derive(Clone, Copy, Default)]
struct FamState {
    end: u32,
    score: f32,
}

const KEPT: u32 = u32::MAX;

#[derive(Default)]
struct SelectScratch {
    fam: Vec<FamState>,
    touched: Vec<u32>,
    /// Per seed passing `h_max`: node range, query span and weight.
    seeds: Vec<(u32, u32, u32, u32, f32)>,
    groups: Vec<(u32, f32)>,
    /// Seed specificity by gene count, for `spec_g_max`.
    spec: Vec<f32>,
    spec_g_max: f32,
}

thread_local! {
    static SELECT: std::cell::RefCell<SelectScratch> = std::cell::RefCell::new(SelectScratch::default());
}

/// [`anchors`] followed by [`candidates`], without materializing the anchors of discarded
/// families: one pass over the seeds computes each family's greedy score in dense per-family
/// state, and anchors are built only for the kept families. `node_fam` is
/// [`node_families`]. Returns the candidates (same order, scores and anchors as
/// `candidates(idx, anchors(idx, seeds, p), p)`) and the number of anchors `anchors` returns.
pub fn select(idx: &Index, node_fam: &[u32], seeds: &[Seed], p: &AnchorParams) -> (Vec<Candidate>, u32) {
    // The streaming greedy score equals the sorted one when anchors arrive in (qs, qe, node)
    // order: seeds sorted by (start, end) with distinct spans, and a label's nodes ascend.
    if !seeds
        .windows(2)
        .all(|w| (w[0].start, w[0].end) < (w[1].start, w[1].end))
    {
        let a = anchors(idx, seeds, p);
        let n = a.len() as u32;
        return (candidates(idx, a, p), n);
    }
    SELECT.with(|sc| select_with(idx, node_fam, seeds, p, &mut sc.borrow_mut()))
}

fn select_with(
    idx: &Index,
    node_fam: &[u32],
    seeds: &[Seed],
    p: &AnchorParams,
    sc: &mut SelectScratch,
) -> (Vec<Candidate>, u32) {
    sc.touched.clear();
    sc.seeds.clear();
    let mut n_anchors = 0u32;
    let mut best = f32::NEG_INFINITY;
    // Label lookups miss in cache twice (directory entry, then labels): pipeline them as
    // `anchors` does, keeping the node ranges of seeds that pass `h_max`.
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
        // Specificity depends on the gene count alone: computed once per count (same
        // arithmetic as `anchors`).
        if sc.spec_g_max != p.g_max {
            sc.spec.clear();
            sc.spec_g_max = p.g_max;
        }
        while sc.spec.len() <= genes {
            let g = sc.spec.len();
            sc.spec.push((1.0 - (g as f32).ln() / p.g_max.ln()).max(0.0));
        }
        let weight = (s.end - s.start).min(p.s_max) as f32 * sc.spec[genes];
        sc.seeds.push((r.start, r.end, s.start, s.end, weight));
        // The family pass below reads the nodes' families: request them now, a whole
        // lookup loop ahead.
        index::prefetch(node_fam.as_ptr().wrapping_add(r.start as usize));
    }
    for i in 0..sc.seeds.len() {
        let (lo, hi, qs, qe, weight) = sc.seeds[i];
        for &f in &node_fam[lo as usize..hi as usize] {
            if f == MASKED_NODE {
                continue;
            }
            n_anchors += 1;
            if f == NO_FAMILY {
                continue;
            }
            if f as usize >= sc.fam.len() {
                sc.fam.resize(f as usize + 1, FamState::default());
            }
            let st = &mut sc.fam[f as usize];
            if st.end == 0 {
                sc.touched.push(f);
            }
            if qs >= st.end {
                st.score += weight;
                st.end = qe;
                // Comparison form (a single max instruction; scores are never NaN).
                best = if st.score > best { st.score } else { best };
            }
        }
    }
    crate::stats::add(&crate::stats::ANCHORS, n_anchors as u64);
    crate::stats::add(&crate::stats::FAMILIES, sc.touched.len() as u64);
    // Rank and cut as `candidates` does; filtering on the score floor first leaves the same
    // sorted prefix with far fewer families to sort.
    sc.groups.clear();
    for &f in &sc.touched {
        let st = std::mem::take(&mut sc.fam[f as usize]);
        if st.score >= p.alpha_fam * best {
            sc.groups.push((f, st.score));
        }
    }
    sc.groups
        .sort_unstable_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
    if let Some(last) = sc.groups.get(p.f_keep.max(1) - 1).map(|g| g.1) {
        let floor = last - p.f_tie * last.abs();
        let n = sc.groups.iter().take_while(|g| g.1 >= floor).count();
        sc.groups.truncate(n.min(p.f_max.max(p.f_keep)));
    }
    crate::stats::add(&crate::stats::KEPT_FAMILIES, sc.groups.len() as u64);
    if sc.groups.is_empty() {
        return (Vec::new(), n_anchors);
    }
    let mut cands: Vec<Candidate> = sc
        .groups
        .iter()
        .enumerate()
        .map(|(i, &(family, score))| {
            sc.fam[family as usize].end = KEPT - i as u32;
            Candidate {
                family,
                score,
                anchors: Vec::new(),
            }
        })
        .collect();
    let max_fam = sc.fam.len();
    let (unitig, offset) = (idx.node_unitig(), idx.node_offset());
    for &(lo, hi, qs, qe, weight) in &sc.seeds {
        for node in lo..hi {
            let f = node_fam[node as usize] as usize;
            if f < max_fam {
                let e = sc.fam[f].end;
                if e != 0 {
                    // Merging and chaining read the node's unitig and offset next.
                    index::prefetch(unitig.as_ptr().wrapping_add(node as usize));
                    index::prefetch(offset.as_ptr().wrapping_add(node as usize));
                    cands[(KEPT - e) as usize]
                        .anchors
                        .push(Anchor { qs, qe, node, weight });
                }
            }
        }
    }
    for g in &sc.groups {
        sc.fam[g.0 as usize] = FamState::default();
    }
    crate::stats::add(
        &crate::stats::KEPT_ANCHORS,
        cands.iter().map(|c| c.anchors.len() as u64).sum(),
    );
    (cands, n_anchors)
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
    MERGE_OPEN.with(|open| merge_collinear_with(idx, c, &mut open.borrow_mut()));
}

thread_local! {
    /// Unitig -> (index of the open run, last anchor's qs, qe, unitig offset), reused.
    static MERGE_OPEN: std::cell::RefCell<rustc_hash::FxHashMap<u32, (u32, u32, u32, u32)>> =
        std::cell::RefCell::new(rustc_hash::FxHashMap::default());
}

fn merge_collinear_with(
    idx: &Index,
    c: &mut Candidate,
    open: &mut rustc_hash::FxHashMap<u32, (u32, u32, u32, u32)>,
) {
    let unitig = idx.node_unitig();
    let offset = idx.node_offset();
    open.clear();
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
            out[i as usize].weight += a.weight * (1.0 - overlap / qlen);
            out[i as usize].qe = out[i as usize].qe.max(a.qe);
            open.insert(u, (i, a.qs, a.qe, off));
            continue;
        }
        open.insert(u, (out.len() as u32, a.qs, a.qe, off));
        out.push(a);
    }
    if !out.is_sorted_by_key(|a| (a.qs, a.qe, a.node)) {
        out.sort_unstable_by_key(|a| (a.qs, a.qe, a.node));
    }
    crate::stats::add(&crate::stats::MERGED_ANCHORS, out.len() as u64);
    c.anchors = out;
}
