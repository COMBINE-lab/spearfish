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
    for s in seeds {
        let r = idx.lookup(s.label);
        let genes = r.len();
        if genes == 0 || genes > p.h_max {
            continue;
        }
        let specificity = 1.0 - (genes as f32).ln() / p.g_max.ln();
        let weight = (s.end - s.start).min(p.s_max) as f32 * specificity.max(0.0);
        for node in r {
            if flags[node as usize] & FLAG_MASKED == 0 {
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
pub fn candidates(idx: &Index, mut anchors: Vec<Anchor>, p: &AnchorParams) -> Vec<Candidate> {
    let fam = idx.unitig_family();
    let unitig = idx.node_unitig();
    // Each anchor's family is looked up once (two dependent random reads), not on every
    // comparison of the sort.
    let mut keyed: Vec<(u32, Anchor)> = anchors
        .drain(..)
        .filter_map(|a| {
            let f = fam[unitig[a.node as usize] as usize];
            (f != NO_FAMILY).then_some((f, a))
        })
        .collect();
    keyed.sort_unstable_by_key(|(f, a)| (*f, a.qs, a.qe, a.node));
    let mut out: Vec<Candidate> = keyed
        .chunk_by(|a, b| a.0 == b.0)
        .map(|g| {
            anchors.clear();
            anchors.extend(g.iter().map(|x| x.1));
            Candidate {
                family: g[0].0,
                score: greedy_score(&anchors),
                anchors: anchors.clone(),
            }
        })
        .collect();
    out.sort_unstable_by(|a, b| b.score.total_cmp(&a.score).then(a.family.cmp(&b.family)));
    if let Some(best) = out.first().map(|c| c.score) {
        out.retain(|c| c.score >= p.alpha_fam * best);
        // Don't cut through a tie: multi-copy families give many genes the same score.
        if let Some(last) = out.get(p.f_keep.max(1) - 1).map(|c| c.score) {
            let floor = last - p.f_tie * last.abs();
            let n = out.iter().take_while(|c| c.score >= floor).count();
            out.truncate(n.min(p.f_max.max(p.f_keep)));
        }
    }
    out
}
