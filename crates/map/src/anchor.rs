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
    /// Keep at most this many genes.
    pub f_keep: usize,
}

impl Default for AnchorParams {
    fn default() -> Self {
        AnchorParams {
            h_max: 64,
            s_max: 100,
            g_max: 64.0,
            alpha_fam: 0.5,
            f_keep: 8,
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
    let family = |a: &Anchor| fam[unitig[a.node as usize] as usize];
    anchors.retain(|a| family(a) != NO_FAMILY);
    anchors.sort_unstable_by_key(|a| (family(a), a.qs, a.qe, a.node));
    let mut out: Vec<Candidate> = anchors
        .chunk_by(|a, b| family(a) == family(b))
        .map(|g| Candidate {
            family: family(&g[0]),
            score: greedy_score(g),
            anchors: g.to_vec(),
        })
        .collect();
    out.sort_unstable_by(|a, b| b.score.total_cmp(&a.score).then(a.family.cmp(&b.family)));
    if let Some(best) = out.first().map(|c| c.score) {
        out.retain(|c| c.score >= p.alpha_fam * best);
        out.truncate(p.f_keep);
    }
    out
}
