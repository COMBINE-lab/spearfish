//! Colinear chaining in transcript coordinates (design doc, "Chaining").
//!
//! `OracleChainer` runs one DP per candidate transcript. It defines the scores the
//! multi-transcript refinement chainer (M6) must reproduce exactly with an unbounded beam
//! (Theorem 2), so the window, validity rule, gain and cost here are the specification:
//!
//! - Anchors of a candidate are indexed in their family order (sorted by query start);
//!   anchor `i'` may extend a chain ending at anchor `i` only if `i' - h <= i < i'`.
//! - In transcript `j`, an anchor occurs at `tpos = occ_pos + node_offset` for each
//!   occurrence of its unitig in `j` (repeats give several).
//! - A transition `i -> i'` is valid if `Δq = qs' - qe` and `Δt = ts' - te` both lie in
//!   `[-o_max, g_max]` and both starts strictly increase.
//! - Gain: `w' * (fraction of i' not overlapping i in the query)`; cost:
//!   `γ = α|Δq - ρΔt| + β log2(1 + max(Δq, Δt, 0))`.
//! - Chains are local: any anchor can start one with score `w`.

use index::Index;

use crate::anchor::Candidate;

#[derive(Clone, Debug)]
pub struct ChainParams {
    pub alpha: f32,
    pub beta: f32,
    pub rho: f32,
    pub o_max: i64,
    pub g_max: i64,
    pub h: usize,
}

impl Default for ChainParams {
    fn default() -> Self {
        ChainParams {
            alpha: 0.15,
            beta: 1.0,
            rho: 1.0,
            o_max: 100,
            g_max: 5000,
            h: 50,
        }
    }
}

/// The best chain of one transcript.
#[derive(Clone, Debug, PartialEq)]
pub struct TxHit {
    pub tx: u32,
    pub score: f32,
    pub qstart: u32,
    pub qend: u32,
    pub tstart: u32,
    pub tend: u32,
    pub anchors: u32,
    /// Query bases covered by the chain's anchors.
    pub anchored: u32,
    /// Largest |Δq - ρΔt| over the chain's transitions.
    pub max_discrepancy: f32,
}

pub trait Chainer {
    /// Best chain per transcript of the candidate, in any order.
    fn chain(&self, idx: &Index, cand: &Candidate) -> Vec<TxHit>;
}

/// An anchor placed in one transcript.
#[derive(Clone, Copy)]
struct Placed {
    /// Index of the anchor in the candidate's family order.
    a: u32,
    ts: u32,
}

pub struct OracleChainer {
    pub p: ChainParams,
}

impl OracleChainer {
    /// Transition score from anchor `x` (at target `tx_s`) to anchor `y` (at `ty_s`), if valid.
    #[inline]
    fn transition(&self, cand: &Candidate, len: &[u32], x: Placed, y: Placed) -> Option<(f32, f32)> {
        let (ax, ay) = (&cand.anchors[x.a as usize], &cand.anchors[y.a as usize]);
        if ay.qs <= ax.qs || y.ts <= x.ts {
            return None;
        }
        let te = x.ts + len[x.a as usize];
        let dq = ay.qs as i64 - ax.qe as i64;
        let dt = y.ts as i64 - te as i64;
        let p = &self.p;
        if dq < -p.o_max || dq > p.g_max || dt < -p.o_max || dt > p.g_max {
            return None;
        }
        let d = (dq as f32 - p.rho * dt as f32).abs();
        let gamma = p.alpha * d + p.beta * (1.0 + dq.max(dt).max(0) as f32).log2();
        let qlen = (ay.qe - ay.qs) as f32;
        let overlap = (ax.qe.saturating_sub(ay.qs)).min(ay.qe - ay.qs) as f32;
        let gain = ay.weight * (1.0 - overlap / qlen);
        Some((gain - gamma, d))
    }
}

impl Chainer for OracleChainer {
    fn chain(&self, idx: &Index, cand: &Candidate) -> Vec<TxHit> {
        let unitig = idx.node_unitig();
        let offset = idx.node_offset();
        let node_len = idx.node_len();
        let len: Vec<u32> = cand.anchors.iter().map(|a| node_len[a.node as usize]).collect();
        // (transcript, anchor, tpos) for every placement.
        let mut placed: Vec<(u32, Placed)> = Vec::new();
        for (i, a) in cand.anchors.iter().enumerate() {
            let (ts, ps) = idx.occurrences(unitig[a.node as usize]);
            let off = offset[a.node as usize];
            for (&t, &p) in ts.iter().zip(ps) {
                placed.push((
                    t,
                    Placed {
                        a: i as u32,
                        ts: p + off,
                    },
                ));
            }
        }
        placed.sort_unstable_by_key(|&(t, x)| (t, x.a, x.ts));

        let mut hits = Vec::new();
        let (mut f, mut back, mut disc): (Vec<f32>, Vec<u32>, Vec<f32>) =
            (Vec::new(), Vec::new(), Vec::new());
        for group in placed.chunk_by(|x, y| x.0 == y.0) {
            let t = group[0].0;
            let el: Vec<Placed> = group.iter().map(|x| x.1).collect();
            f.clear();
            back.clear();
            disc.clear();
            let mut lo = 0;
            for k in 0..el.len() {
                let y = el[k];
                while (el[lo].a as usize) + self.p.h < y.a as usize {
                    lo += 1;
                }
                let (mut best, mut arg, mut bd) = (cand.anchors[y.a as usize].weight, u32::MAX, 0.0f32);
                for m in lo..k {
                    if el[m].a >= y.a {
                        break;
                    }
                    if let Some((delta, d)) = self.transition(cand, &len, el[m], y) {
                        let s = f[m] + delta;
                        if s > best {
                            (best, arg, bd) = (s, m as u32, d);
                        }
                    }
                }
                f.push(best);
                back.push(arg);
                disc.push(bd);
            }
            // Best chain end (ties: earliest element), then trace back.
            let (mut end, mut score) = (0, f32::NEG_INFINITY);
            for (k, &s) in f.iter().enumerate() {
                if s > score {
                    (end, score) = (k, s);
                }
            }
            let (mut k, mut n, mut maxd) = (end, 0u32, 0.0f32);
            let mut spans: Vec<(u32, u32)> = Vec::new();
            loop {
                let a = &cand.anchors[el[k].a as usize];
                spans.push((a.qs, a.qe));
                n += 1;
                maxd = maxd.max(disc[k]);
                if back[k] == u32::MAX {
                    break;
                }
                k = back[k] as usize;
            }
            let first = el[k];
            let last = el[end];
            spans.reverse();
            let mut anchored = 0;
            let mut covered_to = 0;
            for &(s, e) in &spans {
                let from = s.max(covered_to);
                if e > from {
                    anchored += e - from;
                }
                covered_to = covered_to.max(e);
            }
            hits.push(TxHit {
                tx: t,
                score,
                qstart: cand.anchors[first.a as usize].qs,
                qend: cand.anchors[last.a as usize].qe,
                tstart: first.ts,
                tend: last.ts + len[last.a as usize],
                anchors: n,
                anchored,
                max_discrepancy: maxd,
            });
        }
        hits
    }
}
