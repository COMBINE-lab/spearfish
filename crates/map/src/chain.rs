//! Colinear chaining in transcript coordinates (design doc, "Chaining").
//!
//! `OracleChainer` runs one DP per candidate transcript. Transcripts are visited in order of
//! their placement lists (shifted to the first placement), and each DP resumes where its list
//! diverges from the previous transcript's: identical prefixes give identical DP values, so
//! this is exact and shares the work of isoforms that agree over a read's anchors. It defines
//! the scores the multi-transcript refinement chainer (M6, `RefineChainer`) must reproduce
//! exactly with an unbounded beam (Theorem 2), so the window, validity rule, gain and cost here
//! are the specification:
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
            // Tuned on simulated chr21 ONT reads.
            alpha: 1.5,
            beta: 0.0,
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
    /// The chain's anchors in query order: (query start, transcript start, length).
    pub chain: Vec<(u32, u32, u32)>,
    /// The read's reverse complement was chained (set by the mapper).
    pub reverse: bool,
    /// The span was extended past the outermost anchors to the read's ends (set by the mapper).
    pub extended: bool,
}

/// A transcript's hit from its chain, given as (qs, qe, ts, len) per anchor in query order.
pub(crate) fn hit_from_path(
    tx: u32,
    score: f32,
    max_discrepancy: f32,
    path: impl Iterator<Item = (u32, u32, u32, u32)>,
) -> TxHit {
    let (mut chain, mut anchored, mut covered_to) = (Vec::new(), 0, 0);
    let mut qend = 0;
    for (qs, qe, ts, len) in path {
        chain.push((qs, ts, len));
        let from = qs.max(covered_to);
        if qe > from {
            anchored += qe - from;
        }
        covered_to = covered_to.max(qe);
        qend = qe;
    }
    let (first, last) = (chain[0], chain[chain.len() - 1]);
    TxHit {
        tx,
        score,
        qstart: first.0,
        qend,
        tstart: first.1,
        tend: last.1 + last.2,
        anchors: chain.len() as u32,
        anchored,
        max_discrepancy,
        chain,
        reverse: false,
        extended: false,
    }
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

/// One placed anchor, as the transition rule sees it: query span `[qs, qe)`, transcript
/// start `ts`, seed length `len` in the transcript, and weight `w`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Pt {
    pub qs: u32,
    pub qe: u32,
    pub ts: u32,
    pub len: u32,
    pub w: f32,
}

/// The transition rule shared by every chainer: score delta of extending a chain ending at
/// `x` with `y` and its gap discrepancy, or `None` if the transition is invalid.
#[inline]
pub(crate) fn transition(p: &ChainParams, x: Pt, y: Pt) -> Option<(f32, f32)> {
    if y.qs <= x.qs || y.ts <= x.ts {
        return None;
    }
    let dq = y.qs as i64 - x.qe as i64;
    let dt = y.ts as i64 - (x.ts + x.len) as i64;
    if dq < -p.o_max || dq > p.g_max || dt < -p.o_max || dt > p.g_max {
        return None;
    }
    let d = (dq as f32 - p.rho * dt as f32).abs();
    let gamma = if p.beta == 0.0 {
        p.alpha * d
    } else {
        p.alpha * d + p.beta * (1.0 + dq.max(dt).max(0) as f32).log2()
    };
    // Without overlap the gain is y.w exactly (0 / qlen = 0), so skip the division.
    let overlap = x.qe.saturating_sub(y.qs).min(y.qe - y.qs);
    if overlap == 0 {
        return Some((y.w - gamma, d));
    }
    let qlen = (y.qe - y.qs) as f32;
    Some((y.w * (1.0 - overlap as f32 / qlen) - gamma, d))
}

/// Best local colinear chain over points sorted by (qs, ts), with predecessors among the
/// previous `2h` points. Returns (score, chain indices in order, max discrepancy), or `None`
/// for no points.
pub(crate) fn best_chain(p: &ChainParams, pts: &[Pt]) -> Option<(f32, Vec<usize>, f32)> {
    if pts.is_empty() {
        return None;
    }
    let (mut f, mut back, mut disc) = (
        vec![0f32; pts.len()],
        vec![usize::MAX; pts.len()],
        vec![0f32; pts.len()],
    );
    let mut lo = 0;
    for k in 0..pts.len() {
        let y = pts[k];
        while k - lo > 2 * p.h {
            lo += 1;
        }
        let (mut best, mut arg, mut bd) = (y.w, usize::MAX, 0f32);
        for m in lo..k {
            if let Some((delta, d)) = transition(p, pts[m], y) {
                let s = f[m] + delta;
                if s > best {
                    (best, arg, bd) = (s, m, d);
                }
            }
        }
        f[k] = best;
        back[k] = arg;
        disc[k] = bd;
    }
    let mut end = 0;
    for k in 1..pts.len() {
        if f[k] > f[end] {
            end = k;
        }
    }
    let (mut k, mut chain, mut maxd) = (end, Vec::new(), 0f32);
    loop {
        chain.push(k);
        maxd = maxd.max(disc[k]);
        if back[k] == usize::MAX {
            break;
        }
        k = back[k];
    }
    chain.reverse();
    Some((f[end], chain, maxd))
}

/// Bias added to relative target positions so they stay positive.
const REL_BIAS: i64 = 1 << 29;

/// Per-thread buffers for the oracle DP (reused across reads).
#[derive(Default)]
struct DpScratch {
    /// Transcript -> bucket, and per-bucket placements in anchor order.
    bucket_of: rustc_hash::FxHashMap<u32, u32>,
    tx: Vec<u32>,
    buckets: Vec<Vec<Placed>>,
    // Structure-of-arrays for one transcript's placements.
    a: Vec<u32>,
    qs: Vec<i32>,
    qe: Vec<i32>,
    ts: Vec<i32>,
    te: Vec<i32>,
    w: Vec<f32>,
    f: Vec<f32>,
    back: Vec<u32>,
    disc: Vec<f32>,
    /// Transcripts (bucket indices) in placement-list order.
    order: Vec<u32>,
    // Candidate scores and discrepancies of one anchor's predecessors.
    cs: Vec<f32>,
    cd: Vec<f32>,
}

thread_local! {
    static DP: std::cell::RefCell<DpScratch> = std::cell::RefCell::new(DpScratch::default());
}

/// Scores of extending each predecessor `m` in `lo..k` with `y` (`-inf` if invalid) and the
/// gap discrepancies, branch-free so it vectorizes; the arithmetic is exactly `transition`'s
/// (with β = 0), so results are bit-identical.
#[multiversion::multiversion(targets(
    "x86_64+avx512f+avx512bw+avx512vl+avx2+fma",
    "x86_64+avx2+fma",
    "x86_64+avx",
    "aarch64+neon"
))]
#[allow(clippy::too_many_arguments)]
fn pred_scores(
    alpha: f32,
    rho: f32,
    o_max: i32,
    g_max: i32,
    xa: &[u32],
    xqs: &[i32],
    xqe: &[i32],
    xts: &[i32],
    xte: &[i32],
    xf: &[f32],
    (ya, yqs, yqe, yts, yw): (u32, i32, i32, i32, f32),
    cs: &mut [f32],
    cd: &mut [f32],
) {
    let qlen = (yqe - yqs) as f32;
    for m in 0..cs.len() {
        let dq = yqs - xqe[m];
        let dt = yts - xte[m];
        let valid = xa[m] < ya
            && yqs > xqs[m]
            && yts > xts[m]
            && dq >= -o_max
            && dq <= g_max
            && dt >= -o_max
            && dt <= g_max;
        let d = (dq as f32 - rho * dt as f32).abs();
        let gamma = alpha * d;
        let overlap = (xqe[m] - yqs).max(0).min(yqe - yqs) as f32;
        let delta = yw * (1.0 - overlap / qlen) - gamma;
        cs[m] = if valid { xf[m] + delta } else { f32::NEG_INFINITY };
        cd[m] = d;
    }
}

impl Chainer for OracleChainer {
    fn chain(&self, idx: &Index, cand: &Candidate) -> Vec<TxHit> {
        DP.with(|dp| self.chain_with(idx, cand, &mut dp.borrow_mut()))
    }
}

impl OracleChainer {
    fn chain_with(&self, idx: &Index, cand: &Candidate, sc: &mut DpScratch) -> Vec<TxHit> {
        let unitig = idx.node_unitig();
        let offset = idx.node_offset();
        // An anchor's transcript span equals its query span (exact core matches; merged
        // collinear runs).
        let len = |i: usize| cand.anchors[i].qe - cand.anchors[i].qs;
        // Bucket placements by transcript; visiting anchors in order keeps each bucket sorted
        // by (anchor, target position).
        sc.bucket_of.clear();
        sc.tx.clear();
        let mut used = 0;
        for (i, a) in cand.anchors.iter().enumerate() {
            let (ts, ps) = idx.occurrences(unitig[a.node as usize]);
            let off = offset[a.node as usize];
            for (&t, &p) in ts.iter().zip(ps) {
                let b = *sc.bucket_of.entry(t).or_insert_with(|| {
                    if used == sc.buckets.len() {
                        sc.buckets.push(Vec::new());
                    }
                    sc.buckets[used].clear();
                    sc.tx.push(t);
                    used += 1;
                    used as u32 - 1
                });
                sc.buckets[b as usize].push(Placed {
                    a: i as u32,
                    ts: p + off,
                });
            }
        }

        let p = &self.p;
        // Visit transcripts sorted by their placement lists shifted to the first placement:
        // transcripts with a common prefix (relative to their first placement) have identical
        // DP values over it, so each DP resumes where the previous transcript's list diverges.
        sc.order.clear();
        sc.order.extend(0..used as u32);
        {
            let b = &sc.buckets;
            sc.order.sort_unstable_by(|&x, &y| {
                let (u, v) = (&b[x as usize], &b[y as usize]);
                let (u0, v0) = (u[0].ts as i64, v[0].ts as i64);
                for (pu, pv) in u.iter().zip(v) {
                    let c =
                        pu.a.cmp(&pv.a)
                            .then((pu.ts as i64 - u0).cmp(&(pv.ts as i64 - v0)));
                    if c.is_ne() {
                        return c;
                    }
                }
                u.len().cmp(&v.len())
            });
        }
        let mut hits = Vec::with_capacity(used);
        let mut prev: Option<usize> = None;
        for oi in 0..used {
            let bi = sc.order[oi] as usize;
            let t = sc.tx[bi];
            let el = std::mem::take(&mut sc.buckets[bi]);
            let n = el.len();
            // Target positions are stored relative to the first placement (plus a bias that
            // keeps them positive); transitions only use their differences.
            let base = el[0].ts as i64 - REL_BIAS;
            let lcp = prev.map_or(0, |pb| {
                let pe = &sc.buckets[pb];
                let pbase = pe[0].ts as i64 - REL_BIAS;
                el.iter()
                    .zip(pe)
                    .take_while(|(x, y)| x.a == y.a && x.ts as i64 - base == y.ts as i64 - pbase)
                    .count()
            });
            sc.a.truncate(lcp);
            sc.qs.truncate(lcp);
            sc.qe.truncate(lcp);
            sc.ts.truncate(lcp);
            sc.te.truncate(lcp);
            sc.w.truncate(lcp);
            for x in &el[lcp..] {
                let an = &cand.anchors[x.a as usize];
                let rel = (x.ts as i64 - base) as i32;
                sc.a.push(x.a);
                sc.qs.push(an.qs as i32);
                sc.qe.push(an.qe as i32);
                sc.ts.push(rel);
                sc.te.push(rel + len(x.a as usize) as i32);
                sc.w.push(an.weight);
            }
            sc.f.truncate(lcp);
            sc.f.resize(n, 0.0);
            sc.back.truncate(lcp);
            sc.back.resize(n, u32::MAX);
            sc.disc.truncate(lcp);
            sc.disc.resize(n, 0.0);
            let mut lo = 0;
            for k in lcp..n {
                let ya = sc.a[k];
                while (sc.a[lo] as usize) + p.h < ya as usize {
                    lo += 1;
                }
                let y = (ya, sc.qs[k], sc.qe[k], sc.ts[k], sc.w[k]);
                let (mut best, mut arg, mut bd) = (y.4, u32::MAX, 0.0f32);
                if k > lo {
                    if p.beta == 0.0 {
                        sc.cs.resize(k - lo, 0.0);
                        sc.cd.resize(k - lo, 0.0);
                        pred_scores(
                            p.alpha,
                            p.rho,
                            p.o_max as i32,
                            p.g_max as i32,
                            &sc.a[lo..k],
                            &sc.qs[lo..k],
                            &sc.qe[lo..k],
                            &sc.ts[lo..k],
                            &sc.te[lo..k],
                            &sc.f[lo..k],
                            y,
                            &mut sc.cs[..k - lo],
                            &mut sc.cd[..k - lo],
                        );
                        let mx = sc.cs[..k - lo].iter().copied().fold(f32::NEG_INFINITY, f32::max);
                        if mx > best {
                            let j = sc.cs[..k - lo].iter().position(|&v| v == mx).unwrap();
                            (best, arg, bd) = (mx, (lo + j) as u32, sc.cd[j]);
                        }
                    } else {
                        let pt = |m: usize| Pt {
                            qs: sc.qs[m] as u32,
                            qe: sc.qe[m] as u32,
                            ts: sc.ts[m] as u32,
                            len: (sc.te[m] - sc.ts[m]) as u32,
                            w: sc.w[m],
                        };
                        for m in lo..k {
                            if sc.a[m] >= ya {
                                break;
                            }
                            if let Some((delta, d)) = transition(p, pt(m), pt(k)) {
                                let s = sc.f[m] + delta;
                                if s > best {
                                    (best, arg, bd) = (s, m as u32, d);
                                }
                            }
                        }
                    }
                }
                sc.f[k] = best;
                sc.back[k] = arg;
                sc.disc[k] = bd;
            }
            // Best chain end (ties: earliest element), then trace back.
            let (mut end, mut score) = (0, f32::NEG_INFINITY);
            for (k, &s) in sc.f.iter().enumerate() {
                if s > score {
                    (end, score) = (k, s);
                }
            }
            let (mut k, mut maxd) = (end, 0.0f32);
            let mut path: Vec<usize> = Vec::new();
            loop {
                path.push(k);
                maxd = maxd.max(sc.disc[k]);
                if sc.back[k] == u32::MAX {
                    break;
                }
                k = sc.back[k] as usize;
            }
            path.reverse();
            hits.push(hit_from_path(
                t,
                score,
                maxd,
                path.iter().map(|&k| {
                    (
                        sc.qs[k] as u32,
                        sc.qe[k] as u32,
                        (sc.ts[k] as i64 + base) as u32,
                        (sc.te[k] - sc.ts[k]) as u32,
                    )
                }),
            ));
            prev = Some(bi);
            sc.buckets[bi] = el;
        }
        hits
    }
}

/// A chainer chosen at run time.
pub enum AnyChainer {
    Oracle(OracleChainer),
    Refine(crate::RefineChainer),
}

impl Chainer for AnyChainer {
    fn chain(&self, idx: &Index, cand: &Candidate) -> Vec<TxHit> {
        match self {
            AnyChainer::Oracle(c) => c.chain(idx, cand),
            AnyChainer::Refine(c) => c.chain(idx, cand),
        }
    }
}
