//! Colinear chaining in transcript coordinates.
//!
//! `OracleChainer` runs one DP per candidate transcript. Transcripts are visited in order of
//! their placement lists (shifted to the first placement), and each DP resumes where its list
//! diverges from the previous transcript's: identical prefixes give identical DP values, so
//! this is exact and shares the work of isoforms that agree over a read's anchors. It defines
//! the scores the multi-transcript refinement chainer (`RefineChainer`) must reproduce exactly
//! with an unbounded beam, so the window, validity rule, gain and cost here are the
//! specification:
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
use lcp_core::simd::Level;

use crate::anchor::Candidate;

#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
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
    let (mut chain, mut anchored, mut covered_to) = (Vec::with_capacity(path.size_hint().0), 0, 0);
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
#[cfg(test)]
pub(crate) fn best_chain(p: &ChainParams, pts: &[Pt]) -> Option<(f32, Vec<usize>, f32)> {
    if pts.is_empty() {
        return None;
    }
    if p.beta == 0.0 {
        return best_chains(p, &[pts]).pop().unwrap();
    }
    best_chain_scalar(p, pts)
}

/// [`best_chain`] one transition at a time (any β).
fn best_chain_scalar(p: &ChainParams, pts: &[Pt]) -> Option<(f32, Vec<usize>, f32)> {
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

thread_local! {
    static BEST: std::cell::RefCell<DpScratch> = std::cell::RefCell::new(DpScratch::default());
}

/// A chain result: (score, chain indices in order, max discrepancy).
pub(crate) type ChainResult = (f32, Vec<usize>, f32);

/// [`best_chain`] for each of `lists`, sharing work: with β = 0 the lists are visited in
/// order of their points shifted to the first point's transcript start, and each DP resumes
/// where the list diverges from the previous one (transitions depend on target positions
/// only through differences, so equal shifted prefixes have equal DP values). The same
/// results as separate calls.
pub(crate) fn best_chains(p: &ChainParams, lists: &[&[Pt]]) -> Vec<Option<ChainResult>> {
    if p.beta != 0.0 {
        return lists.iter().map(|l| best_chain_scalar(p, l)).collect();
    }
    BEST.with(|sc| best_chains_beta0(p, lists, &mut sc.borrow_mut()))
}

/// A point as the DP sees it: query span, target start relative to the list's first point,
/// length and weight (bits).
fn shifted(l: &[Pt], i: usize) -> (u32, u32, i64, u32, u32) {
    let x = l[i];
    (x.qs, x.qe, x.ts as i64 - l[0].ts as i64, x.len, x.w.to_bits())
}

/// [`best_chains`] with β = 0 through the oracle's vectorized DP: point `k` gets identity `k`
/// and the window `2h`, so `dp_beta0` sees `best_chain`'s predecessors, and its transition
/// arithmetic and tie-break (first maximum above the start score) are `transition`'s.
fn best_chains_beta0(p: &ChainParams, lists: &[&[Pt]], sc: &mut DpScratch) -> Vec<Option<ChainResult>> {
    let mut out: Vec<Option<ChainResult>> = vec![None; lists.len()];
    sc.order.clear();
    sc.order.extend(0..lists.len() as u32);
    if lists.len() > 1 {
        sc.order.sort_unstable_by(|&x, &y| {
            let (u, v) = (lists[x as usize], lists[y as usize]);
            (0..u.len().min(v.len()))
                .map(|i| shifted(u, i).cmp(&shifted(v, i)))
                .find(|o| o.is_ne())
                .unwrap_or(u.len().cmp(&v.len()))
        });
    }
    let wp = ChainParams {
        h: 2 * p.h,
        ..p.clone()
    };
    let level = crate::align::simd_level();
    // The list whose points and DP values the buffers hold.
    let mut prev: Option<&[Pt]> = None;
    let order = std::mem::take(&mut sc.order);
    for &x in &order {
        let pts = lists[x as usize];
        let n = pts.len();
        if n == 0 {
            continue;
        }
        let lcp = prev.map_or(0, |q| {
            (0..q.len().min(n))
                .take_while(|&i| shifted(q, i) == shifted(pts, i))
                .count()
        });
        crate::stats::add(&crate::stats::RESCUE_PTS, n as u64);
        crate::stats::add(&crate::stats::RESCUE_PTS_DP, (n - lcp) as u64);
        for v in [&mut sc.qs, &mut sc.qe, &mut sc.ts, &mut sc.te] {
            if v.len() < n {
                v.resize(n, 0);
            }
        }
        for v in [&mut sc.w, &mut sc.f, &mut sc.disc, &mut sc.cs, &mut sc.cd] {
            if v.len() < n {
                v.resize(n, 0.0);
            }
        }
        if sc.a.len() < n {
            sc.a.resize(n, 0);
            sc.back.resize(n, u32::MAX);
        }
        let ts0 = pts[0].ts as i64;
        for (i, x) in pts.iter().enumerate().skip(lcp) {
            let ts = (x.ts as i64 - ts0) as i32;
            sc.a[i] = i as u32;
            sc.qs[i] = x.qs as i32;
            sc.qe[i] = x.qe as i32;
            sc.ts[i] = ts;
            sc.te[i] = ts + x.len as i32;
            sc.w[i] = x.w;
        }
        dp_beta0(
            level,
            &wp,
            (
                &sc.a[..n],
                &sc.qs[..n],
                &sc.qe[..n],
                &sc.ts[..n],
                &sc.te[..n],
                &sc.w[..n],
            ),
            &mut sc.f[..n],
            &mut sc.back[..n],
            &mut sc.disc[..n],
            lcp,
            &mut sc.cs,
            &mut sc.cd,
        );
        prev = Some(pts);
        let (f, back, disc) = (&sc.f[..n], &sc.back[..n], &sc.disc[..n]);
        let mut end = 0;
        for k in 1..n {
            if f[k] > f[end] {
                end = k;
            }
        }
        let (mut k, mut chain, mut maxd) = (end, Vec::new(), 0f32);
        loop {
            chain.push(k);
            maxd = maxd.max(disc[k]);
            if back[k] == u32::MAX {
                break;
            }
            k = back[k] as usize;
        }
        chain.reverse();
        out[x as usize] = Some((f[end], chain, maxd));
    }
    sc.order = order;
    out
}

/// Bias added to relative target positions so they stay positive.
const REL_BIAS: i64 = 1 << 29;

/// Per-thread buffers for the oracle DP (reused across reads).
#[derive(Default)]
struct DpScratch {
    /// Transcript -> bucket + 1 (all 0 between calls), and per-bucket placements in anchor
    /// order.
    bucket_of: Vec<u32>,
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
    /// Traceback of the current transcript.
    path: Vec<usize>,
}

thread_local! {
    static DP: std::cell::RefCell<DpScratch> = std::cell::RefCell::new(DpScratch::default());
}

/// Scores of extending each predecessor `m` in `lo..k` with `y` (`-inf` if invalid) and the
/// gap discrepancies, branch-free so it vectorizes; the arithmetic is exactly `transition`'s
/// (with β = 0), so results are bit-identical.
#[inline(always)]
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

/// [`pred_scores`] in 4-lane NEON blocks, the rest scalar: the same operations per element
/// (conversions, products, quotient and sums in the same order), so the same bits.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
#[allow(clippy::too_many_arguments)]
unsafe fn pred_scores_neon(
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
    y: (u32, i32, i32, i32, f32),
    cs: &mut [f32],
    cd: &mut [f32],
) {
    use std::arch::aarch64::*;
    let n = cs.len();
    let lens = [
        xa.len(),
        xqs.len(),
        xqe.len(),
        xts.len(),
        xte.len(),
        xf.len(),
        cd.len(),
    ];
    assert!(lens.iter().all(|&l| l >= n));
    let (ya, yqs, yqe, yts, yw) = y;
    let (vya, vyqs, vyts) = (vdupq_n_u32(ya), vdupq_n_s32(yqs), vdupq_n_s32(yts));
    let (lo, hi) = (vdupq_n_s32(-o_max), vdupq_n_s32(g_max));
    let (vlen, vqlen) = (vdupq_n_s32(yqe - yqs), vdupq_n_f32((yqe - yqs) as f32));
    let (valpha, vrho, vyw) = (vdupq_n_f32(alpha), vdupq_n_f32(rho), vdupq_n_f32(yw));
    let (one, ninf) = (vdupq_n_f32(1.0), vdupq_n_f32(f32::NEG_INFINITY));
    let mut m = 0;
    while m + 4 <= n {
        // SAFETY: m + 4 <= n, and every slice holds at least n elements (checked above).
        let (a, qs, qe, ts, te, f) = unsafe {
            (
                vld1q_u32(xa.as_ptr().add(m)),
                vld1q_s32(xqs.as_ptr().add(m)),
                vld1q_s32(xqe.as_ptr().add(m)),
                vld1q_s32(xts.as_ptr().add(m)),
                vld1q_s32(xte.as_ptr().add(m)),
                vld1q_f32(xf.as_ptr().add(m)),
            )
        };
        let dq = vsubq_s32(vyqs, qe);
        let dt = vsubq_s32(vyts, te);
        let valid = vandq_u32(
            vandq_u32(
                vandq_u32(vcltq_u32(a, vya), vcgtq_s32(vyqs, qs)),
                vcgtq_s32(vyts, ts),
            ),
            vandq_u32(
                vandq_u32(vcgeq_s32(dq, lo), vcleq_s32(dq, hi)),
                vandq_u32(vcgeq_s32(dt, lo), vcleq_s32(dt, hi)),
            ),
        );
        let d = vabsq_f32(vsubq_f32(vcvtq_f32_s32(dq), vmulq_f32(vrho, vcvtq_f32_s32(dt))));
        let gamma = vmulq_f32(valpha, d);
        let overlap = vcvtq_f32_s32(vminq_s32(vmaxq_s32(vsubq_s32(qe, vyqs), vdupq_n_s32(0)), vlen));
        let delta = vsubq_f32(vmulq_f32(vyw, vsubq_f32(one, vdivq_f32(overlap, vqlen))), gamma);
        let c = vbslq_f32(valid, vaddq_f32(f, delta), ninf);
        // SAFETY: as for the loads.
        unsafe {
            vst1q_f32(cs.as_mut_ptr().add(m), c);
            vst1q_f32(cd.as_mut_ptr().add(m), d);
        }
        m += 4;
    }
    pred_scores(
        alpha,
        rho,
        o_max,
        g_max,
        &xa[m..n],
        &xqs[m..n],
        &xqe[m..n],
        &xts[m..n],
        &xte[m..n],
        &xf[m..n],
        y,
        &mut cs[m..],
        &mut cd[m..n],
    );
}

/// [`max_score`] of `xs` and the first position equal to it, in one NEON pass: each lane keeps
/// its maximum and where it first reached it (strictly greater values replace it), and the
/// lowest position among the lanes at the overall maximum is the first one. The maximum is
/// [`max_score`]'s up to the sign of a zero, which `==` ignores.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn max_first_neon(xs: &[f32]) -> (f32, usize) {
    use std::arch::aarch64::*;
    let (blocks, rest) = xs.as_chunks::<4>();
    let lane: [u32; 4] = [0, 1, 2, 3];
    // SAFETY: four u32s.
    let mut at = unsafe { vld1q_u32(lane.as_ptr()) };
    let (mut mxv, mut idx) = (vdupq_n_f32(f32::NEG_INFINITY), vdupq_n_u32(0));
    for c in blocks {
        // SAFETY: four f32s.
        let v = unsafe { vld1q_f32(c.as_ptr()) };
        let gt = vcgtq_f32(v, mxv);
        mxv = vbslq_f32(gt, v, mxv);
        idx = vbslq_u32(gt, at, idx);
        at = vaddq_u32(at, vdupq_n_u32(4));
    }
    let mut mx = vmaxvq_f32(mxv);
    let mut j = vminvq_u32(vbslq_u32(
        vceqq_f32(mxv, vdupq_n_f32(mx)),
        idx,
        vdupq_n_u32(u32::MAX),
    )) as usize;
    for (i, &v) in rest.iter().enumerate() {
        if v > mx {
            (mx, j) = (v, 4 * blocks.len() + i);
        }
    }
    (mx, j)
}

/// Maximum of `xs` (never NaN; equal to the in-order maximum up to the sign of a zero).
/// Long windows reduce in 16 independent lanes; a single running maximum is one dependent
/// step per element.
#[inline(always)]
fn max_score(xs: &[f32]) -> f32 {
    // Comparison form: one max instruction per step (scores are never NaN).
    let max = |a: f32, v: f32| if v > a { v } else { a };
    if xs.len() < 32 {
        return xs.iter().fold(f32::NEG_INFINITY, |a, &v| max(a, v));
    }
    let (blocks, rest) = xs.as_chunks::<16>();
    let mut acc = [f32::NEG_INFINITY; 16];
    for c in blocks {
        for l in 0..16 {
            acc[l] = max(acc[l], c[l]);
        }
    }
    rest.iter().chain(&acc).fold(f32::NEG_INFINITY, |a, &v| max(a, v))
}

/// The per-anchor kernels of [`dp_beta0_core`]: the scalar loops (autovectorized for the
/// entry point's level), or NEON.
trait Preds {
    /// [`pred_scores`].
    #[allow(clippy::too_many_arguments)]
    fn scores(
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
        y: (u32, i32, i32, i32, f32),
        cs: &mut [f32],
        cd: &mut [f32],
    );
    /// The first position of the maximum of `xs` (never NaN), if that maximum exceeds `best`.
    fn first_max_above(xs: &[f32], best: f32) -> Option<usize>;
}

/// [`pred_scores`] and [`max_score`].
struct ScalarPreds;

impl Preds for ScalarPreds {
    #[inline(always)]
    fn scores(
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
        y: (u32, i32, i32, i32, f32),
        cs: &mut [f32],
        cd: &mut [f32],
    ) {
        pred_scores(alpha, rho, o_max, g_max, xa, xqs, xqe, xts, xte, xf, y, cs, cd)
    }
    #[inline(always)]
    fn first_max_above(xs: &[f32], best: f32) -> Option<usize> {
        let mx = max_score(xs);
        // The first predecessor at the maximum (its own value is what the in-order maximum
        // would return, whatever the sign of a zero).
        (mx > best).then(|| xs.iter().position(|&v| v == mx).unwrap())
    }
}

/// [`pred_scores_neon`] and [`max_first_neon`] (NEON is part of the aarch64 baseline; LLVM's
/// cost model leaves the scalar loop unvectorized there).
#[cfg(target_arch = "aarch64")]
struct NeonPreds;

#[cfg(target_arch = "aarch64")]
impl Preds for NeonPreds {
    #[inline(always)]
    fn scores(
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
        y: (u32, i32, i32, i32, f32),
        cs: &mut [f32],
        cd: &mut [f32],
    ) {
        // SAFETY: NEON is part of the aarch64 baseline.
        unsafe { pred_scores_neon(alpha, rho, o_max, g_max, xa, xqs, xqe, xts, xte, xf, y, cs, cd) }
    }
    #[inline(always)]
    fn first_max_above(xs: &[f32], best: f32) -> Option<usize> {
        // As the scalar one, the position found in the same pass.
        // SAFETY: NEON is part of the aarch64 baseline.
        let (mx, j) = unsafe { max_first_neon(xs) };
        (mx > best).then_some(j)
    }
}

/// A transcript's placements as structure-of-arrays: identity, query start and end, relative
/// target start and end, weight.
type Soa<'a> = (&'a [u32], &'a [i32], &'a [i32], &'a [i32], &'a [i32], &'a [f32]);

/// One transcript's forward DP with β = 0 over `lcp..n` (`f`, `back`, `disc` hold the shared
/// prefix's values below `lcp`), at `level` (which the CPU must support): compiled per level as
/// a whole so the per-anchor work (predecessor scores, their maximum) stays in one vectorized
/// loop without dispatch. The arithmetic is the same at every level (no FMA contraction), so
/// results are bit-identical.
#[allow(clippy::too_many_arguments)]
fn dp_beta0(
    level: Level,
    p: &ChainParams,
    soa: Soa,
    f: &mut [f32],
    back: &mut [u32],
    disc: &mut [f32],
    lcp: usize,
    cs: &mut [f32],
    cd: &mut [f32],
) {
    debug_assert!(lcp_core::simd::detected().includes(level));
    match level {
        // SAFETY: the CPU has the level's features (callers pass the active level).
        #[cfg(target_arch = "x86_64")]
        Level::Avx512 => unsafe { dp_beta0_avx512(p, soa, f, back, disc, lcp, cs, cd) },
        #[cfg(target_arch = "x86_64")]
        Level::Avx2 => unsafe { dp_beta0_avx2(p, soa, f, back, disc, lcp, cs, cd) },
        // NEON is part of the aarch64 baseline: no entry point of its own.
        #[cfg(target_arch = "aarch64")]
        Level::Neon => dp_beta0_core::<NeonPreds>(p, soa, f, back, disc, lcp, cs, cd),
        _ => dp_beta0_core::<ScalarPreds>(p, soa, f, back, disc, lcp, cs, cd),
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vl,avx2,fma,bmi1,bmi2,popcnt,lzcnt")]
#[allow(clippy::too_many_arguments)]
fn dp_beta0_avx512(
    p: &ChainParams,
    soa: Soa,
    f: &mut [f32],
    back: &mut [u32],
    disc: &mut [f32],
    lcp: usize,
    cs: &mut [f32],
    cd: &mut [f32],
) {
    dp_beta0_core::<ScalarPreds>(p, soa, f, back, disc, lcp, cs, cd)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma,bmi1,bmi2,popcnt,lzcnt")]
#[allow(clippy::too_many_arguments)]
fn dp_beta0_avx2(
    p: &ChainParams,
    soa: Soa,
    f: &mut [f32],
    back: &mut [u32],
    disc: &mut [f32],
    lcp: usize,
    cs: &mut [f32],
    cd: &mut [f32],
) {
    dp_beta0_core::<ScalarPreds>(p, soa, f, back, disc, lcp, cs, cd)
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn dp_beta0_core<K: Preds>(
    p: &ChainParams,
    (a, qs, qe, ts, te, w): Soa,
    f: &mut [f32],
    back: &mut [u32],
    disc: &mut [f32],
    lcp: usize,
    cs: &mut [f32],
    cd: &mut [f32],
) {
    let n = a.len();
    let mut lo = 0;
    for k in lcp..n {
        let ya = a[k];
        while (a[lo] as usize) + p.h < ya as usize {
            lo += 1;
        }
        let y = (ya, qs[k], qe[k], ts[k], w[k]);
        let (mut best, mut arg, mut bd) = (y.4, u32::MAX, 0.0f32);
        if k > lo {
            let m = k - lo;
            K::scores(
                p.alpha,
                p.rho,
                p.o_max as i32,
                p.g_max as i32,
                &a[lo..k],
                &qs[lo..k],
                &qe[lo..k],
                &ts[lo..k],
                &te[lo..k],
                &f[lo..k],
                y,
                &mut cs[..m],
                &mut cd[..m],
            );
            if let Some(j) = K::first_max_above(&cs[..m], best) {
                (best, arg, bd) = (cs[j], (lo + j) as u32, cd[j]);
            }
        }
        f[k] = best;
        back[k] = arg;
        disc[k] = bd;
    }
}

impl Chainer for OracleChainer {
    fn chain(&self, idx: &Index, cand: &Candidate) -> Vec<TxHit> {
        DP.with(|dp| self.chain_with(idx, cand, &mut dp.borrow_mut()))
    }
}

impl OracleChainer {
    fn chain_with(&self, idx: &Index, cand: &Candidate, sc: &mut DpScratch) -> Vec<TxHit> {
        let level = crate::align::simd_level();
        let unitig = idx.node_unitig();
        let offset = idx.node_offset();
        // An anchor's transcript span equals its query span (exact core matches; merged
        // collinear runs).
        let len = |i: usize| cand.anchors[i].qe - cand.anchors[i].qs;
        // Bucket placements by transcript; visiting anchors in order keeps each bucket sorted
        // by (anchor, target position).
        // Transcript -> bucket + 1 (0: none yet), dense; reset below for the transcripts used.
        if sc.bucket_of.len() < idx.n_transcripts() {
            sc.bucket_of.resize(idx.n_transcripts(), 0);
        }
        sc.tx.clear();
        let mut used = 0;
        // Occurrence lookups miss in cache twice (offsets, then lists): issue them ahead.
        const AHEAD: usize = 4;
        for (i, a) in cand.anchors.iter().enumerate() {
            if let Some(x) = cand.anchors.get(i + 2 * AHEAD) {
                idx.prefetch_occ_off(unitig[x.node as usize]);
            }
            if let Some(x) = cand.anchors.get(i + AHEAD) {
                idx.prefetch_occurrences(unitig[x.node as usize]);
            }
            let (ts, ps) = idx.occurrences(unitig[a.node as usize]);
            let off = offset[a.node as usize];
            for (&t, &p) in ts.iter().zip(ps) {
                let mut b = sc.bucket_of[t as usize];
                if b == 0 {
                    if used == sc.buckets.len() {
                        sc.buckets.push(Vec::new());
                    }
                    sc.buckets[used].clear();
                    sc.tx.push(t);
                    used += 1;
                    b = used as u32;
                    sc.bucket_of[t as usize] = b;
                }
                sc.buckets[b as usize - 1].push(Placed {
                    a: i as u32,
                    ts: p + off,
                });
            }
        }
        for &t in &sc.tx {
            sc.bucket_of[t as usize] = 0;
        }

        crate::stats::add(&crate::stats::CHAIN_TX, used as u64);
        crate::stats::add(
            &crate::stats::CHAIN_PLACEMENTS,
            sc.buckets[..used].iter().map(|b| b.len() as u64).sum(),
        );
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
            if p.beta == 0.0 {
                if sc.cs.len() < n {
                    sc.cs.resize(n, 0.0);
                    sc.cd.resize(n, 0.0);
                }
                dp_beta0(
                    level,
                    p,
                    (&sc.a, &sc.qs, &sc.qe, &sc.ts, &sc.te, &sc.w),
                    &mut sc.f,
                    &mut sc.back,
                    &mut sc.disc,
                    lcp,
                    &mut sc.cs,
                    &mut sc.cd,
                );
            }
            let mut lo = 0;
            for k in (lcp..n).filter(|_| p.beta != 0.0) {
                let ya = sc.a[k];
                while (sc.a[lo] as usize) + p.h < ya as usize {
                    lo += 1;
                }
                let (mut best, mut arg, mut bd) = (sc.w[k], u32::MAX, 0.0f32);
                if k > lo {
                    {
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
            let path = &mut sc.path;
            path.clear();
            loop {
                path.push(k);
                maxd = if sc.disc[k] > maxd { sc.disc[k] } else { maxd };
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

#[cfg(test)]
mod best_chain_tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        /// The vectorized β = 0 path returns exactly what the scalar DP does.
        #[test]
        fn beta0_matches_scalar(
            raw in prop::collection::vec((0u32..3000, 0u32..3000, 1u32..40), 0..400),
            h in 1usize..40,
            alpha in prop_oneof![Just(1.5f32), Just(0.25f32)],
        ) {
            let mut pts: Vec<Pt> = raw
                .iter()
                .map(|&(qs, ts, len)| Pt { qs, qe: qs + len, ts, len, w: len as f32 })
                .collect();
            pts.sort_unstable_by_key(|p| (p.qs, p.ts));
            let p = ChainParams { alpha, h, o_max: 30, g_max: 500, ..Default::default() };
            let want = best_chain_scalar(&p, &pts);
            // The vectorized DP at every SIMD level this CPU supports, each forced explicitly.
            for level in lcp_core::simd::supported() {
                let got = crate::align::at_level(level, || best_chain(&p, &pts));
                match (&want, got) {
                    (None, None) => {}
                    (Some(w), Some(g)) => {
                        prop_assert_eq!(w.0.to_bits(), g.0.to_bits(), "at {}", level);
                        prop_assert_eq!(&w.1, &g.1, "at {}", level);
                        prop_assert_eq!(w.2.to_bits(), g.2.to_bits(), "at {}", level);
                    }
                    (w, g) => prop_assert!(false, "{:?} vs {:?} at {}", w, g, level),
                }
            }
        }

        /// The NEON predecessor scores equal the scalar ones bit for bit, including invalid
        /// predecessors, overlaps and gaps at the limits.
        #[cfg(target_arch = "aarch64")]
        #[test]
        fn pred_scores_neon_matches_scalar(
            xs in prop::collection::vec((0u32..60, -50i32..600, 0i32..80, -50i32..600, 0i32..80, -1000.0f32..1000.0), 0..40),
            y in (0u32..60, 0i32..600, 1i32..80, 0i32..600, 0.0f32..100.0),
            alpha in prop_oneof![Just(1.5f32), Just(0.25f32), 0.0f32..4.0],
            rho in prop_oneof![Just(1.0f32), 0.5f32..2.0],
        ) {
            let xa: Vec<u32> = xs.iter().map(|x| x.0).collect();
            let xqs: Vec<i32> = xs.iter().map(|x| x.1).collect();
            let xqe: Vec<i32> = xs.iter().map(|x| x.1 + x.2).collect();
            let xts: Vec<i32> = xs.iter().map(|x| x.3).collect();
            let xte: Vec<i32> = xs.iter().map(|x| x.3 + x.4).collect();
            let xf: Vec<f32> = xs.iter().map(|x| x.5).collect();
            let y = (y.0, y.1, y.1 + y.2, y.3, y.4);
            let n = xs.len();
            let (mut cs, mut cd, mut ws, mut wd) = (vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]);
            pred_scores(alpha, rho, 30, 500, &xa, &xqs, &xqe, &xts, &xte, &xf, y, &mut ws, &mut wd);
            // SAFETY: NEON is part of the aarch64 baseline.
            unsafe { pred_scores_neon(alpha, rho, 30, 500, &xa, &xqs, &xqe, &xts, &xte, &xf, y, &mut cs, &mut cd) };
            let bits = |v: &[f32]| v.iter().map(|x| x.to_bits()).collect::<Vec<u32>>();
            prop_assert_eq!(bits(&cs), bits(&ws));
            prop_assert_eq!(bits(&cd), bits(&wd));
        }

        /// The NEON maximum and first position at it equal the scalar ones, with ties,
        /// signed zeros and unreachable (-inf) predecessors.
        #[cfg(target_arch = "aarch64")]
        #[test]
        fn max_first_neon_matches_scalar(
            xs in prop::collection::vec(prop_oneof![
                Just(f32::NEG_INFINITY), Just(0.0f32), Just(-0.0f32), Just(1.5f32), -100.0f32..100.0,
            ], 0..60),
        ) {
            let want = max_score(&xs);
            // SAFETY: NEON is part of the aarch64 baseline.
            let (mx, j) = unsafe { max_first_neon(&xs) };
            prop_assert!(mx == want, "{} vs {}", mx, want);
            if want > f32::NEG_INFINITY {
                prop_assert_eq!(j, xs.iter().position(|&v| v == want).unwrap());
                prop_assert_eq!(xs[j].to_bits(), xs[xs.iter().position(|&v| v == want).unwrap()].to_bits());
            }
        }

        /// Shared-prefix chaining of several lists equals one call per list: lists are
        /// shifted copies of one point set, cut and extended differently (as isoform windows
        /// are), plus unrelated ones.
        #[test]
        fn shared_lists_match_single(
            raw in prop::collection::vec((0u32..3000, 0u32..3000, 1u32..40), 0..300),
            cuts in prop::collection::vec((0usize..300, 0u32..2000, prop::collection::vec((0u32..3000, 0u32..3000, 1u32..40), 0..40)), 1..12),
            h in 1usize..40,
        ) {
            let mk = |v: &[(u32, u32, u32)], shift: u32| {
                let mut pts: Vec<Pt> = v
                    .iter()
                    .map(|&(qs, ts, len)| Pt { qs, qe: qs + len, ts: ts + shift, len, w: len as f32 })
                    .collect();
                pts.sort_unstable_by_key(|p| (p.qs, p.ts));
                pts
            };
            let base = mk(&raw, 0);
            let lists: Vec<Vec<Pt>> = cuts
                .iter()
                .map(|(cut, shift, extra)| {
                    // A shifted prefix of the base, then other points after it in the query.
                    let mut l: Vec<Pt> = base[..(*cut).min(base.len())]
                        .iter()
                        .map(|p| Pt { ts: p.ts + shift, ..*p })
                        .collect();
                    let q0 = l.last().map_or(0, |p| p.qs + 1);
                    let mut tail = mk(extra, *shift);
                    for p in &mut tail {
                        p.qs += q0;
                        p.qe += q0;
                    }
                    l.extend(tail);
                    l
                })
                .collect();
            let p = ChainParams { h, o_max: 30, g_max: 500, ..Default::default() };
            let refs: Vec<&[Pt]> = lists.iter().map(|l| l.as_slice()).collect();
            for level in lcp_core::simd::supported() {
                let got = crate::align::at_level(level, || best_chains(&p, &refs));
                for (l, g) in lists.iter().zip(&got) {
                    let w = best_chain_scalar(&p, l);
                    match (w, g) {
                        (None, None) => {}
                        (Some(w), Some(g)) => {
                            prop_assert_eq!(w.0.to_bits(), g.0.to_bits(), "at {}", level);
                            prop_assert_eq!(&w.1, &g.1, "at {}", level);
                            prop_assert_eq!(w.2.to_bits(), g.2.to_bits(), "at {}", level);
                        }
                        (w, g) => prop_assert!(false, "{:?} vs {:?} at {}", w, g, level),
                    }
                }
            }
        }
    }
}
