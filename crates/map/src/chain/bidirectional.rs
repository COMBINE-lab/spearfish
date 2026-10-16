//! Exact, certified bidirectional chaining. See `NUMERICS.md` in this directory.
use super::*;
use std::cmp::Ordering;
use std::ops::Range;

/// Prefix/suffix sharing with the original forward oracle as an exact fallback.
/// `OracleChainer` itself remains the forward-only correctness reference.
pub struct BidirectionalChainer {
    pub p: ChainParams,
}

#[derive(Default)]
struct Work {
    a: Vec<u32>,
    qs: Vec<u32>,
    qe: Vec<u32>,
    ts: Vec<u32>,
    len: Vec<u32>,
    w: Vec<f32>,
    split: usize,
    f: Vec<f32>,
    back: Vec<usize>,
    fd: Vec<f32>,
    h: Vec<f32>,
    // Distance to the next placement, so pointers survive suffix-relative reuse.
    next: Vec<usize>,
    hd: Vec<f32>,
}

impl Work {
    fn load(&mut self, cand: &Candidate, el: &[Placed], boundary: usize) -> Option<()> {
        self.a.clear();
        self.qs.clear();
        self.qe.clear();
        self.ts.clear();
        self.len.clear();
        self.w.clear();
        let base = el.first()?.ts as i64 - REL_BIAS;
        for x in el {
            let a = &cand.anchors[x.a as usize];
            let ts = x.ts as i64 - base;
            let len = a.qe - a.qs;
            // The forward SIMD oracle uses biased i32 coordinates. Certify that neither
            // its conversions nor its subtractions overflow, and that scalar Pt agrees.
            if a.qe > i32::MAX as u32
                || ts < 0
                || ts + len as i64 > i32::MAX as i64
                || x.ts.checked_add(len).is_none()
            {
                return None;
            }
            self.a.push(x.a);
            self.qs.push(a.qs);
            self.qe.push(a.qe);
            self.ts.push(ts as u32);
            self.len.push(len);
            self.w.push(a.weight);
        }
        self.split = self.a.partition_point(|&a| (a as usize) < boundary);
        let n = el.len();
        self.f.resize(self.split, 0.0);
        self.back.resize(self.split, usize::MAX);
        self.fd.resize(self.split, 0.0);
        self.h.resize(n, 0.0);
        self.next.resize(n, 0);
        self.hd.resize(n, 0.0);
        Some(())
    }

    fn point(&self, i: usize) -> Pt {
        Pt {
            qs: self.qs[i],
            qe: self.qe[i],
            ts: self.ts[i],
            len: self.len[i],
            w: self.w[i],
        }
    }

    fn side_len(&self, reverse: bool) -> usize {
        if reverse {
            self.a.len() - self.split
        } else {
            self.split
        }
    }

    fn rank(&self, rank: usize, reverse: bool) -> usize {
        if reverse { self.a.len() - 1 - rank } else { rank }
    }

    fn key(&self, rank: usize, reverse: bool) -> (u32, i64) {
        let i = self.rank(rank, reverse);
        let origin = self.rank(0, reverse);
        (self.a[i], self.ts[i] as i64 - self.ts[origin] as i64)
    }

    fn compare(&self, other: &Self, reverse: bool) -> Ordering {
        (0..self.side_len(reverse))
            .map(|r| self.key(r, reverse))
            .cmp((0..other.side_len(reverse)).map(|r| other.key(r, reverse)))
    }

    fn common(&self, other: &Self, reverse: bool) -> usize {
        (0..self.side_len(reverse).min(other.side_len(reverse)))
            .take_while(|&r| self.key(r, reverse) == other.key(r, reverse))
            .count()
    }
}

#[derive(Default)]
struct Scratch {
    plan: Plan,
    failed: Vec<bool>,
    work: Vec<Work>,
    order: Vec<usize>,
    delta: Vec<f32>,
    scores: Vec<f32>,
    disc: Vec<f32>,
    path: Vec<usize>,
    #[cfg(test)]
    reused_suffix: usize,
    #[cfg(test)]
    evaluated_pairs: u128,
}

thread_local! {
    static BIDI: std::cell::RefCell<Scratch> = std::cell::RefCell::new(Scratch::default());
}

// Independent lanes, SoA loads, no allocations. Backward evaluation always evaluates the
// ORIGINAL directed edge fixed -> lane, including the successor's weight/overlap gain.
#[multiversion::multiversion(targets(
    "x86_64+avx512f+avx512bw+avx512vl+avx2+fma",
    "x86_64+avx2+fma",
    "x86_64+avx",
    "aarch64+neon"
))]
fn edge_batch(
    p: &ChainParams,
    w: &Work,
    fixed: usize,
    range: Range<usize>,
    successors: bool,
    delta: &mut [f32],
    disc: &mut [f32],
) {
    let fixed_pt = w.point(fixed);
    if p.beta == 0.0 {
        // Match the original forward SIMD arithmetic. The certificate requires nonempty
        // query spans, so the branch-free gain also agrees with transition at zero overlap.
        for (lane, k) in range.enumerate() {
            let (j, i, x, y) = if successors {
                (fixed, k, fixed_pt, w.point(k))
            } else {
                (k, fixed, w.point(k), fixed_pt)
            };
            let dq = y.qs as i32 - x.qe as i32;
            let dt = y.ts as i32 - (x.ts + x.len) as i32;
            let valid = w.a[j] < w.a[i]
                && (w.a[i] - w.a[j]) as usize <= p.h
                && y.qs > x.qs
                && y.ts > x.ts
                && dq >= -(p.o_max as i32)
                && dq <= p.g_max as i32
                && dt >= -(p.o_max as i32)
                && dt <= p.g_max as i32;
            let d = (dq as f32 - p.rho * dt as f32).abs();
            let overlap = (x.qe as i32 - y.qs as i32).max(0).min((y.qe - y.qs) as i32) as f32;
            let e = y.w * (1.0 - overlap / (y.qe - y.qs) as f32) - p.alpha * d;
            delta[lane] = if valid {
                if e.is_finite() { e } else { f32::NAN }
            } else {
                f32::NEG_INFINITY
            };
            disc[lane] = d;
        }
        return;
    }
    for (lane, k) in range.enumerate() {
        let (j, i, x, y) = if successors {
            (fixed, k, fixed_pt, w.point(k))
        } else {
            (k, fixed, w.point(k), fixed_pt)
        };
        let edge = if w.a[j] < w.a[i] && (w.a[i] - w.a[j]) as usize <= p.h {
            transition(p, x, y)
        } else {
            None
        };
        let (e, d) = edge
            .map(|(e, d)| (if e.is_finite() { e } else { f32::NAN }, d))
            .unwrap_or((f32::NEG_INFINITY, 0.0));
        delta[lane] = e;
        disc[lane] = d;
    }
}

// All weights and increments must be integral. Bounding every intermediate sum by 2^24
// certifies exact f32 arithmetic in BOTH associations, including losing DP alternatives.
#[derive(Clone, Copy)]
struct Certificate {
    terms: f64,
    max: f64,
}
impl Certificate {
    fn candidate(p: &ChainParams, cand: &Candidate) -> Option<Self> {
        if cand.anchors.len() > u32::MAX as usize
            || !p.alpha.is_finite()
            || !p.beta.is_finite()
            || !p.rho.is_finite()
            || !(0..=i32::MAX as i64).contains(&p.o_max)
            || !(0..=i32::MAX as i64).contains(&p.g_max)
        {
            return None;
        }
        let mut cert = Self {
            terms: cand.anchors.len() as f64 + 1.0,
            max: 0.0,
        };
        for a in &cand.anchors {
            if a.qe <= a.qs {
                return None;
            }
            cert.include(a.weight)?;
        }
        Some(cert)
    }
    fn include(&mut self, x: f32) -> Option<()> {
        if !x.is_finite() || x.fract() != 0.0 || (x == 0.0 && x.is_sign_negative()) {
            return None;
        }
        self.max = self.max.max((x as f64).abs());
        (self.terms * self.max <= (1u64 << 24) as f64).then_some(())
    }

    fn batch(&mut self, delta: &[f32]) -> Option<()> {
        // Reduce once per batch instead of a scalar f64 bound/multiply and early-return
        // branch per edge. Accepted magnitudes are <=2^24, so i32 round-tripping is an
        // exact integrality check. Invalid-edge sentinels contribute zero to the bound.
        let mut valid = true;
        let mut magnitude = 0.0f32;
        for &e in delta {
            let x = if e == f32::NEG_INFINITY { 0.0 } else { e };
            valid &= x.is_finite() && x as i32 as f32 == x && x.to_bits() != (-0.0f32).to_bits();
            magnitude = magnitude.max(x.abs());
        }
        if !valid {
            return None;
        }
        self.include(magnitude)
    }
}

/// Transition-count planning over the original buckets with reusable index arrays. Avoid
/// copying SoA/DP state when sharing cannot pay for joins and the extra placement passes.
#[derive(Default)]
struct Plan {
    forward: Vec<usize>,
    backward: Vec<usize>,
    prefix_common: Vec<usize>,
    suffix_common: Vec<usize>,
    split: Vec<usize>,
    certificate: Option<Certificate>,
    forward_pairs: u128,
    bidirectional_pairs: u128,
    placements: u128,
}

fn common_placements(u: &[Placed], v: &[Placed], reverse: bool) -> usize {
    if u.is_empty() || v.is_empty() {
        return 0;
    }
    let (ub, vb) = if reverse {
        (u[u.len() - 1].ts, v[v.len() - 1].ts)
    } else {
        (u[0].ts, v[0].ts)
    };
    let same = |(x, y): (&Placed, &Placed)| x.a == y.a && x.ts as i64 - ub as i64 == y.ts as i64 - vb as i64;
    if reverse {
        u.iter()
            .rev()
            .zip(v.iter().rev())
            .take_while(|&(x, y)| same((x, y)))
            .count()
    } else {
        u.iter().zip(v).take_while(|&(x, y)| same((x, y))).count()
    }
}

impl Plan {
    fn build(
        &mut self,
        p: &ChainParams,
        cand: &Candidate,
        buckets: &[Vec<Placed>],
        order: &[u32],
        cert: Certificate,
    ) -> bool {
        let used = order.len();
        let boundary = cand.anchors.len() / 2;
        self.certificate = Some(cert);
        self.forward.clear();
        self.forward.extend(order.iter().map(|&i| i as usize));
        self.backward.clear();
        self.backward.extend(0..used);
        self.prefix_common.resize(used, 0);
        self.suffix_common.resize(used, 0);
        self.split.clear();
        self.split.extend(
            buckets
                .iter()
                .map(|b| b.partition_point(|x| (x.a as usize) < boundary)),
        );
        self.forward_pairs = 0;
        self.bidirectional_pairs = 0;
        self.placements = 0;
        let mut previous: Option<usize> = None;
        for &bi in &self.forward {
            let el = &buckets[bi];
            let split = self.split[bi];
            let common = previous.map_or(0, |prev| common_placements(el, &buckets[prev], false));
            self.prefix_common[bi] = common.min(split);
            let (mut lo, mut group) = (0, 0);
            for k in 0..el.len() {
                let a = el[k].a as usize;
                while (el[lo].a as usize) < a.saturating_sub(p.h) {
                    lo += 1;
                }
                while el[group].a < el[k].a {
                    group += 1;
                }
                // The original kernel visits same-identity placements and rejects them;
                // the optimized kernel excludes these lanes before evaluating arithmetic.
                if k >= common {
                    self.forward_pairs += (k - lo) as u128;
                }
                if k < split && k >= self.prefix_common[bi] {
                    self.bidirectional_pairs += (group - lo) as u128;
                }
            }
            let mut lo = 0;
            for i in split..el.len() {
                if (el[i].a as usize) - boundary >= p.h {
                    break;
                }
                while lo < split && (el[lo].a as usize) < (el[i].a as usize).saturating_sub(p.h) {
                    lo += 1;
                }
                self.bidirectional_pairs += (split - lo) as u128;
            }
            self.placements += el.len() as u128;
            previous = Some(bi);
        }
        // Even with perfect suffix sharing there is no saving if the prefix/join cost
        // already exhausts the original forward work. Avoid the second sort in this case.
        if self.bidirectional_pairs + self.placements * 4 >= self.forward_pairs {
            return false;
        }
        self.backward.sort_unstable_by(|&x, &y| {
            let (u, v) = (&buckets[x][self.split[x]..], &buckets[y][self.split[y]..]);
            let ub = u.last().map_or(0, |p| p.ts) as i64;
            let vb = v.last().map_or(0, |p| p.ts) as i64;
            u.iter()
                .rev()
                .map(|p| (p.a, p.ts as i64 - ub))
                .cmp(v.iter().rev().map(|p| (p.a, p.ts as i64 - vb)))
        });
        previous = None;
        for &bi in &self.backward {
            let el = &buckets[bi];
            let split = self.split[bi];
            let common = previous.map_or(0, |prev| {
                common_placements(&el[split..], &buckets[prev][self.split[prev]..], true)
            });
            self.suffix_common[bi] = common;
            let (mut hi, mut group) = (el.len(), el.len());
            for k in (split..el.len()).rev() {
                while (el[hi - 1].a - el[k].a) as usize > p.h {
                    hi -= 1;
                }
                while el[group - 1].a > el[k].a {
                    group -= 1;
                }
                if el.len() - 1 - k >= common {
                    self.bidirectional_pairs += (hi - group) as u128;
                }
            }
            previous = Some(bi);
        }
        self.bidirectional_pairs + self.placements * 4 < self.forward_pairs
    }
}

impl Scratch {
    fn batch(
        &mut self,
        p: &ChainParams,
        bi: usize,
        fixed: usize,
        r: Range<usize>,
        successors: bool,
        cert: &mut Certificate,
    ) -> Option<()> {
        #[cfg(test)]
        {
            self.evaluated_pairs += r.len() as u128;
        }
        self.delta.resize(r.len(), 0.0);
        self.disc.resize(r.len(), 0.0);
        edge_batch(
            p,
            &self.work[bi],
            fixed,
            r,
            successors,
            &mut self.delta,
            &mut self.disc,
        );
        cert.batch(&self.delta)
    }

    fn compute_side(
        &mut self,
        p: &ChainParams,
        bi: usize,
        common: usize,
        reverse: bool,
        cert: &mut Certificate,
    ) -> Option<()> {
        // Monotone identity-window cursors replace a binary search at every
        // endpoint. Exclude all other placements of the same identity up front.
        let n = self.work[bi].a.len();
        let (mut lo, mut hi, mut group_start, mut group_end) = (0, n, 0, n);
        for rank in common..self.work[bi].side_len(reverse) {
            let k = self.work[bi].rank(rank, reverse);
            let a = self.work[bi].a[k] as usize;
            let ids = &self.work[bi].a;
            let range = if reverse {
                while (ids[hi - 1] as usize) - a > p.h {
                    hi -= 1;
                }
                while ids[group_end - 1] > ids[k] {
                    group_end -= 1;
                }
                group_end..hi
            } else {
                while (ids[lo] as usize) < a.saturating_sub(p.h) {
                    lo += 1;
                }
                while ids[group_start] < ids[k] {
                    group_start += 1;
                }
                lo..group_start
            };
            self.batch(p, bi, k, range.clone(), reverse, cert)?;
            let w = &mut self.work[bi];
            self.scores.resize(range.len(), 0.0);
            for (lane, j) in range.clone().enumerate() {
                self.scores[lane] = if reverse {
                    self.delta[lane] + w.h[j]
                } else {
                    w.f[j] + self.delta[lane]
                };
            }
            let (mut best, mut arg, mut d) = if reverse {
                (0.0, 0, 0.0)
            } else {
                (w.w[k], usize::MAX, 0.0)
            };
            // Prefix arithmetic and strict-greater selection are the original
            // forward recurrence, including ties. Only suffix ties need fallback:
            // backward selection cannot in general reproduce forward tie order.
            let mut tied = false;
            for (lane, j) in range.enumerate() {
                let score = self.scores[lane];
                if score > best {
                    best = score;
                    arg = if reverse { j - k } else { j };
                    d = self.disc[lane];
                    tied = false;
                } else if score == best {
                    tied = true;
                }
            }
            if reverse && tied {
                return None;
            }
            if reverse {
                w.h[k] = best;
                w.next[k] = arg;
                w.hd[k] = d;
            } else {
                w.f[k] = best;
                w.back[k] = arg;
                w.fd[k] = d;
            }
        }
        Some(())
    }

    fn join(
        &mut self,
        p: &ChainParams,
        bi: usize,
        boundary: usize,
        el: &[Placed],
        transcript: u32,
        cert: &mut Certificate,
    ) -> Option<TxHit> {
        let split = self.work[bi].split;
        let n = self.work[bi].a.len();
        // Represent the winner by (prefix endpoint, first suffix placement, crossing d).
        let mut best = f32::NEG_INFINITY;
        let mut winner = (None, None, 0.0);
        let mut tied = false;
        let mut consider = |score, prefix, suffix, d| {
            if score > best {
                best = score;
                winner = (prefix, suffix, d);
                tied = false;
            } else if score == best {
                tied = true;
            }
        };
        let w = &self.work[bi];
        for j in 0..split {
            consider(w.f[j], Some(j), None, 0.0);
        }
        for i in split..n {
            consider(w.w[i] + w.h[i], None, Some(i), 0.0);
        }
        // Only the last h prefix identities and first h suffix identities can cross.
        // The per-pair identity check in edge_batch handles missing/repeated anchors.
        let first =
            self.work[bi].a[..split].partition_point(|&a| (a as usize) < boundary.saturating_sub(p.h));
        let end = split
            + self.work[bi].a[split..].partition_point(|&a| (a as usize).saturating_sub(boundary) < p.h);
        let mut lo = first;
        for i in split..end {
            let a = self.work[bi].a[i] as usize;
            while lo < split && (self.work[bi].a[lo] as usize) < a.saturating_sub(p.h) {
                lo += 1;
            }
            if lo == split {
                break;
            }
            self.batch(p, bi, i, lo..split, false, cert)?;
            let w = &self.work[bi];
            self.scores.resize(split - lo, 0.0);
            // Keep arithmetic independent across lanes; reduce/trace ties separately.
            for (lane, j) in (lo..split).enumerate() {
                self.scores[lane] = w.f[j] + self.delta[lane] + w.h[i];
            }
            for (lane, j) in (lo..split).enumerate() {
                if self.delta[lane] != f32::NEG_INFINITY {
                    consider(self.scores[lane], Some(j), Some(i), self.disc[lane]);
                }
            }
        }
        if tied {
            return None;
        }
        let w = &self.work[bi];
        let (prefix, suffix, mut maxd) = winner;
        self.path.clear();
        if let Some(mut j) = prefix {
            loop {
                self.path.push(j);
                maxd = maxd.max(w.fd[j]);
                if w.back[j] == usize::MAX {
                    break;
                }
                j = w.back[j];
            }
            self.path.reverse();
        }
        if let Some(mut i) = suffix {
            loop {
                self.path.push(i);
                maxd = maxd.max(w.hd[i]);
                if w.next[i] == 0 {
                    break;
                }
                i += w.next[i];
            }
        }

        Some(hit_from_path(
            transcript,
            best,
            maxd,
            self.path.iter().map(|&k| (w.qs[k], w.qe[k], el[k].ts, w.len[k])),
        ))
    }

    #[cfg(test)]
    fn compute(
        &mut self,
        p: &ChainParams,
        cand: &Candidate,
        buckets: &[Vec<Placed>],
        tx: &[u32],
    ) -> Option<Vec<TxHit>> {
        self.compute_planned(p, cand, buckets, tx, false)
    }

    fn compute_planned(
        &mut self,
        p: &ChainParams,
        cand: &Candidate,
        buckets: &[Vec<Placed>],
        tx: &[u32],
        planned: bool,
    ) -> Option<Vec<TxHit>> {
        let mut cert = if planned {
            self.plan.certificate?
        } else {
            Certificate::candidate(p, cand)?
        };
        #[cfg(test)]
        {
            self.reused_suffix = 0;
            self.evaluated_pairs = 0;
        }
        let used = tx.len();
        self.work.resize_with(used.max(self.work.len()), Work::default);
        let boundary = cand.anchors.len() / 2;
        self.failed.clear();
        self.failed.resize(used, false);
        for (i, el) in buckets.iter().take(used).enumerate() {
            if self.work[i].load(cand, el, boundary).is_none() {
                if !planned {
                    return None;
                }
                self.failed[i] = true;
            }
        }
        // Each ordering has its own normalized signatures. Empty sides sort first.
        for reverse in [false, true] {
            self.order.clear();
            if planned {
                self.order.extend(if reverse {
                    &self.plan.backward
                } else {
                    &self.plan.forward
                });
            } else {
                self.order.extend(0..used);
                self.order
                    .sort_unstable_by(|&x, &y| self.work[x].compare(&self.work[y], reverse));
            }
            let mut previous = None;
            for oi in 0..used {
                let bi = self.order[oi];
                if self.failed[bi] {
                    continue;
                }
                // Cached common lengths describe adjacent planned transcripts. If one
                // failed, compare with the actual last successful transcript instead.
                let common = if planned && previous == oi.checked_sub(1).map(|j| self.order[j]) {
                    if reverse {
                        self.plan.suffix_common[bi]
                    } else {
                        self.plan.prefix_common[bi]
                    }
                } else {
                    previous.map_or(0, |prev: usize| self.work[bi].common(&self.work[prev], reverse))
                };
                if let Some(prev) = previous {
                    for r in 0..common {
                        let k = self.work[bi].rank(r, reverse);
                        let pk = self.work[prev].rank(r, reverse);
                        if reverse {
                            self.work[bi].h[k] = self.work[prev].h[pk];
                            self.work[bi].next[k] = self.work[prev].next[pk];
                            self.work[bi].hd[k] = self.work[prev].hd[pk];
                        } else {
                            self.work[bi].f[k] = self.work[prev].f[pk];
                            self.work[bi].back[k] = self.work[prev].back[pk];
                            self.work[bi].fd[k] = self.work[prev].fd[pk];
                        }
                    }
                }
                #[cfg(test)]
                {
                    if reverse {
                        self.reused_suffix += common;
                    }
                }
                let saved = cert;
                if self.compute_side(p, bi, common, reverse, &mut cert).is_none() {
                    if !planned {
                        return None;
                    }
                    self.failed[bi] = true;
                    cert = saved;
                    continue;
                }
                previous = Some(bi);
            }
        }
        let mut hits = Vec::with_capacity(used);
        for (bi, &transcript) in tx.iter().enumerate() {
            if self.failed[bi] {
                continue;
            }
            let saved = cert;
            if let Some(hit) = self.join(p, bi, boundary, &buckets[bi], transcript, &mut cert) {
                hits.push(hit);
            } else if planned {
                self.failed[bi] = true;
                cert = saved;
            } else {
                return None;
            }
        }
        Some(hits)
    }
}

impl OracleChainer {
    pub(super) fn chain_bidirectional(&self, idx: &Index, cand: &Candidate) -> Vec<TxHit> {
        DP.with(|dp| {
            let sc = &mut *dp.borrow_mut();
            let used = collect_placements(idx, cand, sc);
            // One transcript has nothing to share; tiny candidates cannot amortize setup.
            if used < 2 || cand.anchors.len() < 4 || self.p.h == 0 {
                return self.chain_placed(cand, sc, used);
            }
            let Some(cert) = Certificate::candidate(&self.p, cand) else {
                return self.chain_placed(cand, sc, used);
            };
            sort_placements(sc, used);
            BIDI.with(|bidi| {
                let scratch = &mut *bidi.borrow_mut();
                if !scratch
                    .plan
                    .build(&self.p, cand, &sc.buckets[..used], &sc.order, cert)
                {
                    return None;
                }
                let mut hits = scratch.compute_planned(&self.p, cand, &sc.buckets[..used], &sc.tx, true)?;
                sc.order.retain(|&bi| scratch.failed[bi as usize]);
                let missing = sc.order.len();
                if missing != 0 {
                    hits.extend(self.chain_ordered(cand, sc, missing));
                }
                Some(hits)
            })
            .unwrap_or_else(|| self.chain_ordered(cand, sc, used))
        })
    }
}

impl Chainer for BidirectionalChainer {
    fn chain(&self, idx: &Index, cand: &Candidate) -> Vec<TxHit> {
        OracleChainer { p: self.p.clone() }.chain_bidirectional(idx, cand)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchor::Anchor;

    fn candidate(weights: &[f32]) -> Candidate {
        Candidate {
            family: 0,
            score: 0.0,
            anchors: weights
                .iter()
                .enumerate()
                .map(|(i, &weight)| Anchor {
                    qs: i as u32 * 20,
                    qe: i as u32 * 20 + 10,
                    node: 0,
                    weight,
                })
                .collect(),
        }
    }

    fn placements(items: &[(u32, u32)]) -> Vec<Placed> {
        items.iter().map(|&(a, ts)| Placed { a, ts }).collect()
    }

    fn check(c: &Candidate, p: &ChainParams, buckets: &[Vec<Placed>]) -> (Option<Vec<TxHit>>, Scratch) {
        let tx: Vec<_> = (0..buckets.len() as u32).collect();
        let mut reference = DpScratch {
            buckets: buckets.to_vec(),
            tx: tx.clone(),
            ..Default::default()
        };
        let mut expected = OracleChainer { p: p.clone() }.chain_placed(c, &mut reference, tx.len());
        let mut sc = Scratch::default();
        let fast = sc.compute(p, c, buckets, &tx);
        let mut actual = fast
            .clone()
            .unwrap_or_else(|| OracleChainer { p: p.clone() }.chain_placed(c, &mut reference, tx.len()));
        expected.sort_by_key(|h| h.tx);
        actual.sort_by_key(|h| h.tx);
        assert_eq!(actual, expected);
        for (a, e) in actual.iter().zip(&expected) {
            assert_eq!(a.score.to_bits(), e.score.to_bits());
            assert_eq!(a.max_discrepancy.to_bits(), e.max_discrepancy.to_bits());
        }
        if fast.is_some() {
            // Independently computed, unshared continuation DP using original coordinates.
            // Check EVERY suffix endpoint, including its next pointer and discrepancy.
            for (bi, el) in buckets.iter().enumerate() {
                let w = &sc.work[bi];
                let point = |k: usize| Pt {
                    ts: el[k].ts,
                    ..w.point(k)
                };
                let mut h = vec![0.0; el.len()];
                let mut next = vec![0; el.len()];
                let mut disc = vec![0.0; el.len()];
                for i in (w.split..el.len()).rev() {
                    for k in i + 1..el.len() {
                        if el[k].a <= el[i].a || (el[k].a - el[i].a) as usize > p.h {
                            continue;
                        }
                        if let Some((e, d)) = transition(p, point(i), point(k)) {
                            let s = e + h[k];
                            if s > h[i] {
                                h[i] = s;
                                next[i] = k - i;
                                disc[i] = d;
                            }
                        }
                    }
                }
                assert_eq!(&w.h[w.split..], &h[w.split..]);
                assert_eq!(&w.next[w.split..], &next[w.split..]);
                assert_eq!(&w.hd[w.split..], &disc[w.split..]);
            }
        }
        (fast, sc)
    }

    #[test]
    fn shared_prefixes_and_suffixes_use_independent_orders() {
        let c = candidate(&[40., 44., 52., 60., 68., 76., 84., 92.]);
        let whole: Vec<_> = (0..8).map(|a| (a, a * 20)).collect();
        let mut different_prefix = whole.clone();
        different_prefix[0].1 += 2;
        let mut different_suffix = whole.clone();
        different_suffix[7].1 += 2;
        let mut both = different_prefix.clone();
        both[7].1 += 4;
        let shifted: Vec<_> = whole.iter().map(|&(a, t)| (a, t + 700)).collect();
        let buckets: Vec<_> = [
            whole.clone(),
            different_prefix,
            different_suffix,
            both,
            shifted,
            whole,
        ]
        .iter()
        .map(|v| placements(v))
        .collect();
        let (fast, sc) = check(&c, &ChainParams::default(), &buckets);
        assert!(fast.is_some());
        assert!(
            sc.reused_suffix >= 12,
            "backward sharing must actually be exercised"
        );
    }

    #[test]
    fn missing_and_repeated_placements_share_suffix_without_fallback() {
        let c = candidate(&[40., 44., 52., 60., 68., 76., 84., 92.]);
        let mut first: Vec<_> = (0..8).map(|a| (a, a * 20)).collect();
        first.extend([(1, 2000), (5, 2000)]);
        first.sort_unstable();
        let second: Vec<_> = first.iter().copied().filter(|&(a, _)| a != 0 && a != 2).collect();
        let third: Vec<_> = first
            .iter()
            .copied()
            .filter(|&(a, _)| a != 6)
            .map(|(a, t)| (a, t + 300))
            .collect();
        let (fast, sc) = check(
            &c,
            &ChainParams::default(),
            &[placements(&first), placements(&second), placements(&third)],
        );
        assert!(fast.is_some());
        assert!(sc.reused_suffix >= 5);
    }

    #[test]
    fn exact_overlap_gain_uses_successor_weight() {
        let mut c = candidate(&[112., 80., 48., 32.]);
        for (i, a) in c.anchors.iter_mut().enumerate() {
            a.qs = i as u32 * 8;
            a.qe = a.qs + 16;
        }
        let b = vec![placements(&[(0, 0), (1, 8), (2, 16), (3, 24)])];
        let hit = check(&c, &ChainParams::default(), &b).0.unwrap().remove(0);
        assert_eq!(hit.score, 112. + 40. + 24. + 16.);
        assert_eq!(hit.anchored, 40);
        assert_eq!(hit.anchors, 4);
    }

    #[test]
    fn discrepancy_is_maximum_over_both_halves_and_crossing() {
        let c = candidate(&[40., 60., 80., 100.]);
        let b = vec![placements(&[(0, 0), (1, 22), (2, 48), (3, 80)])];
        let p = ChainParams {
            alpha: 1.,
            ..Default::default()
        };
        let hit = check(&c, &p, &b).0.unwrap().remove(0);
        assert_eq!(hit.score, 260.);
        assert_eq!(hit.max_discrepancy, 12.);
        assert_eq!(hit.anchors, 4);
    }

    #[test]
    fn equal_query_starts_predecessor_ties_preserve_earliest_forward_path() {
        let mut c = candidate(&[10., 20., 20., 40.]);
        c.anchors[2].qs = 20;
        c.anchors[2].qe = 30;
        let b = vec![placements(&[(0, 0), (1, 20), (2, 20), (3, 60)])];
        let p = ChainParams {
            alpha: 0.,
            ..Default::default()
        };
        assert!(check(&c, &p, &b).0.is_none());
        let mut sc = DpScratch {
            buckets: b,
            tx: vec![0],
            ..Default::default()
        };
        let hit = OracleChainer { p }.chain_placed(&c, &mut sc, 1).remove(0);
        assert_eq!(hit.chain, vec![(0, 0, 10), (20, 20, 10), (60, 60, 10)]);
        assert_eq!(hit.score, 70.);
    }

    #[test]
    fn scratch_survives_alternating_candidates_and_fallback() {
        let p = ChainParams::default();
        let c = candidate(&[40., 44., 52., 60., 68., 76.]);
        let b = vec![placements(&[
            (0, 0),
            (1, 20),
            (2, 40),
            (3, 60),
            (4, 80),
            (5, 100),
        ])];
        let tx = vec![0];
        let mut sc = Scratch::default();
        let expected = sc.compute(&p, &c, &b, &tx).unwrap();
        assert!(
            sc.compute(&p, &candidate(&[1.25]), &[placements(&[(0, 0)])], &tx)
                .is_none()
        );
        assert!(sc.compute(&p, &candidate(&[]), &[], &[]).unwrap().is_empty());
        sc.compute(&p, &candidate(&[17.]), &[placements(&[(0, 7)])], &tx)
            .unwrap();
        assert_eq!(sc.compute(&p, &c, &b, &tx).unwrap(), expected);
    }

    #[test]
    fn batched_directed_edges_match_scalar_transition_bits() {
        let mut rng = 0xc0ffee123456789u64;
        let mut next = || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            rng
        };
        for _ in 0..120 {
            let mut c = candidate(&[1.; 10]);
            for a in &mut c.anchors {
                a.qs = (next() % 100) as u32;
                a.qe = a.qs + 1 + (next() % 60) as u32;
                a.weight = (next() % 1000) as f32 / 16.;
            }
            let el: Vec<_> = (0..10)
                .map(|a| Placed {
                    a,
                    ts: (next() % 200) as u32,
                })
                .collect();
            let mut w = Work::default();
            w.load(&c, &el, 5).unwrap();
            let mut delta = vec![0.; 10];
            let mut disc = vec![0.; 10];
            let p = ChainParams {
                h: (next() % 5) as usize,
                beta: if next() % 2 == 0 { 0. } else { 0.5 },
                rho: if next() % 2 == 0 { 1. } else { 1.25 },
                o_max: 30,
                g_max: 70,
                ..Default::default()
            };
            for fixed in 0..10 {
                for successors in [false, true] {
                    edge_batch(&p, &w, fixed, 0..10, successors, &mut delta, &mut disc);
                    for k in 0..10 {
                        let (j, i) = if successors { (fixed, k) } else { (k, fixed) };
                        let expected = if w.a[j] < w.a[i] && (w.a[i] - w.a[j]) as usize <= p.h {
                            transition(&p, w.point(j), w.point(i))
                        } else {
                            None
                        };
                        if let Some((e, d)) = expected {
                            assert_eq!(delta[k].to_bits(), e.to_bits());
                            assert_eq!(disc[k].to_bits(), d.to_bits());
                        } else {
                            assert_eq!(delta[k], f32::NEG_INFINITY);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn compatible_join_uses_lower_scoring_prefix_endpoint_and_crossing_discrepancy() {
        let c = candidate(&[100., 30., 80., 90.]);
        let p = ChainParams {
            alpha: 1.,
            g_max: 50,
            o_max: 5,
            ..Default::default()
        };
        let (fast, _) = check(&c, &p, &[placements(&[(0, 1000), (1, 20), (2, 44), (3, 64)])]);
        let h = &fast.unwrap()[0];
        assert_eq!(h.score, 196.);
        assert_eq!(h.max_discrepancy, 4.);
        assert_eq!(h.chain, vec![(20, 20, 10), (40, 44, 10), (60, 64, 10)]);
        assert_eq!(
            (h.qstart, h.qend, h.tstart, h.tend, h.anchored),
            (20, 70, 20, 74, 30)
        );
    }

    #[test]
    fn crossing_skips_split_neighbors_at_exact_identity_window() {
        let c = candidate(&[100., 1., 2., 3., 200., 4.]);
        let b = vec![placements(&[
            (0, 0),
            (1, 1000),
            (2, 1000),
            (3, 1000),
            (4, 80),
            (5, 1000),
        ])];
        let mut p = ChainParams {
            h: 4,
            g_max: 100,
            ..Default::default()
        };
        let h = check(&c, &p, &b).0.unwrap().remove(0);
        assert_eq!(h.score, 300.);
        assert_eq!(h.chain, vec![(0, 0, 10), (80, 80, 10)]);
        p.h = 3;
        assert_eq!(check(&c, &p, &b).0.unwrap()[0].score, 200.);
    }

    #[test]
    fn empty_sides_tiny_candidates_and_chains_entirely_in_one_half() {
        let c = candidate(&[10., 20., 30., 40.]);
        let p = ChainParams::default();
        let (fast, _) = check(
            &c,
            &p,
            &[placements(&[(0, 0), (1, 20)]), placements(&[(2, 80), (3, 100)])],
        );
        let hits = fast.unwrap();
        assert_eq!(hits[0].score, 30.);
        assert_eq!(hits[1].score, 70.);
        assert!(check(&candidate(&[]), &p, &[]).0.unwrap().is_empty());
        assert_eq!(
            check(&candidate(&[17.]), &p, &[placements(&[(0, 7)])]).0.unwrap()[0].score,
            17.
        );
        let zero_window = ChainParams { h: 0, ..p };
        assert_eq!(
            check(
                &c,
                &zero_window,
                &[placements(&[(0, 0), (1, 20), (2, 40), (3, 60)])]
            )
            .0
            .unwrap()[0]
                .score,
            40.
        );
    }

    #[test]
    fn overlaps_containment_asymmetric_weights_and_repeated_identities() {
        let mut c = candidate(&[60., 30., 90., 120., 150., 180.]);
        for (a, (qs, qe)) in
            c.anchors
                .iter_mut()
                .zip([(0, 30), (10, 20), (10, 40), (25, 55), (40, 70), (60, 90)])
        {
            a.qs = qs;
            a.qe = qe;
        }
        let b = vec![
            placements(&[
                (0, 0),
                (0, 100),
                (1, 10),
                (2, 10),
                (2, 120),
                (3, 25),
                (4, 40),
                (5, 60),
            ]),
            placements(&[(0, 70), (2, 80), (3, 50), (4, 110), (5, 130)]),
        ];
        check(&c, &ChainParams::default(), &b);
        // The directed overlap gain depends on the successor's weight and query length.
        let p = ChainParams {
            alpha: 0.,
            ..Default::default()
        };
        let x = Pt {
            qs: 0,
            qe: 30,
            ts: 0,
            len: 30,
            w: 60.,
        };
        let y = Pt {
            qs: 10,
            qe: 40,
            ts: 10,
            len: 30,
            w: 90.,
        };
        assert_eq!(transition(&p, x, y).unwrap().0, 90. * (1. - 20. / 30.));
    }

    #[test]
    fn validity_gap_boundaries_and_invalid_backward_positions() {
        let c = candidate(&[100., 200., 400., 800.]);
        for h in [1, 2, 4] {
            for gap in [-11, -10, 0, 20, 21] {
                let b = vec![placements(&[
                    (0, 100),
                    (1, 120),
                    (2, (130 + gap) as u32),
                    (3, (150 + gap) as u32),
                ])];
                check(
                    &c,
                    &ChainParams {
                        h,
                        o_max: 10,
                        g_max: 20,
                        ..Default::default()
                    },
                    &b,
                );
            }
        }
        check(
            &c,
            &ChainParams::default(),
            &[placements(&[(0, 100), (1, 80), (2, 60), (3, 40)])],
        );
    }

    #[test]
    fn numerical_and_tie_fallbacks_preserve_forward_hits() {
        let b = vec![placements(&[(0, 0), (1, 20), (2, 40), (3, 60)])];
        let p = ChainParams::default();
        assert!(
            check(&candidate(&[10.25, 20.5, 30.75, 40.125]), &p, &b)
                .0
                .is_none()
        );
        assert!(check(&candidate(&[8_000_000.; 4]), &p, &b).0.is_none());
        // Forward addition loses each +1 at 2^24, while backward association could
        // preserve them. The certificate must reject this before selecting a path.
        assert!(check(&candidate(&[16_777_216., 1., 1., 1.]), &p, &b).0.is_none());
        assert!(check(&candidate(&[-0.0, 20., 30., 40.]), &p, &b).0.is_none());
        let mut shifted = b.clone();
        shifted[0][2].ts += 1;
        assert!(check(&candidate(&[10., 20., 30., 40.]), &p, &shifted).0.is_none());
        // Final endpoint tie: forward chooses the earliest placement.
        let p0 = ChainParams { h: 0, ..p.clone() };
        assert!(check(&candidate(&[10.; 4]), &p0, &b).0.is_none());
        // Continuation stop tie and a local restart/predecessor tie.
        assert!(check(&candidate(&[10., 20., 0., 0.]), &p, &b).0.is_none());
        let fractional = ChainParams {
            beta: 0.5,
            rho: 1.25,
            ..p.clone()
        };
        assert!(
            check(&candidate(&[10., 20., 30., 40.]), &fractional, &b)
                .0
                .is_none()
        );
        // Nondefault rho and beta can still be certified when every actual edge is integral.
        check(
            &candidate(&[10., 20., 30., 40.]),
            &ChainParams {
                beta: 1.,
                rho: 0.5,
                ..p
            },
            &b,
        );
    }

    #[test]
    fn randomized_full_hit_differential_and_independent_continuations() {
        let mut rng = 0x9e3779b97f4a7c15u64;
        let mut next = || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            rng
        };
        let mut certified = 0;
        for case in 0..1600 {
            let n = 1 + next() as usize % 18;
            let mut c = candidate(&(0..n).map(|i| (i * 12 + 12) as f32).collect::<Vec<_>>());
            let simple = case % 4 == 0;
            if !simple {
                let mut qs = 0;
                for a in &mut c.anchors {
                    qs += next() as u32 % 15;
                    a.qs = qs;
                    a.qe = qs + 1 + next() as u32 % 30;
                    a.weight = (next() % 200) as f32 / if case % 7 == 0 { 4. } else { 1. };
                }
            }
            let p = if simple {
                ChainParams::default()
            } else {
                ChainParams {
                    h: next() as usize % 7,
                    alpha: [0., 0.5, 1., 1.5][next() as usize % 4],
                    rho: [0.5, 1., 1.25][next() as usize % 3],
                    beta: if case % 5 == 0 { 0.5 } else { 0. },
                    o_max: (next() % 40) as i64,
                    g_max: (next() % 100) as i64,
                }
            };
            let mut buckets = Vec::new();
            for t in 0..1 + next() % 6 {
                let mut el = Vec::new();
                for i in 0..n {
                    if !simple && next() % 4 == 0 {
                        continue;
                    }
                    for _ in 0..if simple { 1 } else { 1 + next() % 3 } {
                        let ts = if simple {
                            i as u32 * 20 + t as u32 * 100
                        } else {
                            (next() % 300) as u32
                        };
                        el.push(Placed { a: i as u32, ts });
                    }
                }
                el.sort_by_key(|v| (v.a, v.ts));
                if !el.is_empty() {
                    buckets.push(el);
                }
            }
            if check(&c, &p, &buckets).0.is_some() {
                certified += 1;
            }
        }
        assert!(
            certified >= 400,
            "randomized testing must exercise the certified implementation"
        );
    }

    #[test]
    fn randomized_divergent_halves_missing_anchors_and_suffix_multiplicity() {
        let mut rng = 0x12345678abcdef01u64;
        let mut next = || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            rng
        };
        let mut certified = 0;
        for _ in 0..400 {
            let n = 8 + next() as usize % 8;
            let c = candidate(&(0..n).map(|i| 40. + i as f32 * 12.).collect::<Vec<_>>());
            let boundary = n / 2;
            let missing_prefix = next() as usize % boundary;
            let missing_suffix = boundary + next() as usize % (n - boundary);
            let mut buckets = Vec::new();
            for t in 0..3 {
                let mut el = Vec::new();
                for i in 0..n {
                    if t > 0 && (i == missing_prefix || i == missing_suffix) {
                        continue;
                    }
                    let prefix_shift = if i < boundary { 2 * (next() % 4) as u32 } else { 0 };
                    el.push(Placed {
                        a: i as u32,
                        ts: i as u32 * 20 + t * 300 + prefix_shift,
                    });
                }
                // The same repeated suffix placement is present in the last two
                // transcripts. Its pointer ranks must survive differing prefix lengths.
                el.push(Placed {
                    a: (n - 2) as u32,
                    ts: 50_000 + t * 300,
                });
                el.sort_by_key(|p| (p.a, p.ts));
                buckets.push(el);
            }
            if check(&c, &ChainParams::default(), &buckets).0.is_some() {
                certified += 1;
            }
        }
        assert!(certified >= 300);
    }
    fn planned_check(c: &Candidate, p: &ChainParams, buckets: &[Vec<Placed>]) -> (bool, Scratch) {
        let used = buckets.len();
        let mut reference = DpScratch {
            buckets: buckets.to_vec(),
            tx: (0..used as u32).collect(),
            ..Default::default()
        };
        let oracle = OracleChainer { p: p.clone() };
        let mut expected = oracle.chain_placed(c, &mut reference, used);
        let mut scratch = Scratch::default();
        let accepted = scratch.plan.build(
            p,
            c,
            buckets,
            &reference.order,
            Certificate::candidate(p, c).unwrap(),
        );
        let mut actual = if accepted {
            let mut hits = scratch
                .compute_planned(p, c, buckets, &reference.tx, true)
                .unwrap();
            reference.order.retain(|&bi| scratch.failed[bi as usize]);
            let missing = reference.order.len();
            hits.extend(oracle.chain_ordered(c, &mut reference, missing));
            hits
        } else {
            assert!(
                scratch.work.is_empty(),
                "planner rejection must precede DP allocations"
            );
            oracle.chain_ordered(c, &mut reference, used)
        };
        expected.sort_by_key(|h| h.tx);
        actual.sort_by_key(|h| h.tx);
        assert_eq!(actual, expected);
        for (a, b) in actual.iter().zip(&expected) {
            assert_eq!(a.score.to_bits(), b.score.to_bits());
        }
        (accepted, scratch)
    }

    fn shareable(n: usize, ntx: usize) -> (Candidate, Vec<Vec<Placed>>) {
        let c = candidate(&(0..n).map(|i| 40. + i as f32 * 12.).collect::<Vec<_>>());
        let buckets = (0..ntx)
            .map(|t| {
                (0..n)
                    .map(|i| Placed {
                        a: i as u32,
                        ts: i as u32 * 20 + t as u32 * 1000 + if i == 1 { t as u32 * 2 } else { 0 },
                    })
                    .collect()
            })
            .collect();
        (c, buckets)
    }

    #[test]
    fn planning_rejects_unprofitable_work_before_building_dp_arrays() {
        let p = ChainParams::default();
        let (c, mut buckets) = shareable(256, 8);
        let (accepted, sc) = planned_check(&c, &p, &buckets);
        assert!(accepted);
        assert_eq!(sc.plan.bidirectional_pairs, sc.evaluated_pairs);
        assert!(sc.evaluated_pairs + sc.plan.placements * 4 < sc.plan.forward_pairs);
        // Identical complete lists already get maximal reuse from the original forward DP.
        for t in 1..buckets.len() {
            buckets[t] = buckets[0].clone();
        }
        assert!(!planned_check(&c, &p, &buckets).0);
        // No common suffix: separate final anchor identities by omitting different tails.
        let n = buckets.len();
        for (t, b) in buckets.iter_mut().enumerate() {
            b.truncate(c.anchors.len() - n + t);
        }
        assert!(!planned_check(&c, &p, &buckets).0);
    }

    #[test]
    fn selective_fallback_keeps_other_transcripts_and_skips_failed_sharing_sources() {
        let p = ChainParams::default();
        let (c, mut buckets) = shareable(256, 8);
        buckets[2][253].ts += 1; // fractional suffix increment
        buckets[4][1].ts += 1; // fractional prefix increment
        let last = buckets[6].last_mut().unwrap();
        last.ts -= 2;
        let duplicate = Placed {
            a: last.a,
            ts: last.ts + 4,
        };
        buckets[6].push(duplicate); // suffix tie
        let (accepted, sc) = planned_check(&c, &p, &buckets);
        assert!(accepted);
        assert_eq!(
            sc.failed,
            vec![false, false, true, false, true, false, true, false]
        );
    }

    #[test]
    fn prefix_ties_keep_original_predecessor_without_triggering_fallback() {
        let c = candidate(&[0., 20., 30., 40.]);
        let b = vec![placements(&[(0, 0), (1, 20), (2, 40), (3, 60)])];
        let hit = check(&c, &ChainParams::default(), &b).0.unwrap().remove(0);
        assert_eq!(hit.chain, vec![(20, 20, 10), (40, 40, 10), (60, 60, 10)]);
        assert_eq!(hit.score, 90.);
    }

    #[test]
    fn batch_certificate_matches_scalar_certificate_at_numeric_boundaries() {
        let cases = [
            vec![0., 1., 2., -4., f32::NEG_INFINITY],
            vec![-0.],
            vec![f32::NAN],
            vec![f32::INFINITY],
            vec![0.25],
            vec![16_777_216.],
            vec![16_777_215.],
            vec![2_147_483_648.],
            vec![-2_147_483_648.],
            vec![f32::MIN_POSITIVE],
        ];
        for terms in [1., 2., 100., 1000.] {
            for delta in &cases {
                let mut scalar = Certificate { terms, max: 0. };
                let mut batch = scalar;
                let old = delta
                    .iter()
                    .filter(|&&e| e != f32::NEG_INFINITY)
                    .try_for_each(|&e| scalar.include(e));
                let new = batch.batch(delta);
                assert_eq!(old.is_some(), new.is_some());
                if old.is_some() {
                    assert_eq!(scalar.max, batch.max);
                }
            }
        }
    }

    #[test]
    fn randomized_planned_and_selective_fallback_match_all_forward_hits() {
        let mut rng = 0xabcdef0123456789u64;
        let mut next = || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            rng
        };
        let p = ChainParams::default();
        let mut accepted = 0;
        for _ in 0..80 {
            let n = 160 + (next() % 100) as usize;
            let ntx = 4 + (next() % 6) as usize;
            let (c, mut buckets) = shareable(n, ntx);
            for b in &mut buckets {
                if next() % 3 == 0 {
                    let k = 2 + (next() as usize % (n - 2));
                    b[k].ts += 1; // one transcript must not discard every certified result
                }
                if next() % 3 == 0 {
                    b.remove((next() as usize) % (n / 2));
                }
            }
            if planned_check(&c, &p, &buckets).0 {
                accepted += 1;
            }
        }
        assert!(accepted >= 60, "exercise planned sharing and selective fallback");
    }
}
