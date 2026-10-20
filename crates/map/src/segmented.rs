//! Fixed segmented-affine problems, packed deduplication, and candidate-level bounds.
use super::*;
use index::PackedRange;

impl<C: Chainer> Mapper<C> {
    fn raw_ranges<'a, 'r>(
        &'a self,
        tb: &TieBreak,
        raw: Option<&'r RawRead>,
        spans: &[(u32, u32, u32)],
    ) -> Option<(&'r RawRead, Vec<PackedRange<'a>>)> {
        let r = raw.filter(|_| tb.raw_ends)?;
        let texts = spans
            .iter()
            .map(|&(t, lo, hi)| self.idx.tx_range(t, lo, hi, true))
            .collect::<Option<Vec<_>>>()?;
        Some((r, texts))
    }

    pub(super) fn affine_plan<'a>(
        &'a self,
        tb: &TieBreak,
        codes: &[u8],
        raw: Option<&RawRead>,
        hits: &[TxHit],
    ) -> Option<Plan<'a>> {
        #[cfg(test)]
        if self.reference_tiebreak {
            return None;
        }
        let sc = tb.affine.as_ref()?;
        // Capped edit distances and the approximate prefilter are order-dependent.
        // Keep their original scorer, including its sorted-text order.
        if tb.affine_within.is_some() || [sc.mat, sc.mis, sc.open, sc.ext].iter().any(|&x| x < 0) {
            return None;
        }
        let blocks = shared_blocks(hits)?;
        let plen = codes.len() as u32;
        let mut plan = Plan::new(hits.len());
        // Strict interior gaps [previous anchor end, next anchor start).
        for pair in blocks.windows(2) {
            let (left, right) = (&pair[0], &pair[1]);
            let texts = hits
                .iter()
                .enumerate()
                .map(|(i, h)| {
                    self.idx
                        .tx_range(h.tx, left.targets[i].end, right.targets[i].start, false)
                        .unwrap()
                })
                .collect();
            plan.segment(
                tb,
                texts,
                &codes[left.query.end as usize..right.query.start as usize],
                align::Mode::Global,
                false,
                None,
            );
        }
        // Left end: read[0, q0) against transcript windows ending at each candidate's anchor.
        let first = &blocks[0];
        let q0 = first.query.start;
        if q0 > 0 {
            let spans: Vec<(u32, u32, u32)> = hits
                .iter()
                .zip(&first.targets)
                .map(|(h, t)| {
                    (
                        h.tx,
                        t.start.saturating_sub(q0 + tb.pad + q0 / tb.pad_div),
                        t.start,
                    )
                })
                .collect();
            let short: Vec<bool> = spans.iter().map(|&(_, lo, _)| lo == 0).collect();
            match self.raw_ranges(tb, raw, &spans) {
                Some((r, texts)) => plan.end(tb, texts, &r.codes[..r.pos(q0)], true, &short),
                None => plan.end(
                    tb,
                    spans
                        .iter()
                        .map(|&(t, lo, hi)| self.idx.tx_range(t, lo, hi, false).unwrap())
                        .collect(),
                    &codes[..q0 as usize],
                    true,
                    &short,
                ),
            };
        }
        // Right end: read[last block end, plen) against windows from each candidate's block end.
        // Blocks are exact in parsed space only: if the candidates' raw bases differ inside the
        // last block (homopolymer run lengths), the raw read end is aligned through it.
        let last = blocks.last().unwrap();
        let t0 = &last.targets[0];
        let through = raw.is_some()
            && tb.raw_ends
            && hits.iter().zip(&last.targets).skip(1).any(|(h, t)| {
                self.idx
                    .tx_raw_eq(hits[0].tx, t0.start, t0.end, h.tx, t.start, t.end)
                    == Some(false)
            });
        let ql = if through { last.query.start } else { last.query.end };
        let rest = plen - ql;
        let spans: Vec<(u32, u32, u32)> = hits
            .iter()
            .zip(&last.targets)
            .map(|(h, t)| {
                let tlen = self.idx.tx_plen()[h.tx as usize];
                let s = if through { t.start } else { t.end };
                (h.tx, s, (s + rest + tb.pad + rest / tb.pad_div).min(tlen))
            })
            .collect();
        // The read's poly(A) tail, searched for past the last block (raw positions from ql).
        let a_code = lcp_core::encode::CODE[b'A' as usize] & 3;
        // Whether each window's end (its length from the read end's start, in the tail's space)
        // is where the tail starts: within the tolerance of its start, or further but only by
        // annotated A's (an end that includes part of the tail) up to the run's end.
        let ends_at_tail = |texts: &[index::PackedRange<'_>], cut: usize, run_end: usize| -> Vec<bool> {
            let tol = TAIL_TOL + cut / 10;
            texts
                .iter()
                .map(|t| {
                    let n = t.len();
                    n.abs_diff(cut) <= tol
                        || (n > cut && n <= run_end + tol && t.count_other(cut, a_code) * 10 <= n - cut)
                })
                .collect()
        };
        let tail = match raw {
            Some(r) if tb.tail_trim || tb.tail_bonus > 0 => {
                let (r0, r1) = (r.pos(ql), r.pos(last.query.end));
                tail_start(&r.codes[r1..], a_code, tb.tail_run)
                    .map(|(c, e)| (r, r0, r1 - r0 + c, r1 - r0 + e))
            }
            _ => None,
        };
        let short: Vec<bool> = spans
            .iter()
            .map(|&(tx, _, hi)| hi == self.idx.tx_plen()[tx as usize])
            .collect();
        let ends = match self.raw_ranges(tb, raw, &spans) {
            Some((r, texts)) => {
                let rest = &r.codes[r.pos(ql)..];
                let cut = tail.map_or(rest.len(), |(_, _, c, _)| c);
                let run_end = tail.map_or(cut, |(_, _, _, e)| e);
                let at_end = tail.map(|_| ends_at_tail(&texts, cut, run_end));
                let read = if tb.tail_trim { &rest[..cut] } else { rest };
                plan.end(tb, texts, read, false, &short);
                at_end
            }
            None => {
                let texts: Vec<index::PackedRange<'_>> = spans
                    .iter()
                    .map(|&(t, lo, hi)| self.idx.tx_range(t, lo, hi, false).unwrap())
                    .collect();
                // The tail's start and run end in parsed space: the last parsed base starting at
                // or before each.
                let parsed = |x: usize, r: &RawRead| {
                    let (mut lo, mut hi) = (ql, plen);
                    while lo < hi {
                        let mid = (lo + hi).div_ceil(2);
                        if r.pos(mid) <= x { lo = mid } else { hi = mid - 1 }
                    }
                    (lo - ql) as usize
                };
                let cut = tail.map_or((plen - ql) as usize, |(r, r0, c, _)| parsed(r0 + c, r));
                let run_end = tail.map_or(cut, |(r, r0, _, e)| parsed(r0 + e, r));
                let at_end = tail.map(|_| ends_at_tail(&texts, cut, run_end));
                let read = &codes[ql as usize..];
                let read = if tb.tail_trim { &read[..cut] } else { read };
                plan.end(tb, texts, read, false, &short);
                at_end
            }
        };
        // Candidates whose 3' end is where the read's tail starts gain the bonus.
        if let Some(at_end) = ends
            && tb.tail_bonus > 0
        {
            for ((c, &(tx, _, hi)), &ok) in plan.cost.iter_mut().zip(&spans).zip(&at_end) {
                if ok && hi == self.idx.tx_plen()[tx as usize] {
                    *c -= i64::from(tb.tail_bonus);
                }
            }
        }
        // Stay on the original path when the DP's i32 arithmetic cannot be certified.
        plan.jobs
            .iter()
            .all(|j| j.safe_arithmetic(&plan.contexts[j.context], tb))
            .then_some(plan)
    }
}

/// Reporting uses these exact f32 operations, after tie normalization. Restrict this
/// certificate to a single orientation and a finite, nonnegative baseline/cutoff.
/// A completed candidate guarantees that the final best alignment is at least R;
/// hence B - lambda * max(R - U, 0) is an upper bound on the candidate's final score.
/// Converting the integer deficit and then multiplying/subtracting as the reporting
/// path does is monotone, including f32 rounding. Equality must be retained.
#[derive(Clone, Copy)]
pub(super) struct Retention {
    pub baseline: f32,
    pub delta: f32,
    pub lambda: f32,
}

impl Retention {
    fn rejects(self, reference: i64, upper: i64) -> bool {
        if !self.baseline.is_finite()
            || self.baseline < 0.0
            || !self.lambda.is_finite()
            || self.lambda < 0.0
            || !(0.0..=1.0).contains(&self.delta)
        {
            return false;
        }
        let Some(deficit) = reference.checked_sub(upper) else {
            return false;
        };
        let optimistic = self.baseline - self.lambda * deficit.max(0) as f32;
        let floor = self.baseline - self.delta * self.baseline.abs();
        optimistic.is_finite() && floor.is_finite() && optimistic < floor
    }
}

struct Context {
    read: Vec<u8>,
    mode: align::Mode,
    reverse: bool,
    cut: Option<usize>,
}

struct Job<'a> {
    text: PackedRange<'a>,
    context: usize,
    score: Option<i64>,
}

impl Job<'_> {
    fn safe_arithmetic(&self, ctx: &Context, tb: &TieBreak) -> bool {
        let sc = tb.affine.as_ref().unwrap();
        let magnitude = i64::from(sc.mat) + i64::from(sc.mis) + i64::from(sc.open) + i64::from(sc.ext);
        (ctx.read.len() as i64)
            .checked_add(self.text.len() as i64)
            .and_then(|n| n.checked_add(1))
            .and_then(|n| n.checked_mul(magnitude))
            .and_then(|n| n.checked_add(i64::from(tb.end_bonus).abs()))
            .is_some_and(|n| n < i64::from(i32::MAX / 4))
    }

    fn upper(&self, ctx: &Context, tb: &TieBreak) -> Option<i64> {
        // Nonpositive mismatch/gap contributions: each rewarded pair uses one base
        // of each sequence. Bands and X-drop only restrict the visited alignments.
        // Flank suppression chooses one of two such scores, so the same bound applies.
        // No mandatory length penalty: clipped ends may leave either suffix unaligned.
        let sc = tb.affine.as_ref()?;
        i64::from(sc.mat)
            .checked_mul(ctx.read.len().min(self.text.len()) as i64)?
            .checked_add(if ctx.mode == align::Mode::Prefix && tb.clip_ends {
                i64::from(tb.end_bonus.max(0))
            } else {
                0
            })
    }

    fn work(&self, ctx: &Context, tb: &TieBreak) -> u128 {
        let (m, n) = (ctx.read.len(), self.text.len());
        let geom = if ctx.mode == align::Mode::Prefix && tb.clip_ends {
            if tb.end_diagonal {
                m.min(n)
            } else {
                n.max(m + tb.pad as usize + m / tb.pad_div as usize)
            }
        } else {
            n
        };
        // The actual kernel's two band widths, with its window/diagonal geometry.
        // Saturation also handles the unbanded usize::MAX setting.
        let width = tb
            .band
            .saturating_add(m.saturating_sub(geom))
            .saturating_add(tb.band.saturating_add(geom.saturating_sub(m)))
            .saturating_add(1)
            .min(n.saturating_add(1));
        m as u128 * width as u128
    }
}

pub(super) struct Plan<'a> {
    contexts: Vec<Context>,
    jobs: Vec<Job<'a>>,
    required: Vec<Vec<usize>>,
    // Annotation contributions belong to candidates, never to a sequence cache entry.
    cost: Vec<i64>,
    #[cfg(test)]
    evaluated: Vec<usize>,
    #[cfg(test)]
    candidates: Vec<usize>,
    #[cfg(test)]
    decoded: usize,
}

impl<'a> Plan<'a> {
    fn new(n: usize) -> Self {
        Self {
            contexts: Vec::new(),
            jobs: Vec::new(),
            required: vec![Vec::new(); n],
            cost: vec![0; n],
            #[cfg(test)]
            evaluated: Vec::new(),
            #[cfg(test)]
            candidates: Vec::new(),
            #[cfg(test)]
            decoded: 0,
        }
    }

    fn segment(
        &mut self,
        tb: &TieBreak,
        texts: Vec<PackedRange<'a>>,
        read: &[u8],
        mode: align::Mode,
        reverse: bool,
        cut: Option<usize>,
    ) {
        // Decide uniform skips over ALL ORIGINAL candidates, before any filtering.
        // Every text in a segment has the same read, direction, and scoring context.
        if texts.windows(2).all(|w| w[0].sequence_eq(&w[1])) {
            return;
        }
        stats::add(&stats::TIEBREAK_SEGMENTS, 1);
        let context = self.contexts.len();
        self.contexts.push(Context {
            read: if reverse {
                read.iter().rev().copied().collect()
            } else {
                read.to_vec()
            },
            mode,
            reverse,
            cut,
        });
        let start = self.jobs.len();
        for (i, text) in texts.into_iter().enumerate() {
            // Context-local equality prevents reuse across different read intervals,
            // clipping, cuts, band geometry, or direction, even for identical bases.
            let j = (start..self.jobs.len())
                .find(|&j| text.sequence_eq(&self.jobs[j].text))
                .unwrap_or_else(|| {
                    let j = self.jobs.len();
                    self.jobs.push(Job {
                        text,
                        context,
                        score: None,
                    });
                    j
                });
            self.required[i].push(j);
        }
        // Complete exact simple cases immediately, without decoding reference bases.
        for j in start..self.jobs.len() {
            let n = self.jobs[j].text.len();
            if mode == align::Mode::Global && (read.is_empty() || n == 0) {
                let sc = tb.affine.as_ref().unwrap();
                self.jobs[j].score = Some(if n == 0 && read.is_empty() {
                    0
                } else {
                    -i64::from(sc.open) - i64::from(sc.ext) * n.max(read.len()) as i64
                });
            } else if read.is_empty() {
                self.jobs[j].score = Some(if tb.clip_ends {
                    i64::from(tb.end_bonus.max(0))
                } else {
                    0
                });
            } else if n == 0 && self.jobs[j].safe_arithmetic(&self.contexts[context], tb) {
                // Let the existing kernel handle window reachability and end bonuses;
                // there are no reference bases to decode or interior DP cells to visit.
                self.run(j, tb);
            }
        }
    }

    fn end(
        &mut self,
        tb: &TieBreak,
        texts: Vec<PackedRange<'a>>,
        read: &[u8],
        reverse: bool,
        short: &[bool],
    ) {
        let cut = texts
            .iter()
            .zip(short)
            .filter(|&(_, &s)| s)
            .map(|(t, _)| t.len())
            .min()
            .filter(|&c| tb.flank_min > 0 && texts.iter().any(|t| t.len() > c));
        self.segment(tb, texts, read, align::Mode::Prefix, reverse, cut);
    }

    fn run(&mut self, j: usize, tb: &TieBreak) -> i64 {
        self.run_batch(&[j], tb);
        self.jobs[j].score.unwrap()
    }

    /// Score the jobs `js` of one context that have no score yet, in one pass of the sorted
    /// kernels: texts sorted, so rows that read only a prefix shared with the previous text
    /// are not recomputed (equal to separate calls).
    fn run_batch(&mut self, js: &[usize], tb: &TieBreak) {
        let js: Vec<usize> = js
            .iter()
            .copied()
            .filter(|&j| self.jobs[j].score.is_none())
            .collect();
        let Some(&j0) = js.first() else {
            return;
        };
        let ctx = &self.contexts[self.jobs[j0].context];
        debug_assert!(js.iter().all(|&j| self.jobs[j].context == self.jobs[j0].context));
        let sc = tb.affine.as_ref().unwrap();
        let texts: Vec<Vec<u8>> = js
            .iter()
            .map(|&j| {
                let t = &self.jobs[j].text;
                let mut v = if t.is_empty() { Vec::new() } else { t.decode() };
                if ctx.reverse {
                    v.reverse();
                }
                v
            })
            .collect();
        #[cfg(test)]
        {
            self.decoded += js.iter().filter(|&&j| !self.jobs[j].text.is_empty()).count();
            self.evaluated.extend(&js);
        }
        let mut order: Vec<usize> = (0..js.len()).collect();
        order.sort_unstable_by(|&a, &b| texts[a].cmp(&texts[b]));
        let prefix_mode = ctx.mode == align::Mode::Prefix;
        let diag = prefix_mode && tb.clip_ends && tb.end_diagonal;
        if stats::on() {
            // As the original evaluator counts: text bases, those outside the common prefix
            // with the previous text in sorted order, and the banded cells per text.
            let m = ctx.read.len() as u64;
            let mut prev: Option<&Vec<u8>> = None;
            for &k in &order {
                let t = &texts[k];
                let l = prev.map_or(0, |p| t.iter().zip(p).take_while(|(a, b)| a == b).count());
                let n = t.len() as u64;
                let drift = if diag && n >= m { 0 } else { m.abs_diff(n) };
                let width = (tb.band as u64)
                    .saturating_mul(2)
                    .saturating_add(drift)
                    .saturating_add(1);
                let cells = m.saturating_mul(width.min(n + 1));
                stats::add(&stats::AFFINE_CALLS, 1);
                stats::add(&stats::AFFINE_TEXT_BASES, n);
                stats::add(&stats::AFFINE_TEXT_NEW_BASES, (t.len() - l) as u64);
                stats::add(&stats::AFFINE_CELLS, cells);
                if prefix_mode {
                    stats::add(&stats::AFFINE_END_CELLS, cells);
                }
                prev = Some(t);
            }
        }
        let tx: Vec<&[u8]> = order.iter().map(|&k| texts[k].as_slice()).collect();
        let w0 = align::dp_work();
        let xs: Vec<(i32, i32, usize)> = if prefix_mode && tb.clip_ends {
            let window = ctx.read.len() + tb.pad as usize + ctx.read.len() / tb.pad_div as usize;
            align::affine_extend_sorted_cut(
                &ctx.read,
                &tx,
                sc,
                tb.band,
                tb.end_bonus,
                tb.xdrop,
                diag,
                window,
                ctx.cut.unwrap_or(usize::MAX),
            )
        } else {
            align::affine_sorted(&ctx.read, &tx, ctx.mode, sc, tb.band)
                .into_iter()
                .map(|x| (x, x, 0))
                .collect()
        };
        if stats::on() {
            let w1 = align::dp_work();
            let (c, r) = if prefix_mode {
                (&stats::DP_TBEND_CELLS, &stats::DP_TBEND_ROWS)
            } else {
                (&stats::DP_SEG_CELLS, &stats::DP_SEG_ROWS)
            };
            stats::add(c, w1.0 - w0.0);
            stats::add(r, w1.1 - w0.1);
        }
        let cut = ctx.cut;
        for (&k, &(x, near, cols)) in order.iter().zip(&xs) {
            let mut score = i64::from(x);
            if let Some(cut) = cut {
                let gain = score - i64::from(near);
                let past = cols.saturating_sub(cut) as i64;
                if texts[k].len() > cut
                    && gain < tb.flank_min
                    && !(past > 0 && 4 * gain >= 3 * i64::from(sc.mat) * past)
                {
                    score = i64::from(near);
                }
            }
            self.jobs[js[k]].score = Some(score);
        }
    }

    /// Candidate `i`'s score, once all its jobs are scored.
    fn total(&self, i: usize) -> Option<i64> {
        self.required[i]
            .iter()
            .try_fold(-self.cost[i], |sum, &j| Some(sum + self.jobs[j].score?))
    }

    /// An upper bound on candidate `i`'s score: scored jobs as they are, the rest at their
    /// bound. None (no pruning) if a sum would overflow.
    fn upper_bound(&self, i: usize, tb: &TieBreak) -> Option<i64> {
        self.required[i].iter().try_fold(-self.cost[i], |sum, &j| {
            let job = &self.jobs[j];
            sum.checked_add(job.score.or_else(|| job.upper(&self.contexts[job.context], tb))?)
        })
    }

    /// Costs (negated scores) of every candidate, or None for one certified (by `filter`) to
    /// fall below the reporting floor. Every candidate has one job per context (uniform
    /// segments have none); contexts are scored cheapest first, each context's texts still
    /// needed in one sorted batch, so isoforms that agree near the anchor share DP rows. With a
    /// filter, the candidate with the best chain score is aligned first (one text per
    /// context): its score is the reference that certifies the others before each batch.
    pub(super) fn evaluate(
        &mut self,
        tb: &TieBreak,
        hits: &[TxHit],
        filter: Option<Retention>,
    ) -> Vec<Option<i64>> {
        // Input order is already read-salted. Strict improvement preserves that order
        // among equal original chain scores; scores are never changed during evaluation.
        let n = hits.len();
        let first = (0..n)
            .reduce(|a, b| if hits[b].score > hits[a].score { b } else { a })
            .unwrap();
        #[cfg(test)]
        self.candidates.push(first);
        let work = |p: &Self, j: usize| {
            let job = &p.jobs[j];
            if job.score.is_some() {
                0
            } else {
                job.work(&p.contexts[job.context], tb)
            }
        };
        let mut reference = None;
        if filter.is_some() {
            let mut own: Vec<usize> = self.required[first].clone();
            own.sort_by_key(|&j| (work(self, j), j));
            for j in own {
                self.run(j, tb);
            }
            reference = self.total(first);
        }
        let mut alive = vec![true; n];
        let mut order: Vec<usize> = (0..self.contexts.len()).collect();
        let ctx_work: Vec<u128> = order
            .iter()
            .map(|&c| (0..n).map(|i| self.required[i][c]).map(|j| work(self, j)).sum())
            .collect();
        order.sort_by_key(|&c| (ctx_work[c], c));
        for c in order {
            if let (Some(f), Some(r)) = (filter, reference) {
                for i in 0..n {
                    if alive[i] && self.upper_bound(i, tb).is_some_and(|u| f.rejects(r, u)) {
                        alive[i] = false;
                    }
                }
            }
            let mut js: Vec<usize> = (0..n)
                .filter(|&i| alive[i])
                .map(|i| self.required[i][c])
                .collect();
            js.sort_unstable();
            js.dedup();
            self.run_batch(&js, tb);
            if filter.is_some() {
                for i in (0..n).filter(|&i| alive[i]) {
                    if let Some(s) = self.total(i) {
                        reference = Some(reference.map_or(s, |r: i64| r.max(s)));
                    }
                }
            }
        }
        (0..n)
            .map(|i| alive[i].then(|| -self.total(i).unwrap()))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::segmented_tests::{codes, hit, mapper};
    use proptest::prelude::*;

    fn ranges(m: &Mapper<OracleChainer>, lo: u32, hi: u32) -> Vec<PackedRange<'_>> {
        (0..m.idx.n_transcripts() as u32)
            .map(|t| m.idx.tx_range(t, lo, hi, false).unwrap())
            .collect()
    }

    fn filter(baseline: f32, delta: f32) -> Option<Retention> {
        Some(Retention {
            baseline,
            delta,
            lambda: 1.0,
        })
    }

    #[test]
    fn cheap_first_prunes_unneeded_work_but_preserves_shared_jobs() {
        let a = vec![b'A'; 26];
        let mut b = a.clone();
        b[19] = b'C';
        b[20..].fill(b'T');
        let mut c = a.clone();
        c[19] = b'C';
        let mut d = a.clone();
        d[18..20].fill(b'C');
        d[20..].fill(b'G');
        let m = mapper(&[&a, &b, &c, &d], false);
        let tb = TieBreak::default();
        let mut p = Plan::new(4);
        // Deliberately insert the expensive segment first.
        p.segment(
            &tb,
            ranges(&m, 0, 20),
            &codes(&a[..20]),
            align::Mode::Global,
            false,
            None,
        );
        p.segment(
            &tb,
            ranges(&m, 20, 26),
            &codes(&a[20..]),
            align::Mode::Global,
            false,
            None,
        );
        let required = p.required.clone();
        let hits: Vec<_> = (0..4).map(|t| hit(t, 100. - t as f32, &[(0, 0, 1)])).collect();
        let got = p.evaluate(&tb, &hits, filter(100., 0.1));
        assert_eq!(got, [Some(-52), None, Some(-46), None]);
        // The reference candidate's jobs first, cheap before expensive; then the cheap context's
        // remaining texts in one batch; then the expensive context's, for survivors only.
        assert_eq!(p.evaluated[..2], [required[0][1], required[0][0]]);
        let mut rest = p.evaluated[2..].to_vec();
        rest.sort_unstable();
        let mut want = vec![required[1][1], required[2][0], required[3][1]];
        want.sort_unstable();
        assert_eq!(rest, want);
        // Candidate 1's skipped expensive job was still needed by surviving candidate 2.
        assert_eq!(required[1][0], required[2][0]);
        assert!(p.jobs[required[3][0]].score.is_none());
        assert_eq!(p.decoded, 5);
        assert_eq!(p.required, required); // original segmentation is immutable
        assert_eq!(p.jobs.len(), 6);
    }

    #[test]
    fn poor_early_score_can_compensate_and_improve_reference() {
        let a = [vec![b'C'; 20], vec![b'A'; 6]].concat();
        let b = [vec![b'A'; 20], vec![b'T'; 6]].concat();
        let c = [vec![b'G'; 19], vec![b'T'; 6]].concat();
        let m = mapper(&[&a, &b, &c], false);
        let tb = TieBreak::default();
        let mut p = Plan::new(3);
        p.segment(
            &tb,
            vec![
                m.idx.tx_range(0, 0, 20, false).unwrap(),
                m.idx.tx_range(1, 0, 20, false).unwrap(),
                m.idx.tx_range(2, 0, 19, false).unwrap(),
            ],
            &codes(&[b'A'; 20]),
            align::Mode::Global,
            false,
            None,
        );
        p.segment(
            &tb,
            vec![
                m.idx.tx_range(0, 20, 26, false).unwrap(),
                m.idx.tx_range(1, 20, 26, false).unwrap(),
                m.idx.tx_range(2, 19, 25, false).unwrap(),
            ],
            &codes(&[b'A'; 6]),
            align::Mode::Global,
            false,
            None,
        );
        let hits: Vec<_> = (0..3).map(|t| hit(t, 100. - t as f32, &[(0, 0, 1)])).collect();
        let got = p.evaluate(&tb, &hits, filter(100., 0.));
        // Candidate 1 improves R from -68 to 16, after which candidate 2's bound (14) would be
        // rejected; its texts were already scored in the same batches as candidate 1's, so it
        // is kept with its exact cost (below the floor; the retain step drops it).
        assert_eq!(got, [Some(68), Some(-16), Some(106)]);
    }

    #[test]
    fn original_chain_score_then_existing_salted_order() {
        let m = mapper(&[b"AAAA", b"CCCC", b"GGGG"], false);
        let tb = TieBreak::default();
        let mut hits = vec![
            hit(0, 9., &[(0, 0, 1)]),
            hit(1, 11., &[(0, 0, 1)]),
            hit(2, 10., &[(0, 0, 1)]),
        ];
        let mut p = Plan::new(3);
        p.segment(
            &tb,
            ranges(&m, 0, 4),
            &codes(b"AAAA"),
            align::Mode::Global,
            false,
            None,
        );
        p.evaluate(&tb, &hits, None);
        assert_eq!(p.candidates, [1]);
        for h in &mut hits {
            h.score = 10.;
        }
        let mut firsts = std::collections::BTreeSet::new();
        for salt in 0..32 {
            sort_hits(&mut hits, salt);
            let mut p = Plan::new(3);
            p.evaluate(&tb, &hits, None);
            assert_eq!(p.candidates[0], 0);
            firsts.insert(hits[0].tx);
        }
        assert!(firsts.len() > 1); // no transcript-ID preference
    }

    #[test]
    fn inclusive_rounding_and_unavailable_certificates() {
        let f = filter(100., 0.1).unwrap();
        assert!(!f.rejects(20, 10));
        assert!(f.rejects(20, 9));
        assert!(!f.rejects(i64::MAX, i64::MIN));
        assert!(
            !Retention {
                lambda: f32::MAX,
                ..f
            }
            .rejects(100, 0)
        );
        for f in [
            Retention { lambda: -1., ..f },
            Retention { baseline: -1., ..f },
            Retention { delta: 2., ..f },
            Retention { delta: f32::NAN, ..f },
            Retention {
                lambda: f32::INFINITY,
                ..f
            },
        ] {
            assert!(!f.rejects(100, 0));
        }
        let f = filter(16_777_216., 0.).unwrap();
        // Deficits below half an ulp may round back to the baseline.
        assert_eq!(f.rejects(1, 0), (f.baseline - 1.) < f.baseline);
    }

    #[test]
    fn clipped_ends_and_annotation_bonuses_are_in_bounds() {
        let m = mapper(&[b"AAAACCCC", b"AAAA"], false);
        let tb = TieBreak {
            end_bonus: 7,
            pad: 0,
            ..Default::default()
        };
        let mut p = Plan::new(2);
        p.segment(
            &tb,
            vec![
                m.idx.tx_range(0, 0, 8, false).unwrap(),
                m.idx.tx_range(1, 0, 4, false).unwrap(),
            ],
            &codes(b"AAAA"),
            align::Mode::Prefix,
            false,
            None,
        );
        // Length differences are free at clipped ends; both earn the complete-read bonus.
        assert_eq!(p.jobs[1].upper(&p.contexts[0], &tb), Some(15));
        p.cost[1] = -10; // transcript-specific annotation bonus, outside the shared job
        let hits = [hit(0, 100., &[(0, 0, 1)]), hit(1, 99., &[(0, 0, 1)])];
        assert_eq!(p.evaluate(&tb, &hits, filter(100., 0.)), [Some(-15), Some(-25)]);
    }

    #[test]
    fn packed_dedup_skips_uniform_segments_and_keeps_contexts_separate() {
        let m = mapper(&[b"ACGTACGT", b"TTACGTAA", b"ACGTCCCC"], false);
        let tb = TieBreak::default();
        let same = vec![
            m.idx.tx_range(0, 0, 4, false).unwrap(),
            m.idx.tx_range(1, 2, 6, false).unwrap(),
            m.idx.tx_range(2, 0, 4, false).unwrap(),
        ];
        let mut p = Plan::new(3);
        p.segment(&tb, same, &codes(b"ACGT"), align::Mode::Global, false, None);
        assert!(p.jobs.is_empty());
        assert_eq!(p.decoded, 0);
        let texts = vec![
            m.idx.tx_range(0, 0, 4, false).unwrap(),
            m.idx.tx_range(1, 2, 6, false).unwrap(),
            m.idx.tx_range(2, 4, 8, false).unwrap(),
        ];
        p.segment(
            &tb,
            texts.clone(),
            &codes(b"ACGT"),
            align::Mode::Global,
            false,
            None,
        );
        p.segment(
            &tb,
            texts.clone(),
            &codes(b"ACGT"),
            align::Mode::Prefix,
            false,
            None,
        );
        p.segment(&tb, texts, &codes(b"ACGT"), align::Mode::Prefix, true, Some(2));
        assert_eq!(p.jobs.len(), 6); // two sequences in each of three contexts
        let hits: Vec<_> = (0..3).map(|t| hit(t, 100., &[(0, 0, 1)])).collect();
        p.evaluate(&tb, &hits, None);
        assert_eq!(p.decoded, 6);
        assert_ne!(
            p.jobs[0].work(&p.contexts[0], &tb),
            p.jobs[0].work(
                &p.contexts[0],
                &TieBreak {
                    band: 0,
                    ..tb.clone()
                }
            )
        );
    }

    #[test]
    fn unsupported_scoring_uses_original_path() {
        let m = mapper(&[b"ACGTAC", b"ACCCAC"], false);
        let hits = [
            hit(0, 100., &[(0, 0, 2), (4, 4, 2)]),
            hit(1, 99., &[(0, 0, 2), (4, 4, 2)]),
        ];
        for tb in [
            TieBreak {
                affine: None,
                ..Default::default()
            },
            TieBreak {
                affine_within: Some(0),
                ..Default::default()
            },
            TieBreak {
                affine: Some(align::Affine {
                    mat: -1,
                    ..align::Affine::ONT
                }),
                ..Default::default()
            },
        ] {
            assert!(m.affine_plan(&tb, &codes(b"ACGTAC"), None, &hits).is_none());
            assert_eq!(
                m.segmented(&tb, &codes(b"ACGTAC"), None, &hits),
                m.segmented_reference(&tb, &codes(b"ACGTAC"), None, &hits)
            );
        }
        let tb = TieBreak {
            affine: Some(align::Affine {
                mat: i32::MAX,
                ..align::Affine::ONT
            }),
            ..Default::default()
        };
        assert!(m.affine_plan(&tb, &codes(b"ACGTAC"), None, &hits).is_none());
    }

    #[test]
    fn rejection_does_not_reselect_shared_anchors_or_uniform_skips() {
        let q = [b"G".as_slice(), &[b'A'; 6], b"C", &[b'A'; 20], b"T"].concat();
        let mut bad = q.clone();
        bad[1..7].fill(b'T');
        let mut other = q.clone();
        other[27] = b'C';
        let m = mapper(&[&q, &bad, &other], false);
        let detailed = [(0, 0, 1), (7, 7, 1), (18, 18, 1), (28, 28, 1)];
        let hits = [
            hit(0, 100., &detailed),
            hit(1, 99., &[(0, 0, 1), (7, 7, 1), (28, 28, 1)]),
            hit(2, 98., &detailed),
        ];
        let tb = TieBreak::default();
        let mut p = m.affine_plan(&tb, &codes(&q), None, &hits).unwrap();
        assert_eq!(shared_blocks(&hits).unwrap().len(), 3);
        let got = p.evaluate(&tb, &hits, filter(100., 0.1));
        assert_eq!(got, [Some(-52), None, Some(-46)]);
        assert_eq!(
            p.contexts.iter().map(|c| c.read.len()).collect::<Vec<_>>(),
            [6, 20]
        );
        let survivors = [hits[0].clone(), hits[2].clone()];
        assert_eq!(shared_blocks(&survivors).unwrap().len(), 4);
        // Replanning would now skip the cheap gap and split the expensive gap. The
        // evaluated plan retained both original contexts and the original totals.
        let replanned = m.affine_plan(&tb, &codes(&q), None, &survivors).unwrap();
        assert_eq!(
            replanned
                .contexts
                .iter()
                .map(|c| c.read.len())
                .collect::<Vec<_>>(),
            [9]
        );
        let reference = m.segmented_reference(&tb, &codes(&q), None, &hits).unwrap().0;
        assert_eq!((got[0], got[2]), (Some(reference[0]), Some(reference[2])));
    }

    #[test]
    fn exact_boundary_is_retained_and_simple_jobs_do_not_decode() {
        let m = mapper(&[b"AAAA", b"AAAC"], false);
        let tb = TieBreak::default();
        let hits = [hit(0, 100., &[(0, 0, 1)]), hit(1, 99., &[(0, 0, 1)])];
        let mut p = Plan::new(2);
        p.segment(
            &tb,
            ranges(&m, 0, 4),
            &codes(b"AAAA"),
            align::Mode::Global,
            false,
            None,
        );
        assert_eq!(p.evaluate(&tb, &hits, filter(100., 0.06)), [Some(-8), Some(-2)]);
        let mut p = Plan::new(2);
        p.segment(
            &tb,
            vec![
                m.idx.tx_range(0, 0, 0, false).unwrap(),
                m.idx.tx_range(1, 0, 4, false).unwrap(),
            ],
            &[],
            align::Mode::Global,
            false,
            None,
        );
        assert_eq!(p.evaluate(&tb, &hits, None), [Some(0), Some(12)]);
        assert_eq!(p.decoded, 0);
        assert!(p.evaluated.is_empty());
    }

    #[test]
    fn equal_hpc_ranges_can_have_different_raw_bases() {
        let m = mapper(&[b"AACCCGGTTTAA", b"ACGTA"], true);
        let (a, b) = (
            m.idx.tx_range(0, 0, 5, false).unwrap(),
            m.idx.tx_range(1, 0, 5, false).unwrap(),
        );
        assert!(a.sequence_eq(&b));
        let (a, b) = (
            m.idx.tx_range(0, 0, 5, true).unwrap(),
            m.idx.tx_range(1, 0, 5, true).unwrap(),
        );
        assert!(!a.sequence_eq(&b));
        let q = b"AACCCGGTTTAA";
        let (seq, raw) = Mapper::<OracleChainer>::parse(&m.seeder, q);
        let hits = [hit(0, 100., &[(1, 1, 2)]), hit(1, 99., &[(1, 1, 2)])];
        let tb = TieBreak::default();
        assert_eq!(
            m.segmented(&tb, &codes(&seq), raw.as_ref(), &hits),
            m.segmented_reference(&tb, &codes(&seq), raw.as_ref(), &hits)
        );
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(96))]
        #[test]
        fn randomized_segmented_scores_and_reported_metadata(
            gaps in prop::collection::vec(prop::collection::vec(0u8..4, 0..45), 6),
            read_mid in prop::collection::vec(0u8..4, 0..40),
            read_end in prop::collection::vec(0u8..4, 0..30),
            band in prop_oneof![Just(usize::MAX), 0usize..24],
            clip in any::<bool>(), diagonal in any::<bool>(), flank in 0i64..20,
            bonus in -2i32..15, xdrop in prop_oneof![Just(i32::MAX), 0i32..60],
            delta in 0u32..101, lambda in 0u32..5, salt in any::<u64>(),
            scoring in (0i32..4, 0i32..8, 0i32..8, 0i32..5),
        ) {
            let dna = |s: &[u8]| s.iter().map(|&c| b"ACGT"[c as usize]).collect::<Vec<_>>();
            let left = b"GCTAC"; let right = b"GATCG";
            let texts: Vec<Vec<u8>> = (0..3).map(|i| [left.as_slice(), &dna(&gaps[2*i]), right.as_slice(), &dna(&gaps[2*i+1])].concat()).collect();
            let mut m = mapper(&texts.iter().map(Vec::as_slice).collect::<Vec<_>>(), false);
            m.p.delta_rep = delta as f32 / 100.;
            let tb = TieBreak { band, clip_ends: clip, end_diagonal: diagonal, flank_min: flank,
                end_bonus: bonus, xdrop, lambda: lambda as f32 / 2.,
                affine: Some(align::Affine { mat: scoring.0, mis: scoring.1, open: scoring.2, ext: scoring.3 }),
                ..Default::default() };
            let q = [left.as_slice(), &dna(&read_mid), right.as_slice(), &dna(&read_end)].concat();
            let mut hits: Vec<_> = (0..3).map(|i| hit(i as u32, 100. - i as f32,
                &[(0,0,5), (5+read_mid.len() as u32, 5+gaps[2*i].len() as u32,5)])).collect();
            sort_hits(&mut hits, salt);
            let mut plan = m.affine_plan(&tb, &codes(&q), None, &hits).unwrap();
            plan.evaluate(&tb, &hits, None);
            for job in &plan.jobs {
                prop_assert!(job.score.unwrap() <= job.upper(&plan.contexts[job.context], &tb).unwrap());
            }
            prop_assert_eq!(m.segmented(&tb, &codes(&q), None, &hits), m.segmented_reference(&tb, &codes(&q), None, &hits));
            let mut optimized = hits.clone();
            let n = m.tiebreak(&tb, &q, None, &mut optimized, salt, true);
            m.reference_tiebreak = true;
            let nr = m.tiebreak(&tb, &q, None, &mut hits, salt, true);
            let retain = |h: &mut Vec<TxHit>| { let floor = h[0].score - m.p.delta_rep * h[0].score.abs(); h.retain(|x| x.score >= floor); };
            retain(&mut optimized); retain(&mut hits);
            prop_assert_eq!(n, nr);
            prop_assert_eq!(optimized, hits);
        }

        #[test]
        fn randomized_full_mapping_matches_reference(
            body in prop::collection::vec(0u8..4, 80..180),
            variant in prop::collection::vec(0usize..180, 0..12),
            hpc in any::<bool>(), both in any::<bool>(), tail in any::<bool>(),
            flank in 0i64..20, diagonal in any::<bool>(), delta in 0u32..101,
        ) {
            let dna: Vec<_> = body.iter().map(|&c| b"ACGT"[c as usize]).collect();
            let mut other = dna.clone();
            for x in variant { let i = x % other.len(); other[i] = b"ACGT"[((body[i] + 1) % 4) as usize]; }
            let shorter = &dna[..dna.len()-12];
            let mut m = mapper(&[&dna, &other, shorter], hpc);
            m.p.strands = if both { Strands::Both } else { Strands::Forward };
            m.p.delta_rep = delta as f32 / 100.;
            m.p.tiebreak.as_mut().unwrap().flank_min = flank;
            m.p.tiebreak.as_mut().unwrap().end_diagonal = diagonal;
            let q = if tail { [shorter, b"AAAAAAAAAAAAAAAAAAAAAGTCG"].concat() } else { dna };
            let got = m.map(&q);
            m.reference_tiebreak = true;
            prop_assert_eq!(got, m.map(&q));
        }
    }
}
