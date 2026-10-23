//! Transcript quantification by EM over per-read compatibility records (`spearfish em`, a
//! baseline for method development; formerly the `quant` crate).
//!
//! A read compatible with transcripts `j` with scores `s_j` contributes the conditional
//! likelihood `P(r | j) ∝ exp((s_j - s_max) / D)` (optionally divided by an effective
//! length), after dropping transcripts scoring below `min_frac × s_max`. This mirrors
//! oarfish's score-to-probability model, so a comparison isolates the mapper. Reads with the
//! same transcripts and (quantized) weights form one equivalence class.

use rayon::prelude::*;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct QuantParams {
    /// Score denominator D in exp((s - s_max) / D). The default, 2, is `spearfish em`'s.
    /// oarfish weighs spearfish's hits with D = 1 (`spearfish::SCORE_DENOM`, in spearfish's
    /// units) and minimap2 alignments with its own default.
    pub score_denom: f64,
    /// Keep transcripts scoring at least this fraction of the read's best score.
    pub min_frac: f64,
    /// Divide by transcript length (false for full-length long reads).
    pub length_norm: bool,
    pub max_iter: usize,
    /// Stop when the largest relative change among transcripts with at least 0.01 reads is
    /// below this.
    pub tol: f64,
}

impl Default for QuantParams {
    fn default() -> Self {
        QuantParams {
            score_denom: 2.0,
            min_frac: 0.95,
            length_norm: false,
            max_iter: 10_000,
            tol: 1e-6,
        }
    }
}

/// Reads sharing transcripts and weights.
#[derive(Clone, Debug, PartialEq)]
pub struct EqClass {
    pub txs: Vec<u32>,
    pub weights: Vec<f64>,
    pub count: f64,
}

/// Class key: transcripts and quantized weights.
type ClassKey = (Vec<u32>, Vec<u32>);

/// Builds equivalence classes from reads.
#[derive(Default)]
pub struct EqBuilder {
    /// Key -> (weights, read count).
    classes: HashMap<ClassKey, (Vec<f64>, f64)>,
    pub reads: u64,
}

impl EqBuilder {
    /// Add one read's (transcript, score) hits.
    pub fn add(&mut self, hits: &[(u32, f64)], lens: &[u32], p: &QuantParams) {
        let Some(best) = hits.iter().map(|h| h.1).reduce(f64::max) else {
            return;
        };
        let floor = if best > 0.0 { p.min_frac * best } else { best };
        let mut kept: Vec<(u32, f64)> = hits
            .iter()
            .filter(|h| h.1 >= floor)
            .map(|&(t, s)| {
                let w = ((s - best) / p.score_denom).exp();
                (
                    t,
                    if p.length_norm {
                        w / lens[t as usize].max(1) as f64
                    } else {
                        w
                    },
                )
            })
            .collect();
        // A transcript hit on both read orientations keeps its larger weight.
        kept.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(b.1.total_cmp(&a.1)));
        kept.dedup_by_key(|h| h.0);
        let max = kept.iter().map(|h| h.1).fold(0.0, f64::max);
        let txs: Vec<u32> = kept.iter().map(|h| h.0).collect();
        // Weights relative to the largest, quantized to 1e-4 for class identity.
        let q: Vec<u32> = kept.iter().map(|h| (h.1 / max * 1e4).round() as u32).collect();
        let e = self
            .classes
            .entry((txs, q))
            .or_insert_with(|| (kept.iter().map(|h| h.1 / max).collect(), 0.0));
        e.1 += 1.0;
        self.reads += 1;
    }

    pub fn finish(self) -> Vec<EqClass> {
        let mut v: Vec<EqClass> = self
            .classes
            .into_iter()
            .map(|((txs, _), (weights, count))| EqClass { txs, weights, count })
            .collect();
        // Weights compared by total_cmp: the order partial_cmp gives for real runs' weights
        // (finite, non-negative), without a panic on NaN (an infinite score, or D = 0).
        v.sort_unstable_by(|a, b| {
            a.txs.cmp(&b.txs).then_with(|| {
                let mut w = a.weights.iter().zip(&b.weights).map(|(x, y)| x.total_cmp(y));
                w.find(|o| o.is_ne()).unwrap_or(std::cmp::Ordering::Equal)
            })
        });
        v
    }
}

/// Estimated read counts per transcript. Returns (counts, iterations).
///
/// Converges when the largest relative change among transcripts with at least `0.01`
/// estimated reads drops below `tol`. Classes are flattened once into CSR arrays (each class's
/// transcripts and weights), plus the transpose (each transcript's classes and weights). Each
/// iteration is two passes on the current rayon pool, both without synchronization: per class,
/// the scale `count / Σ α_t w_t`; per transcript, `α_t Σ_c w_ct scale_c`. Sums run in class
/// order, so the result doesn't depend on the thread count.
pub fn em(n_tx: usize, classes: &[EqClass], p: &QuantParams) -> (Vec<f64>, usize) {
    let total: f64 = classes.iter().map(|c| c.count).sum();
    if total == 0.0 {
        return (vec![0.0; n_tx], 0);
    }
    // Classes → (transcript, weight).
    let mut off = Vec::with_capacity(classes.len() + 1);
    off.push(0usize);
    let (mut txs, mut ws) = (Vec::new(), Vec::new());
    for c in classes {
        txs.extend(c.txs.iter().copied());
        ws.extend_from_slice(&c.weights);
        off.push(txs.len());
    }
    let counts: Vec<f64> = classes.iter().map(|c| c.count).collect();
    // Transcripts → (class, weight), classes ascending.
    let mut toff = vec![0usize; n_tx + 1];
    for &t in &txs {
        toff[t as usize + 1] += 1;
    }
    for t in 0..n_tx {
        toff[t + 1] += toff[t];
    }
    let mut fill = toff.clone();
    let mut tcls = vec![0u32; txs.len()];
    let mut tws = vec![0f64; txs.len()];
    for c in 0..classes.len() {
        for k in off[c]..off[c + 1] {
            let t = txs[k] as usize;
            tcls[fill[t]] = c as u32;
            tws[fill[t]] = ws[k];
            fill[t] += 1;
        }
    }
    // Start from uniform shares within each class.
    let mut alpha = vec![0.0; n_tx];
    for (c, w) in off.windows(2).zip(&counts) {
        for &t in &txs[c[0]..c[1]] {
            alpha[t as usize] += w / (c[1] - c[0]) as f64;
        }
    }
    let mut next = vec![0.0; n_tx];
    let mut scale = vec![0.0; classes.len()];
    let mut iters = 0;
    while iters < p.max_iter {
        iters += 1;
        scale.par_iter_mut().enumerate().for_each(|(c, s)| {
            let r = off[c]..off[c + 1];
            let denom: f64 = txs[r.clone()]
                .iter()
                .zip(&ws[r])
                .map(|(&t, &w)| alpha[t as usize] * w)
                .sum();
            *s = if denom > 0.0 { counts[c] / denom } else { 0.0 };
        });
        let change = next
            .par_iter_mut()
            .enumerate()
            .map(|(t, n)| {
                let r = toff[t]..toff[t + 1];
                let acc: f64 = tcls[r.clone()]
                    .iter()
                    .zip(&tws[r])
                    .map(|(&c, &w)| w * scale[c as usize])
                    .sum();
                *n = alpha[t] * acc;
                if *n >= 0.01 {
                    (*n - alpha[t]).abs() / *n
                } else {
                    0.0
                }
            })
            .reduce(|| 0.0, f64::max);
        std::mem::swap(&mut alpha, &mut next);
        if change < p.tol {
            break;
        }
    }
    (alpha, iters)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conserves_reads_and_fixes_unique_only_transcripts() {
        // Random classes over transcripts 0..20; transcripts 20..23 only ever get unique reads.
        let p = QuantParams::default();
        let lens = vec![1000; 23];
        let mut b = EqBuilder::default();
        let mut x = 7u64;
        let mut next = |n: u64| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x % n
        };
        for _ in 0..5000 {
            let k = 1 + next(4) as usize;
            let hits: Vec<(u32, f64)> = (0..k)
                .map(|_| (next(20) as u32, 90.0 + next(10) as f64))
                .collect();
            b.add(&hits, &lens, &p);
        }
        for (t, n) in [(20u32, 37), (21, 5), (22, 1)] {
            for _ in 0..n {
                b.add(&[(t, 100.0)], &lens, &p);
            }
        }
        let (c, _) = em(23, &b.finish(), &p);
        assert!((c.iter().sum::<f64>() - 5043.0).abs() < 1e-6);
        assert!(c.iter().all(|&v| v >= 0.0));
        assert!((c[20] - 37.0).abs() < 1e-9 && (c[21] - 5.0).abs() < 1e-9 && (c[22] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn exact_recovery_when_identifiable() {
        // Each transcript emits reads uniquely with probability 1/2 and ambiguously (with
        // its neighbour) otherwise; with equal weights and proportional ambiguity, the MLE is
        // the truth.
        let p = QuantParams::default();
        let lens = [1, 1];
        let mut b = EqBuilder::default();
        for _ in 0..300 {
            b.add(&[(0, 50.0)], &lens, &p);
        }
        for _ in 0..100 {
            b.add(&[(1, 50.0)], &lens, &p);
        }
        for _ in 0..400 {
            b.add(&[(0, 50.0), (1, 50.0)], &lens, &p);
        }
        let (c, _) = em(2, &b.finish(), &p);
        assert!(
            (c[0] - 600.0).abs() < 1e-3 && (c[1] - 200.0).abs() < 1e-3,
            "{c:?}"
        );
    }

    #[test]
    fn score_weights_and_floor() {
        let p = QuantParams::default();
        let lens = [1, 1, 1];
        let mut b = EqBuilder::default();
        // Transcript 2 scores below 95% of the best and is dropped; 1 is down-weighted.
        b.add(&[(0, 100.0), (1, 97.0), (2, 90.0)], &lens, &p);
        let c = b.finish();
        assert_eq!(c[0].txs, vec![0, 1]);
        assert!((c[0].weights[1] - (-3.0f64 / 2.0).exp()).abs() < 1e-12);
    }

    #[test]
    fn infinite_scores_do_not_panic() {
        let p = QuantParams::default();
        let lens = [1, 1];
        let mut b = EqBuilder::default();
        b.add(&[(0, f64::INFINITY), (1, 1.0)], &lens, &p);
        b.add(&[(0, 2.0), (1, 1.0)], &lens, &p);
        b.add(&[(0, 2.0), (1, 2.0)], &lens, &p);
        assert_eq!(b.finish().len(), 3);
    }
}
