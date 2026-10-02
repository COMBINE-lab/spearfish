//! k-mer seed schemes for the seed-scheme ablations: minimizers, syncmers (closed/open, with
//! optional hash subsampling to match a target density), and order-2 randstrobes on syncmer
//! positions (as in strobealign). k <= 32.

use std::collections::VecDeque;

use lcp_core::encode::acgt_segments;
use lcp_core::{Core, Hierarchy, Level};

use crate::SeedScheme;

/// Invertible 64-bit mix (murmur3 fmix64): distinct k-mers get distinct hashes.
#[inline]
pub(crate) fn mix64(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    x ^= x >> 33;
    x
}

const SEED: u64 = 0x5eed_5eed_5eed_5eed;

/// Hashes of all k-mers of `codes` (2-bit), in order: `out[i]` is the k-mer at `i`.
fn kmer_hashes(codes: &[u8], k: usize, salt: u64, out: &mut Vec<u64>) {
    out.clear();
    if codes.len() < k {
        return;
    }
    let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let mut x = 0u64;
    for (i, &c) in codes.iter().enumerate() {
        x = ((x << 2) | c as u64) & mask;
        if i + 1 >= k {
            // Mixing in k and a salt keeps hashes of different k and roles independent.
            out.push(mix64(x ^ SEED ^ salt ^ ((k as u64) << 58)));
        }
    }
}

/// Index of the minimum of each window of `w` consecutive values (leftmost on ties).
fn sliding_argmin(vals: &[u64], w: usize, mut f: impl FnMut(usize, usize)) {
    let mut dq: VecDeque<usize> = VecDeque::new();
    for i in 0..vals.len() {
        while let Some(&b) = dq.back() {
            if vals[b] > vals[i] {
                dq.pop_back();
            } else {
                break;
            }
        }
        dq.push_back(i);
        if i + 1 >= w {
            let start = i + 1 - w;
            while *dq.front().unwrap() < start {
                dq.pop_front();
            }
            f(start, *dq.front().unwrap());
        }
    }
}

/// Run `f(codes)` on each ACGT segment and collect its seeds (spans in segment coordinates,
/// as `(label, start, end)`) into a single-level hierarchy over `seq`.
fn single_level(seq: &[u8], mut f: impl FnMut(&[u8], &mut Vec<(u64, u32, u32)>)) -> Hierarchy {
    let mut out = Vec::new();
    let mut cores = Vec::new();
    for (off, codes) in acgt_segments(seq) {
        out.clear();
        f(&codes, &mut out);
        let off = off as u32;
        cores.extend(out.iter().map(|&(label, s, e)| Core {
            label,
            start: s + off,
            end: e + off,
            child_lo: s + off,
            child_hi: e + off,
        }));
    }
    let ctx_start = cores.iter().map(|c| c.start).collect();
    Hierarchy {
        len: seq.len() as u32,
        levels: vec![Level { cores, ctx_start }],
    }
}

/// (k, w)-minimizers: the leftmost smallest-hash k-mer of each window of w consecutive
/// k-mers. Expected density 2/(w+1).
#[derive(Clone, Debug)]
pub struct Minimizer {
    pub k: usize,
    pub w: usize,
}

impl SeedScheme for Minimizer {
    fn seeds(&self, seq: &[u8]) -> Hierarchy {
        let (k, w) = (self.k, self.w);
        let mut h = Vec::new();
        single_level(seq, |codes, out| {
            kmer_hashes(codes, k, 0, &mut h);
            let mut last = usize::MAX;
            sliding_argmin(&h, w, |_, i| {
                if i != last {
                    out.push((h[i], i as u32, (i + k) as u32));
                    last = i;
                }
            });
        })
    }

    fn name(&self) -> String {
        format!("minimizer:k={},w={}", self.k, self.w)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SyncmerKind {
    /// Smallest s-mer at either end of the k-mer: density ≈ 2/(k-s+1).
    Closed,
    /// Smallest s-mer at offset t: density ≈ 1/(k-s+1).
    Open { t: usize },
}

/// Syncmers, optionally subsampled: a selected k-mer is kept if its hash falls in the lowest
/// `keep` fraction, which scales density by `keep` while staying content-defined.
#[derive(Clone, Debug)]
pub struct Syncmer {
    pub k: usize,
    pub s: usize,
    pub kind: SyncmerKind,
    pub keep: f64,
}

impl Syncmer {
    /// Selected k-mer positions and their hashes in one ACGT segment.
    fn select(&self, codes: &[u8], kh: &mut Vec<u64>, sh: &mut Vec<u64>, out: &mut Vec<(usize, u64)>) {
        out.clear();
        let (k, s) = (self.k, self.s);
        kmer_hashes(codes, k, 0, kh);
        kmer_hashes(codes, s, 0x5ca1_ab1e, sh);
        if kh.is_empty() {
            return;
        }
        let threshold = if self.keep >= 1.0 {
            u64::MAX
        } else {
            (self.keep * u64::MAX as f64) as u64
        };
        sliding_argmin(sh, k - s + 1, |i, m| {
            let off = m - i;
            let hit = match self.kind {
                SyncmerKind::Closed => off == 0 || off == k - s,
                SyncmerKind::Open { t } => off == t,
            };
            if hit && kh[i] <= threshold {
                out.push((i, kh[i]));
            }
        });
    }
}

impl SeedScheme for Syncmer {
    fn seeds(&self, seq: &[u8]) -> Hierarchy {
        let (mut kh, mut sh, mut sel) = (Vec::new(), Vec::new(), Vec::new());
        single_level(seq, |codes, out| {
            self.select(codes, &mut kh, &mut sh, &mut sel);
            out.extend(sel.iter().map(|&(i, h)| (h, i as u32, (i + self.k) as u32)));
        })
    }

    fn name(&self) -> String {
        let kind = match self.kind {
            SyncmerKind::Closed => "closed-syncmer".to_string(),
            SyncmerKind::Open { t } => format!("open-syncmer:t={t},"),
        };
        let kind = kind.trim_end_matches(',');
        let sep = if kind.contains(':') { "," } else { ":" };
        let keep = if self.keep < 1.0 {
            format!(",keep={}", self.keep)
        } else {
            String::new()
        };
        format!("{kind}{sep}k={},s={}{keep}", self.k, self.s)
    }
}

/// Order-2 randstrobes on syncmer positions: each syncmer is linked to the syncmer, among
/// those `w_min..=w_max` syncmers downstream, minimising popcount(h1 ^ h2). The seed spans
/// from the first strobe's start to the second's end.
#[derive(Clone, Debug)]
pub struct Randstrobe {
    pub base: Syncmer,
    pub w_min: usize,
    pub w_max: usize,
}

impl SeedScheme for Randstrobe {
    fn seeds(&self, seq: &[u8]) -> Hierarchy {
        let k = self.base.k;
        let (mut kh, mut sh, mut sel) = (Vec::new(), Vec::new(), Vec::new());
        single_level(seq, |codes, out| {
            self.base.select(codes, &mut kh, &mut sh, &mut sel);
            for j in 0..sel.len() {
                let lo = j + self.w_min;
                if lo >= sel.len() {
                    break;
                }
                let hi = (j + self.w_max).min(sel.len() - 1);
                let (p1, h1) = sel[j];
                let (p2, h2) = (lo..=hi)
                    .map(|m| sel[m])
                    .min_by_key(|&(_, h2)| (h1 ^ h2).count_ones())
                    .unwrap();
                let label = mix64(h1 ^ h2.rotate_left(29) ^ 0x57_0be5);
                out.push((label, p1 as u32, (p2 + k) as u32));
            }
        })
    }

    fn name(&self) -> String {
        format!(
            "randstrobe:{},wmin={},wmax={}",
            self.base.name().split_once(':').unwrap().1,
            self.w_min,
            self.w_max
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x & 3) as usize]
            })
            .collect()
    }

    fn density(s: &dyn SeedScheme, seq: &[u8]) -> f64 {
        s.seeds(seq).level(1).len() as f64 / seq.len() as f64
    }

    #[test]
    fn densities_match_theory() {
        let seq = random(500_000, 3);
        let m = Minimizer { k: 15, w: 10 };
        assert!((density(&m, &seq) - 2.0 / 11.0).abs() < 0.01);
        let c = Syncmer {
            k: 15,
            s: 5,
            kind: SyncmerKind::Closed,
            keep: 1.0,
        };
        assert!((density(&c, &seq) - 2.0 / 11.0).abs() < 0.01);
        let o = Syncmer {
            k: 15,
            s: 5,
            kind: SyncmerKind::Open { t: 3 },
            keep: 1.0,
        };
        assert!((density(&o, &seq) - 1.0 / 11.0).abs() < 0.01);
        let half = Syncmer {
            keep: 0.5,
            ..c.clone()
        };
        assert!((density(&half, &seq) - 1.0 / 11.0).abs() < 0.01);
    }

    #[test]
    fn seeds_are_content_defined() {
        // The same 1 kb embedded in different flanks yields the same seeds in its interior.
        let core = random(1000, 9);
        let a = [random(300, 1), core.clone(), random(300, 2)].concat();
        let b = [random(500, 4), core.clone(), random(200, 5)].concat();
        for s in [
            Box::new(Minimizer { k: 15, w: 10 }) as Box<dyn SeedScheme>,
            Box::new(Syncmer {
                k: 15,
                s: 5,
                kind: SyncmerKind::Closed,
                keep: 1.0,
            }),
            Box::new(Randstrobe {
                base: Syncmer {
                    k: 15,
                    s: 9,
                    kind: SyncmerKind::Open { t: 3 },
                    keep: 1.0,
                },
                w_min: 1,
                w_max: 4,
            }),
        ] {
            let inner = |seq: &[u8], off: u32| -> Vec<(u64, u32)> {
                s.seeds(seq)
                    .level(1)
                    .iter()
                    .filter(|c| c.start >= off + 200 && c.end <= off + 800)
                    .map(|c| (c.label, c.start - off))
                    .collect()
            };
            let (x, y) = (inner(&a, 300), inner(&b, 500));
            assert!(!x.is_empty());
            assert_eq!(x, y, "{}", s.name());
        }
    }

    #[test]
    fn skips_non_acgt() {
        let m = Minimizer { k: 5, w: 3 };
        let seq = b"ACGTACGTNNNNACGTACGTAC";
        for c in m.seeds(seq).level(1) {
            assert!(!seq[c.start as usize..c.end as usize].contains(&b'N'));
        }
    }
}
