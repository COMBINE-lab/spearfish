//! LCP in the lcptools style: level-1 cores from local rules on the 2-bit alphabet, then
//! repeated deterministic coin tossing (DCT) and the same rules on the reduced alphabet.
//!
//! Core *boundaries* reproduce lcptools (BilkentCompGen/lcptools, `lps.c`) exactly; this is
//! checked by the differential test against `lcp-oracle-sys`. Labels differ by design:
//! lcptools labels are lossy summaries, ours are 64-bit hashes of the full content
//! (level 1: the bases; level k: the child labels).

use xxhash_rust::xxh3::xxh3_64_with_seed;

use crate::bits::{BitArena, dct_value, dct_value_words, value_width};
use crate::encode::acgt_segments;
use crate::{Core, Hierarchy, Level, Parser, ParserMeta};

/// Bumped whenever a change alters any parse or label.
pub const LCP_PARSER_VERSION: u32 = 1;

#[derive(Clone, Debug)]
pub struct LcpConfig {
    /// Global label seed, shared by index and reads.
    pub seed: u64,
    /// DCT iterations between levels (lcptools `DCT`, default 1).
    pub dct_iters: u8,
}

impl Default for LcpConfig {
    fn default() -> Self {
        LcpConfig {
            seed: 0x6c63_706d_6170,
            dct_iters: 1,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct LcpParser {
    cfg: LcpConfig,
}

impl LcpParser {
    pub fn new(cfg: LcpConfig) -> Self {
        assert!(cfg.dct_iters >= 1, "at least one DCT iteration is required");
        LcpParser { cfg }
    }

    #[inline]
    fn level_seed(&self, level: u8) -> u64 {
        self.cfg.seed ^ (level as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15)
    }

    /// Parse one run of 2-bit codes. Returns exactly `max_level` levels.
    fn parse_codes(&self, s: &[u8], max_level: u8) -> Vec<Level> {
        let mut levels: Vec<Level> = Vec::with_capacity(max_level as usize);
        if max_level == 0 {
            return levels;
        }

        // Level 1.
        let mut spans = Vec::with_capacity(s.len() / 2);
        parse_rules(s, |lo, hi| spans.push((lo, hi)));
        let seed1 = self.level_seed(1);
        let mut bits = BitArena::with_capacity(spans.len());
        let mut level = Level {
            cores: Vec::with_capacity(spans.len()),
            ctx_start: Vec::with_capacity(spans.len()),
        };
        for &(lo, hi) in &spans {
            level.cores.push(Core {
                label: xxh3_64_with_seed(&s[lo..hi], seed1),
                start: lo as u32,
                end: hi as u32,
                child_lo: lo as u32,
                child_hi: hi as u32,
            });
            level.ctx_start.push(lo as u32);
            for &c in s[lo..hi].iter().rev() {
                bits.push_bits(c as u64, 2);
            }
            bits.finish_item();
        }
        levels.push(level);

        let d = self.cfg.dct_iters as usize;
        let mut label_buf: Vec<u8> = Vec::new();
        for k in 2..=max_level {
            let prev = levels.last().unwrap();
            let n = prev.cores.len();
            if n < d + 1 {
                levels.push(Level::default());
                continue;
            }
            let vals = dct(&bits, d);

            spans.clear();
            parse_rules(&vals[d..], |lo, hi| spans.push((lo + d, hi + d)));

            let seed = self.level_seed(k);
            let mut next_bits = BitArena::with_capacity(spans.len());
            let mut level = Level {
                cores: Vec::with_capacity(spans.len()),
                ctx_start: Vec::with_capacity(spans.len()),
            };
            for &(lo, hi) in &spans {
                let children = &prev.cores[lo..hi];
                label_buf.clear();
                for c in children {
                    label_buf.extend_from_slice(&c.label.to_le_bytes());
                }
                level.cores.push(Core {
                    label: xxh3_64_with_seed(&label_buf, seed),
                    start: children[0].start,
                    end: children.iter().map(|c| c.end).max().unwrap(),
                    child_lo: lo as u32,
                    child_hi: hi as u32,
                });
                // DCT replaced each core's start by that of the core `d` places to its left.
                level.ctx_start.push(prev.ctx_start[lo - d]);
                for &v in vals[lo..hi].iter().rev() {
                    next_bits.push_bits(v, value_width(v));
                }
                next_bits.finish_item();
            }
            levels.push(level);
            bits = next_bits;
        }
        levels
    }
}

impl Parser for LcpParser {
    fn parse(&self, seq: &[u8], max_level: u8) -> Hierarchy {
        let mut h = Hierarchy {
            len: seq.len() as u32,
            levels: vec![Level::default(); max_level as usize],
        };
        for (off, codes) in acgt_segments(seq) {
            let seg = self.parse_codes(&codes, max_level);
            // Child indices of level k point into level k-1, so offset them by the number of
            // level k-1 cores already present.
            let before: Vec<u32> = h.levels.iter().map(|l| l.cores.len() as u32).collect();
            for (k, level) in seg.into_iter().enumerate() {
                let off = off as u32;
                let child_off = if k == 0 { off } else { before[k - 1] };
                let dst = &mut h.levels[k];
                dst.cores.extend(level.cores.into_iter().map(|c| Core {
                    start: c.start + off,
                    end: c.end + off,
                    child_lo: c.child_lo + child_off,
                    child_hi: c.child_hi + child_off,
                    ..c
                }));
                dst.ctx_start.extend(level.ctx_start.into_iter().map(|s| s + off));
            }
        }
        h
    }

    fn meta(&self) -> ParserMeta {
        ParserMeta {
            name: "lcp",
            version: LCP_PARSER_VERSION,
            seed: self.cfg.seed,
            dct_iters: self.cfg.dct_iters,
        }
    }
}

/// DCT values after `d` iterations. Entries `[d..]` are valid; the first `d` cores have no
/// left neighbour at some iteration and are dropped, as in lcptools.
fn dct(bits: &BitArena, d: usize) -> Vec<u64> {
    let n = bits.len();
    let mut vals = vec![0u64; n];
    for (i, v) in vals.iter_mut().enumerate().skip(1) {
        *v = dct_value(bits.get(i - 1), bits.get(i));
    }
    for t in 1..d {
        // Iteration t compresses cores i > t against their left neighbour's previous value.
        for i in (t + 1..n).rev() {
            let (l, r) = (vals[i - 1], vals[i]);
            vals[i] = dct_value_words(l, value_width(l), r, value_width(r));
        }
    }
    vals
}

/// The LCP core rules, shared by all levels (lcptools `parse1`/`parse3`). Calls
/// `emit(lo, hi)` for each core `[lo, hi)` over `s`, in order:
///
/// * RINT: a symbol, a run of length >= 2 of another symbol, and the symbol after the run
///   (skipped if the run reaches the end of `s`);
/// * LMIN: `s[i] > s[i+1] < s[i+2]`;
/// * LMAX: `s[i-1] <= s[i] < s[i+1] > s[i+2] >= s[i+3]`;
/// * SSEQ: fills the gap between consecutive cores, overlapping each by one symbol.
pub(crate) fn parse_rules<T: Ord>(s: &[T], mut emit: impl FnMut(usize, usize)) {
    let n = s.len();
    // End of the previous core; usize::MAX means "none yet" (no SSEQ before the first core).
    let mut prev_end = usize::MAX;
    let mut core = |lo: usize, hi: usize, prev_end: &mut usize| {
        if *prev_end < lo {
            emit(*prev_end - 1, lo + 1);
        }
        emit(lo, hi);
        *prev_end = hi;
    };

    let mut i = 0;
    while i + 2 < n {
        if s[i] == s[i + 1] {
            i += 1;
            continue;
        }
        if s[i + 1] == s[i + 2] {
            let mut t = i + 2;
            while t < n && s[t - 1] == s[t] {
                t += 1;
            }
            if t != n {
                core(i, t + 1, &mut prev_end);
                i += 1;
                continue;
            }
        }
        if s[i] > s[i + 1] && s[i + 1] < s[i + 2] {
            core(i, i + 3, &mut prev_end);
            i += 1;
            continue;
        }
        if i > 0
            && i + 3 < n
            && s[i] < s[i + 1]
            && s[i + 1] > s[i + 2]
            && s[i - 1] <= s[i]
            && s[i + 2] >= s[i + 3]
        {
            core(i, i + 3, &mut prev_end);
        }
        i += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spans(s: &[u8]) -> Vec<(usize, usize)> {
        let codes: Vec<u8> = s.iter().map(|&b| crate::encode::CODE[b as usize]).collect();
        let mut v = Vec::new();
        parse_rules(&codes, |lo, hi| v.push((lo, hi)));
        v
    }

    #[test]
    fn lmin_cores() {
        // A=0 C=1 G=2 T=3. "CAC" is a local minimum at A.
        assert_eq!(spans(b"CACAC"), vec![(0, 3), (2, 5)]);
    }

    #[test]
    fn rint_core_includes_flanks() {
        // C AAAA G: RINT over [0, 6); run reaching the end is not a core.
        assert_eq!(spans(b"CAAAAG"), vec![(0, 6)]);
        assert_eq!(spans(b"CAAAA"), vec![]);
    }

    #[test]
    fn sseq_fills_gaps() {
        // ACGT is strictly increasing: no rule fires until the LMIN "TAC" at 3.
        // Then "ACGTTTC" ... construct: LMIN at 0 ("CAG"), plain ascent, LMIN later.
        let v = spans(b"CAGTTTTGAC");
        // Every gap between consecutive cores is filled.
        for w in v.windows(2) {
            assert!(w[1].0 <= w[0].1, "gap between {:?} and {:?}", w[0], w[1]);
        }
    }

    #[test]
    fn hierarchy_children_are_consistent() {
        let p = LcpParser::default();
        let seq: Vec<u8> = (0..5000u64)
            .map(|i| b"ACGT"[((i.wrapping_mul(0x9e3779b97f4a7c15) >> 61) & 3) as usize])
            .collect();
        let h = p.parse(&seq, 5);
        assert_eq!(h.num_levels(), 5);
        for k in 2..=5u8 {
            let prev = h.level(k - 1);
            for c in h.level(k) {
                let ch = &prev[c.child_lo as usize..c.child_hi as usize];
                assert!(ch.len() >= 2);
                assert_eq!(c.start, ch[0].start);
                assert_eq!(c.end, ch.iter().map(|x| x.end).max().unwrap());
            }
        }
    }

    #[test]
    fn n_splits_parse() {
        let p = LcpParser::default();
        let seq = b"CACACAGTGTCANNNNCACACAGTGTCA";
        let h = p.parse(seq, 1);
        for c in h.level(1) {
            assert!(!seq[c.start as usize..c.end as usize].contains(&b'N'));
        }
        // Same content on both sides of the Ns gives the same labels.
        let left: Vec<u64> = h
            .level(1)
            .iter()
            .filter(|c| c.end <= 12)
            .map(|c| c.label)
            .collect();
        let right: Vec<u64> = h
            .level(1)
            .iter()
            .filter(|c| c.start >= 16)
            .map(|c| c.label)
            .collect();
        assert_eq!(left, right);
    }
}
