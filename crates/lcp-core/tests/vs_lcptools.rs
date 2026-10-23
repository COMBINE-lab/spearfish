//! Differential test: lcp-core's core boundaries equal lcptools' at every level.
//! lcptools reports DCT-extended starts above level 1, which we expose as `ctx_start`.

extern crate spearfish_lcp as lcp_core;
use lcp_core::{LcpParser, Parser};
use lcp_oracle_sys::lcptools_boundaries;
use proptest::prelude::*;

const MAX_LEVEL: u8 = 6;

fn ours(seq: &[u8]) -> Vec<Vec<(u32, u32)>> {
    let h = LcpParser::default().parse(seq, MAX_LEVEL);
    h.levels
        .iter()
        .map(|l| {
            l.cores
                .iter()
                .zip(&l.ctx_start)
                .map(|(c, &s)| (s, c.end))
                .collect()
        })
        .collect()
}

fn check(seq: &[u8]) {
    let a = ours(seq);
    let b = lcptools_boundaries(seq, MAX_LEVEL);
    for k in 0..MAX_LEVEL as usize {
        assert_eq!(a[k], b[k], "level {} differs (len {})", k + 1, seq.len());
    }
}

/// Sequences over a skewed alphabet so runs (RINT cores) and repeats are common.
fn dna(max_len: usize) -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(prop::sample::select(b"AACGTTTG".to_vec()), 0..max_len)
}

/// Tandem repeats: a short unit repeated with a few substitutions.
fn repeats() -> impl Strategy<Value = Vec<u8>> {
    (
        prop::collection::vec(prop::sample::select(b"ACGT".to_vec()), 1..8),
        10usize..400,
        prop::collection::vec((0usize..4000, prop::sample::select(b"ACGT".to_vec())), 0..5),
    )
        .prop_map(|(unit, reps, subs)| {
            let mut s: Vec<u8> = unit.iter().cycle().take(unit.len() * reps).copied().collect();
            for (p, b) in subs {
                if !s.is_empty() {
                    let n = s.len();
                    s[p % n] = b;
                }
            }
            s
        })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    #[test]
    fn random_matches_lcptools(seq in dna(3000)) {
        check(&seq);
    }

    #[test]
    fn repeats_match_lcptools(seq in repeats()) {
        check(&seq);
    }
}

#[test]
fn long_uniform_matches_lcptools() {
    let mut x: u64 = 0x1234_5678;
    let seq: Vec<u8> = (0..200_000)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[(x & 3) as usize]
        })
        .collect();
    check(&seq);
}
