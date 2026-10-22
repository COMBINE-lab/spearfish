//! With an unbounded beam, the refinement chainer's per-transcript scores equal the oracle's
//! exactly, on random multi-isoform genes (alternative exons, repeated exons) and reads drawn
//! from them with errors.

use std::collections::BTreeMap;

use index::{BuildParams, Collapse, Index, TxInput, build};
use map::{AnchorParams, ChainParams, Chainer, MapParams, Mapper, OracleChainer, RefineChainer, anchor};
use proptest::prelude::*;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn seq(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[self.below(4) as usize]).collect()
    }
}

/// Genes of random exons; isoforms are random exon subsets in order, sometimes repeating an
/// exon (tandem duplication) so a unitig occurs twice in one transcript.
fn family(seed: u64, genes: usize) -> Vec<TxInput> {
    let mut rng = Rng(seed | 1);
    let mut txs = Vec::new();
    for g in 0..genes {
        let exons: Vec<Vec<u8>> = (0..3 + rng.below(6) as usize)
            .map(|_| {
                let n = 40 + rng.below(300) as usize;
                rng.seq(n)
            })
            .collect();
        let n_iso = 1 + rng.below(12) as usize;
        for t in 0..n_iso {
            let mut seq = Vec::new();
            for (e, ex) in exons.iter().enumerate() {
                if e == 0 || rng.below(3) != 0 {
                    seq.extend_from_slice(ex);
                    if rng.below(12) == 0 {
                        seq.extend_from_slice(ex);
                    }
                }
            }
            txs.push(TxInput {
                name: format!("g{g}t{t}"),
                gene: format!("g{g}"),
                seq,
            });
        }
    }
    txs
}

fn mutate(rng: &mut Rng, seq: &[u8], rate: u64) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq {
        match rng.below(1000) {
            x if x < rate / 3 => {}
            x if x < 2 * rate / 3 => out.push(b"ACGT"[rng.below(4) as usize]),
            x if x < rate => {
                out.push(b);
                out.push(b"ACGT"[rng.below(4) as usize]);
            }
            _ => out.push(b),
        }
    }
    out
}

fn scores(hits: Vec<map::TxHit>) -> BTreeMap<u32, u32> {
    hits.into_iter().map(|h| (h.tx, h.score.to_bits())).collect()
}

fn check(seed: u64, spec: &str, h: usize, beta: f32) -> Result<(), TestCaseError> {
    let txs = family(seed, 3);
    let (data, _) = build(
        txs.clone(),
        &BuildParams {
            spec: spec.into(),
            hpc: true,
            f_max: 1000,
            allow_collisions: false,
            collapse: Collapse::Gene,
            infer_genes: None,
        },
    )
    .unwrap();
    let idx = std::sync::Arc::new(Index::from_data(&data).unwrap());
    let p = ChainParams {
        h,
        beta,
        ..Default::default()
    };
    let mapper = Mapper::new(idx.clone(), OracleChainer { p: p.clone() }, MapParams::default()).unwrap();
    let oracle = OracleChainer { p: p.clone() };
    let refine = RefineChainer { p, beam: None };
    let mut rng = Rng(seed.wrapping_mul(0x9e3779b97f4a7c15) | 1);
    for _ in 0..8 {
        let t = &txs[rng.below(txs.len() as u64) as usize];
        let a = rng.below(t.seq.len() as u64 / 2) as usize;
        let b = (a + 100 + rng.below(t.seq.len() as u64) as usize).min(t.seq.len());
        let rate = rng.below(80);
        let read = mutate(&mut rng, &t.seq[a..b], rate);
        let seeds = mapper.seeder.seed(&read);
        let anchors = anchor::anchors(&idx, &seeds, &AnchorParams::default());
        for mut c in anchor::candidates(&idx, anchors, &AnchorParams::default()) {
            let (o, r) = (oracle.chain(&idx, &c), refine.chain(&idx, &c));
            prop_assert_eq!(scores(o), scores(r), "seed {} read from {}", seed, t.name);
            // With collinear anchors merged, as the mapper chains candidates by default.
            anchor::merge_collinear(&idx, &mut c);
            let (o, r) = (oracle.chain(&idx, &c), refine.chain(&idx, &c));
            prop_assert_eq!(scores(o), scores(r), "merged: seed {} read from {}", seed, t.name);
        }
    }
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn refine_equals_oracle(seed in any::<u64>(), h in prop_oneof![Just(50usize), 1usize..8]) {
        check(seed, "lcp:K=2", h, 0.0)?;
        check(seed, "lcp:K=3", h, 0.0)?;
        check(seed, "minimizer:k=11,w=5", h, 0.0)?;
    }

    #[test]
    fn refine_equals_oracle_with_beta(seed in any::<u64>()) {
        check(seed, "lcp:K=2", 50, 0.5)?;
    }
}
