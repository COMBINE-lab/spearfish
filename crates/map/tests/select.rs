//! `anchor::select` returns exactly `candidates(anchors(..))`: the same families in the same
//! order, bit-identical scores, the same anchors, and the same anchor count. Genes share a
//! repeat element so that seeds hit several families.

extern crate spearfish_map as map;
use index::{BuildParams, Collapse, Index, TxInput, build};
use map::{AnchorParams, MapParams, Mapper, OracleChainer, anchor};
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

/// Genes of random exons, each isoform a random exon subset; about half the exons carry one
/// of a few shared repeat elements.
fn genes(rng: &mut Rng, n: usize) -> Vec<TxInput> {
    let repeats: Vec<Vec<u8>> = (0..3)
        .map(|_| {
            let n = 60 + rng.below(80) as usize;
            rng.seq(n)
        })
        .collect();
    let mut txs = Vec::new();
    for g in 0..n {
        let exons: Vec<Vec<u8>> = (0..2 + rng.below(5) as usize)
            .map(|_| {
                let n = 40 + rng.below(200) as usize;
                let mut e = rng.seq(n);
                if rng.below(2) == 0 {
                    let r = &repeats[rng.below(repeats.len() as u64) as usize];
                    let at = rng.below(e.len() as u64) as usize;
                    e.splice(at..at, r.iter().copied());
                }
                e
            })
            .collect();
        for t in 0..1 + rng.below(4) as usize {
            let mut seq = Vec::new();
            for (e, ex) in exons.iter().enumerate() {
                if e == 0 || rng.below(3) != 0 {
                    seq.extend_from_slice(ex);
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

fn check(seed: u64, n_genes: usize, alpha_fam: f32, f_keep: usize) -> Result<(), TestCaseError> {
    let mut rng = Rng(seed.wrapping_mul(0x9e3779b97f4a7c15) | 1);
    let txs = genes(&mut rng, n_genes);
    let (data, _) = build(
        txs.clone(),
        &BuildParams {
            spec: "lcp:K=3".into(),
            hpc: true,
            f_max: 1000,
            allow_collisions: false,
            collapse: Collapse::Gene,
            infer_genes: None,
        },
    )
    .unwrap();
    let idx = std::sync::Arc::new(Index::from_data(&data).unwrap());
    let mapper = Mapper::new(
        idx.clone(),
        OracleChainer {
            p: Default::default(),
        },
        MapParams::default(),
    )
    .unwrap();
    let node_fam = anchor::node_families(&idx);
    let p = AnchorParams {
        alpha_fam,
        f_keep,
        f_max: f_keep.max(3),
        g_max: 50.0,
        ..Default::default()
    };
    for _ in 0..6 {
        let t = &txs[rng.below(txs.len() as u64) as usize];
        let read = mutate(&mut rng, &t.seq, 30);
        let seeds = mapper.seeder.seed(&read);
        let a = anchor::anchors(&idx, &seeds, &p);
        let n = a.len() as u32;
        let want = anchor::candidates(&idx, a, &p);
        let (got, got_n) = anchor::select(&idx, &node_fam, &seeds, &p);
        prop_assert_eq!(got_n, n, "anchor count, seed {}", seed);
        prop_assert_eq!(got.len(), want.len(), "candidates, seed {}", seed);
        for (g, w) in got.iter().zip(&want) {
            prop_assert_eq!(g.family, w.family, "seed {}", seed);
            prop_assert_eq!(g.score.to_bits(), w.score.to_bits(), "seed {}", seed);
            prop_assert_eq!(&g.anchors, &w.anchors, "seed {}", seed);
        }
    }
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]
    #[test]
    fn select_matches_anchors_then_candidates(
        seed in any::<u64>(),
        n_genes in 2usize..12,
        alpha_fam in prop_oneof![Just(0.0f32), Just(0.5f32)],
        f_keep in 1usize..6,
    ) {
        check(seed, n_genes, alpha_fam, f_keep)?;
    }
}
