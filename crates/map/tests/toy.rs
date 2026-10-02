//! End-to-end mapping on the toy transcriptome: reads from each transcript, on both strands,
//! must map with their own transcript on top (exon skip, alternative 3' site, paralogs).

use std::path::PathBuf;

use index::{BuildParams, Collapse, Index, TxInput, build};
use lcp_core::encode::revcomp;
use map::{MapParams, Mapper, OracleChainer};

fn toy() -> Vec<TxInput> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/toy/transcripts.fa");
    std::fs::read_to_string(path)
        .unwrap()
        .split('>')
        .skip(1)
        .map(|rec| {
            let mut lines = rec.lines();
            let h: Vec<&str> = lines.next().unwrap().split('|').collect();
            TxInput {
                name: h[0].into(),
                gene: h[1].into(),
                seq: lines.flat_map(|l| l.trim().bytes()).collect(),
            }
        })
        .collect()
}

/// Deterministic substitutions/indels at about `rate` (xorshift).
fn mutate(seq: &[u8], rate: f64, mut x: u64) -> Vec<u8> {
    let mut next = || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        (x >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq {
        let r = next();
        if r < rate / 3.0 {
            continue;
        } else if r < 2.0 * rate / 3.0 {
            out.push(if b == b'A' { b'C' } else { b'A' });
        } else {
            out.push(b);
            if r > 1.0 - rate / 3.0 {
                out.push(b'G');
            }
        }
    }
    out
}

fn check(spec: &str, rate: f64) {
    let txs = toy();
    let (data, _) = build(
        &txs,
        &BuildParams {
            spec: spec.into(),
            hpc: true,
            f_max: 20,
            allow_collisions: false,
            collapse: Collapse::Gene,
        },
    )
    .unwrap();
    let idx = std::sync::Arc::new(Index::from_bytes(&data.to_bytes()).unwrap());
    let mapper = Mapper::new(
        idx.clone(),
        OracleChainer {
            p: Default::default(),
        },
        MapParams::default(),
    )
    .unwrap();
    for (t, tx) in txs.iter().enumerate() {
        for (rev, seed) in [(false, 1u64), (true, 2u64)] {
            let mut read = mutate(&tx.seq, rate, 0x9e37 + seed * 7919 + t as u64);
            if rev {
                read = revcomp(&read);
            }
            let r = mapper.map(&read);
            assert_eq!(r.reverse, rev, "{spec} rate {rate}: {} strand", tx.name);
            let top = r
                .hits
                .first()
                .unwrap_or_else(|| panic!("{spec} rate {rate}: {} unmapped", tx.name));
            assert_eq!(
                idx.tx_name(top.tx),
                tx.name,
                "{spec} rate {rate}: hits {:?}",
                r.hits
            );
            // Strictly better than every other transcript (no ties with the wrong isoform).
            assert!(
                r.hits.iter().skip(1).all(|h| h.score < top.score),
                "{spec} rate {rate}: {:?}",
                r.hits
            );
            assert!(!r.low_coverage || rate > 0.0, "{spec}: {} low coverage", tx.name);
        }
    }
}

#[test]
fn toy_reads_map_to_their_transcript() {
    for spec in ["lcp:K=3", "minimizer:k=19,w=23"] {
        for rate in [0.0, 0.03, 0.06] {
            check(spec, rate);
        }
    }
}
