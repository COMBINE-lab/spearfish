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
        txs.clone(),
        &BuildParams {
            spec: spec.into(),
            hpc: true,
            f_max: 20,
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

/// Antisense genes sharing sequence: a read from the shared region is reported on both
/// orientations with equal scores; a read spanning the unshared part prefers its own.
#[test]
fn antisense_overlap_reports_both_orientations() {
    let mut x = 0x2545f4914f6cdd1du64;
    let seq: Vec<u8> = (0..3000)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[(x >> 62) as usize]
        })
        .collect();
    let txs = vec![
        TxInput {
            name: "sense".into(),
            gene: "A".into(),
            seq: seq.clone(),
        },
        TxInput {
            name: "anti".into(),
            gene: "B".into(),
            seq: revcomp(&seq[1000..2500]),
        },
    ];
    let (data, _) = build(
        txs,
        &BuildParams {
            spec: "lcp:K=3".into(),
            hpc: true,
            f_max: 20,
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
    let r = mapper.map(&mutate(&seq[1200..2300], 0.03, 11));
    let get = |name: &str| r.hits.iter().find(|h| idx.tx_name(h.tx) == name);
    let (s, a) = (get("sense").expect("sense"), get("anti").expect("anti"));
    assert!(!s.reverse && a.reverse, "{:?}", r.hits);
    assert_eq!(s.score, a.score, "{:?}", r.hits);
    let r = mapper.map(&mutate(&seq[200..2300], 0.03, 12));
    assert_eq!(idx.tx_name(r.hits[0].tx), "sense");
    assert!(r.hits.iter().skip(1).all(|h| h.score < r.hits[0].score));
}

/// Reads too noisy for their level-3 chains to cover half of them are rescued at level 2 onto
/// their own transcript; without the rescue they stay low-coverage.
#[test]
fn noisy_reads_are_rescued() {
    let txs = toy();
    let (data, _) = build(
        txs.clone(),
        &BuildParams {
            spec: "lcp:K=3".into(),
            hpc: true,
            f_max: 20,
            allow_collisions: false,
            collapse: Collapse::Gene,
            infer_genes: None,
        },
    )
    .unwrap();
    let idx = std::sync::Arc::new(Index::from_data(&data).unwrap());
    let mapper = |p: MapParams| {
        Mapper::new(
            idx.clone(),
            OracleChainer {
                p: Default::default(),
            },
            p,
        )
        .unwrap()
    };
    let on = mapper(MapParams::default());
    let off = mapper(MapParams {
        rescue: None,
        ..MapParams::default()
    });
    let mut rescued = 0;
    for (t, tx) in txs.iter().enumerate() {
        for (rev, seed) in [(false, 1u64), (true, 2u64)] {
            let mut read = mutate(&tx.seq, 0.2, 0x9e37 + seed * 7919 + t as u64);
            if rev {
                read = revcomp(&read);
            }
            let r = on.map(&read);
            if !r.rescued {
                continue;
            }
            rescued += 1;
            assert!(
                !r.low_coverage && r.reverse == rev,
                "{} rev {rev}: {r:?}",
                tx.name
            );
            assert_eq!(idx.tx_name(r.hits[0].tx), tx.name, "{:?}", r.hits);
            assert!(
                r.hits.iter().skip(1).all(|h| h.score < r.hits[0].score),
                "{:?}",
                r.hits
            );
            let r = off.map(&read);
            assert!(!r.rescued && r.low_coverage, "{} rev {rev}: {r:?}", tx.name);
        }
    }
    // At least half of the 12 reads (7 when written; the others keep half their coverage at
    // level 3, or have no base anchor at all).
    assert!(rescued >= 6, "{rescued} reads rescued");
}
