//! Index construction tests: the round trip (every transcript's seeds are reconstructed
//! exactly from the serialized index) on the toy fixture and on random transcriptomes.

use std::collections::HashSet;
use std::path::PathBuf;

use index::{BuildParams, Collapse, Index, IndexData, InferGenes, TxInput, build};
use proptest::prelude::*;
use seeds::BaseSeeder;

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

/// Build, serialize, reload, and check every transcript's seeds against fresh seeding.
fn check_roundtrip(txs: &[TxInput], spec: &str, hpc: bool, collapse: Collapse) -> Index {
    let p = BuildParams {
        spec: spec.into(),
        hpc,
        f_max: 20,
        allow_collisions: false,
        collapse,
        infer_genes: None,
    };
    let (data, rep) = build(txs.to_vec(), &p).unwrap();
    assert_eq!(rep.inconsistent_lengths, 0, "{spec}");
    let idx = Index::from_bytes(&data.to_bytes()).unwrap();
    assert_eq!(idx.meta("spec"), Some(spec));
    let direct = Index::from_data(&data).unwrap();
    assert_eq!(direct.section_checksums(), idx.section_checksums());
    let seeder = BaseSeeder::new(spec, hpc).unwrap();
    let got = idx.transcript_seeds();
    for (t, tx) in txs.iter().enumerate() {
        // A label's node for transcript t: the one whose unitig t traverses.
        let node = |label: u64| {
            idx.lookup(label)
                .find(|&n| {
                    idx.occurrences(idx.node_unitig()[n as usize])
                        .0
                        .contains(&(t as u32))
                })
                .unwrap()
        };
        let mut want: Vec<(u32, u32)> = seeder
            .seed(&tx.seq)
            .iter()
            .map(|s| (node(s.label), s.start))
            .collect();
        want.sort_unstable_by_key(|&(n, p)| (p, n));
        assert_eq!(got[t], want, "{spec} hpc={hpc}: transcript {}", tx.name);
        assert_eq!(idx.tx_name(t as u32), tx.name);
    }
    // Stored parsed sequences decode back exactly.
    for (t, tx) in txs.iter().enumerate() {
        let parsed = seeder.parsed(&tx.seq).map_or(tx.seq.clone(), |h| h.seq);
        let codes: Vec<u8> = parsed
            .iter()
            .map(|&b| lcp_core::encode::CODE[b as usize])
            .map(|c| if c > 3 { 0 } else { c })
            .collect();
        let n = idx.tx_plen()[t];
        assert_eq!(idx.tx_codes(t as u32, 0, n), codes);
        assert_eq!(idx.tx_range(t as u32, 0, n, false).unwrap().decode(), codes);
    }
    // HPC: raw windows of parsed spans decode to the raw sequence's whole runs.
    for (t, tx) in txs.iter().enumerate() {
        let Some(h) = seeder.parsed(&tx.seq) else {
            assert_eq!(idx.tx_raw_codes(t as u32, 0, 0), None);
            assert!(idx.tx_range(t as u32, 0, 0, true).is_none());
            continue;
        };
        let n = idx.tx_plen()[t];
        let raw: Vec<u8> = tx
            .seq
            .iter()
            .map(|&b| lcp_core::encode::CODE[b as usize] & 3)
            .collect();
        for (lo, hi) in [
            (0, n),
            (0, 0),
            (n, n),
            (n / 3, n / 2),
            (n.saturating_sub(70), n),
            (1.min(n), n),
        ] {
            let (a, b) = h.to_raw(lo, hi);
            let range = idx.tx_range(t as u32, lo, hi, true).unwrap();
            assert_eq!(range.len(), (b - a) as usize);
            assert_eq!(range.decode(), &raw[a as usize..b as usize]);
            assert_eq!(
                idx.tx_raw_codes(t as u32, lo, hi).as_deref(),
                Some(&raw[a as usize..b as usize]),
                "{spec}: transcript {} [{lo}, {hi})",
                tx.name
            );
        }
        // Every parsed position maps to the raw start of its run.
        for q in 0..=n {
            assert_eq!(idx.tx_raw_pos(t as u32, q), Some(h.to_raw(q, q).0), "{spec}: {q}");
        }
    }
    // Colors agree with occurrences.
    for u in 0..idx.n_unitigs() as u32 {
        let (ts, _) = idx.occurrences(u);
        let set: Vec<u32> = ts
            .iter()
            .copied()
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect();
        assert_eq!(idx.color(idx.unitig_color()[u as usize]), &set[..]);
    }
    idx
}

#[test]
fn toy_roundtrip_and_structure() {
    let txs = toy();
    for spec in [
        "lcp:K=3",
        "lcp:K=4",
        "minimizer:k=19,w=23",
        "closed-syncmer:k=15,s=5",
    ] {
        for hpc in [false, true] {
            for collapse in [Collapse::Global, Collapse::Gene] {
                let idx = check_roundtrip(&txs, spec, hpc, collapse);
                let fam = idx.tx_family();
                assert_eq!(fam[0], fam[1], "{spec}: GENEA isoforms");
                assert_eq!(fam[2], fam[3], "{spec}: GENEB isoforms");
                // Paralogs GENEC/GENECP share unitigs only when nodes are global.
                let global = collapse == Collapse::Global;
                assert_eq!(fam[4] == fam[5], global, "{spec} {collapse:?}: paralogs");
                let distinct: HashSet<u32> = fam.iter().copied().collect();
                assert_eq!(distinct.len(), if global { 3 } else { 4 }, "{spec} {collapse:?}");
                // GENEA: a unitig shared by both isoforms, and one specific to the isoform
                // with the cassette exon.
                let colors: HashSet<Vec<u32>> = (0..idx.n_unitigs() as u32)
                    .map(|u| idx.color(idx.unitig_color()[u as usize]).to_vec())
                    .collect();
                assert!(
                    colors.contains(&vec![0, 1]) && colors.contains(&vec![0]),
                    "{spec}"
                );
                let mixed = colors.iter().any(|c| c.contains(&4) && c.contains(&5));
                assert_eq!(mixed, global, "{spec} {collapse:?}");
            }
        }
    }
}

#[test]
fn mismatched_sections_fail_to_load() {
    let p = BuildParams {
        spec: "lcp:K=3".into(),
        hpc: true,
        f_max: 20,
        allow_collisions: false,
        collapse: Collapse::Gene,
        infer_genes: None,
    };
    let (data, _) = build(toy(), &p).unwrap();
    assert!(Index::from_data(&data).unwrap().has_raw());
    let fails = |d: &IndexData, section: &str| match Index::from_data(d) {
        Ok(_) => panic!("{section}: loaded"),
        Err(e) => assert!(e.to_string().contains(section), "{section}: {e}"),
    };
    let mut d = data.clone();
    d.occ_tx.pop();
    fails(&d, "occ_tx");
    let mut d = data.clone();
    d.node_len.push(0);
    fails(&d, "node_len");
    // A damaged raw section is an error, not an index without raw sequence.
    let mut d = data.clone();
    d.tx_raw_ck.1.pop();
    fails(&d, "tx_raw_ck");
}

#[test]
fn file_roundtrip() {
    let txs = toy();
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
    let path = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("toy.sfi");
    data.write(&path).unwrap();
    let idx = Index::open(&path).unwrap();
    assert_eq!(idx.labels(), &data.labels[..]);
    assert_eq!(idx.occ_pos(), &data.occ_pos[..]);
    assert_eq!(idx.gene_name(idx.tx_gene()[4]), "TOYG0000003");
}

/// Random transcriptomes: genes of random exons; isoforms pick exon subsets; some exons
/// carry tandem repeats or are copied into other genes, exercising branching, repeated
/// nodes and shared unitigs.
fn transcriptome() -> impl Strategy<Value = Vec<TxInput>> {
    let exon = prop::collection::vec(prop::sample::select(b"ACGT".to_vec()), 30..300);
    let gene = (
        prop::collection::vec(exon, 2..6),
        prop::collection::vec(prop::collection::vec(any::<bool>(), 6), 1..5),
    );
    prop::collection::vec(gene, 1..5).prop_map(|genes| {
        let mut txs = Vec::new();
        let shared: Vec<u8> = genes[0].0[0].clone();
        for (g, (mut exons, isoforms)) in genes.into_iter().enumerate() {
            // A tandem repeat inside the last exon, and a copy of gene 0's first exon.
            let last = exons.len() - 1;
            let unit: Vec<u8> = exons[last].iter().take(12).copied().collect();
            exons[last].extend(unit.iter().cycle().take(unit.len() * 4));
            if g > 0 {
                exons.push(shared.clone());
            }
            for (i, pick) in isoforms.iter().enumerate() {
                let seq: Vec<u8> = exons
                    .iter()
                    .zip(pick.iter().chain(std::iter::repeat(&true)))
                    .filter(|(_, k)| **k)
                    .flat_map(|(e, _)| e.clone())
                    .collect();
                if !seq.is_empty() {
                    txs.push(TxInput {
                        name: format!("g{g}t{i}"),
                        gene: format!("g{g}"),
                        seq,
                    });
                }
            }
        }
        txs
    })
}

/// Gene inference on the toy transcriptome: isoforms that share exons are grouped, the paralog
/// pair (GENEC, GENECP) stays apart at the default run length, and a run held by more than
/// `max_tx` transcripts links nothing.
#[test]
fn inferred_genes_follow_shared_exons() {
    let mut txs = toy();
    let annotated: Vec<String> = txs.iter().map(|t| t.gene.clone()).collect();
    for t in &mut txs {
        t.gene = t.name.clone();
    }
    let groups = |span, max_tx| {
        let p = BuildParams {
            spec: "lcp:K=3".into(),
            hpc: true,
            f_max: 20,
            allow_collisions: false,
            collapse: Collapse::Gene,
            infer_genes: Some(InferGenes {
                span,
                max_tx,
                attach_span: 4,
                attach_max: 3,
            }),
        };
        build(txs.clone(), &p).unwrap().1.tx_gene
    };
    // Same partition as the annotation: same group iff same annotated gene.
    let g = groups(12, 64);
    for i in 0..txs.len() {
        for j in 0..txs.len() {
            assert_eq!(
                g[i] == g[j],
                annotated[i] == annotated[j],
                "{} {}",
                txs[i].name,
                txs[j].name
            );
        }
    }
    // Runs held by two or more transcripts are all above max_tx = 1: nothing is linked.
    let g = groups(12, 1);
    let distinct: HashSet<u32> = g.iter().copied().collect();
    assert_eq!(distinct.len(), txs.len());
}

#[test]
fn raw_sequence_skipped_where_non_acgt_merges_runs() {
    // "AN" is two parsed runs but packs as "AA" (one run): no raw sequence for that transcript.
    let txs = vec![
        TxInput {
            name: "a".into(),
            gene: "g".into(),
            seq: b"ACGTTGCAANCCGGTTAACGTAGCATGCATCGA".to_vec(),
        },
        TxInput {
            name: "b".into(),
            gene: "g".into(),
            seq: b"ACGTTGCAAGCCGGTTAACGTAGCATGCATCGA".to_vec(),
        },
    ];
    let p = BuildParams {
        spec: "lcp:K=3".into(),
        hpc: true,
        f_max: 20,
        allow_collisions: false,
        collapse: Collapse::Gene,
        infer_genes: None,
    };
    let idx = Index::from_data(&build(txs, &p).unwrap().0).unwrap();
    assert_eq!(idx.tx_raw_codes(0, 0, 3), None);
    assert_eq!(idx.tx_raw_codes(1, 0, 3), Some(vec![0, 1, 2]));
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(150))]

    #[test]
    fn random_roundtrip(txs in transcriptome(), hpc in any::<bool>()) {
        for spec in ["lcp:K=3", "minimizer:k=15,w=10"] {
            check_roundtrip(&txs, spec, hpc, Collapse::Global);
            check_roundtrip(&txs, spec, hpc, Collapse::Gene);
        }
    }
}
