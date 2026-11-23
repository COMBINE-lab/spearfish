//! Library API of the genome mode: build a genome index in memory and map the toy
//! transcripts to it, spliced, on both strands, with and without the annotation.

use std::collections::HashMap;
use std::path::PathBuf;

use lcp_core::encode::revcomp;
use spearfish::genome::{
    BaseLevel, CigarOp, GenomeIndex, GenomeIndexOptions, SpliceOptions, SplicedMapper, TranscriptModel,
    TranscriptStrand, output,
};
use spearfish::{IndexKind, Preset, read_transcripts};

fn toy(file: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/toy")
        .join(file)
}

/// The toy annotation's transcripts as sorted, 0-based, half-open exons.
fn annotated_exons() -> HashMap<String, Vec<(u32, u32)>> {
    let mut m: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    for line in std::fs::read_to_string(toy("annotation.gtf")).unwrap().lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() == 9 && p[2] == "exon" {
            let id = p[8]
                .split("transcript_id \"")
                .nth(1)
                .unwrap()
                .split('"')
                .next()
                .unwrap();
            m.entry(id.to_string())
                .or_default()
                .push((p[3].parse::<u32>().unwrap() - 1, p[4].parse().unwrap()));
        }
    }
    m.values_mut().for_each(|v| v.sort_unstable());
    m
}

fn assert_send_sync<T: Send + Sync>() {}

#[test]
fn map_toy_transcripts_to_the_genome() {
    assert_send_sync::<SplicedMapper>();
    let truth = annotated_exons();
    let txs = read_transcripts(&toy("transcripts.fa")).unwrap();
    let plain = GenomeIndex::build_from_fasta(&toy("genome.fa"), &GenomeIndexOptions::default()).unwrap();
    assert_eq!(plain.n_sequences(), 1);
    assert_eq!(plain.seq_name(0), "toy");
    assert_eq!(plain.total_len(), plain.seq_len(0) as u64);
    let mut annotated = plain.clone();
    // The clone shares the index: an annotation goes into an index nobody else holds.
    assert!(annotated.load_annotation(&toy("annotation.gtf")).is_err());
    drop(annotated);
    let mut annotated =
        GenomeIndex::build_from_fasta(&toy("genome.fa"), &GenomeIndexOptions::default()).unwrap();
    let sum = annotated.load_annotation(&toy("annotation.gtf")).unwrap();
    assert_eq!((sum.format, sum.transcripts), ("gtf", txs.len()));

    for (known, index) in [(false, plain), (true, annotated)] {
        let mapper = SplicedMapper::new(index, &SpliceOptions::preset(Preset::Ont));
        for (name, _, seq) in &txs {
            let exons = &truth[name];
            for rev in [false, true] {
                let read = if rev { revcomp(seq) } else { seq.clone() };
                let r = mapper.map(&read);
                assert_eq!(r.len as usize, read.len());
                let a = &r.alignments[0];
                assert!(a.primary && a.reverse == rev, "{name} rev {rev}");
                // An exact read of an annotated transcript: its exons, base for base, when
                // the annotation is given. Without it the toy genome leaves a junction free
                // to slide: its introns are random sequence, with a GT or an AG only by
                // chance, so the nearest motif a few bases away wins.
                let blocks: Vec<(u32, u32)> = a.blocks().map(|b| (b.target_start, b.target_end)).collect();
                assert_eq!(blocks.len(), exons.len(), "{name} rev {rev}");
                if known {
                    assert_eq!(&blocks, exons, "{name} rev {rev}");
                }
                for (b, e) in blocks.iter().zip(exons) {
                    assert!(
                        b.0.abs_diff(e.0) <= 8 && b.1.abs_diff(e.1) <= 8,
                        "{name} rev {rev}: {b:?} {e:?}"
                    );
                }
                assert_eq!(
                    (a.target_start, a.target_end),
                    (exons[0].0, exons[exons.len() - 1].1)
                );
                assert_eq!((a.query_start, a.query_end), (0, r.len));
                assert_eq!(a.transcript_strand, TranscriptStrand::Forward);
                assert_eq!(a.introns().count(), exons.len() - 1);
                assert!(a.spliced());
                // One M per exon and one N per intron.
                let cigar = a.cigar(r.len);
                if known {
                    assert_eq!(cigar.len(), 2 * exons.len() - 1);
                }
                assert!(matches!(cigar[1], CigarOp::Intron(n) if n == blocks[1].0 - blocks[0].1));
                let mut paf = Vec::new();
                output::paf_lines(&mut paf, name, &r, mapper.index());
                let paf = String::from_utf8(paf).unwrap();
                assert!(
                    paf.starts_with(&format!("{name}\t{}\t0\t{}\t", r.len, r.len)),
                    "{paf}"
                );
                assert!(paf.contains("\ttoy\t") && paf.contains("\tts:A:+"), "{paf}");
            }
        }
    }
}

#[test]
fn direct_rna_reads_map_in_sense_only_and_base_level_alignments_are_exact() {
    let txs = read_transcripts(&toy("transcripts.fa")).unwrap();
    let mut index = GenomeIndex::build_from_fasta(&toy("genome.fa"), &GenomeIndexOptions::default()).unwrap();
    index.load_annotation(&toy("annotation.gtf")).unwrap();
    let mut opts = SpliceOptions::preset(Preset::Drna);
    assert!(opts.sense_reads);
    opts.base_level = BaseLevel::Affine;
    opts.max_hits = 3;
    let mapper = SplicedMapper::new(index, &opts);
    let (name, _, seq) = &txs[0];
    let r = mapper.map(seq);
    let a = &r.alignments[0];
    assert_eq!(a.edits, Some(0), "{name}");
    assert!(
        a.cigar(r.len)
            .iter()
            .all(|op| matches!(op, CigarOp::Match(_) | CigarOp::Intron(_)))
    );
    let (mut sam, mut bam, mut z) = (Vec::new(), Vec::new(), Vec::new());
    output::sam_header(&mut sam, mapper.index(), "test");
    let qual = vec![b'I'; seq.len()];
    output::sam_records(&mut sam, name, seq, Some(&qual), &r, mapper.index());
    output::sam_unmapped(&mut sam, "lost", b"ACGT", Some(b"!!II"));
    let sam = String::from_utf8(sam).unwrap();
    assert!(
        sam.contains("@SQ\tSN:toy\tLN:18090") && sam.contains("\tCL:test\n"),
        "{sam}"
    );
    assert!(sam.contains("\tNM:i:0") && sam.contains("\tXS:A:+"), "{sam}");
    assert!(
        sam.contains(&format!("\t{}\tAS:i:", "I".repeat(seq.len()))),
        "qualities"
    );
    assert!(
        sam.ends_with("lost\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!II\n"),
        "{sam}"
    );
    output::bam_header(&mut bam, mapper.index(), "test");
    output::bam_records(&mut bam, name, seq, Some(&qual), &r, mapper.index());
    output::bam_unmapped(&mut bam, "lost", b"ACGT", None);
    output::bgzf(&bam, 1, &mut z);
    assert_eq!(&z[..4], &[0x1f, 0x8b, 8, 4]);
    // A sense read that maps as its reverse complement is of a reverse-strand transcript:
    // the antisense copy of this read is never placed on the forward-strand transcript.
    let r = mapper.map(&revcomp(seq));
    assert!(
        r.alignments
            .iter()
            .all(|a| a.reverse && a.transcript_strand != TranscriptStrand::Forward),
        "{:?}",
        r.alignments
    );
}

#[test]
fn index_files_say_what_they_are() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("genome-api");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("toy.sfg");
    let (index, bytes, manifest) = GenomeIndex::build_to_file(
        &toy("genome.fa"),
        &GenomeIndexOptions::default(),
        &path,
        vec!["test".into()],
    )
    .unwrap();
    assert!(bytes > 0);
    assert_eq!(manifest.kind, "spearfish-genome-index");
    assert_eq!((manifest.sequences, manifest.bases), (1, 18090));
    assert_eq!(spearfish::index_kind(&path).unwrap(), IndexKind::Genome);
    let opened = GenomeIndex::open(&path).unwrap();
    assert_eq!(opened.manifest().unwrap(), manifest);
    assert_eq!(opened.id(), manifest.index_id);
    assert_eq!(index.n_sequences(), opened.n_sequences());
    assert!(dir.join("toy.sfg.json").exists());
    // Each kind of index refuses the other's file, and says what the file is.
    let e = spearfish::Index::open(&path).err().unwrap().to_string();
    assert!(e.contains("is a genome index"), "{e}");
    let tx = dir.join("toy.sfi");
    spearfish::Index::build_to_file(
        read_transcripts(&toy("transcripts.fa")).unwrap(),
        &Default::default(),
        &tx,
    )
    .unwrap();
    assert_eq!(spearfish::index_kind(&tx).unwrap(), IndexKind::Transcriptome);
    let e = GenomeIndex::open(&tx).err().unwrap().to_string();
    assert!(e.contains("is a transcriptome index"), "{e}");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn annotation_from_models_and_bad_annotations() {
    let mut index = GenomeIndex::build_from_fasta(&toy("genome.fa"), &GenomeIndexOptions::default()).unwrap();
    let models: Vec<TranscriptModel> = annotated_exons()
        .into_values()
        .map(|exons| TranscriptModel {
            seq: 0,
            reverse: false,
            exons,
        })
        .collect();
    let n = models.len();
    let sum = index.set_annotation(models).unwrap().clone();
    assert_eq!(sum.transcripts, n);
    assert_eq!(index.annotation().unwrap().introns, 11);
    // A model off the end of the sequence is an error, not a silent skip.
    let off = TranscriptModel {
        seq: 0,
        reverse: false,
        exons: vec![(18_000, 19_000)],
    };
    assert!(index.set_annotation([off]).is_err());
    // An annotation of another assembly (no sequence name in common) is an error that names
    // both.
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("genome-ann");
    std::fs::create_dir_all(&dir).unwrap();
    let gtf = dir.join("other.gtf");
    std::fs::write(&gtf, "chr1\tx\texon\t10\t50\t.\t+\t.\ttranscript_id \"t1\";\n").unwrap();
    let e = index.load_annotation(&gtf).err().unwrap().to_string();
    assert!(e.contains("chr1") && e.contains("toy"), "{e}");
    std::fs::remove_dir_all(&dir).unwrap();
}
