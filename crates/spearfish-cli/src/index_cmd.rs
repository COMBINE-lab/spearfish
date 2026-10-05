//! `spearfish index`: build the core-graph index from a transcript FASTA, and
//! `spearfish index-stats`: report on a built index.

use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::Result;
use index::{BuildParams, Collapse, Index, TxInput, build};
use paraseq::{ReaderBuilder, Record};

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Transcript FASTA (optionally compressed). GENCODE headers (`tx|gene|...`) give gene
    /// IDs; otherwise use --t2g or --infer-genes, or each transcript is its own gene.
    transcripts: PathBuf,
    /// Transcript-to-gene table: whitespace-separated `transcript gene` per line (further
    /// columns ignored). Overrides gene IDs from the FASTA headers.
    #[arg(long, conflicts_with = "infer_genes")]
    t2g: Option<PathBuf>,
    /// Group transcripts into genes from their sequence, for transcriptomes without a gene
    /// annotation: transcripts sharing a run of --infer-span consecutive seeds are grouped,
    /// ignoring runs held by more than N transcripts (repeats).
    #[arg(long, value_name = "N", num_args = 0..=1, default_missing_value = "64")]
    infer_genes: Option<u32>,
    /// Gene inference: consecutive seeds a shared run must span (12 level-3 seeds is about
    /// 300 bases of shared sequence).
    #[arg(long, default_value_t = 12)]
    infer_span: usize,
    /// Gene inference, attach pass: shorter runs (in seeds; 0 = off) by which small groups join
    /// the one larger group they point to.
    #[arg(long, default_value_t = 4)]
    infer_attach_span: usize,
    /// Gene inference, attach pass: largest group (in transcripts) that may attach.
    #[arg(long, default_value_t = 3)]
    infer_attach_max: usize,
    /// Write the transcript-to-gene grouping used (given or inferred) to this TSV.
    #[arg(long)]
    write_t2g: Option<PathBuf>,
    /// Base seed scheme (its highest level is used).
    #[arg(short, long, default_value = "lcp:K=3")]
    scheme: String,
    /// Homopolymer-compress transcripts (and, at mapping time, reads).
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    hpc: bool,
    /// Graph nodes: one per label (`global`) or per (label, gene) (`gene`).
    #[arg(long, default_value = "gene", value_parser = ["global", "gene"])]
    collapse: String,
    /// Mask unitigs whose color touches more than this many genes.
    #[arg(long, default_value_t = 20)]
    f_max: u32,
    /// Build even if two different seed contents share a label (normally an error).
    #[arg(long)]
    allow_collisions: bool,
    /// Worker threads (0 = all cores).
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
    /// Output index file.
    #[arg(short, long)]
    output: PathBuf,
}

#[derive(clap::Args, Debug)]
pub struct StatsArgs {
    index: PathBuf,
}

pub fn run(args: Args) -> Result<()> {
    let t0 = Instant::now();
    let mut txs = Vec::new();
    let mut reader = ReaderBuilder::path(&args.transcripts).build()?;
    let mut rset = reader.new_record_set();
    while rset.fill(&mut reader)? {
        for rec in rset.iter() {
            let rec = rec?;
            let id = rec.id_str();
            let mut f = id.split('|');
            let name = f
                .next()
                .unwrap()
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            let gene = f.next().map_or_else(|| name.clone(), String::from);
            txs.push(TxInput {
                name,
                gene,
                seq: rec.seq().into_owned(),
            });
        }
    }
    let mut missing = 0usize;
    if let Some(path) = &args.t2g {
        let mut map = std::collections::HashMap::new();
        for line in std::fs::read_to_string(path)?.lines() {
            let mut f = line.split_whitespace();
            if let (Some(t), Some(g)) = (f.next(), f.next()) {
                map.insert(t.to_string(), g.to_string());
            }
        }
        for t in &mut txs {
            match map.get(&t.name) {
                Some(g) => t.gene = g.clone(),
                None => missing += 1,
            }
        }
    }
    let read_s = t0.elapsed().as_secs_f64();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;
    let params = BuildParams {
        spec: args.scheme.clone(),
        hpc: args.hpc,
        f_max: args.f_max,
        allow_collisions: args.allow_collisions,
        collapse: if args.collapse == "gene" {
            Collapse::Gene
        } else {
            Collapse::Global
        },
        infer_genes: args.infer_genes.map(|max_tx| index::InferGenes {
            span: args.infer_span,
            max_tx,
            attach_span: args.infer_attach_span,
            attach_max: args.infer_attach_max,
        }),
    };
    let (data, rep) = pool.install(|| build(&txs, &params))?;
    if let Some(path) = &args.write_t2g {
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
        for (t, g) in txs.iter().zip(&rep.tx_gene) {
            writeln!(w, "{}\t{}", t.name, rep.gene_names[*g as usize])?;
        }
    }
    drop(txs);
    let t = Instant::now();
    let bytes = data.write(&args.output)?;
    let write_s = t.elapsed().as_secs_f64();

    let mut e = std::io::stderr().lock();
    writeln!(
        e,
        "index: {} ({}, hpc={})",
        args.output.display(),
        args.scheme,
        args.hpc
    )?;
    writeln!(
        e,
        "  transcripts {} ({} unseeded), genes {}",
        rep.transcripts, rep.unseeded_transcripts, rep.genes
    )?;
    if let Some(n) = args.infer_genes {
        writeln!(
            e,
            "  genes inferred from shared runs of {} seeds (in at most {n} transcripts); largest group {} transcripts",
            args.infer_span, rep.largest_group_tx
        )?;
    }
    if missing > 0 {
        writeln!(
            e,
            "  warning: {missing} transcripts not in --t2g; their gene IDs come from the FASTA"
        )?;
    }
    if args.infer_genes.is_none()
        && args.collapse == "gene"
        && rep.transcripts > 1
        && rep.genes == rep.transcripts
    {
        writeln!(
            e,
            "  warning: every transcript is its own gene (no gene IDs in the headers): isoforms are \
             not collapsed, so the index is larger and mapping less accurate; give --t2g or --infer-genes"
        )?;
    }
    writeln!(
        e,
        "  seed occurrences {}, nodes {}",
        rep.seed_occurrences, rep.nodes
    )?;
    writeln!(
        e,
        "  unitigs {} ({} cyclic), N50 {} seeds; occurrences {}",
        rep.unitigs, rep.cyclic_unitigs, rep.unitig_n50_seeds, rep.occurrences
    )?;
    writeln!(
        e,
        "  colors {}, masked unitigs {}, families {} (largest {} tx, {} multi-gene)",
        rep.colors, rep.masked_unitigs, rep.families, rep.largest_family_tx, rep.multi_gene_families
    )?;
    writeln!(e, "  label collisions {}", rep.label_collisions)?;
    writeln!(
        e,
        "  unitig heads: transcript start only {}, several predecessors {}, predecessor branches {}",
        rep.heads_start, rep.heads_multi_pred, rep.heads_pred_branches
    )?;
    if rep.inconsistent_lengths > 0 {
        writeln!(
            e,
            "  warning: {} seeds share a label with a different span",
            rep.inconsistent_lengths
        )?;
    }
    let phases: Vec<String> = rep.seconds.iter().map(|(k, s)| format!("{k} {s:.1}s")).collect();
    writeln!(
        e,
        "  time: read {read_s:.1}s, {}, write {write_s:.1}s",
        phases.join(", ")
    )?;
    writeln!(e, "  size {:.1} MB", bytes as f64 / 1e6)?;
    Ok(())
}

pub fn stats(args: StatsArgs) -> Result<()> {
    let t = Instant::now();
    let idx = Index::open(&args.index)?;
    let load = t.elapsed().as_secs_f64();
    let mut out = std::io::stdout().lock();
    writeln!(out, "key\tvalue")?;
    for k in [
        "spec",
        "hpc",
        "collapse",
        "f_max",
        "lcp_parser_version",
        "input_xxh3",
    ] {
        writeln!(out, "{k}\t{}", idx.meta(k).unwrap_or("-"))?;
    }
    writeln!(out, "transcripts\t{}", idx.n_transcripts())?;
    writeln!(out, "nodes\t{}", idx.n_nodes())?;
    writeln!(out, "unitigs\t{}", idx.n_unitigs())?;
    writeln!(out, "occurrences\t{}", idx.occ_tx().len())?;
    writeln!(out, "colors\t{}", idx.color_off().len() - 1)?;
    writeln!(out, "load_seconds\t{load:.4}")?;
    let mut total = 0;
    for (name, bytes) in idx.sizes() {
        total += bytes;
        writeln!(out, "bytes.{name}\t{bytes}")?;
    }
    writeln!(out, "bytes.total\t{total}")?;
    Ok(())
}
