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
    /// IDs; otherwise each transcript is its own gene.
    transcripts: PathBuf,
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
    };
    let (data, rep) = pool.install(|| build(&txs, &params))?;
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
