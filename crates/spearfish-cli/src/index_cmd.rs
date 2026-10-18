//! `spearfish index`: build the core-graph index from a transcript FASTA, and
//! `spearfish index-stats`: report on a built index.

use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::Result;
use clap::ValueEnum;
use spearfish::provenance::{InputInfo, file_xxh3};
use spearfish::{Index, IndexOptions, InferGenesOptions};

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum CollapseArg {
    /// One graph node per (label, gene).
    Gene,
    /// One graph node per label.
    Global,
}

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Transcript FASTA (optionally compressed). GENCODE headers (`tx|gene|...`) give gene
    /// IDs; otherwise use --t2g or --infer-genes, or each transcript is its own gene.
    transcripts: PathBuf,
    /// Output index file (its manifest is also written as <FILE>.json).
    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,
    /// Worker threads (0: every available CPU).
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
    /// Transcript-to-gene table: whitespace-separated `transcript gene` per line (further
    /// columns ignored). Overrides gene IDs from the FASTA headers.
    #[arg(long, conflicts_with = "infer_genes", help_heading = "Genes")]
    t2g: Option<PathBuf>,
    /// Group transcripts into genes from their sequence, for transcriptomes without a gene
    /// annotation: transcripts sharing a run of --infer-span consecutive seeds are grouped,
    /// ignoring runs held by more than N transcripts (repeats).
    #[arg(
        long,
        value_name = "N",
        num_args = 0..=1,
        default_missing_value = "64",
        help_heading = "Genes"
    )]
    infer_genes: Option<u32>,
    /// Gene inference: consecutive seeds a shared run must span (12 level-3 seeds is about
    /// 300 bases of shared sequence).
    #[arg(long, default_value_t = 12, help_heading = "Genes")]
    infer_span: usize,
    /// Gene inference, attach pass: shorter runs (in seeds; 0 = off) by which small groups join
    /// the one larger group they point to.
    #[arg(long, default_value_t = 4, help_heading = "Genes")]
    infer_attach_span: usize,
    /// Gene inference, attach pass: largest group (in transcripts) that may attach.
    #[arg(long, default_value_t = 3, help_heading = "Genes")]
    infer_attach_max: usize,
    /// Write the transcript-to-gene grouping used (given or inferred) to this TSV.
    #[arg(long, value_name = "FILE", help_heading = "Genes")]
    write_t2g: Option<PathBuf>,
    /// Base seed scheme (its highest level is used).
    #[arg(short, long, default_value = "lcp:K=3", help_heading = "Index")]
    scheme: String,
    /// Don't homopolymer-compress transcripts (and reads at mapping time).
    #[arg(long, help_heading = "Index")]
    no_hpc: bool,
    /// Graph nodes: one per (label, gene) or one per label.
    #[arg(long, value_enum, default_value = "gene", help_heading = "Index")]
    collapse: CollapseArg,
    /// Mask unitigs whose color touches more than this many genes.
    #[arg(long, alias = "f-max", default_value_t = 20, help_heading = "Index")]
    mask_genes: u32,
    /// Build even if two different seed contents share a label (normally an error).
    #[arg(long, hide = true)]
    allow_collisions: bool,
}

#[derive(clap::Args, Debug)]
pub struct InspectArgs {
    /// Index built by `spearfish index`.
    index: PathBuf,
    /// Print the manifest and summary as JSON.
    #[arg(long)]
    json: bool,
    /// Also list every section with its size and xxh3-64 checksum (reads the whole index).
    #[arg(long)]
    sections: bool,
}

pub fn run(args: Args) -> Result<()> {
    let t0 = Instant::now();
    let mut txs = spearfish::read_transcripts(&args.transcripts)?;
    let mut input = InputInfo {
        path: args.transcripts.display().to_string(),
        genes_from: "fasta".into(),
        ..Default::default()
    };
    let mut missing = 0usize;
    if let Some(path) = &args.t2g {
        let table =
            std::fs::read_to_string(path).map_err(|e| anyhow::anyhow!("reading {}: {e}", path.display()))?;
        missing = spearfish::apply_t2g(&mut txs, &table);
        input.genes_from = "t2g".into();
        input.t2g = Some((path.display().to_string(), file_xxh3(path)?));
    }
    if args.infer_genes.is_some() {
        input.genes_from = "inferred".into();
    } else if txs.iter().all(|(n, g, _)| n == g) && args.t2g.is_none() {
        input.genes_from = "none".into();
    }
    let read_s = t0.elapsed().as_secs_f64();
    let opts = IndexOptions::builder()
        .scheme(args.scheme.clone())
        .hpc(!args.no_hpc)
        .f_max(args.mask_genes)
        .per_gene(args.collapse == CollapseArg::Gene)
        .allow_collisions(args.allow_collisions)
        .maybe_infer_genes(args.infer_genes.map(|max_tx| {
            InferGenesOptions::builder()
                .max_tx(max_tx)
                .span(args.infer_span)
                .attach_span(args.infer_attach_span)
                .attach_max(args.infer_attach_max)
                .build()
        }))
        .build();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(crate::runinfo::threads(args.threads))
        .build()?;
    let names: Vec<String> = if args.write_t2g.is_some() {
        txs.iter().map(|t| t.0.clone()).collect()
    } else {
        Vec::new()
    };
    let t = Instant::now();
    let (bytes, rep, manifest) = pool
        .install(|| Index::build_to_file_with(txs, &opts, &args.output, input, std::env::args().collect()))?;
    let build_write_s = t.elapsed().as_secs_f64();
    if let Some(path) = &args.write_t2g {
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
        for (t, g) in names.iter().zip(&rep.tx_gene) {
            writeln!(w, "{t}\t{}", rep.gene_names[*g as usize])?;
        }
    }
    let hpc = opts.hpc;
    let collapse_gene = opts.per_gene;

    let mut e = std::io::stderr().lock();
    if log::max_level() < log::LevelFilter::Info {
        return Ok(());
    }
    writeln!(
        e,
        "index: {} ({}, hpc={hpc}, id {})",
        args.output.display(),
        args.scheme,
        manifest.index_id
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
    if args.infer_genes.is_none() && collapse_gene && rep.transcripts > 1 && rep.genes == rep.transcripts {
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
        "  time: read {read_s:.1}s, {} (build and write {build_write_s:.1}s)",
        phases.join(", ")
    )?;
    writeln!(e, "  size {:.1} MB", bytes as f64 / 1e6)?;
    Ok(())
}

pub fn inspect(args: InspectArgs) -> Result<()> {
    let t = Instant::now();
    let idx = Index::open(&args.index)?;
    let load = t.elapsed().as_secs_f64();
    let raw = idx.raw();
    let meta: serde_json::Map<String, serde_json::Value> = raw
        .meta_text()
        .lines()
        .filter_map(|l| l.split_once('='))
        .map(|(k, v)| (k.to_string(), v.into()))
        .collect();
    let sizes = raw.sizes();
    let total: usize = sizes.iter().map(|(_, b)| b).sum();
    let compatible = raw.check_compatible();
    let mut j = serde_json::json!({
        "path": args.index.display().to_string(),
        "index_id": idx.id(),
        "format_version": index::store::FORMAT_VERSION,
        "compatible": compatible.is_ok(),
        "meta": meta,
        "transcripts": idx.n_transcripts(),
        "nodes": raw.n_nodes(),
        "unitigs": raw.n_unitigs(),
        "occurrences": raw.occ_tx().len(),
        "colors": raw.color_off().len() - 1,
        "raw_sequence": idx.has_raw_sequence(),
        "bytes": total,
        "load_seconds": load,
        "manifest": idx.manifest(),
    });
    if args.sections {
        j["sections"] = raw
            .section_checksums()
            .into_iter()
            .map(|(n, b, h)| serde_json::json!({"name": n, "bytes": b, "xxh3": format!("{h:016x}")}))
            .collect();
    }
    let mut out = std::io::stdout().lock();
    if args.json {
        writeln!(out, "{}", serde_json::to_string_pretty(&j)?)?;
        return Ok(());
    }
    writeln!(out, "index\t{}", args.index.display())?;
    writeln!(out, "index_id\t{}", idx.id())?;
    match &compatible {
        Ok(()) => writeln!(out, "compatible\tyes")?,
        Err(e) => writeln!(out, "compatible\tno: {e}")?,
    }
    for (k, v) in &meta {
        writeln!(out, "{k}\t{}", v.as_str().unwrap_or(""))?;
    }
    match idx.manifest() {
        Some(m) => {
            writeln!(out, "built_by\tspearfish {}", m.built_by.short())?;
            writeln!(out, "created_unix\t{}", m.created_unix)?;
            writeln!(out, "input\t{}", m.input.path)?;
            writeln!(out, "genes_from\t{}", m.input.genes_from)?;
            writeln!(out, "genes\t{}", m.input.genes)?;
        }
        None => writeln!(out, "manifest\t- (built before manifests were recorded)")?,
    }
    for k in [
        "transcripts",
        "nodes",
        "unitigs",
        "occurrences",
        "colors",
        "raw_sequence",
        "bytes",
    ] {
        writeln!(out, "{k}\t{}", j[k])?;
    }
    writeln!(out, "load_seconds\t{load:.4}")?;
    if args.sections {
        for (n, b, h) in raw.section_checksums() {
            writeln!(out, "section.{n}\t{b}\t{h:016x}")?;
        }
    } else {
        for (name, bytes) in sizes {
            writeln!(out, "bytes.{name}\t{bytes}")?;
        }
    }
    Ok(())
}
