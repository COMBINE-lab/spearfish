//! `spearfish em`: spearfish's own EM quantification from a mapping TSV (`spearfish map` output, or
//! another mapper's alignments converted to the same columns). Output: transcript, length,
//! estimated reads.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result};
use paraseq::{ReaderBuilder, Record};
use quant::{EqBuilder, QuantParams, em};

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Transcript FASTA (names and lengths; IDs are cut at the first '|' or whitespace, as
    /// `spearfish index` cuts them).
    transcripts: PathBuf,
    /// Mapping TSV.
    mapping: PathBuf,
    /// Score denominator D in exp((s - s_max) / D). This EM's default (2) is what the
    /// mapper's ablations were evaluated with; oarfish uses D = 1 for spearfish's scores
    /// (`spearfish::SCORE_DENOM`) and 3 for minimap2's.
    #[arg(long, default_value_t = 2.0)]
    score_denom: f64,
    /// Keep hits scoring at least this fraction of the read's best.
    #[arg(long, default_value_t = 0.95)]
    min_frac: f64,
    /// Divide likelihoods by transcript length.
    #[arg(long)]
    length_norm: bool,
    /// Also use reads flagged `low_coverage` (best chain spans < half the read).
    #[arg(long)]
    keep_low_coverage: bool,
    /// EM: maximum iterations.
    #[arg(long, default_value_t = 10_000)]
    max_iter: usize,
    /// EM: stop when the largest relative change (transcripts with >= 0.01 reads) is below
    /// this.
    #[arg(long, default_value_t = 1e-6)]
    tol: f64,
    /// Output TSV.
    #[arg(short, long)]
    output: PathBuf,
    /// EM worker threads (0 = all cores).
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
}

pub fn run(args: Args) -> Result<()> {
    let t0 = Instant::now();
    let mut names: Vec<String> = Vec::new();
    let mut lens: Vec<u32> = Vec::new();
    let mut reader = ReaderBuilder::path(&args.transcripts).build()?;
    let mut rset = reader.new_record_set();
    while rset.fill(&mut reader)? {
        for rec in rset.iter() {
            let rec = rec?;
            names.push(spearfish::name_and_gene(rec.id_str()).0);
            lens.push(rec.seq().len() as u32);
        }
    }
    let id: HashMap<&str, u32> = names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i as u32))
        .collect();
    let p = QuantParams {
        score_denom: args.score_denom,
        min_frac: args.min_frac,
        length_norm: args.length_norm,
        max_iter: args.max_iter,
        tol: args.tol,
    };
    let mut b = EqBuilder::default();
    let f =
        std::fs::File::open(&args.mapping).with_context(|| format!("opening {}", args.mapping.display()))?;
    let (mut cur, mut hits): (String, Vec<(u32, f64)>) = (String::new(), Vec::new());
    let (mut unknown, mut unmapped, mut skipped) = (0u64, 0u64, 0u64);
    let mut cur_low = false;
    for line in BufReader::new(f).lines().skip(1) {
        let line = line?;
        let mut it = line.split('\t');
        let read = it.next().unwrap_or("");
        let tx = it.nth(3).unwrap_or("*");
        let score: f64 = it.next().unwrap_or("0").parse().unwrap_or(0.0);
        let flags = line.rsplit('\t').next().unwrap_or("");
        let low = flags.split(',').any(|f| f == "low_coverage");
        if read != cur {
            if !hits.is_empty() {
                b.add(&hits, &lens, &p);
            } else if cur_low {
                skipped += 1;
            }
            hits.clear();
            cur = read.to_string();
            cur_low = false;
        }
        if tx == "*" {
            unmapped += 1;
            continue;
        }
        if low && !args.keep_low_coverage {
            cur_low = true;
            continue;
        }
        match id.get(tx) {
            Some(&t) => hits.push((t, score)),
            None => unknown += 1,
        }
    }
    if !hits.is_empty() {
        b.add(&hits, &lens, &p);
    } else if cur_low {
        skipped += 1;
    }
    let reads = b.reads;
    let classes = b.finish();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;
    let (counts, iters) = pool.install(|| em(names.len(), &classes, &p));
    let mut out = std::io::BufWriter::new(std::fs::File::create(&args.output)?);
    writeln!(out, "transcript\tlength\treads")?;
    for ((n, l), c) in names.iter().zip(&lens).zip(&counts) {
        writeln!(out, "{n}\t{l}\t{c:.4}")?;
    }
    eprintln!(
        "quant: {reads} reads used ({unmapped} unmapped, {skipped} low-coverage skipped, {unknown} hits to unknown transcripts), {} classes, {iters} EM iterations, {:.2} s",
        classes.len(),
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}
