//! `spearfish quant`: EM quantification from a mapping TSV (`spearfish map` output, or
//! minimap2 converted by `bench/scripts/sam_to_maptsv.py`). Output: transcript, length,
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
    /// Transcript FASTA (names and lengths; GENCODE `tx|...` headers are cut at '|').
    transcripts: PathBuf,
    /// Mapping TSV.
    mapping: PathBuf,
    /// Score denominator D in exp((s - s_max) / D). 2 suits spearfish chain scores; use 5
    /// (oarfish's default) for minimap2 alignment scores.
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
    /// Output TSV.
    #[arg(short, long)]
    output: PathBuf,
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
            names.push(rec.id_str().split('|').next().unwrap().to_string());
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
        ..Default::default()
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
        let low = line.ends_with("low_coverage");
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
    let (counts, iters) = em(names.len(), &classes, &p);
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
