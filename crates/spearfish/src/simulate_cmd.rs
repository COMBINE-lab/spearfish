//! `spearfish simulate`: long reads from a transcriptome with known origin, for evaluation.
//!
//! Each read comes from a transcript drawn with probability proportional to abundance ×
//! length-independent weight (one molecule = one read). With probability `--truncate`, the
//! read starts at a uniform position in the first half of the transcript (5' truncation; the
//! 3' end is kept, as in poly(A)-primed cDNA). Errors come from the built-in simulator; half
//! the reads are reverse complemented. Read names: `r<i>|<transcript>|<strand>|<start>|<end>`.
//! `--truth` writes the true read count per transcript.

use std::io::Write;
use std::path::PathBuf;

use anyhow::{Context, Result};
use lcp_core::encode::revcomp;
use paraseq::{ReaderBuilder, Record};

use crate::rng::SplitMix64;
use crate::sim::{profile, simulate};

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Transcript FASTA.
    transcripts: PathBuf,
    /// Number of reads.
    #[arg(short, long, default_value_t = 10_000)]
    n: usize,
    /// Error profile: uniform, ont, hifi.
    #[arg(long, default_value = "ont")]
    profile: String,
    /// Per-base error rate.
    #[arg(long, default_value_t = 0.06)]
    rate: f64,
    /// Probability of 5' truncation.
    #[arg(long, default_value_t = 0.5)]
    truncate: f64,
    /// Skip transcripts shorter than this.
    #[arg(long, default_value_t = 300)]
    min_len: usize,
    /// Abundances: uniform, or log-normal with this sigma.
    #[arg(long, default_value_t = 1.5)]
    sigma: f64,
    #[arg(long, default_value_t = 1)]
    seed: u64,
    /// Output FASTA.
    #[arg(short, long)]
    output: PathBuf,
    /// Output TSV of true counts per transcript.
    #[arg(long)]
    truth: Option<PathBuf>,
}

pub fn run(args: Args) -> Result<()> {
    let prof = profile(&args.profile).with_context(|| format!("unknown profile {}", args.profile))?;
    let mut txs: Vec<(String, Vec<u8>)> = Vec::new();
    let mut reader = ReaderBuilder::path(&args.transcripts).build()?;
    let mut rset = reader.new_record_set();
    while rset.fill(&mut reader)? {
        for rec in rset.iter() {
            let rec = rec?;
            let seq = rec.seq().into_owned();
            if seq.len() >= args.min_len {
                txs.push((rec.id_str().split('|').next().unwrap().to_string(), seq));
            }
        }
    }
    anyhow::ensure!(!txs.is_empty(), "no transcripts of at least {} nt", args.min_len);
    let mut rng = SplitMix64::for_record(args.seed, 0);
    // Log-normal abundances via Box-Muller.
    let weights: Vec<f64> = txs
        .iter()
        .map(|_| {
            let (u1, u2) = (rng.unit().max(1e-12), rng.unit());
            let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            (args.sigma * z).exp()
        })
        .collect();
    let mut cum = Vec::with_capacity(weights.len());
    let mut acc = 0.0;
    for w in &weights {
        acc += w;
        cum.push(acc);
    }
    let mut counts = vec![0u64; txs.len()];
    let mut out = std::io::BufWriter::new(std::fs::File::create(&args.output)?);
    for i in 0..args.n {
        let x = rng.unit() * acc;
        let t = cum.partition_point(|&c| c < x).min(txs.len() - 1);
        counts[t] += 1;
        let (name, seq) = &txs[t];
        let start = if rng.unit() < args.truncate {
            rng.below((seq.len() / 2) as u64) as usize
        } else {
            0
        };
        let read = simulate(&seq[start..], args.rate, &prof, &mut rng);
        let (strand, s) = if rng.unit() < 0.5 {
            ('-', revcomp(&read.seq))
        } else {
            ('+', read.seq)
        };
        writeln!(out, ">r{i}|{name}|{strand}|{start}|{}", seq.len())?;
        out.write_all(&s)?;
        out.write_all(b"\n")?;
    }
    if let Some(p) = &args.truth {
        let mut f = std::io::BufWriter::new(std::fs::File::create(p)?);
        writeln!(f, "transcript\tlength\treads\tweight")?;
        for ((name, seq), (c, w)) in txs.iter().zip(counts.iter().zip(&weights)) {
            writeln!(f, "{name}\t{}\t{c}\t{w:.6}", seq.len())?;
        }
    }
    eprintln!("simulate: {} reads from {} transcripts", args.n, txs.len());
    Ok(())
}
