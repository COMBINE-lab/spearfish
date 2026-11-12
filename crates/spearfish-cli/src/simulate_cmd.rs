//! `spearfish simulate`: long reads from a transcriptome with known origin, for evaluation.
//!
//! Each read comes from a transcript drawn with probability proportional to abundance ×
//! length-independent weight (one molecule = one read). With probability `--truncate`, the
//! read starts at a uniform position in the first half of the transcript (5' truncation; the
//! 3' end is kept, as in poly(A)-primed cDNA). Errors come from the built-in simulator; half
//! the reads are reverse complemented. Read names: `r<i>|<transcript>|<strand>|<start>|<end>`.
//! `--truth` writes the true read count per transcript.
//!
//! Harder sets, all off by default: `--expressed` leaves most transcripts without reads (a
//! sample expresses a fraction of an annotation), `--fragment` makes reads of any part of a
//! transcript, `--rate-sd` spreads the error rate over reads, `--junk` adds random sequence
//! at both read ends (names gain `|<head>|<tail>`, the junk lengths in the read as written),
//! and `--unmappable` adds reads of random sequence, named `r<i>|*|*|0|0` and left out of the
//! truth: a mapper should report none of them.

use std::io::Write;
use std::path::PathBuf;

use anyhow::{Context, Result};
use lcp_core::encode::revcomp;

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
    /// Fraction of transcripts with any abundance.
    #[arg(long, default_value_t = 1.0)]
    expressed: f64,
    /// Probability that a read is a fragment: a uniform start and a uniform length of at
    /// least `--min-fragment` (instead of 5' truncation).
    #[arg(long, default_value_t = 0.0)]
    fragment: f64,
    /// Shortest fragment.
    #[arg(long, default_value_t = 100)]
    min_fragment: usize,
    /// Standard deviation of the per-read error rate around `--rate` (0: every read at `--rate`).
    #[arg(long, default_value_t = 0.0)]
    rate_sd: f64,
    /// Mean length of the random sequence added at each read end (0: none).
    #[arg(long, default_value_t = 0.0)]
    junk: f64,
    /// Fraction of the reads that are random sequence, with the lengths of real reads.
    #[arg(long, default_value_t = 0.0)]
    unmappable: f64,
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
    // Named as the index names them, so the names in reads and truth match its transcripts.
    let txs: Vec<(String, Vec<u8>)> = spearfish::read_transcripts(&args.transcripts)?
        .into_iter()
        .filter(|t| t.2.len() >= args.min_len)
        .map(|(name, _, seq)| (name, seq))
        .collect();
    anyhow::ensure!(!txs.is_empty(), "no transcripts of at least {} nt", args.min_len);
    for (name, x) in [
        ("--expressed", args.expressed),
        ("--fragment", args.fragment),
        ("--unmappable", args.unmappable),
    ] {
        anyhow::ensure!((0.0..=1.0).contains(&x), "{name} must be within 0..=1");
    }
    anyhow::ensure!(args.expressed > 0.0, "--expressed must be above 0");
    anyhow::ensure!(args.unmappable < 1.0, "--unmappable must be under 1");
    anyhow::ensure!(
        args.rate_sd >= 0.0 && args.junk >= 0.0,
        "--rate-sd and --junk must not be negative"
    );
    let mut rng = SplitMix64::for_record(args.seed, 0);
    // A standard normal via Box-Muller.
    let normal = |rng: &mut SplitMix64| {
        let (u1, u2) = (rng.unit().max(1e-12), rng.unit());
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    };
    // Log-normal abundances.
    let mut weights: Vec<f64> = txs
        .iter()
        .map(|_| (args.sigma * normal(&mut rng)).exp())
        .collect();
    // The options below draw from the generator only when set, so a run without them is
    // the same reads as before they existed.
    if args.expressed < 1.0 {
        for w in &mut weights {
            if rng.unit() >= args.expressed {
                *w = 0.0;
            }
        }
        anyhow::ensure!(
            weights.iter().any(|&w| w > 0.0),
            "--expressed left no transcript expressed"
        );
    }
    // Geometric with mean `args.junk`.
    let junk = |rng: &mut SplitMix64| -> Vec<u8> {
        let n = (-rng.unit().max(1e-12).ln() * args.junk) as usize;
        (0..n).map(|_| rng.base()).collect()
    };
    let mut cum = Vec::with_capacity(weights.len());
    let mut acc = 0.0;
    for w in &weights {
        acc += w;
        cum.push(acc);
    }
    let mut counts = vec![0u64; txs.len()];
    let mut out = std::io::BufWriter::new(std::fs::File::create(&args.output)?);
    let mut n_unmappable = 0usize;
    for i in 0..args.n {
        let x = rng.unit() * acc;
        // The first transcript whose cumulative weight reaches `x`: never one of weight zero.
        let t = cum.partition_point(|&c| c < x).min(txs.len() - 1);
        let (name, seq) = &txs[t];
        let (start, end) = if args.fragment > 0.0 && rng.unit() < args.fragment {
            let min = args.min_fragment.clamp(1, seq.len());
            let len = min + rng.below((seq.len() - min + 1) as u64) as usize;
            let start = rng.below((seq.len() - len + 1) as u64) as usize;
            (start, start + len)
        } else if rng.unit() < args.truncate {
            (rng.below((seq.len() / 2) as u64) as usize, seq.len())
        } else {
            (0, seq.len())
        };
        let rate = if args.rate_sd > 0.0 {
            (args.rate + args.rate_sd * normal(&mut rng)).clamp(0.0, 0.5)
        } else {
            args.rate
        };
        let read = simulate(&seq[start..end], rate, &prof, &mut rng);
        let (strand, mut s) = if rng.unit() < 0.5 {
            ('-', revcomp(&read.seq))
        } else {
            ('+', read.seq)
        };
        let mut ends = String::new();
        if args.junk > 0.0 {
            let (head, tail) = (junk(&mut rng), junk(&mut rng));
            ends = format!("|{}|{}", head.len(), tail.len());
            s = [head, s, tail].concat();
        }
        if args.unmappable > 0.0 && rng.unit() < args.unmappable {
            n_unmappable += 1;
            s.iter_mut().for_each(|b| *b = rng.base());
            writeln!(out, ">r{i}|*|*|0|0")?;
        } else {
            counts[t] += 1;
            writeln!(out, ">r{i}|{name}|{strand}|{start}|{end}{ends}")?;
        }
        out.write_all(&s)?;
        out.write_all(b"\n")?;
    }
    out.flush()?;
    if let Some(p) = &args.truth {
        let mut f = std::io::BufWriter::new(std::fs::File::create(p)?);
        writeln!(f, "transcript\tlength\treads\tweight")?;
        for ((name, seq), (c, w)) in txs.iter().zip(counts.iter().zip(&weights)) {
            writeln!(f, "{name}\t{}\t{c}\t{w:.6}", seq.len())?;
        }
        f.flush()?;
    }
    eprintln!(
        "simulate: {} reads from {} of {} transcripts, {n_unmappable} unmappable reads",
        args.n - n_unmappable,
        counts.iter().filter(|&&c| c > 0).count(),
        txs.len()
    );
    Ok(())
}
