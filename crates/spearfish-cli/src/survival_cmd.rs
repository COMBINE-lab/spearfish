//! `spearfish survival`: how many seeds survive sequencing errors, per seed scheme, level,
//! error profile, error rate and HPC setting (milestone M2, the go/no-go gate).
//!
//! Each input sequence is treated as a full-length read: it's mutated by the simulator,
//! both versions are seeded, and a transcript seed survives if the read has a seed with the
//! same label within `--tol` raw bases of the seed's mapped position. Reported per row:
//! density, span, survival, coverage of the transcript by surviving seeds, the fraction of
//! read seeds that match a transcript seed (precision), and gene-level specificity (the
//! fraction of transcript seed occurrences whose label occurs in a single gene).

use std::io::Write;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use paraseq::prelude::*;
use paraseq::{ReaderBuilder, Record};
use rustc_hash::FxHashMap;

use crate::rng::SplitMix64;
use crate::seeding::{SeedSet, UnitSeeds, has_near};
use crate::sim::{ErrorProfile, profile, simulate};

#[derive(clap::Args, Debug)]
pub struct Args {
    /// FASTA/FASTQ transcripts (optionally compressed); `-` or omitted reads stdin.
    input: Option<PathBuf>,
    /// Seed scheme spec (repeatable), e.g. `lcp:K=6` or `minimizer:k=15,w=10`.
    #[arg(short, long = "scheme")]
    schemes: Vec<String>,
    /// File with one seed scheme spec per line (`#` comments allowed).
    #[arg(long)]
    schemes_file: Option<PathBuf>,
    /// Error profiles: uniform, ont, hifi.
    #[arg(long, value_delimiter = ',', default_value = "ont,hifi")]
    profiles: Vec<String>,
    /// Per-base error rates.
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "0.01,0.02,0.04,0.06,0.08,0.1,0.15"
    )]
    rates: Vec<f64>,
    /// HPC settings to evaluate: off, on.
    #[arg(long, value_delimiter = ',', default_value = "off,on")]
    hpc: Vec<String>,
    /// Position tolerance (raw bases) when matching a read seed to a transcript seed.
    #[arg(long, default_value_t = 3)]
    tol: u32,
    /// Skip sequences shorter than this.
    #[arg(long, default_value_t = 300)]
    min_len: usize,
    /// Fraction of sequences to use (deterministic by record index).
    #[arg(long, default_value_t = 1.0)]
    sample: f64,
    /// RNG seed.
    #[arg(long, default_value_t = 1)]
    seed: u64,
    /// Worker threads (0 = all cores).
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
    /// Output TSV (default stdout).
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Clone, Default)]
struct Row {
    seqs: u64,
    bases: u64,
    seeds: u64,
    survived: u64,
    span_sum: u64,
    covered: u64,
    read_seeds: u64,
    read_matched: u64,
}

impl Row {
    fn merge(&mut self, o: &Row) {
        self.seqs += o.seqs;
        self.bases += o.bases;
        self.seeds += o.seeds;
        self.survived += o.survived;
        self.span_sum += o.span_sum;
        self.covered += o.covered;
        self.read_seeds += o.read_seeds;
        self.read_matched += o.read_matched;
    }
}

const MULTI: u64 = u64::MAX;

#[derive(Clone)]
struct Stats {
    rows: Vec<Row>,
    /// Edits and bases per (profile, rate), for the realized error rate.
    edits: Vec<(u64, u64)>,
    /// Per (hpc, unit): label -> (gene, occurrences); gene = MULTI if in several genes.
    labels: Vec<FxHashMap<u64, (u64, u32)>>,
}

impl Stats {
    fn new(cfg: &Cfg) -> Self {
        let nh = cfg.hpc.len() * cfg.set.units.len();
        Stats {
            rows: vec![Row::default(); cfg.profiles.len() * cfg.rates.len() * nh],
            edits: vec![(0, 0); cfg.profiles.len() * cfg.rates.len()],
            labels: vec![FxHashMap::default(); nh],
        }
    }

    fn merge(&mut self, o: Stats) {
        for (a, b) in self.rows.iter_mut().zip(&o.rows) {
            a.merge(b);
        }
        for (a, b) in self.edits.iter_mut().zip(&o.edits) {
            a.0 += b.0;
            a.1 += b.1;
        }
        for (a, b) in self.labels.iter_mut().zip(o.labels) {
            for (label, (gene, n)) in b {
                let e = a.entry(label).or_insert((gene, 0));
                if e.0 != gene {
                    e.0 = MULTI;
                }
                e.1 += n;
            }
        }
    }
}

struct Cfg {
    set: SeedSet,
    profiles: Vec<ErrorProfile>,
    rates: Vec<f64>,
    hpc: Vec<bool>,
    tol: u32,
    min_len: usize,
    sample: f64,
    seed: u64,
}

impl Cfg {
    fn row(&self, prof: usize, rate: usize, hpc: usize, unit: usize) -> usize {
        ((prof * self.rates.len() + rate) * self.hpc.len() + hpc) * self.set.units.len() + unit
    }
}

#[derive(Clone)]
struct Processor {
    cfg: Arc<Cfg>,
    local: Stats,
    shared: Arc<Mutex<Stats>>,
}

impl Processor {
    fn add(&mut self, raw: &[u8], gene: u64, index: u64) {
        let cfg = self.cfg.clone();
        let mut rng = SplitMix64::for_record(cfg.seed, index);
        if raw.len() < cfg.min_len || rng.unit() >= cfg.sample {
            return;
        }
        let orig: Vec<UnitSeeds> = cfg.hpc.iter().map(|&h| cfg.set.seeds(raw, h)).collect();
        for (hi, per_unit) in orig.iter().enumerate() {
            for (ui, seeds) in per_unit.iter().enumerate() {
                let map = &mut self.local.labels[hi * cfg.set.units.len() + ui];
                for &(label, _, _) in seeds {
                    let e = map.entry(label).or_insert((gene, 0));
                    if e.0 != gene {
                        e.0 = MULTI;
                    }
                    e.1 += 1;
                }
            }
        }
        for (pi, prof) in cfg.profiles.iter().enumerate() {
            for (ri, &rate) in cfg.rates.iter().enumerate() {
                let read = simulate(raw, rate, prof, &mut rng);
                let e = &mut self.local.edits[pi * cfg.rates.len() + ri];
                e.0 += read.edits as u64;
                e.1 += raw.len() as u64;
                for (hi, &hpc) in cfg.hpc.iter().enumerate() {
                    let rseeds = cfg.set.seeds(&read.seq, hpc);
                    for ui in 0..cfg.set.units.len() {
                        let row = &mut self.local.rows[cfg.row(pi, ri, hi, ui)];
                        let o = &orig[hi][ui];
                        let r = &rseeds[ui];
                        let mut rsorted: Vec<(u64, u32)> = r.iter().map(|&(l, s, _)| (l, s)).collect();
                        rsorted.sort_unstable();
                        // Transcript seeds at their expected read positions.
                        let mut osorted: Vec<(u64, u32)> =
                            o.iter().map(|&(l, s, _)| (l, read.o2r[s as usize])).collect();
                        osorted.sort_unstable();

                        row.seqs += 1;
                        row.bases += raw.len() as u64;
                        row.seeds += o.len() as u64;
                        let mut cov_end = 0u32;
                        for &(l, s, e) in o {
                            row.span_sum += (e - s) as u64;
                            if has_near(&rsorted, l, read.o2r[s as usize], cfg.tol) {
                                row.survived += 1;
                                // Seeds come in start order, so a sweep gives the union.
                                let from = s.max(cov_end);
                                if e > from {
                                    row.covered += (e - from) as u64;
                                }
                                cov_end = cov_end.max(e);
                            }
                        }
                        row.read_seeds += r.len() as u64;
                        row.read_matched += r
                            .iter()
                            .filter(|&&(l, s, _)| has_near(&osorted, l, s, cfg.tol))
                            .count() as u64;
                    }
                }
            }
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for Processor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        // Genes as the index takes them: GENCODE headers are `tx|gene|...`; otherwise each
        // sequence is its own gene.
        let (_, gene) = spearfish::name_and_gene(&String::from_utf8_lossy(record.id()));
        let gene = xxhash_rust::xxh3::xxh3_64(gene.as_bytes());
        self.add(&record.seq(), gene, record.index());
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        let local = std::mem::replace(&mut self.local, Stats::new(&self.cfg));
        self.shared.lock().unwrap().merge(local);
        Ok(())
    }
}

pub fn run(args: Args) -> Result<()> {
    let set = SeedSet::load(&args.schemes, args.schemes_file.as_deref(), 1)?;
    let profiles = args
        .profiles
        .iter()
        .map(|p| profile(p).with_context(|| format!("unknown error profile {p}")))
        .collect::<Result<Vec<_>>>()?;
    let hpc = args
        .hpc
        .iter()
        .map(|h| match h.as_str() {
            "off" => Ok(false),
            "on" => Ok(true),
            _ => bail!("--hpc takes off and/or on, got {h}"),
        })
        .collect::<Result<Vec<_>>>()?;
    let cfg = Arc::new(Cfg {
        set,
        profiles,
        rates: args.rates.clone(),
        hpc,
        tol: args.tol,
        min_len: args.min_len,
        sample: args.sample,
        seed: args.seed,
    });

    let input = args.input.clone().filter(|p| p.as_os_str() != "-");
    let reader = ReaderBuilder::optional_path(input).build()?;
    let shared = Arc::new(Mutex::new(Stats::new(&cfg)));
    let mut proc = Processor {
        cfg: cfg.clone(),
        local: Stats::new(&cfg),
        shared: shared.clone(),
    };
    let t = Instant::now();
    reader.process_parallel(&mut proc, args.threads)?;
    let stats = std::mem::replace(&mut *shared.lock().unwrap(), Stats::new(&cfg));
    eprintln!("survival: {:.1} s", t.elapsed().as_secs_f64());

    // Gene-level specificity per (hpc, unit).
    let unique: Vec<f64> = stats
        .labels
        .iter()
        .map(|m| {
            let (mut u, mut n) = (0u64, 0u64);
            for &(g, c) in m.values() {
                n += c as u64;
                if g != MULTI {
                    u += c as u64;
                }
            }
            u as f64 / n.max(1) as f64
        })
        .collect();

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(std::io::BufWriter::new(std::fs::File::create(p)?)),
        None => Box::new(std::io::stdout().lock()),
    };
    writeln!(
        out,
        "profile\trate\trealized_rate\thpc\tscheme\tlevel\tseqs\tseeds\tdensity_per_kb\tspan_mean\t\
         survival\tcoverage\tread_seeds\tread_precision\tgene_unique_frac"
    )?;
    for (pi, prof) in cfg.profiles.iter().enumerate() {
        for (ri, rate) in cfg.rates.iter().enumerate() {
            let (ed, eb) = stats.edits[pi * cfg.rates.len() + ri];
            for (hi, hpc) in cfg.hpc.iter().enumerate() {
                for (ui, u) in cfg.set.units.iter().enumerate() {
                    let r = &stats.rows[cfg.row(pi, ri, hi, ui)];
                    let div = |a: u64, b: u64| a as f64 / b.max(1) as f64;
                    writeln!(
                        out,
                        "{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.4}\t{:.4}\t{}\t{:.4}\t{:.4}",
                        prof.name,
                        rate,
                        div(ed, eb),
                        if *hpc { "on" } else { "off" },
                        cfg.set.name(u),
                        u.level + 1,
                        r.seqs,
                        r.seeds,
                        1000.0 * div(r.seeds, r.bases),
                        div(r.span_sum, r.seeds),
                        div(r.survived, r.seeds),
                        div(r.covered, r.bases),
                        r.read_seeds,
                        div(r.read_matched, r.read_seeds),
                        unique[hi * cfg.set.units.len() + ui]
                    )?;
                }
            }
        }
    }
    Ok(())
}
