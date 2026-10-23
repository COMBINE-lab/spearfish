//! `spearfish locality`: measure LCP context radii λ_k and ρ_k by random single edits (an edit
//! changes only cores whose extended window `[s - λ_k, e + ρ_k)` contains it).
//!
//! For each sequence and edit type, apply `--edits` random edits, re-parse, and record for
//! every changed core whether it overlaps the edit or how far left/right of it it lies.
//! Distances are in parsed coordinates (homopolymer-compressed with `--hpc`).

use std::io::Write;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use anyhow::Result;
use lcp_core::hpc::Hpc;
use lcp_core::locality::{Edit, Placement, edit_effects};
use lcp_core::{LcpParser, Parser};
use paraseq::Record;
use paraseq::prelude::*;

use crate::hist::Hist;
use crate::rng::SplitMix64;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// FASTA/FASTQ input (optionally compressed); `-` or omitted reads stdin.
    input: Option<PathBuf>,
    /// Highest LCP level to parse.
    #[arg(short = 'K', long, default_value_t = 7)]
    max_level: u8,
    /// Homopolymer-compress sequences before parsing (edits apply to the compressed sequence).
    #[arg(long)]
    hpc: bool,
    /// Random edits of each type per sequence.
    #[arg(long, default_value_t = 3)]
    edits: u32,
    /// Only use sequences at least this long (after HPC).
    #[arg(long, default_value_t = 200)]
    min_len: usize,
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

const KINDS: [&str; 3] = ["sub", "ins", "del"];
const CAP: usize = 1 << 16;

#[derive(Clone)]
struct Cell {
    edits: u64,
    /// Changed cores per edit, and those overlapping the edit.
    changed: Hist,
    overlap: Hist,
    lambda: Hist,
    rho: Hist,
}

impl Cell {
    fn new() -> Self {
        Cell {
            edits: 0,
            changed: Hist::new(CAP),
            overlap: Hist::new(CAP),
            lambda: Hist::new(CAP),
            rho: Hist::new(CAP),
        }
    }

    fn merge(&mut self, o: &Cell) {
        self.edits += o.edits;
        self.changed.merge(&o.changed);
        self.overlap.merge(&o.overlap);
        self.lambda.merge(&o.lambda);
        self.rho.merge(&o.rho);
    }
}

#[derive(Clone)]
struct Stats {
    /// `cells[kind][level]`.
    cells: Vec<Vec<Cell>>,
    core_len: Vec<Hist>,
    seqs: u64,
}

impl Stats {
    fn new(k: u8) -> Self {
        Stats {
            cells: (0..KINDS.len())
                .map(|_| (0..k).map(|_| Cell::new()).collect())
                .collect(),
            core_len: (0..k).map(|_| Hist::new(CAP)).collect(),
            seqs: 0,
        }
    }

    fn merge(&mut self, o: &Stats) {
        for (a, b) in self.cells.iter_mut().flatten().zip(o.cells.iter().flatten()) {
            a.merge(b);
        }
        for (a, b) in self.core_len.iter_mut().zip(&o.core_len) {
            a.merge(b);
        }
        self.seqs += o.seqs;
    }
}

#[derive(Clone)]
struct Processor {
    parser: LcpParser,
    args: Arc<Args>,
    local: Stats,
    shared: Arc<Mutex<Stats>>,
}

impl Processor {
    fn add(&mut self, raw: &[u8], index: u64) {
        let k = self.args.max_level;
        let seq = if self.args.hpc {
            Hpc::new(raw).seq
        } else {
            raw.to_vec()
        };
        if seq.len() < self.args.min_len {
            return;
        }
        self.local.seqs += 1;
        let orig = self.parser.parse(&seq, k);
        for (l, level) in orig.levels.iter().enumerate() {
            for c in &level.cores {
                self.local.core_len[l].add(c.len() as u64);
            }
        }
        let mut rng = SplitMix64::for_record(self.args.seed, index);
        let n = seq.len() as u64;
        for (kind, cells) in self.local.cells.iter_mut().enumerate() {
            for _ in 0..self.args.edits {
                let (edit, base) = match kind {
                    0 => {
                        let p = rng.below(n) as u32;
                        (Edit::Sub(p), rng.other_base(seq[p as usize]))
                    }
                    1 => (Edit::Ins(rng.below(n + 1) as u32), rng.base()),
                    _ => (Edit::Del(rng.below(n) as u32), b'A'),
                };
                let mutated = self.parser.parse(&edit.apply(&seq, base), k);
                let mut changed = vec![0u64; k as usize];
                let mut overlap = vec![0u64; k as usize];
                edit_effects(&orig, &mutated, edit, |l, place| {
                    changed[l] += 1;
                    match place {
                        Placement::Overlap => overlap[l] += 1,
                        Placement::Right(d) => cells[l].lambda.add(d as u64),
                        Placement::Left(d) => cells[l].rho.add(d as u64),
                    }
                });
                for l in 0..k as usize {
                    cells[l].edits += 1;
                    cells[l].changed.add(changed[l]);
                    cells[l].overlap.add(overlap[l]);
                }
            }
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for Processor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        self.add(&record.seq(), record.index());
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        self.shared.lock().unwrap().merge(&self.local);
        self.local = Stats::new(self.args.max_level);
        Ok(())
    }
}

pub fn run(args: Args) -> Result<()> {
    let input = args.input.clone().filter(|p| p.as_os_str() != "-");
    let reader = spearfish::open_fastx(input.as_deref(), None)?;
    let k = args.max_level;
    let threads = args.threads;
    let args = Arc::new(args);
    let shared = Arc::new(Mutex::new(Stats::new(k)));
    let mut proc = Processor {
        parser: LcpParser::default(),
        args: args.clone(),
        local: Stats::new(k),
        shared: shared.clone(),
    };
    let t = Instant::now();
    reader.process_parallel(&mut proc, threads)?;
    let stats = shared.lock().unwrap().clone();
    eprintln!(
        "locality: {} sequences, {:.1} s",
        stats.seqs,
        t.elapsed().as_secs_f64()
    );

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(std::io::BufWriter::new(std::fs::File::create(p)?)),
        None => Box::new(std::io::stdout().lock()),
    };
    writeln!(
        out,
        "edit\tlevel\thpc\tedits\tcore_len_mean\tchanged_mean\tchanged_p99\toverlap_mean\t\
         lambda_n\tlambda_p50\tlambda_p90\tlambda_p99\tlambda_max\trho_n\trho_p50\trho_p90\trho_p99\trho_max"
    )?;
    for (kind, cells) in stats.cells.iter().enumerate() {
        for (l, c) in cells.iter().enumerate() {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{:.2}\t{:.3}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                KINDS[kind],
                l + 1,
                args.hpc,
                c.edits,
                stats.core_len[l].mean(),
                c.changed.mean(),
                c.changed.quantile(0.99),
                c.overlap.mean(),
                c.lambda.n(),
                c.lambda.quantile(0.5),
                c.lambda.quantile(0.9),
                c.lambda.quantile(0.99),
                c.lambda.max(),
                c.rho.n(),
                c.rho.quantile(0.5),
                c.rho.quantile(0.9),
                c.rho.quantile(0.99),
                c.rho.max()
            )?;
        }
    }
    out.flush()?;
    Ok(())
}
