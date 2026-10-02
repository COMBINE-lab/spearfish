//! `spearfish lcp-stats`: per-level census of LCP cores over a FASTA/FASTQ file.
//!
//! Lengths, spacings and uncovered ends are reported in raw coordinates, also with `--hpc`
//! (spans are mapped back through the HPC offset map).

use std::io::Write;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use anyhow::Result;
use lcp_core::hpc::Hpc;
use lcp_core::{LcpParser, Parser};
use paraseq::prelude::*;
use paraseq::{ReaderBuilder, Record};
use rustc_hash::FxHashMap;

use crate::hist::Hist;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// FASTA/FASTQ input (optionally compressed); `-` or omitted reads stdin.
    input: Option<PathBuf>,
    /// Highest LCP level to parse.
    #[arg(short = 'K', long, default_value_t = 7)]
    max_level: u8,
    /// Homopolymer-compress sequences before parsing.
    #[arg(long)]
    hpc: bool,
    /// Worker threads (0 = all cores).
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
    /// Output TSV (default stdout).
    #[arg(short, long)]
    output: Option<PathBuf>,
}

const LEN_CAP: usize = 1 << 16;

#[derive(Clone)]
struct LevelStats {
    labels: FxHashMap<u64, u32>,
    len: Hist,
    spacing: Hist,
    /// Bases before the first core / after the last core, per sequence that has a core.
    head: Hist,
    tail: Hist,
    /// Sequences with no core at this level.
    empty: u64,
}

impl LevelStats {
    fn new() -> Self {
        LevelStats {
            labels: FxHashMap::default(),
            len: Hist::new(LEN_CAP),
            spacing: Hist::new(LEN_CAP),
            head: Hist::new(LEN_CAP),
            tail: Hist::new(LEN_CAP),
            empty: 0,
        }
    }

    fn merge(&mut self, o: LevelStats) {
        if self.labels.len() < o.labels.len() {
            let mine = std::mem::replace(&mut self.labels, o.labels);
            for (k, v) in mine {
                *self.labels.entry(k).or_default() += v;
            }
        } else {
            for (k, v) in o.labels {
                *self.labels.entry(k).or_default() += v;
            }
        }
        self.len.merge(&o.len);
        self.spacing.merge(&o.spacing);
        self.head.merge(&o.head);
        self.tail.merge(&o.tail);
        self.empty += o.empty;
    }
}

#[derive(Clone)]
struct Stats {
    levels: Vec<LevelStats>,
    seqs: u64,
    bases: u64,
}

impl Stats {
    fn new(k: u8) -> Self {
        Stats {
            levels: (0..k).map(|_| LevelStats::new()).collect(),
            seqs: 0,
            bases: 0,
        }
    }

    fn merge(&mut self, o: Stats) {
        for (a, b) in self.levels.iter_mut().zip(o.levels) {
            a.merge(b);
        }
        self.seqs += o.seqs;
        self.bases += o.bases;
    }
}

#[derive(Clone)]
struct Processor {
    parser: LcpParser,
    max_level: u8,
    hpc: bool,
    local: Stats,
    shared: Arc<Mutex<Stats>>,
}

impl Processor {
    fn add(&mut self, raw: &[u8]) {
        self.local.seqs += 1;
        self.local.bases += raw.len() as u64;
        let hpc = self.hpc.then(|| Hpc::new(raw));
        let seq: &[u8] = hpc.as_ref().map_or(raw, |h| &h.seq);
        let to_raw = |s: u32, e: u32| hpc.as_ref().map_or((s, e), |h| h.to_raw(s, e));
        let h = self.parser.parse(seq, self.max_level);
        for (level, st) in h.levels.iter().zip(self.local.levels.iter_mut()) {
            let cores = &level.cores;
            if cores.is_empty() {
                st.empty += 1;
                continue;
            }
            let mut prev_start = None;
            for c in cores {
                *st.labels.entry(c.label).or_default() += 1;
                let (s, e) = to_raw(c.start, c.end);
                st.len.add((e - s) as u64);
                if let Some(p) = prev_start {
                    st.spacing.add((s - p) as u64);
                }
                prev_start = Some(s);
            }
            let first = to_raw(cores[0].start, cores[0].end).0;
            let last_end = cores.iter().map(|c| to_raw(c.start, c.end).1).max().unwrap();
            st.head.add(first as u64);
            st.tail.add((raw.len() as u32 - last_end) as u64);
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for Processor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        self.add(&record.seq());
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        let local = std::mem::replace(&mut self.local, Stats::new(self.max_level));
        self.shared.lock().unwrap().merge(local);
        Ok(())
    }
}

pub fn run(args: Args) -> Result<()> {
    let input = args.input.filter(|p| p.as_os_str() != "-");
    let reader = ReaderBuilder::optional_path(input).build()?;
    let shared = Arc::new(Mutex::new(Stats::new(args.max_level)));
    let mut proc = Processor {
        parser: LcpParser::default(),
        max_level: args.max_level,
        hpc: args.hpc,
        local: Stats::new(args.max_level),
        shared: shared.clone(),
    };
    let t = Instant::now();
    reader.process_parallel(&mut proc, args.threads)?;
    let secs = t.elapsed().as_secs_f64();
    let stats = std::mem::replace(&mut *shared.lock().unwrap(), Stats::new(0));
    eprintln!(
        "lcp-stats: {} sequences, {:.1} Mbp, {:.2} s ({:.1} Mbp/s total)",
        stats.seqs,
        stats.bases as f64 / 1e6,
        secs,
        stats.bases as f64 / 1e6 / secs
    );

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(std::io::BufWriter::new(std::fs::File::create(p)?)),
        None => Box::new(std::io::stdout().lock()),
    };
    writeln!(
        out,
        "level\thpc\toccurrences\tdistinct\tsingleton_frac\tlen_mean\tlen_p10\tlen_p50\tlen_p90\tlen_max\t\
         spacing_mean\tspacing_p50\tspacing_p90\thead_p50\thead_p90\ttail_p50\ttail_p90\tseqs_without_cores"
    )?;
    for (k, st) in stats.levels.iter().enumerate() {
        let singletons = st.labels.values().filter(|&&c| c == 1).count();
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{:.4}\t{:.2}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            k + 1,
            args.hpc,
            st.len.n(),
            st.labels.len(),
            singletons as f64 / st.labels.len().max(1) as f64,
            st.len.mean(),
            st.len.quantile(0.1),
            st.len.quantile(0.5),
            st.len.quantile(0.9),
            st.len.max(),
            st.spacing.mean(),
            st.spacing.quantile(0.5),
            st.spacing.quantile(0.9),
            st.head.quantile(0.5),
            st.head.quantile(0.9),
            st.tail.quantile(0.5),
            st.tail.quantile(0.9),
            st.empty
        )?;
    }
    Ok(())
}
