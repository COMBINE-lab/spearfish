//! `spearfish map`: map long reads to the indexed transcriptome.
//!
//! Output TSV, one line per reported (read, transcript), best first, in input order:
//! read, read_len, strand, rank, transcript, score, qstart, qend, tstart, tend, anchors,
//! anchored, max_disc, flags. Coordinates are parsed (HPC) coordinates for now. Unmapped
//! reads get one line with transcript `*`.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use anyhow::Result;
use index::Index;
use map::{MapParams, Mapper, OracleChainer};
use paraseq::prelude::*;
use paraseq::{ReaderBuilder, Record};

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Index built by `spearfish index`.
    index: PathBuf,
    /// Reads (FASTA/FASTQ, optionally compressed); `-` reads stdin.
    reads: PathBuf,
    /// Worker threads (0 = all cores).
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
    /// Output TSV (default stdout).
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// Write records in completion order instead of input order.
    #[arg(long)]
    unordered: bool,
    /// Reads per work batch (smaller batches balance load better when read costs vary).
    #[arg(long, default_value_t = 32)]
    batch_size: usize,
    /// Chaining: cost per base of gap discrepancy |Δq - ρΔt|.
    #[arg(long, default_value_t = 1.5)]
    alpha: f32,
    /// Chaining: cost per log2 of gap length.
    #[arg(long, default_value_t = 0.0)]
    beta: f32,
    /// Chaining: predecessor window (anchors).
    #[arg(long, default_value_t = 50)]
    window: usize,
    /// Disable the hierarchical rescue of low-coverage reads.
    #[arg(long)]
    no_rescue: bool,
    /// Rescue seed scheme (applied to parsed read and candidate transcripts).
    #[arg(long, default_value = "lcp:K=2")]
    rescue_scheme: String,
    /// Rescue: candidate genes per strand.
    #[arg(long, default_value_t = 16)]
    rescue_genes: usize,
    /// Rescue: transcripts per gene (those with the most base anchors).
    #[arg(long, default_value_t = 128)]
    rescue_tx: usize,
    /// Disable the base-level tie-break among near-best transcripts.
    #[arg(long)]
    no_tiebreak: bool,
    /// Tie-break: score charged per extra edit.
    #[arg(long, default_value_t = 6.0)]
    kappa: f32,
    /// Tie-break: align whole reads instead of only where tied transcripts differ.
    #[arg(long)]
    whole_read_tiebreak: bool,
    /// Tie-break: absolute score margin defining the tie set.
    #[arg(long, default_value_t = 30.0)]
    tie_margin: f32,
    /// Candidate genes kept per read (best-scoring first).
    #[arg(long, default_value_t = 8)]
    f_keep: usize,
    /// Ignore seed labels found in more than this many genes.
    #[arg(long, default_value_t = 1000)]
    h_max: usize,
    /// Gene count at which a seed's weight reaches zero (specificity weighting; the
    /// default effectively disables it).
    #[arg(long, default_value_t = 1e9)]
    g_max: f32,
    /// Report transcripts within this fraction of the best score.
    #[arg(long, default_value_t = 0.3)]
    delta_rep: f32,
}

/// Output sink with a reorder buffer: batches finish in any order, are keyed by their first
/// record index, and are written as soon as every earlier record has been written. Workers
/// never wait for their turn; whoever holds the lock flushes what is ready.
struct Sink {
    w: Box<dyn Write + Send>,
    ordered: bool,
    /// Index of the next record to write (ordered mode).
    next: u64,
    /// Finished batches not yet writable: first record -> (records, bytes).
    pending: BTreeMap<u64, (u64, Vec<u8>)>,
}

impl Sink {
    fn push(&mut self, first: u64, n: u64, buf: Vec<u8>) -> std::io::Result<()> {
        if !self.ordered {
            return self.w.write_all(&buf);
        }
        self.pending.insert(first, (n, buf));
        while let Some(entry) = self.pending.first_entry() {
            if *entry.key() != self.next {
                break;
            }
            let (n, buf) = entry.remove();
            self.w.write_all(&buf)?;
            self.next += n;
        }
        Ok(())
    }
}

type Out = Arc<Mutex<Sink>>;

#[derive(Clone)]
struct Processor {
    mapper: Arc<Mapper<OracleChainer>>,
    buf: Vec<u8>,
    out: Out,
    stats: Arc<Mutex<(u64, u64, u64)>>,
    local: (u64, u64, u64),
    /// First record index and record count of the current batch.
    batch: Option<(u64, u64)>,
}

impl<Rf: Record> ParallelProcessor<Rf> for Processor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        let b = self.batch.get_or_insert((record.index(), 0));
        b.1 += 1;
        let seq = record.seq();
        let r = self.mapper.map(&seq);
        let name = record.id_str().split_whitespace().next().unwrap_or("");
        let strand = if r.reverse { '-' } else { '+' };
        self.local.0 += 1;
        self.local.2 += seq.len() as u64;
        if r.hits.is_empty() {
            let _ = writeln!(
                self.buf,
                "{name}\t{}\t*\t0\t*\t0\t0\t0\t0\t0\t0\t0\t0\tunmapped",
                r.len
            );
            return Ok(());
        }
        self.local.1 += 1;
        let flags = match (r.low_coverage, r.rescued) {
            (true, true) => "rescued,low_coverage",
            (true, false) => "low_coverage",
            (false, true) => "rescued",
            (false, false) => "-",
        };
        for (rank, h) in r.hits.iter().enumerate() {
            let _ = writeln!(
                self.buf,
                "{name}\t{}\t{strand}\t{rank}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.0}\t{flags}",
                r.len,
                self.mapper.idx.tx_name(h.tx),
                h.score,
                h.qstart,
                h.qend,
                h.tstart,
                h.tend,
                h.anchors,
                h.anchored,
                h.max_discrepancy
            );
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        if let Some((first, n)) = self.batch.take() {
            let buf = std::mem::take(&mut self.buf);
            self.out.lock().unwrap().push(first, n, buf)?;
        }
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        let mut s = self.stats.lock().unwrap();
        s.0 += self.local.0;
        s.1 += self.local.1;
        s.2 += self.local.2;
        self.local = (0, 0, 0);
        Ok(())
    }
}

pub fn run(args: Args) -> Result<()> {
    let t0 = Instant::now();
    let idx = Arc::new(Index::open(&args.index)?);
    let chain = map::ChainParams {
        alpha: args.alpha,
        beta: args.beta,
        h: args.window,
        ..Default::default()
    };
    let tiebreak = (!args.no_tiebreak).then(|| map::TieBreak {
        kappa: args.kappa,
        margin_abs: args.tie_margin,
        segmented: !args.whole_read_tiebreak,
        ..Default::default()
    });
    let anchor = map::AnchorParams {
        f_keep: args.f_keep,
        h_max: args.h_max,
        g_max: args.g_max,
        ..Default::default()
    };
    let rescue = (!args.no_rescue).then(|| map::RescueParams {
        spec: args.rescue_scheme.clone(),
        max_genes: args.rescue_genes,
        max_tx_per_gene: args.rescue_tx,
        ..Default::default()
    });
    let params = MapParams {
        anchor,
        rescue,
        delta_rep: args.delta_rep,
        chain: chain.clone(),
        tiebreak,
        ..Default::default()
    };
    let mapper = Arc::new(Mapper::new(idx, OracleChainer { p: chain }, params)?);
    let load = t0.elapsed().as_secs_f64();
    let sink: Box<dyn Write + Send> = match &args.output {
        Some(p) => Box::new(std::io::BufWriter::new(std::fs::File::create(p)?)),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };
    let out: Out = Arc::new(Mutex::new(Sink {
        w: sink,
        ordered: !args.unordered,
        next: 0,
        pending: BTreeMap::new(),
    }));
    out.lock().unwrap().w.write_all(
        b"read\tread_len\tstrand\trank\ttranscript\tscore\tqstart\tqend\ttstart\ttend\tanchors\tanchored\tmax_disc\tflags\n",
    )?;
    let stats = Arc::new(Mutex::new((0u64, 0u64, 0u64)));
    let mut proc = Processor {
        mapper,
        buf: Vec::new(),
        out: out.clone(),
        stats: stats.clone(),
        local: (0, 0, 0),
        batch: None,
    };
    let input = Some(args.reads.clone()).filter(|p| p.as_os_str() != "-");
    let reader = ReaderBuilder::optional_path(input)
        .batch_size(args.batch_size)
        .build()?;
    let t1 = Instant::now();
    reader.process_parallel(&mut proc, args.threads)?;
    {
        let mut o = out.lock().unwrap();
        anyhow::ensure!(
            o.pending.is_empty(),
            "internal error: {} batches left unwritten",
            o.pending.len()
        );
        o.w.flush()?;
    }
    let secs = t1.elapsed().as_secs_f64();
    let (n, mapped, bases) = *stats.lock().unwrap();
    eprintln!(
        "map: {n} reads ({:.1}% mapped), {:.1} Mbp in {secs:.2} s ({:.0} reads/s, {:.1} Mbp/s); index load {load:.2} s",
        100.0 * mapped as f64 / n.max(1) as f64,
        bases as f64 / 1e6,
        n as f64 / secs,
        bases as f64 / 1e6 / secs
    );
    Ok(())
}
