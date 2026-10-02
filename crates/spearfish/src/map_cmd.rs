//! `spearfish map`: map long reads to the indexed transcriptome.
//!
//! Output TSV, one line per reported (read, transcript), best first. A read's lines are
//! contiguous; reads appear in completion order, not input order:
//! read, read_len, strand, rank, transcript, score, qstart, qend, tstart, tend, anchors,
//! anchored, max_disc, flags. Coordinates are parsed (HPC) coordinates for now. Unmapped
//! reads get one line with transcript `*`.

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

type Out = Arc<Mutex<Box<dyn Write + Send>>>;

#[derive(Clone)]
struct Processor {
    mapper: Arc<Mapper<OracleChainer>>,
    buf: Vec<u8>,
    out: Out,
    stats: Arc<Mutex<(u64, u64, u64)>>,
    local: (u64, u64, u64),
}

impl<Rf: Record> ParallelProcessor<Rf> for Processor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
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

    /// Write the batch's records. All mappings of a read are contiguous; reads appear in
    /// completion order, not input order.
    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        self.out.lock().unwrap().write_all(&self.buf)?;
        self.buf.clear();
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
    let out: Out = Arc::new(Mutex::new(sink));
    out.lock().unwrap().write_all(
        b"read\tread_len\tstrand\trank\ttranscript\tscore\tqstart\tqend\ttstart\ttend\tanchors\tanchored\tmax_disc\tflags\n",
    )?;
    let stats = Arc::new(Mutex::new((0u64, 0u64, 0u64)));
    let mut proc = Processor {
        mapper,
        buf: Vec::new(),
        out: out.clone(),
        stats: stats.clone(),
        local: (0, 0, 0),
    };
    let input = Some(args.reads.clone()).filter(|p| p.as_os_str() != "-");
    let reader = ReaderBuilder::optional_path(input)
        .batch_size(args.batch_size)
        .build()?;
    let t1 = Instant::now();
    reader.process_parallel(&mut proc, args.threads)?;
    out.lock().unwrap().flush()?;
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
