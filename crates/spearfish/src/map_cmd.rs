//! `spearfish map`: map long reads to the indexed transcriptome.
//!
//! Output TSV, one line per reported (read, transcript), best first, in input order:
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
        let flags = if r.low_coverage { "low_coverage" } else { "-" };
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

    fn requires_ordering(&self) -> bool {
        true
    }
}

pub fn run(args: Args) -> Result<()> {
    let t0 = Instant::now();
    let idx = Arc::new(Index::open(&args.index)?);
    let mapper = Arc::new(Mapper::new(
        idx,
        OracleChainer {
            p: Default::default(),
        },
        MapParams::default(),
    )?);
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
    let reader = ReaderBuilder::optional_path(input).build()?;
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
