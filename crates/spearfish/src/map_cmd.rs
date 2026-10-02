//! `spearfish map`: map long reads to the indexed transcriptome.
//!
//! Output TSV, one line per reported (read, transcript), best first. A read's lines are
//! contiguous; reads appear in completion order, not input order. The strand is per line: a
//! read whose orientations nearly tie (antisense transcripts) is reported on both.
//! read, read_len, strand, rank, transcript, score, qstart, qend, tstart, tend, anchors,
//! anchored, max_disc, flags. Coordinates are parsed (HPC) coordinates for now. Unmapped
//! reads get one line with transcript `*`.

use std::io::Write;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use anyhow::Result;
use index::Index;
use map::{AnyChainer, MapParams, Mapper, OracleChainer, RefineChainer};
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
    /// Chainer: `oracle` (per-transcript DP with shared prefixes) or `refine` (multi-transcript
    /// state DP with lazy color refinement). Both give identical scores with no beam.
    #[arg(long, default_value = "oracle")]
    chainer: String,
    /// Refine chainer: keep at most this many states per anchor (default: unbounded, exact).
    #[arg(long)]
    beam: Option<usize>,
    /// Disable the hierarchical rescue of low-coverage reads.
    #[arg(long)]
    no_rescue: bool,
    /// Rescue seed scheme (applied to parsed read and candidate transcripts).
    #[arg(long, default_value = "lcp:K=2")]
    rescue_scheme: String,
    /// Rescue: candidate genes per strand.
    #[arg(long, default_value_t = 16)]
    rescue_genes: usize,
    /// Rescue: re-seed whole transcripts instead of the window the read can occupy.
    #[arg(long)]
    rescue_full: bool,
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
    /// Read type: `ont` (cDNA, both orientations), `drna` (ONT direct RNA: sense
    /// orientation only), or `hifi` (PacBio HiFi). Sets the tie-break's affine scores
    /// (minimap2 map-ont / map-hifi) and the orientations mapped.
    #[arg(long, default_value = "ont")]
    preset: String,
    /// Orientations to map, overriding the preset: `both`, `fwd` (as given) or `rev`.
    #[arg(long)]
    strands: Option<String>,
    /// Tie-break: score segments by unit-cost edit distance instead of affine gaps.
    #[arg(long)]
    unit_tiebreak: bool,
    /// Tie-break: affine-score only segment texts within this many edits of the best
    /// (unit-cost prefilter); default from the preset (drna: 4, others: all).
    #[arg(long)]
    affine_within: Option<u32>,
    /// Tie-break: transcript bases beyond the read's projected ends (plus 1/`tie_pad_div` of
    /// the unanchored end).
    /// Default from the preset (drna: 20 + 1/10, others: 50 + 1/6).
    #[arg(long)]
    tie_pad: Option<u32>,
    #[arg(long)]
    tie_pad_div: Option<u32>,
    /// Tie-break: force whole read ends through the alignment instead of clipping what
    /// doesn't align (adapters, junk).
    #[arg(long)]
    no_tie_clip: bool,
    /// Tie-break: bonus for aligning through to the read's end (minimap2's --end-bonus).
    #[arg(long, default_value_t = 5)]
    end_bonus: i32,
    /// Tie-break: score charged per affine-score point below the best.
    #[arg(long, default_value_t = 1.0)]
    lambda: f32,
    /// Tie-break: align at most this many near-best transcripts.
    #[arg(long, default_value_t = 256)]
    tie_cands: usize,
    /// Candidate genes kept per read (best-scoring first).
    #[arg(long, default_value_t = 8)]
    f_keep: usize,
    /// Also keep genes within this relative score of the last kept one (ties)...
    #[arg(long, default_value_t = 0.01)]
    f_tie: f32,
    /// ...up to this many genes in total.
    #[arg(long, default_value_t = 64)]
    f_max: usize,
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
    mapper: Arc<Mapper<AnyChainer>>,
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
            let strand = if h.reverse { '-' } else { '+' };
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
    let (affine, preset_strands, preset_within, preset_pad) = match args.preset.as_str() {
        "ont" => (map::align::Affine::ONT, "both", None, (50, 6)),
        "drna" => (map::align::Affine::ONT, "fwd", Some(4), (20, 10)),
        "hifi" => (map::align::Affine::HIFI, "both", None, (50, 6)),
        x => anyhow::bail!("unknown preset {x:?} (ont, drna, hifi)"),
    };
    let strands = match args.strands.as_deref().unwrap_or(preset_strands) {
        "both" => map::Strands::Both,
        "fwd" => map::Strands::Forward,
        "rev" => map::Strands::Reverse,
        x => anyhow::bail!("unknown strands {x:?} (both, fwd, rev)"),
    };
    let tiebreak = (!args.no_tiebreak).then(|| map::TieBreak {
        kappa: args.kappa,
        margin_abs: args.tie_margin,
        affine: (!args.unit_tiebreak).then_some(affine),
        lambda: args.lambda,
        clip_ends: !args.no_tie_clip,
        end_bonus: args.end_bonus,
        pad: args.tie_pad.unwrap_or(preset_pad.0),
        pad_div: args.tie_pad_div.unwrap_or(preset_pad.1).max(1),
        affine_within: args.affine_within.or(preset_within),
        max_cands: args.tie_cands,
        segmented: !args.whole_read_tiebreak,
        ..Default::default()
    });
    let anchor = map::AnchorParams {
        f_keep: args.f_keep,
        f_tie: args.f_tie,
        f_max: args.f_max,
        h_max: args.h_max,
        g_max: args.g_max,
        ..Default::default()
    };
    let rescue = (!args.no_rescue).then(|| map::RescueParams {
        spec: args.rescue_scheme.clone(),
        max_genes: args.rescue_genes,
        max_tx_per_gene: args.rescue_tx,
        window: if args.rescue_full {
            None
        } else {
            map::RescueParams::default().window
        },
        ..Default::default()
    });
    let params = MapParams {
        anchor,
        rescue,
        delta_rep: args.delta_rep,
        chain: chain.clone(),
        tiebreak,
        strands,
        ..Default::default()
    };
    let chainer = match args.chainer.as_str() {
        "oracle" => AnyChainer::Oracle(OracleChainer { p: chain }),
        "refine" => AnyChainer::Refine(RefineChainer {
            p: chain,
            beam: args.beam,
        }),
        c => anyhow::bail!("unknown chainer {c:?} (oracle, refine)"),
    };
    let mapper = Arc::new(Mapper::new(idx, chainer, params)?);
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
    if map::stats::on() {
        for (k, v) in map::stats::snapshot() {
            eprintln!("stats\t{k}\t{v}");
        }
    }
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
