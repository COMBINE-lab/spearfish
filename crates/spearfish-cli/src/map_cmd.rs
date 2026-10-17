//! `spearfish map`: map long reads to the indexed transcriptome.
//!
//! Output TSV, one line per reported (read, transcript), best first. A read's lines are
//! contiguous; reads appear in completion order, not input order. The strand is per line: a
//! read whose orientations nearly tie (antisense transcripts) is reported on both.
//! read, read_len, strand, rank, transcript, score, qstart, qend, tstart, tend, anchors,
//! anchored, max_disc, flags. qstart/qend are raw read coordinates on the read's forward
//! strand and tstart/tend raw transcript coordinates of the alignment, extended past the
//! outermost anchors for hits within `--extend-within` of the best score (`--extend-ends`);
//! other hits' spans run between their outermost anchors and carry the flag `anchor_span`.
//! anchored counts compressed bases. Unmapped reads get one line with transcript `*`.
//!
//! With `--rad`, the mappings are also written as a RAD file of `bulk_long_read` records (see
//! the `long-read-rad` crate), one per mapped read, in chunks compressed per `--rad-codec`.
//! oarfish quantifies these directly; the file's tags carry the preset, the score scale, and
//! whether spans were extended to the read ends.

use std::io::Write;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use anyhow::Result;
use paraseq::prelude::*;
use paraseq::{ReaderBuilder, Record};
use spearfish::{ChainerKind, Index, MapOptions, Mapper, Preset};

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Index built by `spearfish index`.
    index: PathBuf,
    /// Reads (FASTA/FASTQ, optionally compressed); `-` reads stdin.
    reads: PathBuf,
    /// Worker threads (0 = all cores).
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
    /// Output TSV (default stdout, unless `--rad` is given).
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// Write the mappings as a RAD file (`bulk_long_read` records, which oarfish reads).
    #[arg(long)]
    rad: Option<PathBuf>,
    /// RAD chunk compression: `lz4`, `zstd` or `none`.
    #[arg(long, default_value = "lz4")]
    rad_codec: String,
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
    /// Rescue: re-seed the other orientation only if its base-level best is within this
    /// margin of the chosen one's (negative: never; default: always).
    #[arg(long, allow_hyphen_values = true)]
    rescue_other_margin: Option<f32>,
    /// Rescue: window slack beyond the projected read, as PAD,DIV (PAD + read/DIV bases).
    #[arg(long)]
    rescue_window: Option<String>,
    /// Rescue: candidate genes must score this fraction of the best gene's base score.
    #[arg(long)]
    rescue_alpha: Option<f32>,
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
    /// Default 20 + 1/10.
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
    /// Tie-break: when tied hits share no anchor, group them by shared anchors instead of
    /// aligning every hit's whole read (fewer whole-read alignments; default for hifi).
    #[arg(long)]
    tie_group: bool,
    /// Near-tied read orientations: compare by whole-read edit distance instead of clipped
    /// affine scores.
    #[arg(long)]
    orient_edit: bool,
    /// Tie-break: align read ends in raw (uncompressed) space, needing an index with raw
    /// transcripts (default: the preset's; on except for drna).
    #[arg(long)]
    raw_ends: Option<bool>,
    /// Tie-break: A's appended to transcript windows that reach the 3' end, standing for the
    /// mRNA's poly(A) tail (0: off).
    #[arg(long)]
    tail_a: Option<u32>,
    /// Chain every anchor separately (no merging of collinear anchors within a unitig).
    #[arg(long)]
    no_merge: bool,
    /// Extend reported spans past the outermost anchors to where the read's alignment ends
    /// (default: the preset's; on for ont). Off: spans run from the first to the last
    /// anchor.
    #[arg(long)]
    extend_ends: Option<bool>,
    /// Extend only hits scoring within this fraction of the read's best.
    #[arg(long, default_value_t = 0.05)]
    extend_within: f32,
    /// Tie-break: affine DP band (cells off the diagonals spanning each segment).
    #[arg(long)]
    tie_band: Option<usize>,
    /// Tie-break: X-drop for read-end extensions (default: the preset's).
    #[arg(long)]
    xdrop: Option<i32>,
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
    mapper: Arc<Mapper>,
    buf: Vec<u8>,
    out: Option<Out>,
    rad: Option<long_read_rad::Sink>,
    rec: long_read_rad::ReadRecord,
    stats: Arc<Mutex<(u64, u64, u64)>>,
    local: (u64, u64, u64),
}

impl<Rf: Record> ParallelProcessor<Rf> for Processor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        let seq = record.seq();
        let r = self.mapper.map(&seq);
        self.local.0 += 1;
        self.local.2 += seq.len() as u64;
        if !r.hits.is_empty() {
            self.local.1 += 1;
        }
        if let Some(rad) = &mut self.rad
            && !r.hits.is_empty()
        {
            use long_read_rad::{Aln, aln_flags, read_flags};
            let rec = &mut self.rec;
            rec.clear();
            rec.read_len = r.len;
            rec.flags = if r.low_coverage {
                read_flags::LOW_COVERAGE
            } else {
                0
            } | if r.rescued { read_flags::RESCUED } else { 0 };
            for h in &r.hits {
                rec.push(Aln {
                    tid: h.tx,
                    fw: !h.reverse,
                    ref_start: h.target_start,
                    ref_end: h.target_end,
                    query_start: h.query_start,
                    query_end: h.query_end,
                    score: h.score,
                    flags: if h.extended { 0 } else { aln_flags::ANCHOR_SPAN },
                });
            }
            rad.push(rec).map_err(std::io::Error::other)?;
        }
        if self.out.is_none() {
            return Ok(());
        }
        let name = record.id_str().split_whitespace().next().unwrap_or("");
        if r.hits.is_empty() {
            let _ = writeln!(
                self.buf,
                "{name}\t{}\t*\t0\t*\t0\t0\t0\t0\t0\t0\t0\t0\tunmapped",
                r.len
            );
            return Ok(());
        }
        let flags = match (r.low_coverage, r.rescued) {
            (true, true) => "rescued,low_coverage",
            (true, false) => "low_coverage",
            (false, true) => "rescued",
            (false, false) => "-",
        };
        for (rank, h) in r.hits.iter().enumerate() {
            let strand = if h.reverse { '-' } else { '+' };
            let flags = match (h.extended, flags) {
                (true, f) => f.to_string(),
                (false, "-") => "anchor_span".to_string(),
                (false, f) => format!("{f},anchor_span"),
            };
            let _ = writeln!(
                self.buf,
                "{name}\t{}\t{strand}\t{rank}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.0}\t{flags}",
                r.len,
                self.mapper.index().tx_name(h.tx),
                h.score,
                h.query_start,
                h.query_end,
                h.target_start,
                h.target_end,
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
        if let Some(out) = &self.out {
            out.lock().unwrap().write_all(&self.buf)?;
            self.buf.clear();
        }
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        if let Some(rad) = &mut self.rad {
            rad.flush().map_err(std::io::Error::other)?;
        }
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
    let mapper = Arc::new(Mapper::new(idx, &options(&args)?)?);
    let load = t0.elapsed().as_secs_f64();
    let sink: Option<Box<dyn Write + Send>> = match (&args.output, &args.rad) {
        (Some(p), _) => Some(Box::new(std::io::BufWriter::new(std::fs::File::create(p)?))),
        (None, None) => Some(Box::new(std::io::BufWriter::new(std::io::stdout()))),
        (None, Some(_)) => None,
    };
    let out: Option<Out> = sink.map(|s| Arc::new(Mutex::new(s)));
    if let Some(out) = &out {
        out.lock().unwrap().write_all(
            b"read\tread_len\tstrand\trank\ttranscript\tscore\tqstart\tqend\ttstart\ttend\tanchors\tanchored\tmax_disc\tflags\n",
        )?;
    }
    let rad = match &args.rad {
        Some(p) => Some(rad_writer(p, &args, &mapper)?),
        None => None,
    };
    let stats = Arc::new(Mutex::new((0u64, 0u64, 0u64)));
    let mut proc = Processor {
        mapper,
        buf: Vec::new(),
        out: out.clone(),
        rad: rad.as_ref().map(|w| w.sink()),
        rec: long_read_rad::ReadRecord::default(),
        stats: stats.clone(),
        local: (0, 0, 0),
    };
    let input = Some(args.reads.clone()).filter(|p| p.as_os_str() != "-");
    let reader = ReaderBuilder::optional_path(input)
        .batch_size(args.batch_size)
        .build()?;
    let t1 = Instant::now();
    reader.process_parallel(&mut proc, args.threads)?;
    if let Some(out) = &out {
        out.lock().unwrap().flush()?;
    }
    drop(proc);
    let secs = t1.elapsed().as_secs_f64();
    if map::stats::on() {
        for (k, v) in map::stats::snapshot() {
            eprintln!("stats\t{k}\t{v}");
        }
    }
    let (n, mapped, bases) = *stats.lock().unwrap();
    if let Some(w) = rad {
        w.finish(long_read_rad::Counts {
            reads: n,
            mapped,
            bases,
        })?;
    }
    eprintln!(
        "map: {n} reads ({:.1}% mapped), {:.1} Mbp in {secs:.2} s ({:.0} reads/s, {:.1} Mbp/s); index load {load:.2} s",
        100.0 * mapped as f64 / n.max(1) as f64,
        bases as f64 / 1e6,
        n as f64 / secs,
        bases as f64 / 1e6 / secs
    );
    Ok(())
}

/// A RAD writer over the index's transcripts, with the tags oarfish uses to read the scores
/// and spans.
fn rad_writer(path: &std::path::Path, args: &Args, mapper: &Mapper) -> Result<long_read_rad::RadWriter> {
    let codec = match args.rad_codec.as_str() {
        "lz4" => long_read_rad::ChunkCodec::Lz4,
        "zstd" => long_read_rad::ChunkCodec::Zstd,
        "none" => long_read_rad::ChunkCodec::None,
        c => anyhow::bail!("unknown --rad-codec {c:?} (lz4, zstd, none)"),
    };
    let idx = mapper.index();
    let refs: Vec<(&str, u32)> = (0..idx.n_transcripts() as u32)
        .map(|t| (idx.tx_name(t), idx.tx_len(t)))
        .collect();
    let opts = options(args)?;
    let meta = long_read_rad::Meta {
        mapper: "spearfish".into(),
        mapper_version: env!("CARGO_PKG_VERSION").into(),
        cmdline: std::env::args().collect::<Vec<_>>().join(" "),
        params: format!("{:?}", opts.params),
        preset: args.preset.clone(),
        seq_tech: match args.preset.parse::<Preset>()? {
            Preset::Ont => "ont-cdna",
            Preset::Drna => "ont-drna",
            Preset::Hifi => "pac-bio-hifi",
        }
        .into(),
        // spearfish weighs a hit scoring d below the best by exp(-d / 2).
        score_denom: 2.0,
        ends_extended: opts.params.ends.is_some(),
        strands: match opts.params.strands {
            map::Strands::Both => "both",
            map::Strands::Forward => "fwd",
            map::Strands::Reverse => "rev",
        }
        .into(),
    };
    long_read_rad::RadWriter::create(path, &refs, &meta, codec)
}

/// Mapping options: the preset's, then any explicit flags.
fn options(args: &Args) -> Result<MapOptions> {
    let preset: Preset = args.preset.parse()?;
    let mut o = MapOptions::preset(preset);
    let p = &mut o.params;
    p.chain.alpha = args.alpha;
    p.chain.beta = args.beta;
    p.chain.h = args.window;
    if let Some(s) = &args.strands {
        p.strands = match s.as_str() {
            "both" => map::Strands::Both,
            "fwd" => map::Strands::Forward,
            "rev" => map::Strands::Reverse,
            x => anyhow::bail!("unknown strands {x:?} (both, fwd, rev)"),
        };
    }
    if args.no_tiebreak {
        p.tiebreak = None;
    } else if let Some(tb) = p.tiebreak.as_mut() {
        tb.kappa = args.kappa;
        tb.margin_abs = args.tie_margin;
        if args.unit_tiebreak {
            tb.affine = None;
        }
        tb.lambda = args.lambda;
        tb.clip_ends = !args.no_tie_clip;
        tb.end_bonus = args.end_bonus;
        tb.group_fallback |= args.tie_group;
        tb.orient_clipped = !args.orient_edit;
        if let Some(r) = args.raw_ends {
            tb.raw_ends = r;
        }
        if let Some(n) = args.tail_a {
            tb.tail_a = n;
        }
        if let Some(b) = args.tie_band {
            tb.band = b;
        }
        if let Some(x) = args.xdrop {
            tb.xdrop = x;
        }
        if let Some(x) = args.tie_pad {
            tb.pad = x;
        }
        if let Some(x) = args.tie_pad_div {
            tb.pad_div = x.max(1);
        }
        if args.affine_within.is_some() {
            tb.affine_within = args.affine_within;
        }
        tb.max_cands = args.tie_cands;
        tb.segmented = !args.whole_read_tiebreak;
    }
    p.anchor.f_keep = args.f_keep;
    p.anchor.f_tie = args.f_tie;
    p.anchor.f_max = args.f_max;
    p.anchor.h_max = args.h_max;
    p.anchor.g_max = args.g_max;
    if args.no_rescue {
        p.rescue = None;
    } else if let Some(r) = p.rescue.as_mut() {
        r.spec = args.rescue_scheme.clone();
        r.max_genes = args.rescue_genes;
        r.max_tx_per_gene = args.rescue_tx;
        if args.rescue_full {
            r.window = None;
        }
        if let Some(w) = &args.rescue_window {
            let (a, b) = w
                .split_once(',')
                .ok_or_else(|| anyhow::anyhow!("--rescue-window PAD,DIV"))?;
            r.window = Some((a.parse()?, b.parse()?));
        }
        if let Some(a) = args.rescue_alpha {
            r.alpha = a;
        }
        if args.rescue_other_margin.is_some() {
            r.other_margin = args.rescue_other_margin;
        }
    }
    p.delta_rep = args.delta_rep;
    p.merge_collinear = !args.no_merge;
    match args.extend_ends {
        Some(false) => o.params.ends = None,
        Some(true) => o.extend_ends(preset),
        None => {}
    }
    if let Some(e) = &mut o.params.ends {
        e.within = args.extend_within;
    }
    o.chainer = match args.chainer.as_str() {
        "oracle" => ChainerKind::Oracle,
        "refine" => ChainerKind::Refine { beam: args.beam },
        c => anyhow::bail!("unknown chainer {c:?} (oracle, refine)"),
    };
    Ok(o)
}
