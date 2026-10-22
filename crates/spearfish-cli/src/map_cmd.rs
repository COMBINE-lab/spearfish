//! `spearfish map`: map long reads to the indexed transcriptome.
//!
//! Outputs, any combination (default: TSV on stdout):
//!
//! - **TSV** (`-o`): one line per reported (read, transcript), best first. A read's lines are
//!   contiguous; reads appear in completion order, not input order. The strand is per line: a
//!   read whose orientations nearly tie (antisense transcripts) is reported on both. Columns:
//!   read, read_len, strand, rank, transcript, score, qstart, qend, tstart, tend, anchors,
//!   anchored, max_disc, flags. qstart/qend are raw read coordinates on the read's forward
//!   strand and tstart/tend raw transcript coordinates of the alignment, extended past the
//!   outermost anchors for hits within `--extend-within` of the best score (`--extend-ends`);
//!   other hits' spans run between their outermost anchors and carry the flag `anchor_span`.
//!   anchored counts compressed bases. Unmapped reads get one line with transcript `*`.
//! - **PAF** (`--paf`): the same hits as PAF lines (no CIGAR; unmapped reads omitted; see
//!   [`paf_line`]).
//! - **RAD** (`--rad`): `bulk_long_read` records (see the `long-read-rad` crate), one per
//!   mapped read, in chunks compressed per `--rad-codec`. oarfish quantifies these directly;
//!   the file's tags carry the preset, the score scale, the index identifier, and whether
//!   spans were extended to the read ends.
//!
//! Every run also writes a run record, `<first output>.meta.json` (or `--meta`): the build,
//! command line, index and its manifest, the settings changed, every parameter used, and run
//! statistics.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::ValueEnum;
use paraseq::prelude::*;
use paraseq::{ReaderBuilder, Record};
use serde::Serialize;
use spearfish::{Index, MapOptions, MapSettings, MappedRead, Mapper};

use crate::map_args::MapArgs;
use crate::runinfo::{self, IndexInfo, Invocation, Usage};

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum Codec {
    Lz4,
    Zstd,
    None,
}

#[derive(clap::Args, Debug)]
#[command(after_long_help = preset_table())]
pub struct Args {
    /// Index built by `spearfish index`.
    #[arg(required_unless_present = "print_params")]
    index: Option<PathBuf>,
    /// Reads (FASTA/FASTQ, optionally gzip/bzip2/xz/zstd-compressed); `-` reads stdin.
    #[arg(required_unless_present = "print_params")]
    reads: Option<PathBuf>,
    /// Worker threads (0: every available CPU).
    #[arg(short, long, default_value_t = 0, help_heading = "Input and output")]
    threads: usize,
    /// TSV output (`-`: stdout; the default when no output is given).
    #[arg(short, long, value_name = "FILE", help_heading = "Input and output")]
    output: Option<PathBuf>,
    /// PAF output (`-`: stdout).
    #[arg(long, value_name = "FILE", help_heading = "Input and output")]
    paf: Option<PathBuf>,
    /// RAD output (`bulk_long_read` records, which oarfish quantifies with `-a`).
    #[arg(long, value_name = "FILE", help_heading = "Input and output")]
    rad: Option<PathBuf>,
    /// RAD chunk compression.
    #[arg(
        long,
        value_enum,
        default_value = "lz4",
        hide_short_help = true,
        help_heading = "Input and output"
    )]
    rad_codec: Codec,
    /// Run record (JSON) [default: <first output file>.meta.json; none if no output is a regular
    /// file, e.g. all go to stdout or /dev/null].
    #[arg(long, value_name = "FILE", help_heading = "Input and output")]
    meta: Option<PathBuf>,
    /// Print the complete parameters the settings resolve to (JSON) and exit.
    #[arg(long, help_heading = "Input and output")]
    print_params: bool,
    /// Reads per work batch (smaller batches balance load when read costs vary).
    #[arg(
        long,
        default_value_t = 32,
        hide_short_help = true,
        help_heading = "Input and output"
    )]
    batch_size: usize,
    #[command(flatten)]
    settings: MapArgs,
}

/// The settings that differ between presets, for `--help`.
fn preset_table() -> String {
    let opts: Vec<serde_json::Value> = spearfish::Preset::ALL
        .iter()
        .map(|&p| serde_json::to_value(MapOptions::preset(p)).unwrap())
        .collect();
    let mut rows = Vec::new();
    fn walk(path: &str, vals: &[&serde_json::Value], rows: &mut Vec<(String, Vec<String>)>) {
        if vals.iter().all(|v| v.is_object()) {
            let mut keys: Vec<&String> = vals.iter().flat_map(|v| v.as_object().unwrap().keys()).collect();
            keys.sort();
            keys.dedup();
            for k in keys {
                let sub: Vec<&serde_json::Value> = vals
                    .iter()
                    .map(|v| v.get(k).unwrap_or(&serde_json::Value::Null))
                    .collect();
                let p = if path.is_empty() {
                    k.clone()
                } else {
                    format!("{path}.{k}")
                };
                walk(&p, &sub, rows);
            }
        } else if vals.windows(2).any(|w| w[0] != w[1]) && path != "preset" {
            let show = |v: &serde_json::Value| match v {
                serde_json::Value::Null => "off".to_string(),
                serde_json::Value::Object(_) => "on".to_string(),
                v => v.to_string().trim_matches('"').to_string(),
            };
            rows.push((path.to_string(), vals.iter().map(|v| show(v)).collect()));
        }
    }
    walk("", &opts.iter().collect::<Vec<_>>(), &mut rows);
    let w = rows.iter().map(|r| r.0.len()).max().unwrap_or(0);
    let mut s = format!(
        "Preset differences (all other defaults are shared):\n  {:w$}  ont  drna  hifi\n",
        ""
    );
    for (k, v) in rows {
        s += &format!("  {k:w$}  {}\n", v.join("  "));
    }
    s + "\n`spearfish map --print-params [--preset P] [settings]` prints every parameter."
}

/// Per-thread and total counts.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct Counts {
    pub reads: u64,
    pub mapped: u64,
    pub unmapped: u64,
    /// Mapped reads with more than one reported hit.
    pub multi_mapped: u64,
    pub rescued: u64,
    pub low_coverage: u64,
    pub hits: u64,
    pub bases: u64,
}

impl Counts {
    fn add_read(&mut self, r: &MappedRead) {
        self.reads += 1;
        self.bases += r.len as u64;
        if r.hits.is_empty() {
            self.unmapped += 1;
        } else {
            self.mapped += 1;
            self.multi_mapped += (r.hits.len() > 1) as u64;
        }
        self.rescued += r.rescued as u64;
        self.low_coverage += r.low_coverage as u64;
        self.hits += r.hits.len() as u64;
    }

    fn add(&mut self, o: &Counts) {
        self.reads += o.reads;
        self.mapped += o.mapped;
        self.unmapped += o.unmapped;
        self.multi_mapped += o.multi_mapped;
        self.rescued += o.rescued;
        self.low_coverage += o.low_coverage;
        self.hits += o.hits;
        self.bases += o.bases;
    }
}

type Out = Arc<Mutex<Box<dyn Write + Send>>>;

#[derive(Clone)]
struct Processor {
    mapper: Arc<Mapper>,
    tsv: Option<(Out, Vec<u8>)>,
    paf: Option<(Out, Vec<u8>)>,
    rad: Option<long_read_rad::Sink>,
    rec: long_read_rad::ReadRecord,
    total: Arc<Mutex<Counts>>,
    progress: Arc<AtomicU64>,
    local: Counts,
    batch_reads: u64,
}

/// TSV lines of one read.
fn tsv_lines(buf: &mut Vec<u8>, name: &str, r: &MappedRead, idx: &Index) {
    if r.hits.is_empty() {
        let _ = writeln!(
            buf,
            "{name}\t{}\t*\t0\t*\t0\t0\t0\t0\t0\t0\t0\t0\tunmapped",
            r.len
        );
        return;
    }
    let read_flags = read_flags(r);
    for (rank, h) in r.hits.iter().enumerate() {
        let strand = if h.reverse { '-' } else { '+' };
        let _ = writeln!(
            buf,
            "{name}\t{}\t{strand}\t{rank}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.0}\t{}",
            r.len,
            idx.tx_name(h.tx),
            h.score,
            h.query_start,
            h.query_end,
            h.target_start,
            h.target_end,
            h.anchors,
            h.anchored,
            h.max_discrepancy,
            hit_flags(read_flags, h.extended)
        );
    }
}

fn read_flags(r: &MappedRead) -> &'static str {
    match (r.low_coverage, r.rescued) {
        (true, true) => "rescued,low_coverage",
        (true, false) => "low_coverage",
        (false, true) => "rescued",
        (false, false) => "-",
    }
}

fn hit_flags(read_flags: &str, extended: bool) -> std::borrow::Cow<'_, str> {
    match (extended, read_flags) {
        (true, f) => f.into(),
        (false, "-") => "anchor_span".into(),
        (false, f) => format!("{f},anchor_span").into(),
    }
}

/// PAF lines of one read (none if unmapped). Columns as minimap2 writes them without `-c`:
/// residue matches are the read bases the hit's seeds cover (homopolymer-compressed), the
/// block length the longer of the two spans, and the mapping quality 255 (not computed: tied
/// isoforms are reported as hits instead). Tags: `tp:A` P for the best hit, S for the others;
/// `AS:i` the score rounded; `sc:f` the score; `cm:i` anchors in the chain; `rk:i` the rank;
/// `fl:Z` the TSV flags, when any.
pub fn paf_line(buf: &mut Vec<u8>, name: &str, r: &MappedRead, idx: &Index) {
    let read_flags = read_flags(r);
    for (rank, h) in r.hits.iter().enumerate() {
        let strand = if h.reverse { '-' } else { '+' };
        let block = (h.query_end - h.query_start).max(h.target_end - h.target_start);
        let _ = write!(
            buf,
            "{name}\t{}\t{}\t{}\t{strand}\t{}\t{}\t{}\t{}\t{}\t{block}\t255\ttp:A:{}\tAS:i:{:.0}\tsc:f:{:.2}\tcm:i:{}\trk:i:{rank}",
            r.len,
            h.query_start,
            h.query_end,
            idx.tx_name(h.tx),
            idx.tx_len(h.tx),
            h.target_start,
            h.target_end,
            h.anchored.min(block),
            if rank == 0 { 'P' } else { 'S' },
            h.score,
            h.score,
            h.anchors,
        );
        let f = hit_flags(read_flags, h.extended);
        if f != "-" {
            let _ = write!(buf, "\tfl:Z:{f}");
        }
        buf.push(b'\n');
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for Processor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        let seq = record.seq();
        let r = self.mapper.map(&seq);
        self.local.add_read(&r);
        self.batch_reads += 1;
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
        if self.tsv.is_none() && self.paf.is_none() {
            return Ok(());
        }
        let name = record.id_str().split_whitespace().next().unwrap_or("");
        let idx = self.mapper.index();
        if let Some((_, buf)) = &mut self.tsv {
            tsv_lines(buf, name, &r, idx);
        }
        if let Some((_, buf)) = &mut self.paf {
            paf_line(buf, name, &r, idx);
        }
        Ok(())
    }

    /// Write the batch's records. All mappings of a read are contiguous; reads appear in
    /// completion order, not input order.
    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        for (out, buf) in [&mut self.tsv, &mut self.paf].into_iter().flatten() {
            out.lock().unwrap().write_all(buf)?;
            buf.clear();
        }
        self.progress.fetch_add(self.batch_reads, Ordering::Relaxed);
        self.batch_reads = 0;
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        if let Some(rad) = &mut self.rad {
            rad.flush().map_err(std::io::Error::other)?;
        }
        self.total.lock().unwrap().add(&self.local);
        self.local = Counts::default();
        Ok(())
    }
}

/// Open an output: a file, or stdout for `-`.
fn sink(p: &Path) -> Result<Out> {
    let w: Box<dyn Write + Send> = if p.as_os_str() == "-" {
        Box::new(std::io::BufWriter::new(std::io::stdout()))
    } else {
        let f = std::fs::File::create(p).with_context(|| format!("creating {}", p.display()))?;
        Box::new(std::io::BufWriter::new(f))
    };
    Ok(Arc::new(Mutex::new(w)))
}

/// Draws `mapped N reads (R reads/s)` on stderr every second until dropped.
struct Progress {
    stop: Arc<AtomicBool>,
    handle: Option<std::thread::JoinHandle<()>>,
}

impl Progress {
    fn start(count: Arc<AtomicU64>) -> Self {
        let stop = Arc::new(AtomicBool::new(false));
        let handle = crate::logging::progress_enabled().then(|| {
            let stop = stop.clone();
            std::thread::spawn(move || {
                let t0 = Instant::now();
                let mut drawn = false;
                while !stop.load(Ordering::Relaxed) {
                    std::thread::sleep(std::time::Duration::from_millis(250));
                    let s = t0.elapsed().as_secs_f64();
                    if s < 1.0 {
                        continue;
                    }
                    let n = count.load(Ordering::Relaxed);
                    eprint!("\r\x1b[2Kmapped {n} reads ({:.0} reads/s)", n as f64 / s);
                    drawn = true;
                }
                if drawn {
                    eprint!("\r\x1b[2K");
                }
            })
        });
        Progress { stop, handle }
    }
}

impl Drop for Progress {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

#[derive(Serialize)]
struct Outputs {
    #[serde(skip_serializing_if = "Option::is_none")]
    tsv: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    paf: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    rad: Option<String>,
}

#[derive(Serialize)]
struct RunStats {
    #[serde(flatten)]
    counts: Counts,
    mean_read_length: f64,
    mean_hits_per_mapped_read: f64,
    index_load_seconds: f64,
    mapping_seconds: f64,
    wall_seconds: f64,
    reads_per_second: f64,
    mbp_per_second: f64,
    #[serde(flatten)]
    usage: Usage,
    /// Internal counters (`SPEARFISH_STATS=1`).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    internal: Vec<(String, u64)>,
}

/// Whether `p` names a regular file to write (or one to be created), so a run record can sit
/// next to it: not `-` (stdout), and not a device, pipe or other special file such as
/// `/dev/null`.
fn is_file_output(p: &Path) -> bool {
    p.as_os_str() != "-" && std::fs::metadata(p).map_or(true, |m| m.is_file())
}

/// The run record (`<output>.meta.json`).
#[derive(Serialize)]
struct RunRecord<'a> {
    kind: &'static str,
    #[serde(flatten)]
    invocation: &'a Invocation,
    index: IndexInfo,
    reads: String,
    outputs: Outputs,
    /// The settings changed from the preset.
    settings: &'a MapSettings,
    /// Every parameter used.
    parameters: &'a MapOptions,
    score_denom: f32,
    finished_unix: u64,
    stats: RunStats,
}

/// One mapping run.
pub struct Job {
    pub index: PathBuf,
    pub reads: PathBuf,
    pub threads: usize,
    pub output: Option<PathBuf>,
    pub paf: Option<PathBuf>,
    pub rad: Option<PathBuf>,
    pub rad_codec: Codec,
    pub meta: Option<PathBuf>,
    pub batch_size: usize,
    pub settings: MapSettings,
}

pub fn run(args: Args) -> Result<()> {
    let settings = args.settings.settings();
    if args.print_params {
        println!("{}", serde_json::to_string_pretty(&settings.resolve()?)?);
        return Ok(());
    }
    execute(Job {
        index: args.index.unwrap(),
        reads: args.reads.unwrap(),
        threads: args.threads,
        output: args.output,
        paf: args.paf,
        rad: args.rad,
        rad_codec: args.rad_codec,
        meta: args.meta,
        batch_size: args.batch_size,
        settings,
    })
    .map(|_| ())
}

/// Run a mapping job; returns the read counts.
pub fn execute(args: Job) -> Result<Counts> {
    let settings = args.settings.clone();
    let opts = settings.resolve()?;
    let (index_path, reads) = (args.index.clone(), args.reads.clone());
    let stdout_outputs = [&args.output, &args.paf]
        .iter()
        .filter(|o| o.as_ref().is_some_and(|p| p.as_os_str() == "-"))
        .count();
    if stdout_outputs > 1 {
        bail!("only one output can go to stdout");
    }
    let tsv_path = match (&args.output, &args.paf, &args.rad) {
        (None, None, None) => Some(PathBuf::from("-")),
        (o, _, _) => o.clone(),
    };
    let meta_path = args.meta.clone().or_else(|| {
        [&tsv_path, &args.paf, &args.rad]
            .into_iter()
            .flatten()
            .find(|p| is_file_output(p))
            .map(|p| {
                let mut s = p.as_os_str().to_owned();
                s.push(".meta.json");
                PathBuf::from(s)
            })
    });
    let threads = runinfo::threads(args.threads);
    let invocation = Invocation::new(threads);

    let input = Some(reads.clone()).filter(|p| p.as_os_str() != "-");
    if let Some(p) = &input {
        std::fs::metadata(p).with_context(|| format!("reads {}", p.display()))?;
    }
    let t0 = Instant::now();
    let idx = Arc::new(Index::open(&index_path).with_context(|| format!("index {}", index_path.display()))?);
    let mapper = Arc::new(Mapper::new(idx.clone(), &opts)?);
    let load = t0.elapsed().as_secs_f64();
    let index_info = IndexInfo::of(&index_path, &idx);
    if index_info.manifest.is_none() {
        log::debug!("the index has no manifest (built before manifests were recorded)");
    }

    let tsv = match &tsv_path {
        Some(p) => {
            let out = sink(p)?;
            out.lock().unwrap().write_all(
                b"read\tread_len\tstrand\trank\ttranscript\tscore\tqstart\tqend\ttstart\ttend\tanchors\tanchored\tmax_disc\tflags\n",
            )?;
            Some((out, Vec::new()))
        }
        None => None,
    };
    let paf = match &args.paf {
        Some(p) => Some((sink(p)?, Vec::new())),
        None => None,
    };
    let rad = match &args.rad {
        Some(p) => Some(rad_writer(p, args.rad_codec, &mapper, &settings)?),
        None => None,
    };
    let total = Arc::new(Mutex::new(Counts::default()));
    let progress = Arc::new(AtomicU64::new(0));
    let mut proc = Processor {
        mapper,
        tsv: tsv.clone(),
        paf: paf.clone(),
        rad: rad.as_ref().map(|w| w.sink()),
        rec: long_read_rad::ReadRecord::default(),
        total: total.clone(),
        progress: progress.clone(),
        local: Counts::default(),
        batch_reads: 0,
    };
    let reader = ReaderBuilder::optional_path(input)
        .batch_size(args.batch_size)
        .build()
        .with_context(|| format!("reading {}", reads.display()))?;
    let t1 = Instant::now();
    {
        let _p = Progress::start(progress);
        reader.process_parallel(&mut proc, threads)?;
    }
    for (out, _) in [&tsv, &paf].into_iter().flatten() {
        out.lock().unwrap().flush()?;
    }
    drop(proc);
    let secs = t1.elapsed().as_secs_f64();
    let internal = if map::stats::on() {
        let v: Vec<(String, u64)> = map::stats::snapshot()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect();
        for (k, v) in &v {
            log::info!("stats\t{k}\t{v}");
        }
        v
    } else {
        Vec::new()
    };
    let c = *total.lock().unwrap();
    if let Some(w) = rad {
        w.finish(long_read_rad::Counts {
            reads: c.reads,
            mapped: c.mapped,
            bases: c.bases,
        })?;
    }
    log::info!(
        "map: {} reads ({:.1}% mapped, {:.1}% rescued), {:.1} Mbp in {secs:.2} s ({:.0} reads/s, {:.1} Mbp/s); index load {load:.2} s",
        c.reads,
        100.0 * c.mapped as f64 / c.reads.max(1) as f64,
        100.0 * c.rescued as f64 / c.reads.max(1) as f64,
        c.bases as f64 / 1e6,
        c.reads as f64 / secs,
        c.bases as f64 / 1e6 / secs
    );
    if let Some(m) = &meta_path {
        // Absolute paths, so the record still names the files when read from elsewhere.
        let abs = |p: &Path| {
            if p.as_os_str() == "-" {
                "-".to_string()
            } else {
                spearfish::provenance::absolute(p)
            }
        };
        let show = |p: &Option<PathBuf>| p.as_deref().map(abs);
        let record = RunRecord {
            kind: "spearfish-map-run",
            invocation: &invocation,
            index: index_info,
            reads: abs(&reads),
            outputs: Outputs {
                tsv: show(&tsv_path),
                paf: show(&args.paf),
                rad: show(&args.rad),
            },
            settings: &settings,
            parameters: &opts,
            score_denom: spearfish::SCORE_DENOM,
            finished_unix: runinfo::unix_now(),
            stats: RunStats {
                counts: c,
                mean_read_length: c.bases as f64 / c.reads.max(1) as f64,
                mean_hits_per_mapped_read: c.hits as f64 / c.mapped.max(1) as f64,
                index_load_seconds: load,
                mapping_seconds: secs,
                wall_seconds: t0.elapsed().as_secs_f64(),
                reads_per_second: c.reads as f64 / secs,
                mbp_per_second: c.bases as f64 / 1e6 / secs,
                usage: Usage::now(),
                internal,
            },
        };
        runinfo::write_json(m, &record)?;
    }
    Ok(c)
}

/// A RAD writer over the index's transcripts, with the tags oarfish uses to read the scores
/// and spans.
fn rad_writer(
    path: &Path,
    codec: Codec,
    mapper: &Mapper,
    settings: &MapSettings,
) -> Result<long_read_rad::RadWriter> {
    let codec = match codec {
        Codec::Lz4 => long_read_rad::ChunkCodec::Lz4,
        Codec::Zstd => long_read_rad::ChunkCodec::Zstd,
        Codec::None => long_read_rad::ChunkCodec::None,
    };
    let idx = mapper.index();
    let refs: Vec<(&str, u32)> = (0..idx.n_transcripts() as u32)
        .map(|t| (idx.tx_name(t), idx.tx_len(t)))
        .collect();
    let opts = mapper.options();
    let meta = long_read_rad::Meta {
        mapper: "spearfish".into(),
        mapper_version: spearfish::provenance::BuildInfo::current().short(),
        cmdline: std::env::args().collect::<Vec<_>>().join(" "),
        params: serde_json::to_string(&serde_json::json!({
            "settings": settings,
            "parameters": opts,
        }))?,
        preset: opts.preset.name().into(),
        seq_tech: opts.preset.seq_tech().into(),
        score_denom: mapper.score_denom(),
        ends_extended: opts.params.ends.is_some(),
        strands: opts.params.strands.name().into(),
        index_id: idx.id(),
    };
    long_read_rad::RadWriter::create(path, &refs, &meta, codec)
}
