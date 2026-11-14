//! `spearfish map`: map long reads to the indexed transcriptome.
//!
//! Outputs, any combination (default: TSV on stdout):
//!
//! - **TSV** (`-o`): one line per reported (read, transcript), best first; reads appear in
//!   completion order, not input order. Unmapped reads get one line with transcript `*`.
//! - **PAF** (`--paf`): the same hits as PAF lines (no CIGAR; unmapped reads omitted).
//! - **RAD** (`--rad`): `bulk_long_read` records, one per mapped read, in chunks compressed per
//!   `--rad-codec`, which oarfish quantifies directly.
//!
//! [`spearfish::output`] writes all three and documents their columns and tags.
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
use paraseq::Record;
use paraseq::prelude::*;
use serde::Serialize;
use spearfish::{Index, MapOptions, MapSettings, MappedRead, Mapper, output};

use crate::map_args::MapArgs;
use crate::runinfo::{self, IndexInfo, Invocation, Usage};

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
pub enum Codec {
    #[default]
    Lz4,
    Zstd,
    None,
}

/// Reads per work batch, unless `--batch-size` says otherwise.
pub const DEFAULT_BATCH_SIZE: usize = 32;

/// What `--bump-alloc` asks for: the arenas with this configuration, or no arenas.
#[derive(Clone, Copy, Debug)]
struct BumpAlloc(Option<nudgy::Config>);

fn parse_bump_alloc(s: &str) -> Result<BumpAlloc, nudgy::ParseConfigError> {
    if s.eq_ignore_ascii_case("off") || s.eq_ignore_ascii_case("none") {
        return Ok(BumpAlloc(None));
    }
    // Arenas that serve blocks of no size are no arenas.
    let config: nudgy::Config = s.parse()?;
    Ok(BumpAlloc((config.max_block() > 0).then_some(config)))
}

#[derive(clap::Args, Debug)]
#[command(after_long_help = preset_table())]
pub struct Args {
    /// Index built by `spearfish index`.
    #[arg(required_unless_present = "print_params")]
    index: Option<PathBuf>,
    /// Reads (FASTA/FASTQ, optionally gzip/bzip2/zstd-compressed); `-` reads stdin.
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
    /// RAD output, a file and not `-` (`bulk_long_read` records, which oarfish quantifies with
    /// `-a`).
    #[arg(long, value_name = "FILE", help_heading = "Input and output")]
    rad: Option<PathBuf>,
    /// RAD chunk compression.
    #[arg(
        long,
        value_enum,
        default_value_t = Codec::default(),
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
        default_value_t = DEFAULT_BATCH_SIZE,
        hide_short_help = true,
        help_heading = "Input and output"
    )]
    batch_size: usize,
    /// Take the mapper's small allocations from per-thread arenas [default: on, where the
    /// system has them]: less CPU time for a few MB of memory a thread (an arena is at most
    /// 4 MB), and the same output. A size is the largest allocation the arenas serve
    /// (`--bump-alloc=1KB`; 8KB without one, at most 32KB); `off` or `none` is the same as
    /// --no-bump-alloc.
    #[arg(
        long,
        value_name = "SIZE",
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "8KB",
        value_parser = parse_bump_alloc,
        overrides_with = "no_bump_alloc",
        hide_short_help = true,
        help_heading = "Input and output"
    )]
    bump_alloc: Option<BumpAlloc>,
    /// Take every allocation from the general allocator (no per-thread arenas): a lower peak
    /// memory, for more CPU time.
    #[arg(
        long,
        overrides_with = "bump_alloc",
        hide_short_help = true,
        help_heading = "Input and output"
    )]
    no_bump_alloc: bool,
    #[command(flatten)]
    settings: MapArgs,
}

/// The settings that differ between presets, for `--help` (the `-fast` presets, which turn
/// off what came after 0.1.0, are described in a line).
fn preset_table() -> String {
    use spearfish::Preset;
    let opts: Vec<serde_json::Value> = [Preset::Ont, Preset::Drna, Preset::Hifi]
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
    s + "\nont-fast, drna-fast and hifi-fast map as 0.1.0 did: faster on noisy reads and in less\n\
         memory, less accurate, and every read with a chain is reported.\n\
         \n`spearfish map --print-params [--preset P] [settings]` prints every parameter."
}

/// Per-thread and total counts.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct Counts {
    pub reads: u64,
    pub mapped: u64,
    pub unmapped: u64,
    /// Mapped reads with more than one reported hit.
    pub multi_mapped: u64,
    /// Mapped reads whose best hit is on the read's forward strand.
    pub forward: u64,
    pub rescued: u64,
    pub low_coverage: u64,
    /// Reads whose best chain scored under the minimum (unmapped unless kept).
    pub low_score: u64,
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
            self.forward += !r.hits[0].reverse as u64;
        }
        self.rescued += r.rescued as u64;
        self.low_coverage += r.low_coverage as u64;
        self.low_score += r.low_score as u64;
        self.hits += r.hits.len() as u64;
    }

    fn add(&mut self, o: &Counts) {
        self.reads += o.reads;
        self.mapped += o.mapped;
        self.unmapped += o.unmapped;
        self.multi_mapped += o.multi_mapped;
        self.forward += o.forward;
        self.rescued += o.rescued;
        self.low_coverage += o.low_coverage;
        self.low_score += o.low_score;
        self.hits += o.hits;
        self.bases += o.bases;
    }
}

/// Fewest mapped reads from which [`oriented`] concludes anything.
const ORIENTED_MIN_READS: u64 = 1000;
/// Share of the mapped reads on one strand from which a library counts as oriented. Reads of
/// both strands split about evenly; among oriented reads, only wrongly placed ones (and reads
/// of sequence shared antisense) map to the other strand.
const ORIENTED_SHARE: f64 = 0.99;

/// The `--strands` value that suits reads mapped in both orientations nearly all of which
/// mapped in one: such reads were oriented before mapping. A file of reads does not say so.
fn oriented(c: &Counts, both_strands: bool) -> Option<&'static str> {
    if !both_strands || c.mapped < ORIENTED_MIN_READS {
        return None;
    }
    let fwd = c.forward as f64 / c.mapped as f64;
    if fwd >= ORIENTED_SHARE {
        Some("fwd")
    } else if fwd <= 1.0 - ORIENTED_SHARE {
        Some("rev")
    } else {
        None
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

impl<Rf: Record> ParallelProcessor<Rf> for Processor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        let seq = record.seq();
        let r = self.mapper.map(&seq);
        self.local.add_read(&r);
        self.batch_reads += 1;
        if let Some(rad) = &mut self.rad
            && output::rad_record(&mut self.rec, &r)
        {
            rad.push(&self.rec).map_err(std::io::Error::other)?;
        }
        if self.tsv.is_none() && self.paf.is_none() {
            return Ok(());
        }
        // Lossy: a header that is not UTF-8 gets U+FFFD instead of failing the run.
        let id = String::from_utf8_lossy(record.id());
        let name = id.split_whitespace().next().unwrap_or("");
        let idx = self.mapper.index();
        if let Some((_, buf)) = &mut self.tsv {
            output::tsv_lines(buf, name, &r, idx);
        }
        if let Some((_, buf)) = &mut self.paf {
            output::paf_lines(buf, name, &r, idx);
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
    /// The per-thread arenas (`--bump-alloc`), when they were on.
    #[serde(skip_serializing_if = "Option::is_none")]
    bump_alloc: Option<ArenaStats>,
    /// Internal counters (`SPEARFISH_STATS=1`).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    internal: Vec<(String, u64)>,
}

/// What the per-thread arenas did, from `nudgy::stats()`.
#[derive(Serialize)]
struct ArenaStats {
    /// The largest allocation, in bytes, that they served.
    max_block: usize,
    /// The most threads that had an arena at once.
    arenas_peak: usize,
    /// Times a thread moved on to an empty chunk of its arena.
    chunk_refills: u64,
    /// Times a thread found no arena free, and times its arena had no empty chunk.
    no_arena: u64,
    arena_full: u64,
}

impl ArenaStats {
    fn now() -> Option<ArenaStats> {
        let s = nudgy::stats();
        s.enabled.then_some(ArenaStats {
            max_block: s.max_block,
            arenas_peak: s.arenas_peak,
            chunk_refills: s.chunk_refills,
            no_arena: s.no_arena,
            arena_full: s.arena_full,
        })
    }
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
    match (args.no_bump_alloc, args.bump_alloc) {
        (true, _) | (_, Some(BumpAlloc(None))) => {}
        // The default: the arenas where the system has them, and none, silently, elsewhere.
        (false, None) => {
            nudgy::enable();
        }
        (false, Some(BumpAlloc(Some(config)))) => {
            if !nudgy::enable_with(config) {
                log::warn!("--bump-alloc is not available on this system; mapping without it");
            }
        }
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
    let both_strands = opts.params.strands == spearfish::Strands::Both;
    let (index_path, reads) = (args.index.clone(), args.reads.clone());
    let stdout_outputs = [&args.output, &args.paf]
        .iter()
        .filter(|o| o.as_ref().is_some_and(|p| p.as_os_str() == "-"))
        .count();
    if stdout_outputs > 1 {
        bail!("only one output can go to stdout");
    }
    if args.rad.as_deref().is_some_and(|p| p.as_os_str() == "-") {
        bail!(
            "RAD output needs a file (its read counts are written when it is finished); `-` (stdout) is not supported"
        );
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
            .map(|p| runinfo::with_suffix(p, ".meta.json"))
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
            out.lock()
                .unwrap()
                .write_all(output::tsv_header(mapper.quant_scores()).as_bytes())?;
            Some((out, Vec::new()))
        }
        None => None,
    };
    let paf = match &args.paf {
        Some(p) => Some((sink(p)?, Vec::new())),
        None => None,
    };
    let rad = match &args.rad {
        Some(p) => Some(rad_writer(p, args.rad_codec, &mapper, &settings, &invocation)?),
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
    let reader = spearfish::open_fastx(input.as_deref(), Some(args.batch_size))
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
        "map: {} reads ({:.1}% mapped, {:.1}% rescued, {:.1}% low score), {:.1} Mbp in {secs:.2} s ({:.0} reads/s, {:.1} Mbp/s); index load {load:.2} s",
        c.reads,
        100.0 * c.mapped as f64 / c.reads.max(1) as f64,
        100.0 * c.rescued as f64 / c.reads.max(1) as f64,
        100.0 * c.low_score as f64 / c.reads.max(1) as f64,
        c.bases as f64 / 1e6,
        c.reads as f64 / secs,
        c.bases as f64 / 1e6 / secs
    );
    let arenas = nudgy::stats();
    if arenas.enabled {
        log::debug!(
            "bump-alloc: blocks of up to {} bytes, {} arenas at most, {} chunks started, no arena free {} times, arena full {} times (at most {} allocations left to the general allocator)",
            arenas.max_block,
            arenas.arenas_peak,
            arenas.chunk_refills,
            arenas.no_arena,
            arenas.arena_full,
            arenas.fallback_allocations()
        );
    }
    if let Some(strand) = oriented(&c, both_strands) {
        log::info!(
            "map: {:.1}% of the mapped reads map forward: the reads look oriented, and `--strands {strand}` maps only that orientation, which is faster",
            100.0 * c.forward as f64 / c.mapped as f64
        );
    }
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
                bump_alloc: ArenaStats::now(),
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
    invocation: &Invocation,
) -> Result<long_read_rad::RadWriter> {
    let codec = match codec {
        Codec::Lz4 => long_read_rad::ChunkCodec::Lz4,
        Codec::Zstd => long_read_rad::ChunkCodec::Zstd,
        Codec::None => long_read_rad::ChunkCodec::None,
    };
    let cmdline = invocation.command_line.join(" ");
    output::rad_writer(path, mapper, Some(settings), &cmdline, codec)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn oriented_needs_both_orientations_and_one_sided_reads() {
        let counts = |mapped, forward| Counts {
            mapped,
            forward,
            ..Counts::default()
        };
        assert_eq!(oriented(&counts(10_000, 9_950), true), Some("fwd"));
        assert_eq!(oriented(&counts(10_000, 60), true), Some("rev"));
        // Reads of both strands, too few reads, and a run that mapped one orientation anyway.
        assert_eq!(oriented(&counts(10_000, 5_100), true), None);
        assert_eq!(oriented(&counts(10_000, 9_800), true), None);
        assert_eq!(oriented(&counts(900, 900), true), None);
        assert_eq!(oriented(&counts(10_000, 10_000), false), None);
    }

    #[test]
    fn bump_alloc_takes_a_size_or_off() {
        let limit = |s: &str| parse_bump_alloc(s).map(|b| b.0.map(|c| c.max_block()));
        assert_eq!(limit("off"), Ok(None));
        assert_eq!(limit("none"), Ok(None));
        assert_eq!(limit("None"), Ok(None));
        assert_eq!(limit("0"), Ok(None));
        assert_eq!(limit("1KB"), Ok(Some(1024)));
        assert_eq!(limit("8KB"), Ok(Some(nudgy::DEFAULT_MAX_BLOCK)));
        assert_eq!(limit("32KB"), Ok(Some(nudgy::LARGEST_MAX_BLOCK)));
        assert!(limit("33KB").is_err());
        assert!(limit("on").is_err());
    }
}
