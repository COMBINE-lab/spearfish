//! `spearfish map-spliced`: map long reads, spliced, to a genome index; and the genome side of
//! `spearfish index` (`--genome`) and `inspect`.
//!
//! Spliced alignments are written as PAF (`--paf`, the default on stdout), SAM (`--sam`) or
//! BAM (`--bam`), with introns as `N` in the CIGAR. Reads appear in completion order; a
//! read's records are adjacent.
//!
//! Every run also writes a run record, `<first output>.meta.json` (or `--meta`).

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::ValueEnum;
use paraseq::Record;
use paraseq::prelude::*;
use serde::Serialize;
use spearfish::Preset;
use spearfish::genome::{
    BaseLevel, GenomeIndex, GenomeIndexOptions, GenomeManifest, SpliceOptions, SplicedMapper, SplicedRead,
    Workspace, output,
};

use crate::runinfo::{self, Invocation, Usage};

const HEADING: &str = "Mapping";
const IO: &str = "Input and output";

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum BaseLevelArg {
    /// Exon blocks only (fastest; enough for quantification and exon structure).
    #[default]
    Blocks,
    /// Each block aligned base by base with the fewest edits.
    Edits,
    /// Base by base with affine gap costs, as minimap2 places indels.
    Affine,
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
pub enum PresetArg {
    /// ONT cDNA (both orientations).
    #[default]
    Ont,
    /// ONT direct RNA (reads are the transcripts' sense strand).
    Drna,
    /// PacBio HiFi, Kinnex and Iso-Seq (as ont).
    Hifi,
}

impl From<PresetArg> for Preset {
    fn from(p: PresetArg) -> Self {
        match p {
            PresetArg::Ont => Preset::Ont,
            PresetArg::Drna => Preset::Drna,
            PresetArg::Hifi => Preset::Hifi,
        }
    }
}

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Genome index built by `spearfish index --genome`.
    index: PathBuf,
    /// Reads (FASTA/FASTQ, optionally gzip/bzip2/zstd-compressed); `-` reads stdin.
    reads: PathBuf,
    /// Worker threads (0: every available CPU).
    #[arg(short, long, default_value_t = 0, help_heading = IO)]
    threads: usize,
    /// PAF output (`-`: stdout; the default when no output is given).
    #[arg(long, value_name = "FILE", help_heading = IO)]
    paf: Option<PathBuf>,
    /// Run record (JSON) [default: <first output file>.meta.json; none if no output is a regular
    /// file].
    #[arg(long, value_name = "FILE", help_heading = IO)]
    meta: Option<PathBuf>,
    /// Reads per work batch.
    #[arg(long, default_value_t = crate::map_cmd::DEFAULT_BATCH_SIZE, hide_short_help = true, help_heading = IO)]
    batch_size: usize,
    /// Take the mapper's small allocations from per-thread arenas [default: on]; see
    /// `spearfish map --help`.
    #[arg(
        long,
        value_name = "SIZE",
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "8KB",
        value_parser = crate::map_cmd::parse_bump_alloc,
        overrides_with = "no_bump_alloc",
        hide_short_help = true,
        help_heading = IO
    )]
    bump_alloc: Option<crate::map_cmd::BumpAlloc>,
    /// Take every allocation from the general allocator.
    #[arg(long, overrides_with = "bump_alloc", hide_short_help = true, help_heading = IO)]
    no_bump_alloc: bool,
    /// Read type.
    #[arg(long, value_enum, default_value = "ont", help_heading = HEADING)]
    preset: PresetArg,
    #[command(flatten)]
    genome: GenomeArgs,
}

/// `spearfish map-spliced`.
pub fn run(args: Args) -> Result<()> {
    if spearfish::index_kind(&args.index).with_context(|| format!("index {}", args.index.display()))?
        == spearfish::IndexKind::Transcriptome
    {
        bail!(
            "{} is a transcriptome index: map reads to it with `spearfish map`, or build a genome index with `spearfish index --genome`",
            args.index.display()
        );
    }
    crate::map_cmd::enable_arenas(args.no_bump_alloc, args.bump_alloc);
    execute(Job {
        index: args.index,
        reads: args.reads,
        threads: args.threads,
        paf: args.paf,
        meta: args.meta,
        batch_size: args.batch_size,
        preset: args.preset.into(),
        args: args.genome,
    })
    .map(|_| ())
}

/// The mapping and output options of `spearfish map-spliced`.
#[derive(clap::Args, Debug, Default)]
pub struct GenomeArgs {
    /// Known transcripts (GTF, GFF3 or BED12, optionally compressed). Optional: their splice
    /// sites, junctions and exons are preferred where a read fits them, and junctions not
    /// among them are still found. Sequence names must match the genome's.
    #[arg(short = 'a', long, alias = "junctions", value_name = "FILE", help_heading = HEADING)]
    pub annotation: Option<PathBuf>,
    /// SAM output (`-`: stdout).
    #[arg(long, value_name = "FILE", help_heading = IO)]
    pub sam: Option<PathBuf>,
    /// BAM output (`-`: stdout), unsorted, a read's records together.
    #[arg(long, value_name = "FILE", help_heading = IO)]
    pub bam: Option<PathBuf>,
    /// BAM compression level, 0 (none) to 9.
    #[arg(
        long,
        value_name = "N",
        default_value_t = 1,
        value_parser = clap::value_parser!(u8).range(0..=9),
        hide_short_help = true,
        help_heading = IO
    )]
    pub bam_level: u8,
    /// Report up to N secondary alignments per read.
    #[arg(long, value_name = "N", default_value_t = 0, help_heading = HEADING)]
    pub secondary: usize,
    /// How much of an alignment to compute and report.
    #[arg(long, value_enum, value_name = "LEVEL", default_value = "blocks", help_heading = HEADING)]
    pub base_level: BaseLevelArg,
    /// Longest intron.
    #[arg(long, value_name = "BASES", default_value_t = 200_000, help_heading = HEADING)]
    pub max_intron: u32,
    /// A reference skip shorter than this between two seed matches is a deletion, not an
    /// intron.
    #[arg(long, value_name = "BASES", default_value_t = 40, hide_short_help = true, help_heading = HEADING)]
    pub min_intron: u32,
    /// Aligned hits scoring below this are dropped.
    #[arg(long, value_name = "SCORE", default_value_t = 100, hide_short_help = true, help_heading = HEADING)]
    pub min_splice_score: i32,
    /// Seeds with more occurrences than this are not used (unless the read gets no hit).
    #[arg(long, value_name = "N", default_value_t = 64, hide_short_help = true, help_heading = HEADING)]
    pub seed_cap: usize,
    /// Loci chained per read.
    #[arg(long, value_name = "N", default_value_t = 6, hide_short_help = true, help_heading = HEADING)]
    pub max_loci: usize,
    /// Use pairs of shorter seeds only for reads the main seeds anchor weakly (faster on
    /// accurate reads, loses a little on noisy ones).
    #[arg(long, hide_short_help = true, help_heading = HEADING)]
    pub no_pair_seeds: bool,
}

impl GenomeArgs {
    /// The library options these arguments select under `preset`.
    pub fn options(&self, preset: Preset) -> SpliceOptions {
        let mut o = SpliceOptions::preset(preset);
        o.max_intron = self.max_intron;
        o.min_intron = self.min_intron;
        o.max_hits = 1 + self.secondary;
        o.base_level = match self.base_level {
            BaseLevelArg::Blocks => BaseLevel::Blocks,
            BaseLevelArg::Edits => BaseLevel::Edits,
            BaseLevelArg::Affine => BaseLevel::Affine,
        };
        o.min_score = self.min_splice_score;
        o.seed_cap = self.seed_cap;
        o.max_loci = self.max_loci;
        o.pair_seeds = !self.no_pair_seeds;
        o
    }
}

/// `spearfish index --genome`.
pub fn index(fasta: &Path, out: &Path, threads: usize, opts: &GenomeIndexOptions) -> Result<()> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(runinfo::threads(threads))
        .build()?;
    let (index, bytes, manifest) =
        pool.install(|| GenomeIndex::build_to_file(fasta, opts, out, runinfo::command_line()))?;
    log::info!(
        "index: {} (genome, hpc={}, id {})",
        out.display(),
        opts.hpc,
        manifest.index_id
    );
    log::info!(
        "  sequences {}, bases {}, build {:.1}s, size {:.1} MB",
        manifest.sequences,
        manifest.bases,
        manifest.build_seconds,
        bytes as f64 / 1e6
    );
    for (name, size) in index.sizes() {
        log::debug!("  section {name}\t{:.1} MB", size as f64 / 1e6);
    }
    Ok(())
}

/// `spearfish inspect` on a genome index.
pub fn inspect(path: &Path, json: bool) -> Result<()> {
    let t = Instant::now();
    let idx = GenomeIndex::open(path)?;
    let load = t.elapsed().as_secs_f64();
    let sizes = idx.sizes();
    let total: usize = sizes.iter().map(|(_, b)| b).sum();
    let manifest = idx.manifest();
    let mut out = std::io::stdout().lock();
    if json {
        let j = serde_json::json!({
            "path": path.display().to_string(),
            "kind": "genome",
            "index_id": idx.id(),
            "format_version": index::store::FORMAT_VERSION,
            "sequences": idx.n_sequences(),
            "bases": idx.total_len(),
            "hpc": idx.hpc(),
            "bytes": total,
            "load_seconds": load,
            "manifest": manifest,
        });
        writeln!(out, "{}", serde_json::to_string_pretty(&j)?)?;
        return Ok(());
    }
    writeln!(out, "index\t{}", path.display())?;
    writeln!(out, "kind\tgenome")?;
    writeln!(out, "index_id\t{}", idx.id())?;
    match &manifest {
        Some(m) => {
            writeln!(out, "built_by\tspearfish {}", m.built_by.short())?;
            writeln!(out, "created_unix\t{}", m.created_unix)?;
            writeln!(out, "input\t{}", m.input)?;
            writeln!(out, "max_occ\t{}", m.options.max_occ)?;
        }
        None => writeln!(out, "manifest\t- (written without one)")?,
    }
    writeln!(out, "sequences\t{}", idx.n_sequences())?;
    writeln!(out, "bases\t{}", idx.total_len())?;
    writeln!(out, "hpc\t{}", idx.hpc())?;
    writeln!(out, "bytes\t{total}")?;
    writeln!(out, "load_seconds\t{load:.4}")?;
    for (name, bytes) in sizes {
        writeln!(out, "bytes.{name}\t{bytes}")?;
    }
    Ok(())
}

/// One spliced mapping run.
pub struct Job {
    pub index: PathBuf,
    pub reads: PathBuf,
    pub threads: usize,
    pub paf: Option<PathBuf>,
    pub meta: Option<PathBuf>,
    pub batch_size: usize,
    pub preset: Preset,
    pub args: GenomeArgs,
}

/// Per-thread and total counts.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct Counts {
    pub reads: u64,
    pub mapped: u64,
    pub unmapped: u64,
    /// Mapped reads whose primary alignment has an intron.
    pub spliced: u64,
    /// Introns in primary alignments.
    pub junctions: u64,
    /// Mapped reads with more than one reported alignment.
    pub multi_mapped: u64,
    /// Mapped reads whose primary alignment is on the read's forward strand.
    pub forward: u64,
    pub alignments: u64,
    pub bases: u64,
}

impl Counts {
    fn add_read(&mut self, r: &SplicedRead) {
        self.reads += 1;
        self.bases += r.len as u64;
        match r.alignments.first() {
            None => self.unmapped += 1,
            Some(a) => {
                self.mapped += 1;
                self.spliced += a.spliced() as u64;
                self.junctions += a.introns().count() as u64;
                self.multi_mapped += (r.alignments.len() > 1) as u64;
                self.forward += !a.reverse as u64;
            }
        }
        self.alignments += r.alignments.len() as u64;
    }

    fn add(&mut self, o: &Counts) {
        self.reads += o.reads;
        self.mapped += o.mapped;
        self.unmapped += o.unmapped;
        self.spliced += o.spliced;
        self.junctions += o.junctions;
        self.multi_mapped += o.multi_mapped;
        self.forward += o.forward;
        self.alignments += o.alignments;
        self.bases += o.bases;
    }
}

type Out = crate::map_cmd::Out;

struct Processor {
    mapper: Arc<SplicedMapper>,
    work: Workspace,
    paf: Option<(Out, Vec<u8>)>,
    sam: Option<(Out, Vec<u8>)>,
    /// BAM: the sink, the batch's records, and their BGZF blocks.
    bam: Option<(Out, Vec<u8>, Vec<u8>)>,
    bam_level: u8,
    total: Arc<Mutex<(Counts, spearfish::genome::InternalStats)>>,
    progress: Arc<AtomicU64>,
    local: Counts,
    batch_reads: u64,
}

impl Clone for Processor {
    fn clone(&self) -> Self {
        Processor {
            mapper: self.mapper.clone(),
            work: Workspace::default(),
            paf: self.paf.clone(),
            sam: self.sam.clone(),
            bam: self.bam.clone(),
            bam_level: self.bam_level,
            total: self.total.clone(),
            progress: self.progress.clone(),
            local: Counts::default(),
            batch_reads: 0,
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for Processor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        let seq = record.seq();
        let r = self.mapper.map_with(&seq, &mut self.work);
        self.local.add_read(&r);
        self.batch_reads += 1;
        if r.alignments.is_empty() {
            return Ok(());
        }
        // Lossy: a header that is not UTF-8 gets U+FFFD instead of failing the run.
        let id = String::from_utf8_lossy(record.id());
        let name = id.split_whitespace().next().unwrap_or("");
        let idx = self.mapper.index();
        if let Some((_, buf)) = &mut self.paf {
            output::paf_lines(buf, name, &r, idx);
        }
        if let Some((_, buf)) = &mut self.sam {
            output::sam_records(buf, name, &seq, &r, idx);
        }
        if let Some((_, raw, _)) = &mut self.bam {
            output::bam_records(raw, name, &seq, &r, idx);
        }
        Ok(())
    }

    /// Write the batch's records. For BAM they are compressed here, by the thread that mapped
    /// them, into complete BGZF blocks, so the shared writer only appends bytes.
    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        for (out, buf) in [&mut self.paf, &mut self.sam].into_iter().flatten() {
            if !buf.is_empty() {
                out.lock().unwrap().write_all(buf)?;
                buf.clear();
            }
        }
        if let Some((out, raw, z)) = &mut self.bam
            && !raw.is_empty()
        {
            z.clear();
            output::bgzf(raw, self.bam_level, z);
            raw.clear();
            out.lock().unwrap().write_all(z)?;
        }
        self.progress.fetch_add(self.batch_reads, Ordering::Relaxed);
        self.batch_reads = 0;
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        let mut t = self.total.lock().unwrap();
        t.0.add(&self.local);
        t.1.add(self.work.stats());
        self.local = Counts::default();
        Ok(())
    }
}

#[derive(Serialize)]
struct Outputs {
    #[serde(skip_serializing_if = "Option::is_none")]
    paf: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    sam: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    bam: Option<String>,
}

#[derive(Serialize)]
struct IndexInfo {
    path: String,
    kind: &'static str,
    index_id: String,
    sequences: usize,
    bases: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    manifest: Option<GenomeManifest>,
}

#[derive(Serialize)]
struct AnnotationInfo {
    path: String,
    format: &'static str,
    transcripts: usize,
    introns: usize,
    exons: usize,
    /// Records on sequences the index does not hold.
    skipped_records: usize,
}

#[derive(Serialize)]
struct RunStats {
    #[serde(flatten)]
    counts: Counts,
    mean_read_length: f64,
    index_load_seconds: f64,
    mapping_seconds: f64,
    wall_seconds: f64,
    reads_per_second: f64,
    mbp_per_second: f64,
    #[serde(flatten)]
    usage: Usage,
    /// Internal counters (`SPEARFISH_STATS=1`).
    #[serde(skip_serializing_if = "Option::is_none")]
    internal: Option<String>,
}

#[derive(Serialize)]
struct RunRecord<'a> {
    kind: &'static str,
    #[serde(flatten)]
    invocation: &'a Invocation,
    index: IndexInfo,
    #[serde(skip_serializing_if = "Option::is_none")]
    annotation: Option<AnnotationInfo>,
    reads: String,
    outputs: Outputs,
    preset: &'static str,
    /// Every parameter used.
    parameters: &'a SpliceOptions,
    bam_deflate: &'static str,
    finished_unix: u64,
    stats: RunStats,
}

/// Run a spliced mapping job; returns the read counts.
pub fn execute(job: Job) -> Result<Counts> {
    let a = &job.args;
    let opts = a.options(job.preset);
    // PAF on stdout when no output is named.
    let paf_path = match (&job.paf, &a.sam, &a.bam) {
        (None, None, None) => Some(PathBuf::from("-")),
        (p, _, _) => p.clone(),
    };
    let outs = [&paf_path, &a.sam, &a.bam];
    if outs
        .iter()
        .filter(|o| o.as_ref().is_some_and(|p| p.as_os_str() == "-"))
        .count()
        > 1
    {
        bail!("only one output can go to stdout");
    }
    let meta_path = job.meta.clone().or_else(|| {
        outs.into_iter()
            .flatten()
            .find(|p| crate::map_cmd::is_file_output(p))
            .map(|p| runinfo::with_suffix(p, ".meta.json"))
    });
    let threads = runinfo::threads(job.threads);
    let invocation = Invocation::new(threads);
    let command_line = invocation.command_line.join(" ");

    let input = Some(job.reads.clone()).filter(|p| p.as_os_str() != "-");
    if let Some(p) = &input {
        std::fs::metadata(p).with_context(|| format!("reads {}", p.display()))?;
    }
    let t0 = Instant::now();
    let mut idx = GenomeIndex::open(&job.index).with_context(|| format!("index {}", job.index.display()))?;
    let annotation = match &a.annotation {
        Some(path) => {
            let s = idx.load_annotation(path)?;
            log::info!(
                "annotation: {} transcripts, {} introns, {} exons ({})",
                s.transcripts,
                s.introns,
                s.exons,
                s.format
            );
            if s.skipped_records > 0 {
                log::warn!(
                    "annotation: {} records are on sequences the index does not hold ({}{})",
                    s.skipped_records,
                    s.skipped_sequences.join(", "),
                    if s.skipped_sequences.len() >= 5 {
                        ", ..."
                    } else {
                        ""
                    }
                );
            }
            Some(AnnotationInfo {
                path: spearfish::provenance::absolute(path),
                format: s.format,
                transcripts: s.transcripts,
                introns: s.introns,
                exons: s.exons,
                skipped_records: s.skipped_records,
            })
        }
        None => None,
    };
    let index_info = IndexInfo {
        path: spearfish::provenance::absolute(&job.index),
        kind: "genome",
        index_id: idx.id(),
        sequences: idx.n_sequences(),
        bases: idx.total_len(),
        manifest: idx.manifest(),
    };
    let mapper = Arc::new(SplicedMapper::new(idx, &opts));
    let load = t0.elapsed().as_secs_f64();

    let open = |p: &Option<PathBuf>| p.as_deref().map(crate::map_cmd::sink).transpose();
    let paf = open(&paf_path)?.map(|o| (o, Vec::new()));
    let sam = open(&a.sam)?.map(|o| (o, Vec::new()));
    let bam = open(&a.bam)?.map(|o| (o, Vec::new(), Vec::new()));
    if let Some((out, _)) = &sam {
        let mut h = Vec::new();
        output::sam_header(&mut h, mapper.index(), &command_line);
        out.lock().unwrap().write_all(&h)?;
    }
    if let Some((out, _, _)) = &bam {
        let (mut h, mut z) = (Vec::new(), Vec::new());
        output::bam_header(&mut h, mapper.index(), &command_line);
        output::bgzf(&h, a.bam_level, &mut z);
        out.lock().unwrap().write_all(&z)?;
    }
    let total = Arc::new(Mutex::new((Counts::default(), Default::default())));
    let progress = Arc::new(AtomicU64::new(0));
    let mut proc = Processor {
        mapper,
        work: Workspace::default(),
        paf: paf.clone(),
        sam: sam.clone(),
        bam: bam.clone(),
        bam_level: a.bam_level,
        total: total.clone(),
        progress: progress.clone(),
        local: Counts::default(),
        batch_reads: 0,
    };
    let reader = spearfish::open_fastx(input.as_deref(), Some(job.batch_size))
        .with_context(|| format!("reading {}", job.reads.display()))?;
    let t1 = Instant::now();
    {
        let _p = crate::map_cmd::Progress::start(progress);
        reader.process_parallel(&mut proc, threads)?;
    }
    if let Some((out, _, _)) = &bam {
        out.lock().unwrap().write_all(&output::BGZF_EOF)?;
    }
    for out in [
        paf.as_ref().map(|x| &x.0),
        sam.as_ref().map(|x| &x.0),
        bam.as_ref().map(|x| &x.0),
    ]
    .into_iter()
    .flatten()
    {
        out.lock().unwrap().flush()?;
    }
    let secs = t1.elapsed().as_secs_f64();
    let (c, internal) = *total.lock().unwrap();
    log::info!(
        "map: {} reads ({:.1}% mapped, {:.1}% of them spliced), {:.1} Mbp in {secs:.2} s ({:.0} reads/s, {:.1} Mbp/s); index load {load:.2} s",
        c.reads,
        100.0 * c.mapped as f64 / c.reads.max(1) as f64,
        100.0 * c.spliced as f64 / c.mapped.max(1) as f64,
        c.bases as f64 / 1e6,
        c.reads as f64 / secs,
        c.bases as f64 / 1e6 / secs
    );
    let internal = std::env::var_os("SPEARFISH_STATS").map(|_| {
        let s = format!("{internal:?}");
        log::info!("stats\t{s}");
        s
    });
    if let Some(m) = &meta_path {
        let abs = |p: &Path| {
            if p.as_os_str() == "-" {
                "-".to_string()
            } else {
                spearfish::provenance::absolute(p)
            }
        };
        let show = |p: &Option<PathBuf>| p.as_deref().map(abs);
        let record = RunRecord {
            kind: "spearfish-genome-map-run",
            invocation: &invocation,
            index: index_info,
            annotation,
            reads: abs(&job.reads),
            outputs: Outputs {
                paf: show(&paf_path),
                sam: show(&a.sam),
                bam: show(&a.bam),
            },
            preset: job.preset.name(),
            parameters: &opts,
            bam_deflate: output::deflate_engine(),
            finished_unix: runinfo::unix_now(),
            stats: RunStats {
                counts: c,
                mean_read_length: c.bases as f64 / c.reads.max(1) as f64,
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

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[derive(Parser)]
    struct Cli {
        #[command(flatten)]
        a: GenomeArgs,
    }

    /// `spearfish map-spliced` without options maps as the library does with a
    /// preset's defaults.
    #[test]
    fn cli_defaults_are_library_defaults() {
        let a = Cli::try_parse_from(["x"]).unwrap().a;
        for p in [Preset::Ont, Preset::Drna, Preset::Hifi] {
            assert_eq!(a.options(p), SpliceOptions::preset(p));
        }
        assert_eq!(SpliceOptions::preset(Preset::Ont), SpliceOptions::default());
        let a = Cli::try_parse_from(["x", "--secondary", "4", "--base-level", "affine"])
            .unwrap()
            .a;
        assert_eq!(a.options(Preset::Ont).max_hits, 5);
        assert_eq!(a.options(Preset::Ont).base_level, BaseLevel::Affine);
    }
}
