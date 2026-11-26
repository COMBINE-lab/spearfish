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
use std::sync::Arc;
use std::sync::atomic::AtomicU64;
use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::ValueEnum;
use serde::Serialize;
use spearfish::Preset;
use spearfish::genome::{
    BaseLevel, GenomeIndex, GenomeIndexOptions, GenomeManifest, SpliceOptions, SplicedMapper, output,
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
    /// Known transcripts (GTF, GFF3 or BED12) or known introns (BED, one intron per line with the
    /// strand in column 6, as `minimap2 --write-junc` writes), optionally compressed. Optional: their splice
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
    /// Leave reads that do not map out of SAM and BAM (PAF never has them).
    #[arg(long, hide_short_help = true, help_heading = IO)]
    pub no_unmapped: bool,
    /// Leave base qualities out of SAM and BAM (smaller files, written faster; FASTA reads
    /// have none).
    #[arg(long, hide_short_help = true, help_heading = IO)]
    pub no_qualities: bool,
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

/// Read counts of a run.
pub type Counts = spearfish::genome::FileCounts;

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
    intron_records: usize,
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
            if s.intron_records > 0 {
                log::info!("annotation: {} lines of introns alone", s.intron_records);
            }
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
                intron_records: s.intron_records,
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

    let open = |p: &Option<PathBuf>| -> Result<Option<Box<dyn Write + Send>>> {
        Ok(match p {
            Some(p) if p.as_os_str() == "-" => Some(Box::new(std::io::stdout())),
            Some(p) => Some(Box::new(
                std::fs::File::create(p).with_context(|| format!("creating {}", p.display()))?,
            )),
            None => None,
        })
    };
    let outputs = spearfish::genome::FileOutputs {
        paf: open(&paf_path)?,
        sam: open(&a.sam)?,
        bam: open(&a.bam)?,
        bam_level: a.bam_level,
        unmapped: !a.no_unmapped,
        qualities: !a.no_qualities,
        command_line,
    };
    let progress = Arc::new(AtomicU64::new(0));
    let t1 = Instant::now();
    let (c, internal) = {
        let _p = crate::map_cmd::Progress::start(progress.clone());
        spearfish::genome::map_file(
            &mapper,
            input.as_deref(),
            outputs,
            threads,
            job.batch_size,
            Some(progress),
        )
        .with_context(|| format!("mapping {}", job.reads.display()))?
    };
    let secs = t1.elapsed().as_secs_f64();
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
