//! Spliced mapping of long RNA-seq reads to a genome (feature `genome`, on by default).
//!
//! A read is placed as a series of exon blocks separated by introns. The mapper anchors exons
//! with sparse seeds and then settles, from the few read bases beside each one, where every
//! junction is; it aligns no read base across an intron. An annotation is optional: without
//! one, junctions are placed by the read and the splice-site motifs alone.
//!
//! ```no_run
//! use spearfish::Preset;
//! use spearfish::genome::{GenomeIndex, SpliceOptions, SplicedMapper};
//!
//! let mut index = GenomeIndex::open("genome.sfg".as_ref())?;
//! index.load_annotation("genes.gtf".as_ref())?; // optional: GTF, GFF3 or BED12
//! let mapper = SplicedMapper::new(index, &SpliceOptions::preset(Preset::Ont));
//! let read = mapper.map(b"ACGT...");
//! if let Some(best) = read.alignments.first() {
//!     let name = mapper.index().seq_name(best.seq);
//!     println!("{name}:{}-{} mapq {}", best.target_start, best.target_end, best.mapq);
//!     for b in best.blocks() {
//!         println!("  exon {}-{}", b.target_start, b.target_end);
//!     }
//! }
//! # Ok::<(), anyhow::Error>(())
//! ```
//!
//! [`SplicedMapper`] is `Send + Sync`: share one across threads. Coordinates are 0-based and
//! half-open, on the named sequence. [`output`] writes alignments as PAF, SAM and BAM, as
//! `spearfish map-spliced` does.
//!
//! A genome index and a transcriptome index ([`crate::Index`]) are different files in the same
//! container; [`crate::index_kind`] tells which one a path holds.

use std::cell::RefCell;
use std::io::BufRead;
use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use splice::base::{Aligner, Mode};
use splice::index::BuildOptions;
use splice::junctions::Junctions;
use splice::map::{MapOptions, Scratch, SplicedHit};

use crate::Preset;
use crate::provenance::{self, BuildInfo};

pub use splice::junctions::{Model as TranscriptModel, Summary as AnnotationSummary};
pub use splice::map::Stats as InternalStats;

/// How a genome index is built.
///
/// ```
/// let opts = spearfish::genome::GenomeIndexOptions::builder().max_occ(500).build();
/// assert!(!opts.hpc);
/// ```
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, bon::Builder)]
#[serde(default)]
#[non_exhaustive]
pub struct GenomeIndexOptions {
    /// Seeds that occur more often than this in the genome (repeats) are left out of the seed
    /// table.
    #[builder(default = 1000)]
    pub max_occ: u32,
    /// Seed the homopolymer-compressed genome: longer, sparser seeds in a smaller index. Off
    /// by default: seeds of the sequence as it is find more loci and more exons, in less time.
    #[builder(default)]
    pub hpc: bool,
}

impl Default for GenomeIndexOptions {
    fn default() -> Self {
        GenomeIndexOptions::builder().build()
    }
}

/// What built a genome index and from what: stored in the index and written next to it as
/// `<index>.json`. Fields missing from a manifest read back take neutral defaults.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct GenomeManifest {
    /// Always `spearfish-genome-index`.
    pub kind: String,
    /// Container format version.
    pub format_version: u32,
    /// Identifier of the contents: a hash of the sequence names and lengths, the input's hash
    /// and the build options.
    pub index_id: String,
    /// The build that wrote the index.
    pub built_by: BuildInfo,
    /// Seconds since the Unix epoch when the build finished.
    pub created_unix: u64,
    /// Command line, if built by the `spearfish` binary.
    pub command_line: Vec<String>,
    /// Directory the command ran in.
    pub working_directory: String,
    /// The build options.
    pub options: GenomeIndexOptions,
    /// LCP parser version and seed, which the seeds depend on.
    pub lcp_parser_version: u32,
    /// See `lcp_parser_version`.
    pub lcp_seed: String,
    /// Genome FASTA, as an absolute path where it can be resolved (empty when built from a
    /// reader), and the xxh3-64 of its bytes as stored.
    pub input: String,
    /// See `input`.
    pub input_xxh3: String,
    /// Sequences and bases in the index.
    pub sequences: usize,
    /// See `sequences`.
    pub bases: u64,
    /// Seconds the build took (reading the FASTA included).
    pub build_seconds: f64,
}

/// Where a genome index was built from, for its manifest.
#[derive(Clone, Debug, Default)]
pub struct GenomeInput {
    /// The FASTA's path, as recorded.
    pub path: String,
    /// xxh3-64 of the FASTA's bytes, as 16 hex digits.
    pub xxh3: String,
    /// The command line that asked for the build.
    pub command_line: Vec<String>,
}

/// A buffered reader over `path`: gzip, bzip2 and zstd are decompressed (feature `fasta`).
fn open_text(path: &Path) -> Result<Box<dyn BufRead>> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    #[cfg(feature = "fasta")]
    {
        let (stream, _) = niffler::get_reader(Box::new(std::io::BufReader::new(f))).map_err(|e| match e {
            niffler::Error::FeatureDisabled => anyhow::anyhow!(
                "{}: xz-compressed input is not supported: decompress it, or recompress it with gzip or zstd",
                path.display()
            ),
            e => anyhow::Error::from(e).context(format!("reading {}", path.display())),
        })?;
        Ok(Box::new(std::io::BufReader::with_capacity(1 << 22, stream)))
    }
    #[cfg(not(feature = "fasta"))]
    Ok(Box::new(std::io::BufReader::with_capacity(1 << 22, f)))
}

/// A genome index (memory-mapped when opened from a file), with the annotation loaded for
/// this run, if any. Cloning shares the index.
#[derive(Clone)]
pub struct GenomeIndex {
    inner: Arc<splice::index::GenomeIndex>,
    annotation: Option<AnnotationSummary>,
}

impl GenomeIndex {
    /// Open an index written by `spearfish index --genome` (memory-mapped).
    pub fn open(path: &Path) -> Result<Self> {
        if !splice::index::GenomeIndex::is_genome_index(path)? {
            bail!(
                "{} is a transcriptome index: open it with spearfish::Index, or build a genome index with `spearfish index --genome`",
                path.display()
            );
        }
        Ok(Self::wrap(splice::index::GenomeIndex::open(path)?))
    }

    fn wrap(inner: splice::index::GenomeIndex) -> Self {
        GenomeIndex {
            inner: Arc::new(inner),
            annotation: None,
        }
    }

    /// Build an index in memory from a genome FASTA (optionally gzip, bzip2 or zstd
    /// compressed, with feature `fasta`). Uses rayon's current thread pool.
    pub fn build_from_fasta(path: &Path, opts: &GenomeIndexOptions) -> Result<Self> {
        Self::build_from_reader(open_text(path)?, opts).with_context(|| format!("reading {}", path.display()))
    }

    /// Build an index in memory from FASTA text.
    pub fn build_from_reader(fasta: impl BufRead, opts: &GenomeIndexOptions) -> Result<Self> {
        let o = BuildOptions {
            max_occ: opts.max_occ,
            hpc: opts.hpc,
        };
        Ok(Self::wrap(splice::index::GenomeIndex::build_from_reader(
            fasta, &o,
        )?))
    }

    /// Build an index from a genome FASTA and write it to `out`, with its manifest inside and
    /// next to it (`<out>.json`). Returns the index, the bytes written and the manifest.
    pub fn build_to_file(
        fasta: &Path,
        opts: &GenomeIndexOptions,
        out: &Path,
        command_line: Vec<String>,
    ) -> Result<(Self, u64, GenomeManifest)> {
        let t = std::time::Instant::now();
        let input = GenomeInput {
            path: provenance::absolute(fasta),
            xxh3: provenance::file_xxh3(fasta).with_context(|| format!("reading {}", fasta.display()))?,
            command_line,
        };
        let o = BuildOptions {
            max_occ: opts.max_occ,
            hpc: opts.hpc,
        };
        // The index is written table by table as it is built, which takes much less memory
        // than building it whole and writing it then.
        let mut manifest = None;
        let (inner, bytes) =
            splice::index::GenomeIndex::build_to_file(open_text(fasta)?, &o, out, |names, offsets| {
                let m = Self::manifest_of(names, offsets, opts, input, t.elapsed().as_secs_f64());
                let json = serde_json::to_vec_pretty(&m)?;
                manifest = Some(m);
                Ok(Some(json))
            })
            .with_context(|| format!("building {} from {}", out.display(), fasta.display()))?;
        let manifest = manifest.expect("made with the index");
        let side = provenance::sidecar(out);
        std::fs::write(&side, serde_json::to_vec_pretty(&manifest)?)
            .with_context(|| format!("writing {}", side.display()))?;
        Ok((Self::wrap(inner), bytes, manifest))
    }

    /// The manifest of this index as just built with `opts` from `input`.
    pub fn new_manifest(
        &self,
        opts: &GenomeIndexOptions,
        input: GenomeInput,
        build_seconds: f64,
    ) -> GenomeManifest {
        Self::manifest_of(
            self.inner.names(),
            self.inner.offsets(),
            opts,
            input,
            build_seconds,
        )
    }

    /// The manifest of an index of sequences `names` that start at `offsets` (the total
    /// length last).
    fn manifest_of(
        names: &[String],
        offsets: &[u32],
        opts: &GenomeIndexOptions,
        input: GenomeInput,
        build_seconds: f64,
    ) -> GenomeManifest {
        let mut h = xxhash_rust::xxh3::Xxh3::new();
        for (name, w) in names.iter().zip(offsets.windows(2)) {
            h.update(name.as_bytes());
            h.update(&(w[1] - w[0]).to_le_bytes());
        }
        h.update(input.xxh3.as_bytes());
        h.update(&opts.max_occ.to_le_bytes());
        h.update(&[opts.hpc as u8]);
        GenomeManifest {
            kind: "spearfish-genome-index".into(),
            format_version: index::store::FORMAT_VERSION,
            index_id: format!("{:016x}", h.digest()),
            built_by: BuildInfo::current(),
            created_unix: std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |d| d.as_secs()),
            command_line: input.command_line,
            working_directory: std::env::current_dir()
                .map(|d| d.display().to_string())
                .unwrap_or_default(),
            options: opts.clone(),
            lcp_parser_version: lcp_core::LCP_PARSER_VERSION,
            lcp_seed: format!("{:#x}", lcp_core::LcpConfig::default().seed),
            input: input.path,
            input_xxh3: input.xxh3,
            sequences: names.len(),
            bases: offsets.last().copied().unwrap_or(0) as u64,
            build_seconds,
        }
    }

    /// Write the index to `path`, with `manifest` inside it and as `<path>.json` when given.
    /// Returns the bytes written.
    pub fn write(&self, path: &Path, manifest: Option<&GenomeManifest>) -> Result<u64> {
        let json = manifest.map(serde_json::to_vec_pretty).transpose()?;
        let bytes = self
            .inner
            .write_with_manifest(path, json.as_deref())
            .with_context(|| format!("writing {}", path.display()))?;
        if let Some(j) = &json {
            let side = provenance::sidecar(path);
            std::fs::write(&side, j).with_context(|| format!("writing {}", side.display()))?;
        }
        Ok(bytes)
    }

    /// The manifest stored in the index, if it has one that can be read.
    pub fn manifest(&self) -> Option<GenomeManifest> {
        serde_json::from_slice(self.inner.manifest()?)
            .map_err(|e| log::debug!("the genome index manifest can't be read: {e}"))
            .ok()
    }

    /// The manifest's identifier of the contents, or `unknown`.
    pub fn id(&self) -> String {
        self.manifest().map_or_else(|| "unknown".into(), |m| m.index_id)
    }

    /// Number of reference sequences.
    pub fn n_sequences(&self) -> usize {
        self.inner.names().len()
    }

    /// Name of sequence `i` (file order).
    pub fn seq_name(&self, i: u32) -> &str {
        &self.inner.names()[i as usize]
    }

    /// Length of sequence `i` in bases.
    pub fn seq_len(&self, i: u32) -> u32 {
        self.inner.seq_len(i as usize)
    }

    /// Bases in all sequences.
    pub fn total_len(&self) -> u64 {
        self.inner.total_len() as u64
    }

    /// Whether the seeds are of the homopolymer-compressed genome.
    pub fn hpc(&self) -> bool {
        self.inner.hpc()
    }

    /// The index's sections and their sizes in bytes.
    pub fn sizes(&self) -> Vec<(String, usize)> {
        self.inner.sizes()
    }

    /// The underlying index (the `spearfish-splice` crate's; its interface is not stable).
    pub fn raw(&self) -> &Arc<splice::index::GenomeIndex> {
        &self.inner
    }

    /// Load known transcript models for this run from a GTF, GFF3 or BED12 file (optionally
    /// compressed; the format is told from the content). Known splice sites then cost nothing
    /// where others pay for their motif, and known exons and junctions are candidates where
    /// read bases are left unexplained. Nothing is written to the index. Replaces an
    /// annotation loaded before; must be called before the index is given to a mapper.
    ///
    /// An error when no transcript in the file is on a sequence of the index (sequence names
    /// must match: `chr1` is not `1`).
    pub fn load_annotation(&mut self, path: &Path) -> Result<&AnnotationSummary> {
        let (j, sum) = Junctions::from_reader(open_text(path)?, &self.inner)
            .with_context(|| format!("reading {}", path.display()))?;
        if sum.transcripts == 0 {
            let names = sum.skipped_sequences.join(", ");
            bail!(
                "{}: no transcript is on a sequence of the index{}",
                path.display(),
                if names.is_empty() {
                    " (the file has no exon records)".to_string()
                } else {
                    format!(
                        " (it names {names}, ...; the index has {}, ...)",
                        self.seq_name(0)
                    )
                }
            );
        }
        self.set(j, sum)
    }

    /// As [`GenomeIndex::load_annotation`], from transcript models the caller already holds.
    pub fn set_annotation(
        &mut self,
        models: impl IntoIterator<Item = TranscriptModel>,
    ) -> Result<&AnnotationSummary> {
        let mut n = 0;
        let nseq = self.n_sequences();
        let checked: Vec<TranscriptModel> = models.into_iter().inspect(|_| n += 1).collect();
        for m in &checked {
            if m.seq >= nseq
                || m.exons
                    .iter()
                    .any(|&(s, e)| s >= e || e > self.seq_len(m.seq as u32))
            {
                bail!(
                    "a transcript model is not on the index's sequences (sequence {})",
                    m.seq
                );
            }
        }
        let j = Junctions::from_models(checked, &self.inner);
        let sum = AnnotationSummary {
            format: "models",
            transcripts: n,
            introns: j.len(),
            ..Default::default()
        };
        self.set(j, sum)
    }

    fn set(&mut self, j: Junctions, sum: AnnotationSummary) -> Result<&AnnotationSummary> {
        let Some(inner) = Arc::get_mut(&mut self.inner) else {
            bail!("the annotation must be loaded before the index is shared with a mapper or cloned");
        };
        inner.junctions = Some(j);
        Ok(self.annotation.insert(sum))
    }

    /// What the loaded annotation held, if one was loaded.
    pub fn annotation(&self) -> Option<&AnnotationSummary> {
        self.annotation.as_ref()
    }
}

/// How aligned blocks are reported.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
#[non_exhaustive]
pub enum BaseLevel {
    /// Exon blocks only: each block is one `M` run in a CIGAR (with one `I` or `D` when the
    /// read and genome lengths of the block differ). Enough for quantification and for
    /// tools that read exon structure; the fastest.
    #[default]
    Blocks,
    /// Each block aligned base by base with the fewest edits.
    Edits,
    /// As `Edits`, then every stretch with an indel aligned again with affine gap costs (as
    /// minimap2's long-read scoring), which places indels where a base-level consumer expects
    /// them.
    Affine,
}

/// Every setting of the spliced mapper. [`SpliceOptions::preset`] gives a preset's values;
/// the builder starts from the `ont` preset.
///
/// ```
/// use spearfish::{Preset, genome::SpliceOptions};
/// let mut opts = SpliceOptions::preset(Preset::Drna);
/// assert!(opts.sense_reads);
/// opts.max_hits = 5;
/// let opts = SpliceOptions::builder().max_intron(500_000).build();
/// assert!(!opts.sense_reads);
/// ```
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, bon::Builder)]
#[serde(default)]
#[non_exhaustive]
pub struct SpliceOptions {
    /// The reads are the transcripts' sense strand (direct RNA): a read mapped forward is of
    /// a forward-strand transcript, and splice sites are sought on that strand only.
    #[builder(default)]
    pub sense_reads: bool,
    /// Longest intron.
    #[builder(default = 200_000)]
    pub max_intron: u32,
    /// A reference skip of at least this many bases between two anchors is a junction; a
    /// shorter one is a deletion.
    #[builder(default = 40)]
    pub min_intron: u32,
    /// Alignments reported per read: the primary and up to `max_hits - 1` secondary ones.
    #[builder(default = 1)]
    pub max_hits: usize,
    /// How aligned blocks are reported (base-level alignment is of the primary only).
    #[builder(default)]
    pub base_level: BaseLevel,
    /// Aligned hits that score below this are dropped (a lone, well-anchored hit is not
    /// aligned at all and always kept).
    #[builder(default = 100)]
    pub min_score: i32,
    /// Seeds with more occurrences than this are not used, unless the read gets no hit
    /// without them, when the limit is `seed_cap_retry`.
    #[builder(default = 64)]
    pub seed_cap: usize,
    /// See `seed_cap`.
    #[builder(default = 1000)]
    pub seed_cap_retry: usize,
    /// Loci chained per read.
    #[builder(default = 6)]
    pub max_loci: usize,
    /// Seed every read with pairs of shorter (level-2) seeds as well, which reads with many
    /// errors need; otherwise only reads the main seeds anchor weakly get them.
    #[builder(default = true)]
    pub pair_seeds: bool,
    /// Junction cost per read base (in 64ths) from which unexplained read bases are searched
    /// for a missed exon.
    #[builder(default = 16)]
    pub exon_search_trigger: u32,
}

impl Default for SpliceOptions {
    fn default() -> Self {
        SpliceOptions::builder().build()
    }
}

impl SpliceOptions {
    /// The settings of a preset: `drna` (and `drna-fast`) reads are the transcripts' sense
    /// strand; the other presets differ only for a transcriptome.
    pub fn preset(p: Preset) -> Self {
        SpliceOptions {
            sense_reads: p.name().starts_with("drna"),
            ..Default::default()
        }
    }

    fn map_options(&self) -> MapOptions {
        MapOptions {
            cap: self.seed_cap,
            cap_retry: self.seed_cap_retry,
            max_intron: self.max_intron,
            min_intron: self.min_intron,
            max_loci: self.max_loci,
            forward_only: self.sense_reads,
            pairs: self.pair_seeds,
            trigger_x64: self.exon_search_trigger,
            min_score: self.min_score,
            keep_seeds: self.base_level != BaseLevel::Blocks,
        }
    }
}

/// The strand of the transcript a read came from, told from its splice sites.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum TranscriptStrand {
    /// The genome's forward strand.
    Forward,
    /// The genome's reverse strand.
    Reverse,
    /// No evidence: the alignment has no junction, or its junctions do not agree.
    Unknown,
}

/// One aligned block of an alignment: read `[query_start, query_end)` against the sequence's
/// `[target_start, target_end)`. Read coordinates are in the orientation the read maps in
/// (of the reverse complement when [`SplicedAlignment::reverse`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Block {
    /// Start in the read, in the mapped orientation.
    pub query_start: u32,
    /// End in the read (exclusive).
    pub query_end: u32,
    /// Start on the reference sequence.
    pub target_start: u32,
    /// End on the reference sequence (exclusive).
    pub target_end: u32,
}

/// One CIGAR operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    /// Aligned bases, matching or not (`M`).
    Match(u32),
    /// Read bases absent from the genome (`I`).
    Insertion(u32),
    /// Genome bases absent from the read (`D`).
    Deletion(u32),
    /// An intron (`N`).
    Intron(u32),
    /// Unaligned read bases at an end (`S`).
    SoftClip(u32),
}

impl CigarOp {
    fn from_bam(op: u32) -> CigarOp {
        let n = op >> 4;
        match op & 15 {
            0 => CigarOp::Match(n),
            1 => CigarOp::Insertion(n),
            2 => CigarOp::Deletion(n),
            3 => CigarOp::Intron(n),
            _ => CigarOp::SoftClip(n),
        }
    }

    /// The operation's length and its SAM letter.
    pub fn len_and_letter(self) -> (u32, char) {
        match self {
            CigarOp::Match(n) => (n, 'M'),
            CigarOp::Insertion(n) => (n, 'I'),
            CigarOp::Deletion(n) => (n, 'D'),
            CigarOp::Intron(n) => (n, 'N'),
            CigarOp::SoftClip(n) => (n, 'S'),
        }
    }
}

/// One alignment of a read: a locus, a strand, and exon blocks.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct SplicedAlignment {
    /// The reference sequence (index order; [`GenomeIndex::seq_name`]).
    pub seq: u32,
    /// Start of the first block on the sequence.
    pub target_start: u32,
    /// End of the last block on the sequence (exclusive).
    pub target_end: u32,
    /// Aligned span on the read as sequenced (its forward strand).
    pub query_start: u32,
    /// End of the aligned span on the read (exclusive).
    pub query_end: u32,
    /// The read maps as its reverse complement.
    pub reverse: bool,
    /// The transcript's strand, from the splice sites.
    pub transcript_strand: TranscriptStrand,
    /// Mapping quality, 0 to 60: 60 for a read with one hit, less as the second-best
    /// approaches the best. The same for every alignment of a read.
    pub mapq: u8,
    /// Score; higher is better. Comparable between the alignments of one read.
    pub score: i32,
    /// The best alignment of the read (the others are secondary).
    pub primary: bool,
    /// Seed matches in the chain.
    pub anchors: u32,
    /// Edits in the base-level alignment, when one was asked for and made.
    pub edits: Option<u32>,
    hit: SplicedHit,
    /// Start of `seq` on the index's concatenated axis.
    offset: u32,
    /// The base-level CIGAR as BAM encodes it (empty when there is none).
    cigar: Vec<u32>,
}

impl SplicedAlignment {
    /// The aligned blocks, in genome order.
    pub fn blocks(&self) -> impl ExactSizeIterator<Item = Block> + '_ {
        self.hit.blocks.iter().map(|b| Block {
            query_start: b.qs,
            query_end: b.qe,
            target_start: b.ts - self.offset,
            target_end: b.te - self.offset,
        })
    }

    /// The introns as `(start, end)` on the sequence, in genome order.
    pub fn introns(&self) -> impl Iterator<Item = (u32, u32)> + '_ {
        self.hit
            .blocks
            .windows(2)
            .map(|w| (w[0].te - self.offset, w[1].ts - self.offset))
    }

    /// Whether the alignment has an intron.
    pub fn spliced(&self) -> bool {
        self.hit.blocks.len() > 1
    }

    /// The CIGAR for a read of `read_len` bases: the base-level one when it was made
    /// ([`SpliceOptions::base_level`]), otherwise one `M` run per block.
    pub fn cigar(&self, read_len: u32) -> Vec<CigarOp> {
        if !self.cigar.is_empty() {
            return self.cigar.iter().map(|&op| CigarOp::from_bam(op)).collect();
        }
        let mut ops = Vec::with_capacity(2 * self.hit.blocks.len() + 2);
        let (first, last) = (self.hit.blocks[0], self.hit.blocks[self.hit.blocks.len() - 1]);
        if first.qs > 0 {
            ops.push(CigarOp::SoftClip(first.qs));
        }
        for (i, b) in self.hit.blocks.iter().enumerate() {
            if i > 0 {
                ops.push(CigarOp::Intron(b.ts - self.hit.blocks[i - 1].te));
            }
            let (lq, lt) = (b.qe - b.qs, b.te - b.ts);
            if lq.min(lt) > 0 {
                ops.push(CigarOp::Match(lq.min(lt)));
            }
            if lq > lt {
                ops.push(CigarOp::Insertion(lq - lt));
            } else if lt > lq {
                ops.push(CigarOp::Deletion(lt - lq));
            }
        }
        if last.qe < read_len {
            ops.push(CigarOp::SoftClip(read_len - last.qe));
        }
        ops
    }
}

/// All alignments of one read, best first; none when the read does not map.
#[derive(Clone, Debug, Default)]
#[non_exhaustive]
pub struct SplicedRead {
    /// Read length in bases.
    pub len: u32,
    /// The reported alignments: the primary, then the secondary ones.
    pub alignments: Vec<SplicedAlignment>,
}

/// A thread's working memory for [`SplicedMapper::map_with`]; reuse one per thread.
#[derive(Default)]
pub struct Workspace {
    scratch: Scratch,
    aligner: Aligner,
    ops: Vec<u32>,
}

impl Workspace {
    /// Internal counters of the reads mapped with this workspace (their meaning is not
    /// stable).
    pub fn stats(&self) -> &InternalStats {
        &self.scratch.stats
    }
}

thread_local! {
    static WORKSPACE: RefCell<Workspace> = RefCell::new(Workspace::default());
}

/// A spliced mapper over one genome index; `Send + Sync`.
pub struct SplicedMapper {
    index: GenomeIndex,
    inner: splice::map::Mapper,
    opts: SpliceOptions,
}

impl SplicedMapper {
    /// A mapper over `index` (with the annotation loaded into it, if any) with settings
    /// `opts`.
    pub fn new(index: GenomeIndex, opts: &SpliceOptions) -> Self {
        let inner = splice::map::Mapper::new(index.inner.clone(), opts.map_options());
        SplicedMapper {
            index,
            inner,
            opts: opts.clone(),
        }
    }

    /// The index mapped against.
    pub fn index(&self) -> &GenomeIndex {
        &self.index
    }

    /// The settings the mapper runs with.
    pub fn options(&self) -> &SpliceOptions {
        &self.opts
    }

    /// Map one read (ASCII bases; `U` is read as `T`), with working memory kept per thread.
    pub fn map(&self, seq: &[u8]) -> SplicedRead {
        WORKSPACE.with_borrow_mut(|w| self.map_with(seq, w))
    }

    /// Map one read with the caller's working memory.
    pub fn map_with(&self, seq: &[u8], w: &mut Workspace) -> SplicedRead {
        let ix = &*self.index.inner;
        let hits = self.inner.map(seq, &mut w.scratch);
        let mapq = splice::map::mapq(&hits);
        let len = seq.len() as u32;
        w.aligner.mode = match self.opts.base_level {
            BaseLevel::Affine => Mode::Affine,
            _ => Mode::Edits,
        };
        let base_level = self.opts.base_level != BaseLevel::Blocks;
        let alignments = hits
            .into_iter()
            .take(self.opts.max_hits.max(1))
            .map(|hit| {
                let (first, last) = (hit.blocks[0], hit.blocks[hit.blocks.len() - 1]);
                let (rid, start) = ix.locate(first.ts);
                let offset = first.ts - start;
                // (secondary alignments keep their blocks)
                let (cigar, edits) = if base_level && hit.primary {
                    let e = w.aligner.cigar(ix, seq, &hit, &mut w.ops);
                    (w.ops.clone(), Some(e))
                } else {
                    (Vec::new(), None)
                };
                let (query_start, query_end) = if hit.reverse {
                    (len - last.qe, len - first.qs)
                } else {
                    (first.qs, last.qe)
                };
                SplicedAlignment {
                    seq: rid as u32,
                    target_start: start,
                    target_end: last.te - offset,
                    query_start,
                    query_end,
                    reverse: hit.reverse,
                    transcript_strand: match hit.tx_strand {
                        b'+' => TranscriptStrand::Forward,
                        b'-' => TranscriptStrand::Reverse,
                        _ => TranscriptStrand::Unknown,
                    },
                    mapq,
                    score: hit.score,
                    primary: hit.primary,
                    anchors: hit.anchors,
                    edits,
                    hit,
                    offset,
                    cigar,
                }
            })
            .collect();
        SplicedRead { len, alignments }
    }
}

// ─── Mapping a file ─────────────────────────────────────────────────────────────────────────

/// Where [`map_file`] writes, and what. Any combination of the three formats; each sink
/// gets the whole file of its format, header included.
#[cfg(feature = "fasta")]
#[derive(Default)]
pub struct FileOutputs {
    /// PAF (mapped reads only).
    pub paf: Option<Box<dyn std::io::Write + Send>>,
    /// SAM.
    pub sam: Option<Box<dyn std::io::Write + Send>>,
    /// BAM, unsorted.
    pub bam: Option<Box<dyn std::io::Write + Send>>,
    /// BAM compression level, 0 to 9.
    pub bam_level: u8,
    /// Write a record for each read that does not map (SAM and BAM).
    pub unmapped: bool,
    /// Write base qualities when the reads have them (SAM and BAM).
    pub qualities: bool,
    /// The command line for the `@PG` header line (may be empty).
    pub command_line: String,
}

#[cfg(feature = "fasta")]
use std::io::Write as _;
#[cfg(feature = "fasta")]
use std::sync::Mutex;
#[cfg(feature = "fasta")]
use std::sync::atomic::{AtomicU64, Ordering};

/// Read counts of a file mapped by [`map_file`].
#[cfg(feature = "fasta")]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
#[non_exhaustive]
pub struct FileCounts {
    /// Reads in the file.
    pub reads: u64,
    /// Reads with at least one alignment.
    pub mapped: u64,
    /// Reads with none.
    pub unmapped: u64,
    /// Mapped reads whose primary alignment has an intron.
    pub spliced: u64,
    /// Introns in primary alignments.
    pub junctions: u64,
    /// Mapped reads with more than one reported alignment.
    pub multi_mapped: u64,
    /// Mapped reads whose primary alignment is on the read's forward strand.
    pub forward: u64,
    /// Alignments reported.
    pub alignments: u64,
    /// Read bases.
    pub bases: u64,
}

#[cfg(feature = "fasta")]
impl FileCounts {
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

    fn add(&mut self, o: &FileCounts) {
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

#[cfg(feature = "fasta")]
type Out = Arc<Mutex<Box<dyn std::io::Write + Send>>>;

#[cfg(feature = "fasta")]
struct Processor {
    mapper: Arc<SplicedMapper>,
    work: Workspace,
    paf: Option<(Out, Vec<u8>)>,
    sam: Option<(Out, Vec<u8>)>,
    /// BAM: the sink, the batch's records, and their BGZF blocks.
    bam: Option<(Out, Vec<u8>, Vec<u8>)>,
    bam_level: u8,
    /// Write unmapped reads, and base qualities, to SAM and BAM.
    unmapped: bool,
    qualities: bool,
    total: Arc<Mutex<(FileCounts, InternalStats)>>,
    progress: Option<Arc<AtomicU64>>,
    local: FileCounts,
    batch_reads: u64,
}

#[cfg(feature = "fasta")]
impl Clone for Processor {
    fn clone(&self) -> Self {
        Processor {
            mapper: self.mapper.clone(),
            work: Workspace::default(),
            paf: self.paf.clone(),
            sam: self.sam.clone(),
            bam: self.bam.clone(),
            bam_level: self.bam_level,
            unmapped: self.unmapped,
            qualities: self.qualities,
            total: self.total.clone(),
            progress: self.progress.clone(),
            local: FileCounts::default(),
            batch_reads: 0,
        }
    }
}

#[cfg(feature = "fasta")]
impl<Rf: paraseq::Record> paraseq::prelude::ParallelProcessor<Rf> for Processor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        let seq = record.seq();
        let r = self.mapper.map_with(&seq, &mut self.work);
        self.local.add_read(&r);
        self.batch_reads += 1;
        // A read without an alignment is written only to SAM and BAM, and only when asked.
        if r.alignments.is_empty() && (!self.unmapped || (self.sam.is_none() && self.bam.is_none())) {
            return Ok(());
        }
        // Lossy: a header that is not UTF-8 gets U+FFFD instead of failing the run.
        let id = String::from_utf8_lossy(record.id());
        let name = id.split_whitespace().next().unwrap_or("");
        let idx = self.mapper.index();
        let qual = record.qual().filter(|_| self.qualities);
        if r.alignments.is_empty() {
            if let Some((_, buf)) = &mut self.sam {
                output::sam_unmapped(buf, name, &seq, qual);
            }
            if let Some((_, raw, _)) = &mut self.bam {
                output::bam_unmapped(raw, name, &seq, qual);
            }
            return Ok(());
        }
        if let Some((_, buf)) = &mut self.paf {
            output::paf_lines(buf, name, &r, idx);
        }
        if let Some((_, buf)) = &mut self.sam {
            output::sam_records(buf, name, &seq, qual, &r, idx);
        }
        if let Some((_, raw, _)) = &mut self.bam {
            output::bam_records(raw, name, &seq, qual, &r, idx);
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
        if let Some(p) = &self.progress {
            p.fetch_add(self.batch_reads, Ordering::Relaxed);
        }
        self.batch_reads = 0;
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        let mut t = self.total.lock().unwrap();
        t.0.add(&self.local);
        t.1.add(self.work.stats());
        self.local = FileCounts::default();
        Ok(())
    }
}

/// Map a FASTA/FASTQ file (optionally gzip, bzip2 or zstd compressed; `None` reads stdin) on
/// `threads` threads and write the alignments to `out`. Reads are written as they finish,
/// not in input order; the records of a read are adjacent. BAM records are compressed by the
/// thread that mapped them. `progress`, when given, counts the reads done as they finish.
/// Returns the read counts and the mapper's internal counters.
#[cfg(feature = "fasta")]
pub fn map_file(
    mapper: &Arc<SplicedMapper>,
    reads: Option<&Path>,
    out: FileOutputs,
    threads: usize,
    batch_size: usize,
    progress: Option<Arc<AtomicU64>>,
) -> Result<(FileCounts, InternalStats)> {
    use paraseq::prelude::*;
    let share = |w: Option<Box<dyn std::io::Write + Send>>| -> Option<Out> {
        w.map(|w| {
            let w: Box<dyn std::io::Write + Send> = Box::new(std::io::BufWriter::with_capacity(1 << 20, w));
            Arc::new(Mutex::new(w))
        })
    };
    let (paf, sam, bam) = (share(out.paf), share(out.sam), share(out.bam));
    if let Some(o) = &sam {
        let mut h = Vec::new();
        output::sam_header(&mut h, mapper.index(), &out.command_line);
        o.lock().unwrap().write_all(&h)?;
    }
    if let Some(o) = &bam {
        let (mut h, mut z) = (Vec::new(), Vec::new());
        output::bam_header(&mut h, mapper.index(), &out.command_line);
        output::bgzf(&h, out.bam_level, &mut z);
        o.lock().unwrap().write_all(&z)?;
    }
    let total = Arc::new(Mutex::new((FileCounts::default(), InternalStats::default())));
    let mut proc = Processor {
        mapper: mapper.clone(),
        work: Workspace::default(),
        paf: paf.clone().map(|o| (o, Vec::new())),
        sam: sam.clone().map(|o| (o, Vec::new())),
        bam: bam.clone().map(|o| (o, Vec::new(), Vec::new())),
        bam_level: out.bam_level,
        unmapped: out.unmapped,
        qualities: out.qualities,
        total: total.clone(),
        progress,
        local: FileCounts::default(),
        batch_reads: 0,
    };
    let reader = crate::open_fastx(reads, Some(batch_size.max(1)))?;
    reader.process_parallel(&mut proc, threads.max(1))?;
    if let Some(o) = &bam {
        o.lock().unwrap().write_all(&output::BGZF_EOF)?;
    }
    for o in [&paf, &sam, &bam].into_iter().flatten() {
        o.lock().unwrap().flush()?;
    }
    let t = *total.lock().unwrap();
    Ok(t)
}

/// Alignments as `spearfish map-spliced` writes them: PAF, SAM and BAM.
pub mod output {
    use super::{GenomeIndex, SplicedAlignment, SplicedRead};

    fn base(a: &SplicedAlignment) -> Option<(&[u32], u32)> {
        a.edits.map(|e| (&a.cigar[..], e))
    }

    /// PAF lines of one read's alignments (none if unmapped). The `cg:Z` tag holds the CIGAR
    /// with introns as `N`: per block (`M`, and `I` or `D` for a length difference), or base
    /// by base when the alignment has one, then with `NM:i`. Other tags: `tp:A` P or S,
    /// `AS:i` the score, `cm:i` the seed matches chained, `ts:A` the transcript strand (`+`,
    /// `-` or `.`), and `jo:Z` one letter per junction for how it was found (for diagnosis).
    pub fn paf_lines(buf: &mut Vec<u8>, name: &str, r: &SplicedRead, idx: &GenomeIndex) {
        for a in &r.alignments {
            splice::map::paf_line(&idx.inner, name, r.len, &a.hit, a.mapq, base(a), buf);
        }
    }

    /// The SAM header for the index's sequences (unsorted, grouped by read). The `@PG` line
    /// carries this library's version and `command_line`.
    pub fn sam_header(buf: &mut Vec<u8>, idx: &GenomeIndex, command_line: &str) {
        splice::map::sam_header(&idx.inner, Some((crate::provenance::VERSION, command_line)), buf);
    }

    /// SAM records of one read's alignments (none if unmapped; see [`sam_unmapped`]). `seq`
    /// is the read as sequenced and `qual` its base qualities (Phred+33, as in FASTQ), if it
    /// has them; the primary record carries both, secondary records neither. Tags: `AS:i`,
    /// `cm:i`, `XS:A` (the transcript strand, when known) and `NM:i` with a base-level
    /// alignment.
    pub fn sam_records(
        buf: &mut Vec<u8>,
        name: &str,
        seq: &[u8],
        qual: Option<&[u8]>,
        r: &SplicedRead,
        idx: &GenomeIndex,
    ) {
        for a in &r.alignments {
            splice::map::sam_line(&idx.inner, name, seq, qual, &a.hit, a.mapq, base(a), buf);
        }
    }

    /// The SAM record of a read that did not map (flag 4), with its sequence and qualities.
    pub fn sam_unmapped(buf: &mut Vec<u8>, name: &str, seq: &[u8], qual: Option<&[u8]>) {
        splice::map::sam_unmapped(name, seq, qual, buf);
    }

    /// The uncompressed BAM header (pass it through [`bgzf`]).
    pub fn bam_header(buf: &mut Vec<u8>, idx: &GenomeIndex, command_line: &str) {
        splice::bam::header(&idx.inner, Some((crate::provenance::VERSION, command_line)), buf);
    }

    /// Uncompressed BAM records of one read's alignments, as [`sam_records`]. Collect the
    /// records of some reads and pass them through [`bgzf`].
    pub fn bam_records(
        buf: &mut Vec<u8>,
        name: &str,
        seq: &[u8],
        qual: Option<&[u8]>,
        r: &SplicedRead,
        idx: &GenomeIndex,
    ) {
        let mut ops = Vec::new();
        for a in &r.alignments {
            if a.edits.is_some() {
                ops.clone_from(&a.cigar);
            }
            splice::bam::record(
                &idx.inner, name, seq, qual, &a.hit, a.mapq, &mut ops, a.edits, buf,
            );
        }
    }

    /// The uncompressed BAM record of a read that did not map, as [`sam_unmapped`].
    pub fn bam_unmapped(buf: &mut Vec<u8>, name: &str, seq: &[u8], qual: Option<&[u8]>) {
        splice::bam::unmapped(name, seq, qual, buf);
    }

    /// Compress `data` into complete BGZF blocks appended to `out`, at DEFLATE level `level`
    /// (0 to 9). Blocks are independent, so threads can each compress their own records and
    /// a writer need only append the results. End the file with [`BGZF_EOF`].
    pub fn bgzf(data: &[u8], level: u8, out: &mut Vec<u8>) {
        splice::bam::bgzf(data, splice::bam::Level(level), out);
    }

    /// The empty block that ends a BGZF (BAM) file.
    pub const BGZF_EOF: [u8; 28] = splice::bam::EOF;

    /// The DEFLATE engine BAM output was built with.
    pub fn deflate_engine() -> &'static str {
        splice::bam::engine()
    }
}
