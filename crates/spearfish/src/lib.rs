//! spearfish as a library: build or load an index, map long reads to the transcriptome,
//! and quantify.
//!
//! ```no_run
//! use std::sync::Arc;
//! use spearfish::{Index, MapSettings, Preset};
//!
//! let index = Arc::new(Index::build_from_fasta("transcripts.fa".as_ref(), &Default::default())?);
//! let mapper = MapSettings::builder(Preset::Ont).mapper(index.clone())?;
//! let read = mapper.map(b"ACGT...");
//! for hit in &read.hits {
//!     println!("{} {} {}", index.tx_name(hit.tx), hit.score, hit.reverse);
//! }
//! # Ok::<(), anyhow::Error>(())
//! ```
//!
//! [`Mapper`] is `Send + Sync`; share one across threads (scratch buffers are per thread).
//! Mapping works internally in homopolymer-compressed coordinates; [`Hit`] reports query
//! coordinates in raw read bases and transcript coordinates in raw transcript bases, both exact
//! (transcript coordinates fall back to scaling from compressed lengths for an index without
//! raw sequence). Spans run from where the read's alignment starts to where it ends, extended
//! past the outermost seeds (see [`Ends`]).
//!
//! Settings: [`MapSettings`] holds a [`Preset`] and the settings changed from it, and resolves
//! to the complete [`MapOptions`] (see [`settings`]). Index construction takes
//! [`IndexOptions`], also with a builder. Score differences between a read's hits convert to
//! probabilities as `exp(-d / SCORE_DENOM)` ([`SCORE_DENOM`]).

#![warn(missing_docs)]

pub mod provenance;
pub mod settings;

use std::path::Path;
use std::sync::Arc;

use anyhow::{Result, bail};
use index::{BuildParams, Collapse, TxInput};
use lcp_core::hpc::Hpc;
use map::{AnyChainer, OracleChainer, RefineChainer};
use serde::{Deserialize, Serialize};

pub use index::{BuildReport, InferGenes};
pub use map::align::Affine;
pub use map::{AnchorParams, ChainParams, Ends, MapParams, RescueParams, Strands, TieBreak};
pub use quant;
pub use settings::{MapSettings, RescueWindow};

/// Scale of score differences: a hit scoring `d` below its read's best is `exp(-d /
/// SCORE_DENOM)` times as likely. spearfish's scores are on minimap2's alignment-score scale;
/// D = 1 was the sharpest weighting tried with oarfish's EM on simulated reads and gave the
/// most accurate estimates (minimap2 alignments do best with a wider D, e.g. oarfish's 3).
/// RAD output records it, and oarfish's spearfish backend uses it.
pub const SCORE_DENOM: f32 = 1.0;

/// Read technology / protocol presets.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Preset {
    /// ONT cDNA: both orientations, minimap2 map-ont scores.
    #[default]
    Ont,
    /// ONT direct RNA: sense orientation only, cheaper tie-break.
    Drna,
    /// PacBio HiFi: both orientations, minimap2 map-hifi scores.
    Hifi,
}

impl std::str::FromStr for Preset {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self> {
        Ok(match s {
            "ont" => Preset::Ont,
            "drna" => Preset::Drna,
            "hifi" => Preset::Hifi,
            x => bail!("unknown preset {x:?} (ont, drna, hifi)"),
        })
    }
}

impl Preset {
    /// Every preset.
    pub const ALL: [Preset; 3] = [Preset::Ont, Preset::Drna, Preset::Hifi];

    /// The preset's name (as `--preset` takes it).
    pub fn name(self) -> &'static str {
        match self {
            Preset::Ont => "ont",
            Preset::Drna => "drna",
            Preset::Hifi => "hifi",
        }
    }

    /// The sequencing technology as oarfish's `--seq-tech` names it.
    pub fn seq_tech(self) -> &'static str {
        match self {
            Preset::Ont => "ont-cdna",
            Preset::Drna => "ont-drna",
            Preset::Hifi => "pac-bio-hifi",
        }
    }
}

/// Which chainer to run (both give identical scores with an unbounded beam).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ChainerKind {
    /// Per-transcript DP with prefix sharing (fastest).
    #[default]
    Oracle,
    /// Multi-transcript state DP with lazy color refinement; `beam` bounds states per anchor.
    Refine {
        /// States kept per anchor (`None`: unbounded, exact).
        beam: Option<usize>,
    },
}

/// The complete mapping options: every value a mapper runs with. [`MapSettings::resolve`]
/// derives them from a preset and the settings changed; [`MapOptions::preset`] gives a
/// preset's.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MapOptions {
    /// The preset these options started from.
    pub preset: Preset,
    /// Seeding, chaining, tie-break, rescue and reporting parameters.
    pub params: MapParams,
    /// The chainer.
    pub chainer: ChainerKind,
}

impl MapOptions {
    /// Turn on span extension (if it is off), scored as preset `p` scores. Prefer
    /// [`MapOptions::set_extend_ends`], which uses the options' own preset.
    pub fn extend_ends(&mut self, p: Preset) {
        if self.params.ends.is_none() {
            let sc = match p {
                Preset::Hifi => Affine::HIFI,
                Preset::Ont | Preset::Drna => Affine::ONT,
            };
            self.params.ends = Some(map::Ends {
                sc,
                ..Default::default()
            });
        }
    }

    /// Turn span extension on (scored as the preset scores; kept as is if already on) or off.
    pub fn set_extend_ends(&mut self, on: bool) {
        if on {
            self.extend_ends(self.preset);
        } else {
            self.params.ends = None;
        }
    }

    /// The scale of score differences ([`SCORE_DENOM`]).
    pub fn score_denom(&self) -> f32 {
        SCORE_DENOM
    }

    /// A preset's options.
    pub fn preset(p: Preset) -> Self {
        let mut params = MapParams::default();
        let tb = params.tiebreak.as_mut().expect("tie-break on by default");
        match p {
            Preset::Ont => {}
            Preset::Drna => {
                params.strands = Strands::Forward;
                tb.affine_within = Some(4);
                // No template-switch Gs at direct RNA read ends: raw ends cost 12% for nothing.
                tb.raw_ends = false;
                // Direct RNA reads end a few bases into their tail (basecalled poly(A) is short
                // and noisy): a 10-base run marks it, where cDNA's longer tails need 15 to stay
                // clear of A-runs near annotated 3' ends.
                tb.tail_run = 10;
                // Direct RNA spans extended to the read ends change no spike-in estimate under
                // the default model and cost 31% more instructions. Positional models
                // turn them back on.
                params.ends = None;
            }
            Preset::Hifi => {
                tb.affine = Some(Affine::HIFI);
                // Accurate reads: tied hits without a common anchor are grouped by shared
                // anchors rather than aligned whole, and level-2 rescue finds nothing level-3
                // seeds missed. Half the mapping work at identical accuracy.
                tb.group_fallback = true;
                tb.end_diagonal = true;
                params.rescue = None;
                // HiFi anchors already end within a few bases of the read's ends: extending
                // spans changes no spike-in estimate. Positional models turn it on.
                params.ends = None;
            }
        }
        MapOptions {
            preset: p,
            params,
            chainer: ChainerKind::Oracle,
        }
    }

    /// Mutable access to the tie-break options (None if disabled).
    pub fn tiebreak_mut(&mut self) -> Option<&mut TieBreak> {
        self.params.tiebreak.as_mut()
    }

    /// Mutable access to the chaining parameters.
    pub fn chain_mut(&mut self) -> &mut ChainParams {
        &mut self.params.chain
    }
}

impl Default for MapOptions {
    fn default() -> Self {
        Self::preset(Preset::Ont)
    }
}

/// Index construction options.
///
/// ```
/// let opts = spearfish::IndexOptions::builder().scheme("lcp:K=3").f_max(20).build();
/// assert!(opts.hpc);
/// ```
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, bon::Builder)]
pub struct IndexOptions {
    /// Base seed scheme (its highest level is used), e.g. `lcp:K=3` or `minimizer:k=19,w=23`.
    #[builder(into, default = "lcp:K=3".to_string())]
    pub scheme: String,
    /// Homopolymer-compress transcripts (and reads at mapping time).
    #[builder(default = true)]
    pub hpc: bool,
    /// Mask unitigs whose color touches more than this many genes.
    #[builder(default = 20)]
    pub f_max: u32,
    /// Graph nodes per (label, gene) (true) or per label (false).
    #[builder(default = true)]
    pub per_gene: bool,
    /// Group transcripts into genes from their sequence (for transcriptomes without gene
    /// IDs) instead of using the given genes.
    pub infer_genes: Option<InferGenesOptions>,
    /// Build even if two different seed contents share a label (normally an error).
    #[builder(default)]
    pub allow_collisions: bool,
}

impl Default for IndexOptions {
    fn default() -> Self {
        IndexOptions::builder().build()
    }
}

/// Gene inference from shared sequence (see [`InferGenes`]): transcripts sharing a run of
/// `span` consecutive seeds held by at most `max_tx` transcripts are grouped; then groups of
/// at most `attach_max` transcripts join the one larger group they share runs of
/// `attach_span` seeds with (0: no attach pass).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, bon::Builder)]
pub struct InferGenesOptions {
    /// Ignore runs held by more than this many transcripts (repeats).
    #[builder(default = 64)]
    pub max_tx: u32,
    /// Consecutive seeds a shared run must span.
    #[builder(default = 12)]
    pub span: usize,
    /// Attach pass: run length in seeds (0: off).
    #[builder(default = 4)]
    pub attach_span: usize,
    /// Attach pass: largest group that may attach.
    #[builder(default = 3)]
    pub attach_max: usize,
}

impl Default for InferGenesOptions {
    fn default() -> Self {
        InferGenesOptions::builder().build()
    }
}

impl IndexOptions {
    /// The index crate's build parameters.
    pub fn build_params(&self) -> BuildParams {
        BuildParams {
            spec: self.scheme.clone(),
            hpc: self.hpc,
            f_max: self.f_max,
            allow_collisions: self.allow_collisions,
            collapse: if self.per_gene {
                Collapse::Gene
            } else {
                Collapse::Global
            },
            infer_genes: self.infer_genes.map(|g| InferGenes {
                span: g.span,
                max_tx: g.max_tx,
                attach_span: g.attach_span,
                attach_max: g.attach_max,
            }),
        }
    }
}

/// A transcriptome index (memory-mapped when opened from a file).
pub struct Index {
    inner: Arc<index::Index>,
}

impl Index {
    /// Open an index written by `spearfish index` (memory-mapped).
    pub fn open(path: &Path) -> Result<Self> {
        Ok(Index {
            inner: Arc::new(index::Index::open(path)?),
        })
    }

    /// Build an index in memory from transcripts given as (name, gene, sequence).
    pub fn build(txs: Vec<(String, String, Vec<u8>)>, opts: &IndexOptions) -> Result<(Self, BuildReport)> {
        let (data, rep, _) = Self::build_data(txs, opts, Default::default(), Vec::new())?;
        Ok((
            Index {
                inner: Arc::new(index::Index::from_bytes(&data.to_bytes())?),
            },
            rep,
        ))
    }

    /// Build an index and write it to `path`, with its manifest next to it as
    /// `<path>.json`; returns bytes written.
    pub fn build_to_file(
        txs: Vec<(String, String, Vec<u8>)>,
        opts: &IndexOptions,
        path: &Path,
    ) -> Result<(u64, BuildReport)> {
        let (n, rep, _) = Self::build_to_file_with(txs, opts, path, Default::default(), Vec::new())?;
        Ok((n, rep))
    }

    /// [`Index::build_to_file`], recording where the input came from and the command line in
    /// the manifest, which is also returned.
    pub fn build_to_file_with(
        txs: Vec<(String, String, Vec<u8>)>,
        opts: &IndexOptions,
        path: &Path,
        input: provenance::InputInfo,
        command_line: Vec<String>,
    ) -> Result<(u64, BuildReport, provenance::IndexManifest)> {
        let (data, rep, manifest) = Self::build_data(txs, opts, input, command_line)?;
        let n = data.write(path)?;
        std::fs::write(
            provenance::sidecar(path),
            serde_json::to_string_pretty(&manifest)? + "\n",
        )?;
        Ok((n, rep, manifest))
    }

    fn build_data(
        txs: Vec<(String, String, Vec<u8>)>,
        opts: &IndexOptions,
        input: provenance::InputInfo,
        command_line: Vec<String>,
    ) -> Result<(index::IndexData, BuildReport, provenance::IndexManifest)> {
        let txs: Vec<TxInput> = txs
            .into_iter()
            .map(|(name, gene, seq)| TxInput { name, gene, seq })
            .collect();
        let (mut data, rep) = index::build(&txs, &opts.build_params())?;
        let mut manifest = provenance::IndexManifest::new(&data, opts, &rep, input);
        manifest.command_line = command_line;
        manifest.working_directory = std::env::current_dir()
            .map(|d| d.display().to_string())
            .unwrap_or_default();
        data.manifest = Some(serde_json::to_string(&manifest)?);
        Ok((data, rep, manifest))
    }

    /// Build an index in memory from a (optionally compressed) transcript FASTA. GENCODE
    /// headers (`tx|gene|...`) give gene IDs; otherwise each transcript is its own gene.
    #[cfg(feature = "fasta")]
    pub fn build_from_fasta(path: &Path, opts: &IndexOptions) -> Result<Self> {
        Ok(Self::build(read_transcripts(path)?, opts)?.0)
    }

    /// Number of transcripts.
    pub fn n_transcripts(&self) -> usize {
        self.inner.n_transcripts()
    }

    /// Transcript name (the FASTA ID up to the first `|` or whitespace).
    pub fn tx_name(&self, t: u32) -> &str {
        self.inner.tx_name(t)
    }

    /// Transcript length in raw bases.
    pub fn tx_len(&self, t: u32) -> u32 {
        self.inner.tx_len()[t as usize]
    }

    /// The underlying index (for advanced use).
    pub fn raw(&self) -> &index::Index {
        &self.inner
    }

    /// The index holds raw (uncompressed) transcript sequence.
    pub fn has_raw_sequence(&self) -> bool {
        self.inner.has_raw()
    }

    /// What built the index and from what (see [`provenance::IndexManifest`]); `None` for
    /// an index written before manifests were recorded.
    pub fn manifest(&self) -> Option<provenance::IndexManifest> {
        provenance::IndexManifest::of(&self.inner)
    }

    /// A short identifier of the index's contents: the hash of its input transcripts and
    /// build parameters (stable across rebuilds with the same inputs).
    pub fn id(&self) -> String {
        provenance::index_id(&self.inner)
    }
}

/// Name and gene from a FASTA ID: GENCODE-style `tx|gene|...`, else the ID (up to
/// whitespace) for both.
pub fn name_and_gene(id: &str) -> (String, String) {
    let mut f = id.split('|');
    let name = f
        .next()
        .unwrap_or("")
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string();
    let gene = f.next().map_or_else(|| name.clone(), String::from);
    (name, gene)
}

/// Read transcripts from a (optionally compressed) FASTA as (name, gene, sequence), names
/// and genes as [`name_and_gene`] splits the IDs.
#[cfg(feature = "fasta")]
pub fn read_transcripts(path: &Path) -> Result<Vec<(String, String, Vec<u8>)>> {
    use paraseq::{ReaderBuilder, Record};
    let mut txs = Vec::new();
    let mut reader = ReaderBuilder::path(path).build()?;
    let mut rset = reader.new_record_set();
    while rset.fill(&mut reader)? {
        for rec in rset.iter() {
            let rec = rec?;
            let (name, gene) = name_and_gene(rec.id_str());
            txs.push((name, gene, rec.seq().into_owned()));
        }
    }
    Ok(txs)
}

/// Replace genes from a transcript-to-gene table: whitespace-separated `transcript gene` per
/// line (further columns ignored). Returns how many transcripts the table lacks (they keep
/// their genes).
pub fn apply_t2g(txs: &mut [(String, String, Vec<u8>)], table: &str) -> usize {
    let map: std::collections::HashMap<&str, &str> = table
        .lines()
        .filter_map(|l| {
            let mut f = l.split_whitespace();
            Some((f.next()?, f.next()?))
        })
        .collect();
    let mut missing = 0;
    for (name, gene, _) in txs.iter_mut() {
        match map.get(name.as_str()) {
            Some(g) => *gene = g.to_string(),
            None => missing += 1,
        }
    }
    missing
}

/// One reported transcript for a read.
#[derive(Clone, Debug, PartialEq)]
pub struct Hit {
    /// Transcript id (index order, as in the FASTA the index was built from).
    pub tx: u32,
    /// The read's reverse complement matched the transcript.
    pub reverse: bool,
    /// Chain score after the tie-break; higher is better. Differences are on minimap2's
    /// alignment-score scale; see [`SCORE_DENOM`] for converting them to probabilities.
    pub score: f32,
    /// Query span on the read's forward strand, raw bases.
    pub query_start: u32,
    /// End of the query span (exclusive).
    pub query_end: u32,
    /// Transcript span, raw bases (exact if the index holds raw transcript sequence, else
    /// scaled from compressed coordinates).
    pub target_start: u32,
    /// End of the transcript span (exclusive).
    pub target_end: u32,
    /// Anchors in the chain.
    pub anchors: u32,
    /// Read bases the anchors cover (compressed).
    pub anchored: u32,
    /// Largest gap discrepancy |Δq − ρΔt| between consecutive anchors.
    pub max_discrepancy: f32,
    /// The span reaches the read's alignment ends; otherwise it runs from the first to the
    /// last anchor (hits scoring well below the best, or extension off).
    pub extended: bool,
}

/// All hits of one read, best first.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MappedRead {
    /// Read length (raw bases).
    pub len: u32,
    /// The reported hits, best first.
    pub hits: Vec<Hit>,
    /// The best chain spans less than half the read (spearfish quant skips such reads).
    pub low_coverage: bool,
    /// Hits came from the lower-level rescue.
    pub rescued: bool,
}

/// A mapper over one index; `Send + Sync`.
pub struct Mapper {
    index: Arc<Index>,
    inner: map::Mapper<AnyChainer>,
    hpc: bool,
    opts: MapOptions,
}

impl Mapper {
    /// A mapper over `index` with the complete options `opts` (see also
    /// [`MapSettings::mapper`]). Fails if the index can't serve the options.
    pub fn new(index: Arc<Index>, opts: &MapOptions) -> Result<Self> {
        if let Some(tb) = &opts.params.tiebreak
            && tb.raw_ends
            && index.raw().meta("hpc") == Some("true")
            && !index.has_raw_sequence()
        {
            log::warn!(
                "raw read ends need an index with raw transcript sequence; aligning read ends in compressed bases"
            );
        }
        let chain = opts.params.chain.clone();
        let chainer = match opts.chainer {
            ChainerKind::Oracle => AnyChainer::Oracle(OracleChainer { p: chain }),
            ChainerKind::Refine { beam } => AnyChainer::Refine(RefineChainer { p: chain, beam }),
        };
        let inner = map::Mapper::new(index.inner.clone(), chainer, opts.params.clone())?;
        let hpc = inner.seeder.hpc;
        Ok(Mapper {
            index,
            inner,
            hpc,
            opts: opts.clone(),
        })
    }

    /// The index mapped against.
    pub fn index(&self) -> &Arc<Index> {
        &self.index
    }

    /// The options the mapper runs with.
    pub fn options(&self) -> &MapOptions {
        &self.opts
    }

    /// The scale of score differences ([`SCORE_DENOM`]).
    pub fn score_denom(&self) -> f32 {
        SCORE_DENOM
    }

    /// Map one read (ASCII bases; `U` is read as `T`).
    pub fn map(&self, seq: &[u8]) -> MappedRead {
        let r = self.inner.map(seq);
        let len = seq.len() as u32;
        let plen = r.plen;
        // The read's runs as the mapper compressed it, with U read as T (a T next to a U is
        // one run).
        let hpc = self.hpc.then(|| {
            let dna: Vec<u8> = seq
                .iter()
                .map(|&b| if b == b'U' || b == b'u' { b'T' } else { b })
                .collect();
            Hpc::new(&dna)
        });
        // Parsed query span in the hit's orientation -> raw span on the forward strand.
        let q_raw = |s: u32, e: u32, rev: bool| -> (u32, u32) {
            let (fs, fe) = if rev { (plen - e, plen - s) } else { (s, e) };
            match &hpc {
                Some(h) => h.to_raw(fs, fe),
                None => (fs, fe),
            }
        };
        let idx = self.index.raw();
        let hits = r
            .hits
            .iter()
            .map(|h| {
                let (qs, qe) = q_raw(h.qstart, h.qend, h.reverse);
                let (tl, tpl) = (
                    idx.tx_len()[h.tx as usize] as u64,
                    idx.tx_plen()[h.tx as usize].max(1) as u64,
                );
                let scale = |x: u32| idx.tx_raw_pos(h.tx, x).unwrap_or(((x as u64 * tl) / tpl) as u32);
                Hit {
                    tx: h.tx,
                    reverse: h.reverse,
                    score: h.score,
                    query_start: qs,
                    query_end: qe,
                    target_start: scale(h.tstart),
                    target_end: scale(h.tend).max(scale(h.tstart)),
                    anchors: h.anchors,
                    anchored: h.anchored,
                    max_discrepancy: h.max_discrepancy,
                    extended: h.extended,
                }
            })
            .collect();
        MappedRead {
            len,
            hits,
            low_coverage: r.low_coverage,
            rescued: r.rescued,
        }
    }
}
