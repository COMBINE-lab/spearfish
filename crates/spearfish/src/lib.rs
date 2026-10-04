//! spearfish as a library: build or load an index, map long reads to the transcriptome,
//! and quantify.
//!
//! ```no_run
//! use std::sync::Arc;
//! use spearfish::{Index, MapOptions, Mapper, Preset};
//!
//! let index = Arc::new(Index::build_from_fasta("transcripts.fa".as_ref(), &Default::default())?);
//! let mapper = Mapper::new(index.clone(), &MapOptions::preset(Preset::Ont))?;
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

use std::path::Path;
use std::sync::Arc;

use anyhow::{Result, bail};
use index::{BuildParams, Collapse, TxInput};
use lcp_core::hpc::Hpc;
use map::{AnyChainer, ChainParams, MapParams, OracleChainer, RefineChainer, Strands, TieBreak};

pub use index::BuildReport;
pub use map::Ends;
pub use map::align::Affine;
pub use quant;

/// Read technology / protocol presets.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Preset {
    /// ONT cDNA: both orientations, minimap2 map-ont scores.
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

/// Which chainer to run (both give identical scores with an unbounded beam).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChainerKind {
    /// Per-transcript DP with prefix sharing (fastest).
    Oracle,
    /// Multi-transcript state DP with lazy color refinement; `beam` bounds states per anchor.
    Refine { beam: Option<usize> },
}

/// All mapping options. Start from [`MapOptions::preset`] and adjust fields.
#[derive(Clone, Debug)]
pub struct MapOptions {
    pub params: MapParams,
    pub chainer: ChainerKind,
}

impl MapOptions {
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
                // Direct RNA spans extended to the read ends change no spike-in estimate under
                // the default model and cost 31% more instructions (D27). Positional models
                // turn them back on.
                params.ends = None;
            }
            Preset::Hifi => {
                tb.affine = Some(Affine::HIFI);
                // Accurate reads: tied hits without a common anchor are grouped by shared
                // anchors rather than aligned whole, and level-2 rescue finds nothing level-3
                // seeds missed. Half the mapping work at identical accuracy (D26).
                tb.group_fallback = true;
                params.rescue = None;
                if let Some(e) = &mut params.ends {
                    e.sc = Affine::HIFI;
                }
            }
        }
        MapOptions {
            params,
            chainer: ChainerKind::Oracle,
        }
    }

    /// Mutable access to the tie-break options (None if disabled).
    pub fn tiebreak_mut(&mut self) -> Option<&mut TieBreak> {
        self.params.tiebreak.as_mut()
    }

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
#[derive(Clone, Debug)]
pub struct IndexOptions {
    /// Base seed scheme (its highest level is used), e.g. `lcp:K=3` or `minimizer:k=19,w=23`.
    pub scheme: String,
    /// Homopolymer-compress transcripts (and reads at mapping time).
    pub hpc: bool,
    /// Mask unitigs whose color touches more than this many genes.
    pub f_max: u32,
    /// Graph nodes per (label, gene) (true) or per label (false).
    pub per_gene: bool,
}

impl Default for IndexOptions {
    fn default() -> Self {
        IndexOptions {
            scheme: "lcp:K=3".into(),
            hpc: true,
            f_max: 20,
            per_gene: true,
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
        let (data, rep) = Self::build_data(txs, opts)?;
        Ok((
            Index {
                inner: Arc::new(index::Index::from_bytes(&data.to_bytes())?),
            },
            rep,
        ))
    }

    /// Build an index and write it to `path` (as `spearfish index`); returns bytes written.
    pub fn build_to_file(
        txs: Vec<(String, String, Vec<u8>)>,
        opts: &IndexOptions,
        path: &Path,
    ) -> Result<(u64, BuildReport)> {
        let (data, rep) = Self::build_data(txs, opts)?;
        Ok((data.write(path)?, rep))
    }

    fn build_data(
        txs: Vec<(String, String, Vec<u8>)>,
        opts: &IndexOptions,
    ) -> Result<(index::IndexData, BuildReport)> {
        let txs: Vec<TxInput> = txs
            .into_iter()
            .map(|(name, gene, seq)| TxInput { name, gene, seq })
            .collect();
        let params = BuildParams {
            spec: opts.scheme.clone(),
            hpc: opts.hpc,
            f_max: opts.f_max,
            allow_collisions: false,
            collapse: if opts.per_gene {
                Collapse::Gene
            } else {
                Collapse::Global
            },
        };
        index::build(&txs, &params)
    }

    /// Build an index in memory from a (optionally compressed) transcript FASTA. GENCODE
    /// headers (`tx|gene|...`) give gene IDs; otherwise each transcript is its own gene.
    #[cfg(feature = "fasta")]
    pub fn build_from_fasta(path: &Path, opts: &IndexOptions) -> Result<Self> {
        Ok(Self::build(read_transcripts(path)?, opts)?.0)
    }

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

/// Read transcripts from a FASTA as (name, gene, sequence).
#[cfg(feature = "fasta")]
pub fn read_transcripts(path: &Path) -> Result<Vec<(String, String, Vec<u8>)>> {
    use paraseq::{ReaderBuilder, Record};
    let mut txs = Vec::new();
    let mut reader = ReaderBuilder::path(path).build()?;
    let mut rset = reader.new_record_set();
    while rset.fill(&mut reader)? {
        for rec in rset.iter() {
            let rec = rec?;
            let id = rec.id_str();
            let mut f = id.split('|');
            let name = f
                .next()
                .unwrap()
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            let gene = f.next().map_or_else(|| name.clone(), String::from);
            txs.push((name, gene, rec.seq().into_owned()));
        }
    }
    Ok(txs)
}

/// One reported transcript for a read.
#[derive(Clone, Debug, PartialEq)]
pub struct Hit {
    /// Transcript id (index order, as in the FASTA the index was built from).
    pub tx: u32,
    /// The read's reverse complement matched the transcript.
    pub reverse: bool,
    /// Chain score after the tie-break; higher is better. Differences are on minimap2's
    /// alignment-score scale (λ = 1 per point); spearfish's EM uses exp(Δ / 2).
    pub score: f32,
    /// Query span on the read's forward strand, raw bases.
    pub query_start: u32,
    pub query_end: u32,
    /// Transcript span, raw bases (exact if the index holds raw transcript sequence, else
    /// scaled from compressed coordinates).
    pub target_start: u32,
    pub target_end: u32,
    /// Anchors in the chain, the read bases they cover (compressed), and the largest gap
    /// discrepancy |Δq − ρΔt| between consecutive anchors.
    pub anchors: u32,
    pub anchored: u32,
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
}

impl Mapper {
    pub fn new(index: Arc<Index>, opts: &MapOptions) -> Result<Self> {
        let chain = opts.params.chain.clone();
        let chainer = match opts.chainer {
            ChainerKind::Oracle => AnyChainer::Oracle(OracleChainer { p: chain }),
            ChainerKind::Refine { beam } => AnyChainer::Refine(RefineChainer { p: chain, beam }),
        };
        let inner = map::Mapper::new(index.inner.clone(), chainer, opts.params.clone())?;
        let hpc = inner.seeder.hpc;
        Ok(Mapper { index, inner, hpc })
    }

    pub fn index(&self) -> &Arc<Index> {
        &self.index
    }

    /// Map one read (ASCII bases; `U` is read as `T`).
    pub fn map(&self, seq: &[u8]) -> MappedRead {
        let r = self.inner.map(seq);
        let len = seq.len() as u32;
        let plen = r.plen;
        let hpc = self.hpc.then(|| Hpc::new(seq));
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
