//! Read mapping: seeding both strands, anchor lookup, candidate-gene selection, chaining, and
//! result extraction.

pub mod align;
pub mod anchor;
pub mod chain;
pub mod refine;
pub mod rescue;
pub mod stats;

use std::sync::Arc;

use anyhow::{Context, Result};
use index::Index;
use lcp_core::encode::{codes, revcomp};
use seeds::BaseSeeder;

pub use anchor::{Anchor, AnchorParams, Candidate};
pub use chain::{AnyChainer, ChainParams, Chainer, OracleChainer, TxHit};
pub use refine::RefineChainer;
pub use rescue::{RescueParams, Rescuer};

/// Targeted base-level tie-break among near-best transcripts.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TieBreak {
    /// Transcripts within max(margin_abs, margin_rel × best) of the best chain score form
    /// the tie set.
    pub margin_abs: f32,
    pub margin_rel: f32,
    /// Score charged per extra edit relative to the tie set's best alignment.
    pub kappa: f32,
    /// Transcript bases added beyond the read's projected ends: `pad` + 1/`pad_div` of the
    /// unanchored read end.
    pub pad: u32,
    pub pad_div: u32,
    /// Align at most this many transcripts per read.
    pub max_cands: usize,
    /// Band: distances beyond the best candidate's + cap are not computed exactly.
    pub cap: u32,
    /// Align only where tied transcripts differ (false: whole read per candidate).
    pub segmented: bool,
    /// Segmented mode: score segments with affine gaps (None: unit-cost edit distance with
    /// `kappa` per edit). Scores are charged `lambda` per point below the best.
    pub affine: Option<align::Affine>,
    pub lambda: f32,
    /// Affine band (cells off the diagonals spanning the segment lengths).
    pub band: usize,
    /// Affine-score only segment texts within this many unit-cost edits of the best (None:
    /// all); farther texts are charged one mismatch per extra edit.
    pub affine_within: Option<u32>,
    /// Read ends: clip what doesn't align (adapters, junk; as minimap2's extension) instead
    /// of forcing every read base into the alignment.
    pub clip_ends: bool,
    /// With `clip_ends`: bonus for aligning through to the read's end.
    pub end_bonus: i32,
    /// When tied hits share no anchor, group them by shared anchors (segmented within groups,
    /// group representatives placed by clipped alignments) instead of aligning every hit's
    /// whole read.
    pub group_fallback: bool,
    /// Whole-read distances within a tie set cover only the read between the tied hits'
    /// outermost anchors (the read orientations are still compared whole).
    pub whole_trim: bool,
    /// With `clip_ends`: stop extending a read end once a row falls this far below the best
    /// score so far (X-drop; `i32::MAX`: never).
    pub xdrop: i32,
    /// With `clip_ends`: keep read-end bands on the diagonal from the shared anchor instead of
    /// widening them over the window's slack. Accurate reads barely drift from it, and long
    /// read ends then cost `2 * band + 1` cells per row instead of `slack` more.
    pub end_diagonal: bool,
    /// Near-tied read orientations: compare their best hits by clipped affine alignment
    /// scores (read ends past the transcript are clipped) instead of whole-read edit distance,
    /// which charges a short transcript for the read's adapters.
    pub orient_clipped: bool,
    /// HPC indexes with raw transcripts: align read ends in raw (uncompressed) space, where
    /// homopolymers at the read ends (e.g. a template-switch oligo's Gs) can't merge into the
    /// transcript's first or last run.
    pub raw_ends: bool,
    /// Right read end: a read whose raw right end shows a poly(A) tail (a `tail_run`-base window
    /// of A's with at most one other base) is aligned only up to the tail, so its A's can't
    /// align into an isoform's A-rich continuation, and each candidate whose transcript ends
    /// within `TAIL_TOL` (+ 1/10 of the read end) of the tail's start gains `tail_bonus` points:
    /// the read's 3' end is the molecule's. Cleavage sites and annotated ends vary by a few
    /// bases, so ends that close are not told apart (0: no bonus; the tail is still trimmed).
    pub tail_bonus: i32,
    /// Trim a detected poly(A) tail from the read's right end before aligning it.
    pub tail_trim: bool,
    /// Read ends: points a candidate gains past the shortest annotated end among the windows
    /// (by aligning read bases into its longer flank) count only from this many on, or from a
    /// short nearly exact extension; less is what primers and adapters score against flanks by
    /// chance (0: all count, the default: 10 raises top-1 on reads with primers but merges
    /// isoforms whose ends differ by a few bases).
    pub flank_min: i64,
    /// Shortest poly(A) tail recognized: a window this long with at most one base not A,
    /// followed by at most `TAIL_AFTER` read bases. Tails are searched for in the raw read,
    /// which without HPC is the parsed one.
    pub tail_run: usize,
    /// Near-tied read orientations: points for the orientation whose raw read ends in a poly(A)
    /// tail (the molecule's 3' end) over one that doesn't. Where a transcript's sequence lies
    /// antisense inside another's (antisense isoforms, a pseudogene in a host gene's intron),
    /// both orientations align equally well and the tail is the evidence of strand (0: off).
    pub orient_tail: i32,
    /// Both read orientations are tie-broken when the other's best chain scores within
    /// max(`margin_abs`, this × best) of the chosen one's. Wider than `margin_rel`: LCP parsing
    /// is not strand-symmetric, so the two orientations of a read draw different seeds, and
    /// where transcripts share sequence antisense the true orientation's chain can trail by
    /// seed luck alone. Only the orientation test uses it, not the tie sets within one.
    pub orient_margin_rel: f32,
    /// Report both read orientations with equal alignment evidence, whatever their chain
    /// scores, when their best hits are one sequence on opposite strands: the transcript
    /// bases under the read are reverse complements of each other. Chains and banded
    /// alignments of such a pair differ only by seed luck and heuristics, yet one read in a
    /// few would favor the antisense transcript by a decisive margin and abundance estimation
    /// would follow those reads. With equal scores it follows the transcripts' other reads
    /// (and `orient_tail`).
    pub orient_same: bool,
    /// Near-tied orientations whose clipped alignments leave a lead of at most this many
    /// points are compared again by exact local alignment: every diagonal within
    /// [`ORIENT_EXACT_BAND`] (and a twentieth of the read) of each chain's, so the chain only
    /// places the search and forces none of its anchors into the alignment (0: off). The
    /// clipped alignment follows the chain, and a chain with a stray anchor or a gap the band
    /// cannot cross is scored below what the read aligns at.
    #[serde(default)]
    pub orient_exact: f32,
    /// Compare the read between shared blocks in raw bases, cut at the blocks' midpoints
    /// (with `raw_ends` and an index that stores raw positions). Shared blocks are exact in
    /// parsed space only: candidates may differ in run lengths inside them, and a difference
    /// next to a run changes the parsed alignment around it.
    #[serde(default)]
    pub raw_mid: bool,
    /// A candidate with its own A-run where the read's poly(A) tail starts explains the tail
    /// as well as one that ends there: it gains `tail_bonus` too, and what the read's bases
    /// after the run gain against its continuation.
    #[serde(default)]
    pub tail_templated: bool,
}

impl TieBreak {
    /// Transcript bases searched for an unanchored read end of `x` bases: `x + pad + x / pad_div`.
    fn window(&self, x: u32) -> u32 {
        x + self.pad + x / self.pad_div
    }
}

impl Default for TieBreak {
    fn default() -> Self {
        TieBreak {
            margin_abs: 30.0,
            margin_rel: 0.05,
            kappa: 6.0,
            pad: 20,
            pad_div: 10,
            max_cands: 256,
            cap: 16,
            segmented: true,
            affine: Some(align::Affine::default()),
            lambda: 1.0,
            band: 16,
            affine_within: None,
            clip_ends: true,
            end_bonus: 5,
            xdrop: 50,
            end_diagonal: false,
            group_fallback: true,
            whole_trim: true,
            orient_clipped: true,
            raw_ends: true,
            tail_bonus: 10,
            tail_trim: true,
            tail_run: 15,
            orient_tail: 30,
            orient_margin_rel: 0.35,
            orient_same: true,
            orient_exact: 150.0,
            raw_mid: true,
            tail_templated: true,
            flank_min: 0,
        }
    }
}

#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct MapParams {
    pub anchor: AnchorParams,
    pub chain: ChainParams,
    /// Re-rank near-ties by edit distance to the read (None: chain scores only).
    pub tiebreak: Option<TieBreak>,
    /// Re-seed low-coverage reads at a lower LCP level within candidate genes (None: off).
    pub rescue: Option<RescueParams>,
    /// Report transcripts scoring within this fraction of the best.
    pub delta_rep: f32,
    /// Report transcripts scoring within this many points of the best as well: a fraction of
    /// a low best score is less than the chains of one read's equally good placements differ
    /// by.
    #[serde(default)]
    pub rep_abs: f32,
    /// Flag reads whose best chain spans less than this fraction of the read (cf. oarfish's
    /// minimum aligned fraction, 0.5).
    pub pi_cov: f32,
    /// Read orientations to map (direct RNA reads are sense-strand).
    pub strands: Strands,
    /// Merge collinear anchors within a unitig before chaining (see
    /// [`anchor::merge_collinear`]).
    pub merge_collinear: bool,
    /// Extend reported spans from the outermost anchors to where the read's alignment ends
    /// (None: spans end at the outermost anchors).
    pub ends: Option<Ends>,
    /// Reject reads whose best chain scores no better than chance (None: every read with a
    /// chain is reported).
    #[serde(default)]
    pub min_score: Option<MinScore>,
}

/// The lowest best-chain score a read is reported at. A read that matches nothing still
/// chains a chance seed or two, and the best chance chain grows with the logarithm of the
/// read's length (an extreme value over the places a chain can start), so the threshold does
/// too: `at_1kb + per_ln_len · ln(len / 1000)` for a read of `len` raw bases.
///
/// The defaults sit above what random sequence and reads foreign to the reference reach
/// against a human transcriptome under the default index scheme (about one random read in
/// 10^4 of up to a kilobase scores higher, a few in 10^3 at several kilobases), and below
/// nearly all reads with a true placement of 150 bases or more at ONT error rates. Chain
/// scores depend on the index scheme: recalibrate for another.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct MinScore {
    /// Threshold for a 1,000-base read.
    pub at_1kb: f32,
    /// Added per e-fold of read length.
    pub per_ln_len: f32,
    /// A chain of this many base-level anchors or more passes at any score (0: none does).
    /// Separate seeds that line up on one transcript are far rarer by chance than their
    /// summed weight suggests; the denser seeds of the hierarchical rescue are not.
    #[serde(default)]
    pub anchors: u32,
    /// Keep the hits of a read under the threshold; it is flagged `low_score` either way.
    pub keep: bool,
    /// A weak best hit from the rescue must also align (None: its chain score decides).
    #[serde(default)]
    pub verify: Option<Verify>,
}

/// Alignment check of a best hit that comes from the hierarchical rescue and either spans
/// less than `pi_cov` of the read or scores within `margin` of the threshold. Rescue seeds are
/// dense enough that chance matches chain past the threshold on low-complexity and foreign
/// reads; their alignments do not follow. The hit passes if its clipped affine alignment
/// score (see [`TieBreak::affine`]) reaches `base + per_base` × its parsed read span; without
/// an affine tie-break nothing is checked.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Verify {
    pub base: f32,
    pub per_base: f32,
    pub margin: f32,
}

impl Default for Verify {
    fn default() -> Self {
        Verify {
            base: 40.0,
            per_base: 1.0,
            margin: 30.0,
        }
    }
}

impl MinScore {
    /// The threshold for a read of `len` raw bases.
    pub fn at(&self, len: usize) -> f32 {
        self.at_1kb + self.per_ln_len * (len.max(1) as f32 / 1000.0).ln()
    }

    /// Whether `h`, the best chain of a read of `len` raw bases, is more than a chance match
    /// (`rescued`: it comes from the hierarchical rescue).
    pub fn passes(&self, h: &TxHit, len: usize, rescued: bool) -> bool {
        h.score >= self.at(len) || (!rescued && self.anchors > 0 && h.anchors >= self.anchors)
    }
}

impl Default for MinScore {
    fn default() -> Self {
        MinScore {
            at_1kb: 36.0,
            per_ln_len: 4.0,
            anchors: 0,
            keep: false,
            verify: Some(Verify::default()),
        }
    }
}

/// X-drop extension of each reported hit from its first and last anchors toward the read's
/// ends, as minimap2 extends its chains. The anchors stop short of the read's ends by up to a
/// seed spacing; where the read really starts and ends on the transcript is what positional
/// models (truncation, 3' bias) need, and adapters or poly(A) tails are clipped by the X-drop.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Ends {
    pub sc: align::Affine,
    pub band: usize,
    pub end_bonus: i32,
    pub xdrop: i32,
    /// Transcript flank searched beyond the read flank's length: `pad + flank / pad_div`.
    pub pad: u32,
    pub pad_div: u32,
    /// Extend only hits scoring within this fraction of the read's best (oarfish keeps
    /// alignments within 5% of the best by default); the others keep their anchor spans.
    pub within: f32,
}

/// End extensions of one read, keyed by (left side, reverse, read position): the transcript
/// flank prefix each extension read, whether that was the whole flank, and the read and
/// transcript bases it covers.
type EndCache = rustc_hash::FxHashMap<(bool, bool, u32), Vec<(Vec<u8>, bool, (u32, u32))>>;

impl Default for Ends {
    fn default() -> Self {
        Ends {
            sc: align::Affine::ONT,
            band: 16,
            end_bonus: 5,
            xdrop: 50,
            pad: 20,
            pad_div: 10,
            within: 0.05,
        }
    }
}

/// Which orientations of a read to map against the transcripts.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Strands {
    #[serde(rename = "both")]
    Both,
    /// The read as given (sense: direct RNA).
    #[serde(rename = "fwd")]
    Forward,
    /// Only its reverse complement.
    #[serde(rename = "rev")]
    Reverse,
}

impl Strands {
    /// The name used on the command line and in output headers.
    pub fn name(self) -> &'static str {
        match self {
            Strands::Both => "both",
            Strands::Forward => "fwd",
            Strands::Reverse => "rev",
        }
    }
}

impl Default for MapParams {
    fn default() -> Self {
        MapParams {
            anchor: AnchorParams::default(),
            chain: ChainParams::default(),
            tiebreak: Some(TieBreak::default()),
            rescue: Some(RescueParams::default()),
            delta_rep: 0.3,
            rep_abs: 30.0,
            pi_cov: 0.5,
            strands: Strands::Both,
            merge_collinear: true,
            ends: Some(Ends::default()),
            min_score: Some(MinScore::default()),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadResult {
    /// Read length (raw) and parsed length.
    pub len: u32,
    pub plen: u32,
    /// True if the reverse complement mapped best.
    pub reverse: bool,
    /// Reported transcripts, best first. Near-tied orientations are both reported; each
    /// hit records its own.
    pub hits: Vec<TxHit>,
    /// The best chain spans less than `pi_cov` of the (parsed) read.
    pub low_coverage: bool,
    /// The hits come from the hierarchical rescue.
    pub rescued: bool,
    /// The best chain fails [`MapParams::min_score`] ([`MinScore::passes`]): it is no better
    /// than a chance match. The read has no hits unless [`MinScore::keep`] is set.
    pub low_score: bool,
    /// With HPC, the read's homopolymer runs, forward orientation: parsed query positions to
    /// raw ones (`to_raw`; its `seq` is not kept).
    pub hpc: Option<lcp_core::hpc::Hpc>,
}

pub struct Mapper<C: Chainer> {
    pub idx: Arc<Index>,
    pub seeder: BaseSeeder,
    pub chainer: C,
    pub p: MapParams,
    rescuer: Option<Rescuer>,
    /// [`anchor::node_families`] of `idx`.
    node_fam: Vec<u32>,
    /// Per transcript, its [`AntiSeed`]s by position (made on first use).
    anti: Vec<std::sync::OnceLock<Box<[AntiSeed]>>>,
}

/// A base-level seed of a transcript's reverse complement that occurs in another transcript:
/// the two share that sequence on opposite strands.
#[derive(Clone, Copy)]
struct AntiSeed {
    /// The seed's first base in the transcript (forward coordinates).
    x: u32,
    len: u32,
    /// The index nodes with the seed's label.
    nodes: (u32, u32),
}

/// Transcript -> (weight, lowest and highest diagonal), while a hit's antisense mates are
/// gathered.
type Mates = rustc_hash::FxHashMap<u32, (u32, i64, i64)>;

thread_local! {
    static MATES: std::cell::RefCell<Mates> = std::cell::RefCell::new(Mates::default());
}

/// A read orientation's raw bases (2-bit codes) and HPC run map, for raw-space read ends.
struct RawRead {
    codes: Vec<u8>,
    hpc: lcp_core::hpc::Hpc,
}

impl RawRead {
    /// Raw position where parsed base `p` starts (the raw length at the parsed length).
    fn pos(&self, p: u32) -> usize {
        self.hpc.to_raw(p, p).0 as usize
    }
}

/// A shared exact block, in half-open parsed query/target coordinates.
#[derive(Debug, PartialEq)]
struct SharedBlock {
    query: std::ops::Range<u32>,
    targets: Vec<std::ops::Range<u32>>,
}

/// Match complete query intervals; TxHit does not retain seed/unitig identity.
/// Overlapping intervals can cover one exact block only if their query-to-target
/// translations agree in every chain; an overlapping anchor whose placements disagree, or
/// one whose targets run backwards in some chain, is dropped. `None` if no (query start,
/// length) is in every chain. Adjacent intervals remain separate: a target insertion may lie
/// between them.
fn shared_blocks(hits: &[TxHit]) -> Option<Vec<SharedBlock>> {
    let first = hits.first()?;
    // Chains are already ordered by query start. Keep one monotone cursor per
    // candidate, rather than rebuilding placement hash maps for every tie set.
    let mut cursors = vec![0usize; hits.len()];
    let mut targets = Vec::with_capacity(hits.len());
    let mut lengths = Vec::new();
    let mut blocks: Vec<SharedBlock> = Vec::new();
    let mut group = 0;
    while group < first.chain.len() {
        let q = first.chain[group].0;
        let mut next = group + 1;
        while next < first.chain.len() && first.chain[next].0 == q {
            next += 1;
        }
        lengths.clear();
        lengths.extend(first.chain[group..next].iter().map(|a| a.2));
        // Only equal-start groups need ordering; normal chains have one anchor
        // per query start. Distinct lengths still denote distinct anchors.
        lengths.sort_unstable();
        lengths.dedup();
        group = next;
        for (h, cursor) in hits.iter().zip(&mut cursors) {
            while *cursor < h.chain.len() && h.chain[*cursor].0 < q {
                *cursor += 1;
            }
        }
        for &len in &lengths {
            targets.clear();
            for (h, &cursor) in hits.iter().zip(&cursors) {
                let placement = h.chain[cursor..]
                    .iter()
                    .take_while(|a| a.0 == q)
                    .find(|a| a.2 == len);
                let Some(&(_, t, _)) = placement else { break };
                targets.push(t..t.checked_add(len)?);
            }
            if targets.len() != hits.len() {
                continue;
            }
            let end = q.checked_add(len)?;
            if let Some(b) = blocks.last_mut() {
                if q < b.query.end {
                    // An overlapping anchor on another diagonal in some candidate: drop it and
                    // keep the block (the gap after the block is aligned instead).
                    let d = q - b.query.start;
                    if b.targets
                        .iter()
                        .zip(&targets)
                        .any(|(old, new)| old.start + d != new.start)
                    {
                        continue;
                    }
                    for (old, new) in b.targets.iter_mut().zip(&targets) {
                        old.end = old.end.max(new.end);
                    }
                    b.query.end = b.query.end.max(end);
                    continue;
                }
                // Target ranges that run backwards in some candidate cannot be aligned
                // globally: drop the anchor.
                if b.targets.iter().zip(&targets).any(|(a, b)| a.end > b.start) {
                    continue;
                }
            }
            blocks.push(SharedBlock {
                query: q..end,
                targets: targets.clone(),
            });
        }
    }
    (!blocks.is_empty()).then_some(blocks)
}

/// One orientation of a read: hits (sorted), the parsed sequence and its base seeds.
struct Strand {
    hits: Vec<TxHit>,
    seq: Vec<u8>,
    /// `seq` as 2-bit codes, made on first use ([`strand_codes`]).
    codes: std::cell::OnceCell<Vec<u8>>,
    raw: Option<RawRead>,
    seeds: Vec<seeds::Seed>,
}

/// A strand's `codes`, made from its `seq` on first use (a function of the two fields rather
/// than a method, so the strand's hits can be borrowed mutably meanwhile).
fn strand_codes<'a>(codes: &'a std::cell::OnceCell<Vec<u8>>, seq: &[u8]) -> &'a [u8] {
    codes.get_or_init(|| lcp_core::encode::codes(seq))
}

impl<C: Chainer> Mapper<C> {
    /// A mapper seeding exactly as the index was built.
    pub fn new(idx: Arc<Index>, chainer: C, p: MapParams) -> Result<Self> {
        p.chain.check()?;
        idx.check_compatible()?;
        let spec = idx.meta("spec").context("index lacks its seed spec")?;
        let hpc = idx.meta("hpc").context("index lacks its HPC flag")? == "true";
        let seeder = BaseSeeder::new(spec, hpc)?;
        let mut rescuer = p.rescue.clone().map(Rescuer::new).transpose()?;
        if let Some(r) = &mut rescuer {
            r.build_tables(&idx);
        }
        let node_fam = anchor::node_families(&idx);
        Ok(Mapper {
            node_fam,
            anti: (0..idx.n_transcripts()).map(|_| Default::default()).collect(),
            idx,
            seeder,
            chainer,
            p,
            rescuer,
        })
    }

    /// Re-rank the near-tie set of `hits` (sorted, best first) by edit distance between the
    /// parsed read (`codes` gives its 2-bit codes, asked for only if there is a tie) and each
    /// transcript. Returns the tie-set size. `salt` orders equal scores (see [`sort_hits`]).
    fn tiebreak<'c>(
        &self,
        tb: &TieBreak,
        codes: impl FnOnce() -> &'c [u8],
        raw: Option<&RawRead>,
        hits: &mut [TxHit],
        salt: u64,
    ) -> u32 {
        let Some(best) = hits.first().map(|h| h.score) else {
            return 0;
        };
        let margin = tb.margin_abs.max(tb.margin_rel * best.abs());
        let n = hits
            .iter()
            .take(tb.max_cands)
            .take_while(|h| h.score >= best - margin)
            .count();
        if n < 2 {
            return 0;
        }
        stats::add(&stats::TIEBREAK_CALLS, 1);
        stats::add(&stats::TIEBREAK_TIED, n as u64);
        let codes = codes();
        // Costs (edits, or negated affine scores) and the score charged per unit.
        let seg = if tb.segmented {
            self.segmented(tb, codes, raw, &hits[..n])
        } else {
            None
        };
        // Alignment deficit relative to the tie set's best alignment, per hit.
        let penalty: Vec<f32> = match seg {
            Some((cost, unit)) => {
                let cmin = *cost.iter().min().unwrap();
                cost.iter().map(|&c| unit * (c - cmin) as f32).collect()
            }
            None => {
                stats::add(&stats::TIEBREAK_WHOLE_READ, 1);
                if tb.segmented && tb.group_fallback {
                    self.grouped(tb, codes, raw, &hits[..n])
                } else {
                    let d = self.whole_read(tb, codes, &hits[..n]);
                    let dmin = *d.iter().min().unwrap();
                    d.iter().map(|&x| tb.kappa * (x - dmin) as f32).collect()
                }
            }
        };
        for (h, &p) in hits[..n].iter_mut().zip(&penalty) {
            h.score = best - p;
        }
        sort_hits(hits, salt);
        n as u32
    }

    /// Penalties for a tie set whose hits share no common anchor (typically paralogs): the
    /// hits are grouped greedily so each group shares an anchor, scored by [`Self::segmented`]
    /// within groups, and the groups are placed relative to each other by the whole-read edit
    /// distance of each group's best hit (`kappa` per edit), as the two read orientations are.
    /// Aligns one whole-read region per group instead of one per hit.
    fn grouped(&self, tb: &TieBreak, codes: &[u8], raw: Option<&RawRead>, hits: &[TxHit]) -> Vec<f32> {
        type AnchorGroup = (rustc_hash::FxHashSet<(u32, u32)>, Vec<usize>);
        let mut groups: Vec<AnchorGroup> = Vec::new();
        for (i, h) in hits.iter().enumerate() {
            let qs: rustc_hash::FxHashSet<(u32, u32)> = h.chain.iter().map(|a| (a.0, a.2)).collect();
            match groups.iter_mut().find(|g| g.0.iter().any(|q| qs.contains(q))) {
                Some(g) => {
                    g.0.retain(|q| qs.contains(q));
                    g.1.push(i);
                }
                None => groups.push((qs, vec![i])),
            }
        }
        let mut penalty = vec![0f32; hits.len()];
        let mut reps = Vec::with_capacity(groups.len());
        for (_, members) in &groups {
            let mut rep = members[0];
            if members.len() > 1 {
                let sub: Vec<TxHit> = members.iter().map(|&i| hits[i].clone()).collect();
                // Members share an anchor (the group's set is the intersection of theirs), so
                // `shared_blocks` always finds a block.
                let (cost, unit) = self
                    .segmented(tb, codes, raw, &sub)
                    .expect("group members share an anchor");
                let cmin = *cost.iter().min().unwrap();
                for (k, &i) in members.iter().enumerate() {
                    penalty[i] = unit * (cost[k] - cmin) as f32;
                }
                rep = members[cost.iter().position(|&c| c == cmin).unwrap()];
            }
            reps.push(rep);
        }
        if groups.len() > 1 {
            // Costs of each group's representative, and the score charged per unit: clipped
            // affine alignments where available (as for read orientations), else whole-read
            // edit distances.
            let (cost, unit): (Vec<i64>, f32) = match &tb.affine {
                Some(sc) if tb.orient_clipped => (
                    reps.iter()
                        .map(|&i| -i64::from(self.clipped(tb, sc, codes, raw, &hits[i])))
                        .collect(),
                    tb.lambda,
                ),
                _ => {
                    let rep_hits: Vec<TxHit> = reps.iter().map(|&i| hits[i].clone()).collect();
                    let d = self.whole_read(tb, codes, &rep_hits);
                    (d.into_iter().map(i64::from).collect(), tb.kappa)
                }
            };
            let cmin = *cost.iter().min().unwrap();
            for ((_, members), &cg) in groups.iter().zip(&cost) {
                let off = unit * (cg - cmin) as f32;
                for &i in members {
                    penalty[i] += off;
                }
            }
        }
        penalty
    }

    /// Distances of `read` against several texts with a shared band: the first is exact,
    /// later ones are capped at (best so far + `cap`) + 1, as one [`align::distance_with`]
    /// call per text in order would give them. `distinct`: no text repeats, so the texts are
    /// not hashed to find repeats.
    ///
    /// Identical texts (isoforms agreeing over the region) are aligned once: the band only
    /// tightens (k = best so far + cap), so a repeat's result is the first occurrence's
    /// re-capped at the current k. Distinct texts are aligned in batches of eight
    /// ([`align::batch_distances`]) with the threshold in force at each batch's first text:
    /// later texts' thresholds are no larger (the best only falls), so a batch result is exact
    /// wherever the sequential pass's would be, and re-capping it at each text's own threshold,
    /// in order, gives the sequential results exactly.
    fn banded<T: AsRef<[u8]>>(
        read: &[u8],
        texts: &[T],
        mode: align::Mode,
        cap: u32,
        distinct: bool,
    ) -> Vec<u32> {
        let (pat, m) = (align::Pattern::new(read), read.len());
        // Each text's distinct-text index, and the distinct texts in first-occurrence order.
        let mut uniq: Vec<&[u8]> = Vec::with_capacity(texts.len());
        let mut at: Vec<usize> = Vec::new();
        let id: Vec<u32> = if distinct {
            uniq.extend(texts.iter().map(|t| t.as_ref()));
            at.extend(0..texts.len());
            (0..texts.len() as u32).collect()
        } else {
            let mut first: rustc_hash::FxHashMap<&[u8], u32> = rustc_hash::FxHashMap::default();
            texts
                .iter()
                .enumerate()
                .map(|(i, t)| {
                    *first.entry(t.as_ref()).or_insert_with(|| {
                        uniq.push(t.as_ref());
                        at.push(i);
                        uniq.len() as u32 - 1
                    })
                })
                .collect()
        };
        // Batch results (exact when within the batch's threshold), then the sequential values.
        let mut raw: Vec<Option<u32>> = Vec::with_capacity(uniq.len());
        let mut val = vec![0u32; uniq.len()];
        let mut best = u32::MAX;
        let mut out = Vec::with_capacity(texts.len());
        for (i, &u) in id.iter().enumerate() {
            let u = u as usize;
            let k = best.saturating_add(cap);
            if at[u] != i {
                // A repeat: the first occurrence's value re-capped at the current threshold.
                let d = val[u];
                out.push(if d <= k { d } else { k.saturating_add(1) });
                continue;
            }
            if u == raw.len() {
                let batch = &uniq[u..(u + 8).min(uniq.len())];
                for t in batch {
                    stats::add(&stats::MYERS_CALLS, 1);
                    stats::add(&stats::MYERS_CELLS, (m as u64).div_ceil(64) * t.len() as u64);
                }
                raw.extend(align::batch_distances(&pat, batch, mode, k));
            }
            let d = match raw[u] {
                Some(d) if d <= k => d,
                _ => k.saturating_add(1),
            };
            best = best.min(d);
            val[u] = d;
            out.push(d);
        }
        out
    }

    /// Whole read against each candidate's covered transcript region (projected to the read
    /// ends, padded), infix mode.
    fn whole_read(&self, tb: &TieBreak, codes: &[u8], hits: &[TxHit]) -> Vec<u32> {
        // Only the read between the outermost anchors of any hit: primers and adapters lie
        // outside every hit's chain, and aligning them whole charges each hit for them unless
        // its transcript has a flank there to absorb them.
        let (lo, hi) = if tb.whole_trim {
            (
                hits.iter().map(|h| h.qstart).min().unwrap_or(0),
                hits.iter().map(|h| h.qend).max().unwrap_or(codes.len() as u32),
            )
        } else {
            (0, codes.len() as u32)
        };
        let texts: Vec<Vec<u8>> = hits.iter().map(|h| self.hit_region(tb, lo, hi, h)).collect();
        Self::banded(
            &codes[lo as usize..hi as usize],
            &texts,
            align::Mode::Infix,
            tb.cap,
            false,
        )
    }

    /// Whether the best hits of the two read orientations, `a` and `b`, are one sequence on
    /// opposite strands (antisense transcripts sharing sequence): each transcript's bases
    /// under the read span of both chains occur reverse-complemented in the other's padded
    /// region. Compared in raw bases where the index has them.
    fn same_antisense(&self, tb: &TieBreak, plen: u32, a: &TxHit, b: &TxHit) -> bool {
        // A hit's transcript interval projected to the read span of both chains, with `pad`
        // more bases on each side.
        let region = |h: &TxHit, o: &TxHit, pad: u32| {
            let tlen = self.idx.tx_plen()[h.tx as usize];
            let lo = h.qstart.min(plen.saturating_sub(o.qend));
            let hi = h.qend.max(plen.saturating_sub(o.qstart));
            (
                h.tstart.saturating_sub(h.qstart - lo + pad),
                (h.tend + (hi - h.qend) + pad).min(tlen),
            )
        };
        let raw = self.idx.has_raw();
        let codes = |h: &TxHit, (lo, hi): (u32, u32)| match raw {
            true => self.idx.tx_raw_codes(h.tx, lo, hi),
            false => Some(self.idx.tx_codes(h.tx, lo, hi)),
        };
        let within = |h: &TxHit, o: &TxHit| {
            let (lo, hi) = region(h, o, 0);
            // Indels in the read move the other chain's projection.
            let (wlo, whi) = region(o, h, tb.window(hi - lo) - (hi - lo));
            if hi - lo > whi - wlo {
                return false;
            }
            let (Some(mut needle), Some(text)) = (codes(h, (lo, hi)), codes(o, (wlo, whi))) else {
                return false;
            };
            needle.reverse();
            needle.iter_mut().for_each(|c| *c = 3 - *c);
            !needle.is_empty() && text.windows(needle.len()).any(|w| w == &needle[..])
        };
        within(a, b) && within(b, a)
    }

    /// A hit's transcript region projected to the ends of the read span `[lo, hi)`, padded.
    fn hit_region(&self, tb: &TieBreak, lo: u32, hi: u32, h: &TxHit) -> Vec<u8> {
        let tlen = self.idx.tx_plen()[h.tx as usize];
        let (left, right) = (h.qstart - lo, hi - h.qend);
        let lo = h.tstart.saturating_sub(tb.window(left));
        let hi = (h.tend + tb.window(right)).min(tlen);
        self.idx.tx_codes(h.tx, lo, hi)
    }

    /// Affine score of the read against one hit with both read ends clipped: the chain is cut
    /// at anchors at least two bands apart and each piece aligned globally, the best run of
    /// consecutive pieces is kept, and each read end is extended from that run into the padded
    /// transcript flank. The band so follows the chain (one band around the whole span's
    /// diagonals loses a long read whose indels drift further than the band from them), and
    /// a stray anchor far out in an unalignable read end is clipped with that end instead of
    /// being aligned through. With raw read ends the pieces are aligned in raw space too: hits
    /// compared by this score (read orientations, anchor groups) have different outermost
    /// anchors, and a score that counts part of the read in parsed bases and part in raw ones
    /// would favor whichever hit leaves more of the read to its raw ends.
    fn clipped(
        &self,
        tb: &TieBreak,
        sc: &align::Affine,
        codes: &[u8],
        raw: Option<&RawRead>,
        h: &TxHit,
    ) -> i32 {
        let plen = codes.len() as u32;
        let tlen = self.idx.tx_plen()[h.tx as usize];
        let w0 = align::dp_work();
        let mut cuts = vec![(h.qstart, h.tstart)];
        for &(q, t, _) in &h.chain {
            let &(pq, pt) = cuts.last().unwrap();
            if q >= pq + 2 * tb.band as u32 && t > pt && q < h.qend && t < h.tend {
                cuts.push((q, t));
            }
        }
        cuts.push((h.qend, h.tend));
        // Best run of consecutive pieces: (score, first cut, last cut).
        let (mut mid, mut from, mut to) = (i32::MIN, 0, 0);
        let (mut run, mut start) = (0i32, 0);
        for (i, w) in cuts.windows(2).enumerate() {
            let ((q0, t0), (q1, t1)) = (w[0], w[1]);
            let x = match self.raw_texts(tb, raw, &[(h.tx, t0, t1)]) {
                Some((r, mut t)) => align::affine(
                    &r.codes[r.pos(q0)..r.pos(q1)],
                    &t.pop().unwrap(),
                    align::Mode::Global,
                    sc,
                    tb.band,
                ),
                None => align::affine(
                    &codes[q0 as usize..q1 as usize],
                    &self.idx.tx_codes(h.tx, t0, t1),
                    align::Mode::Global,
                    sc,
                    tb.band,
                ),
            };
            if run <= 0 {
                (run, start) = (0, i);
            }
            run = run.saturating_add(x);
            if run > mid {
                (mid, from, to) = (run, start, i + 1);
            }
        }
        let ((qs, ts), (qe, te)) = (cuts[from], cuts[to]);
        let w1 = align::dp_work();
        let left = if qs > 0 {
            let span = [(h.tx, ts.saturating_sub(tb.window(qs)), ts)];
            let (read, text) = match self.raw_texts(tb, raw, &span) {
                Some((r, mut t)) => (&r.codes[..r.pos(qs)], t.pop().unwrap()),
                None => (&codes[..qs as usize], self.idx.tx_codes(h.tx, span[0].1, ts)),
            };
            let read: Vec<u8> = read.iter().rev().copied().collect();
            let text: Vec<u8> = text.into_iter().rev().collect();
            align::affine_extend(&read, &text, sc, tb.band, tb.end_bonus, tb.xdrop)
        } else {
            tb.end_bonus
        };
        let rest = plen - qe;
        let right = if rest > 0 {
            let span = [(h.tx, te, (te + tb.window(rest)).min(tlen))];
            let (read, text) = match self.raw_texts(tb, raw, &span) {
                Some((r, mut t)) => (&r.codes[r.pos(qe)..], t.pop().unwrap()),
                None => (&codes[qe as usize..], self.idx.tx_codes(h.tx, te, span[0].2)),
            };
            align::affine_extend(read, &text, sc, tb.band, tb.end_bonus, tb.xdrop)
        } else {
            tb.end_bonus
        };
        if stats::on() {
            let w2 = align::dp_work();
            stats::add(&stats::DP_CLIP_MID_CELLS, w1.0 - w0.0);
            stats::add(&stats::DP_CLIP_MID_ROWS, w1.1 - w0.1);
            stats::add(&stats::DP_CLIP_END_CELLS, w2.0 - w1.0);
            stats::add(&stats::DP_CLIP_END_ROWS, w2.1 - w1.1);
        }
        mid + left + right
    }

    /// Local affine score of the whole read (ends free to clip, a bonus per read end reached)
    /// against one hit's transcript, over every diagonal within `w` of the chain's: the chain
    /// only places the search, it forces nothing into the alignment.
    fn local_exact(
        &self,
        tb: &TieBreak,
        sc: &align::Affine,
        codes: &[u8],
        raw: Option<&RawRead>,
        h: &TxHit,
        w: u32,
    ) -> i32 {
        let plen = codes.len() as u32;
        let tlen = self.idx.tx_plen()[h.tx as usize];
        let wlo = h.tstart.saturating_sub(tb.window(h.qstart) + w);
        let whi = (h.tend + tb.window(plen - h.qend) + w).min(tlen);
        let ends = [(h.qstart, h.tstart), (h.qend, h.tend)];
        let pts = ends.into_iter().chain(
            h.chain
                .iter()
                .map(|&(q, t, _)| (q, t))
                .filter(|&(q, t)| q <= plen && t >= wlo && t <= whi),
        );
        let rawr = raw.filter(|_| tb.raw_ends && self.idx.tx_raw_pos(h.tx, 0).is_some());
        let (read, text, d): (&[u8], Vec<u8>, Vec<i64>) = match rawr {
            Some(r) => {
                let base = self.idx.tx_raw_pos(h.tx, wlo).unwrap();
                (
                    &r.codes[..],
                    self.idx.tx_raw_codes(h.tx, wlo, whi).unwrap(),
                    pts.map(|(q, t)| {
                        i64::from(self.idx.tx_raw_pos(h.tx, t).unwrap() - base) - r.pos(q) as i64
                    })
                    .collect(),
                )
            }
            None => (
                codes,
                self.idx.tx_codes(h.tx, wlo, whi),
                pts.map(|(q, t)| i64::from(t - wlo) - i64::from(q)).collect(),
            ),
        };
        let (dlo, dhi) = (
            d.iter().min().unwrap() - i64::from(w),
            d.iter().max().unwrap() + i64::from(w),
        );
        align::local_band(read, &text, sc, dlo, dhi, tb.end_bonus)
    }

    /// Whole-read infix distance against one hit's region, capped: `k + 1` if above `k`.
    fn whole_read_k(&self, tb: &TieBreak, codes: &[u8], h: &TxHit, k: u32) -> u32 {
        let text = self.hit_region(tb, 0, codes.len() as u32, h);
        stats::add(&stats::MYERS_CALLS, 1);
        stats::add(
            &stats::MYERS_CELLS,
            (codes.len() as u64).div_ceil(64) * text.len() as u64,
        );
        align::distance(codes, &text, align::Mode::Infix, k).unwrap_or(k.saturating_add(1))
    }

    /// [`Self::raw_texts`], as packed ranges of the index (nothing decoded).
    fn raw_ranges<'a, 'r>(
        &'a self,
        tb: &TieBreak,
        raw: Option<&'r RawRead>,
        spans: &[(u32, u32, u32)],
    ) -> Option<(&'r RawRead, Vec<index::PackedRange<'a>>)> {
        let r = raw.filter(|_| tb.raw_ends)?;
        let texts = spans
            .iter()
            .map(|&(t, lo, hi)| self.idx.tx_range(t, lo, hi, true))
            .collect::<Option<Vec<_>>>()?;
        Some((r, texts))
    }

    /// Raw transcript texts for parsed spans `(transcript, lo, hi)`, if raw read ends are on
    /// and the index has every transcript's raw sequence.
    fn raw_texts<'r>(
        &self,
        tb: &TieBreak,
        raw: Option<&'r RawRead>,
        spans: &[(u32, u32, u32)],
    ) -> Option<(&'r RawRead, Vec<Vec<u8>>)> {
        let r = raw.filter(|_| tb.raw_ends)?;
        let texts = spans
            .iter()
            .map(|&(t, lo, hi)| self.idx.tx_raw_codes(t, lo, hi))
            .collect::<Option<Vec<_>>>()?;
        Some((r, texts))
    }

    /// Distances through anchors shared by every candidate: segments between consecutive
    /// shared anchors are aligned globally, but only where the candidates' transcript
    /// sequences differ (identical segments cost every candidate the same); the read ends are
    /// aligned against windows anchored at the outermost shared anchors. Anchor bases are
    /// excluded from every DP slice. `None` if no anchor (query start, length) is shared by
    /// every candidate.
    fn segmented(
        &self,
        tb: &TieBreak,
        codes: &[u8],
        raw: Option<&RawRead>,
        hits: &[TxHit],
    ) -> Option<(Vec<i64>, f32)> {
        let blocks = shared_blocks(hits)?;
        let plen = codes.len() as u32;
        let mut cost = vec![0i64; hits.len()];
        // Each candidate's cost for one segment, and with a `cut`, its cost against its text cut
        // at that many bases from the segment's start and the text bases its best alignment
        // covers, from the same pass (None if every candidate's text is the same: equal costs).
        let seg = |texts: Vec<index::PackedRange<'_>>,
                   read: &[u8],
                   mode: align::Mode,
                   rev: bool,
                   cut: Option<usize>|
         -> Option<(Vec<i64>, Vec<i64>, Vec<usize>)> {
            if texts.windows(2).all(|w| w[0].sequence_eq(&w[1])) {
                return None; // same cost for every candidate
            }
            let mut delta = vec![0i64; texts.len()];
            let mut delta_cut = vec![0i64; texts.len()];
            let mut cols = vec![0usize; texts.len()];
            stats::add(&stats::TIEBREAK_SEGMENTS, 1);
            // An insertion-only global gap depends only on its length. Avoid
            // copying DP inputs, sorting distinct strings, and allocating DP rows.
            if mode == align::Mode::Global
                && read.is_empty()
                && tb.affine_within.is_none()
                && let Some(sc) = &tb.affine
            {
                for (c, text) in delta.iter_mut().zip(&texts) {
                    if !text.is_empty() {
                        *c += i64::from(sc.open) + i64::from(sc.ext) * text.len() as i64;
                    }
                }
                return Some((delta.clone(), delta, cols));
            }
            // Tied isoforms often share a segment's sequence: compare the packed texts and
            // decode each distinct one once (`rep0`: each candidate's distinct text).
            let mut distinct: Vec<index::PackedRange<'_>> = Vec::new();
            let rep0: Vec<usize> = texts
                .iter()
                .map(|t| {
                    distinct.iter().position(|d| d.sequence_eq(t)).unwrap_or_else(|| {
                        distinct.push(*t);
                        distinct.len() - 1
                    })
                })
                .collect();
            let decoded: Vec<Vec<u8>> = distinct
                .iter()
                .map(|t| {
                    let mut v = if t.is_empty() { Vec::new() } else { t.decode() };
                    if rev {
                        v.reverse();
                    }
                    v
                })
                .collect();
            let rread: Vec<u8>;
            let read = if rev {
                rread = read.iter().rev().copied().collect();
                &rread[..]
            } else {
                read
            };
            match &tb.affine {
                Some(sc) => {
                    // Each distinct text is scored once, in sorted order (`rep`: each
                    // candidate's place in it). With `affine_within`, a capped unit-cost (Myers)
                    // pass comes first and only texts within that many edits of the best get an
                    // affine score; the rest are charged one mismatch (mat + mis) per extra edit.
                    let mut uniq: Vec<usize> = (0..decoded.len()).collect();
                    uniq.sort_unstable_by(|&a, &b| decoded[a].cmp(&decoded[b]));
                    let mut place = vec![0usize; decoded.len()];
                    for (k, &u) in uniq.iter().enumerate() {
                        place[u] = k;
                    }
                    let rep: Vec<usize> = rep0.iter().map(|&d| place[d]).collect();
                    let unit = match tb.affine_within {
                        Some(w) if uniq.len() > 1 => {
                            // The distinct texts: none repeats.
                            let ut: Vec<&[u8]> = uniq.iter().map(|&u| decoded[u].as_slice()).collect();
                            Some((Self::banded(read, &ut, mode, w, true), w))
                        }
                        _ => None,
                    };
                    let dmin = unit.as_ref().map_or(0, |(d, _)| *d.iter().min().unwrap());
                    if stats::on() {
                        // Text bases, and those outside the common prefix with the previous
                        // distinct text in sorted order (what prefix sharing would compute).
                        let mut prev: Option<&Vec<u8>> = None;
                        for &u in &uniq {
                            let t = &decoded[u];
                            let l = prev.map_or(0, |p| t.iter().zip(p).take_while(|(a, b)| a == b).count());
                            stats::add(&stats::AFFINE_TEXT_BASES, t.len() as u64);
                            stats::add(&stats::AFFINE_TEXT_NEW_BASES, (t.len() - l) as u64);
                            prev = Some(t);
                        }
                    }
                    let mut score = vec![None::<i64>; uniq.len()];
                    let mut score_cut = vec![None::<i64>; uniq.len()];
                    let mut score_cols = vec![0usize; uniq.len()];
                    let mut best: Option<(u32, i64)> = None; // (unit distance, affine score)
                    // Distinct texts to score, in sorted order, so isoforms that agree near
                    // the anchor share the rows that read only those bases.
                    let todo: Vec<usize> = (0..uniq.len())
                        .filter(|&k| unit.as_ref().is_none_or(|(d, w)| d[k] <= dmin + w))
                        .collect();
                    let diag = mode == align::Mode::Prefix && tb.clip_ends && tb.end_diagonal;
                    for &k in &todo {
                        let (m, n) = (read.len() as u64, decoded[uniq[k]].len() as u64);
                        stats::add(&stats::AFFINE_CALLS, 1);
                        let drift = if diag && n >= m { 0 } else { m.abs_diff(n) };
                        let cells = m * (2 * tb.band as u64 + drift + 1).min(n + 1);
                        stats::add(&stats::AFFINE_CELLS, cells);
                        if mode == align::Mode::Prefix {
                            stats::add(&stats::AFFINE_END_CELLS, cells);
                        }
                    }
                    let tx: Vec<&[u8]> = todo.iter().map(|&k| decoded[uniq[k]].as_slice()).collect();
                    let w0 = align::dp_work();
                    let xs = if mode == align::Mode::Prefix && tb.clip_ends {
                        // The band is laid out for the window the read end would get against a
                        // transcript that doesn't end inside it (`tb.window`, as the windows
                        // below are cut), or for the longest text: texts are windows sized in
                        // parsed bases, which in raw bases (raw read ends) or against a
                        // tail-trimmed read end can be longer. Every text then gets the same band,
                        // and sorted texts share rows.
                        let window = (tb.window(read.len() as u32) as usize)
                            .max(tx.iter().map(|t| t.len()).max().unwrap_or(0));
                        align::affine_extend_sorted_cut(
                            read,
                            &tx,
                            sc,
                            tb.band,
                            tb.end_bonus,
                            tb.xdrop,
                            diag,
                            window,
                            cut.unwrap_or(usize::MAX),
                        )
                    } else {
                        align::affine_sorted(read, &tx, mode, sc, tb.band)
                            .into_iter()
                            .map(|x| (x, x, 0))
                            .collect()
                    };
                    if stats::on() {
                        let w1 = align::dp_work();
                        let (c, r) = if mode == align::Mode::Prefix {
                            (&stats::DP_TBEND_CELLS, &stats::DP_TBEND_ROWS)
                        } else {
                            (&stats::DP_SEG_CELLS, &stats::DP_SEG_ROWS)
                        };
                        stats::add(c, w1.0 - w0.0);
                        stats::add(r, w1.1 - w0.1);
                    }
                    for (&k, &(x, xc, cc)) in todo.iter().zip(&xs) {
                        let x = i64::from(x);
                        score[k] = Some(x);
                        score_cut[k] = Some(i64::from(xc));
                        score_cols[k] = cc;
                        let du = unit.as_ref().map_or(0, |(d, _)| d[k]);
                        if best.is_none_or(|(bd, bx)| (du, -x) < (bd, -bx)) {
                            best = Some((du, x));
                        }
                    }
                    let per_edit = i64::from(sc.mat + sc.mis);
                    for (i, c) in delta.iter_mut().enumerate() {
                        let k = rep[i];
                        let x = score[k].unwrap_or_else(|| {
                            let (bd, bx) = best.unwrap();
                            let d = unit.as_ref().unwrap().0[k];
                            bx - per_edit * i64::from(d - bd)
                        });
                        *c -= x;
                        delta_cut[i] -= score_cut[k].unwrap_or(x);
                        cols[i] = score_cols[k];
                    }
                }
                None => {
                    let texts: Vec<&[u8]> = rep0.iter().map(|&d| decoded[d].as_slice()).collect();
                    for (c, x) in delta
                        .iter_mut()
                        .zip(Self::banded(read, &texts, mode, tb.cap, false))
                    {
                        *c += i64::from(x);
                    }
                    delta_cut.copy_from_slice(&delta);
                }
            }
            Some((delta, delta_cut, cols))
        };
        // A read end against each candidate's window; `short` marks windows that stop at their
        // transcript's end. Past the shortest such end the read has no counterpart in that
        // candidate, and what a longer candidate gains there (read bases aligned by chance into
        // its flank: primers, adapters) counts only if it reaches `flank_min` points.
        let end = |texts: Vec<index::PackedRange<'_>>,
                   read: &[u8],
                   rev: bool,
                   short: &[bool]|
         -> Option<Vec<i64>> {
            let cut = texts
                .iter()
                .zip(short)
                .filter(|&(_, &s)| s)
                .map(|(t, _)| t.len())
                .min();
            let Some(cut) = cut.filter(|&c| tb.flank_min > 0 && texts.iter().any(|t| t.len() > c)) else {
                return seg(texts, read, align::Mode::Prefix, rev, None).map(|(d, _, _)| d);
            };
            // The cut counts from the read end's start (a left-end window is aligned backwards,
            // so it keeps the bases next to the block).
            let long: Vec<bool> = texts.iter().map(|t| t.len() > cut).collect();
            let (mut out, near, cols) = seg(texts, read, align::Mode::Prefix, rev, Some(cut))?;
            // Gains past the cut count if they reach `flank_min`, or if they come from a short,
            // nearly exact extension (at least 3/4 of a match per text base past the cut):
            // isoforms whose ends differ by a few bases are told apart, a primer's chance
            // matches against a long flank are not.
            let mat = tb.affine.as_ref().map_or(1, |sc| i64::from(sc.mat));
            for (((o, &n), &l), &c) in out.iter_mut().zip(&near).zip(&long).zip(&cols) {
                let gain = n - *o;
                let past = c.saturating_sub(cut) as i64;
                if l && gain < tb.flank_min && !(past > 0 && 4 * gain >= 3 * mat * past) {
                    *o = n;
                }
            }
            Some(out)
        };
        let mut add = |d: Option<Vec<i64>>| {
            if let Some(d) = d {
                for (c, x) in cost.iter_mut().zip(d) {
                    *c += x;
                }
            }
        };
        // Interior segments in raw bases (`TieBreak::raw_mid`), from each shared block's
        // midpoint to the next one's: next to a block's edge the candidates differ, and a
        // boundary there forces the read's bases around it onto one side.
        let raw_mid = raw.filter(|_| {
            tb.raw_ends && tb.raw_mid && hits.iter().all(|h| self.idx.tx_raw_pos(h.tx, 0).is_some())
        });
        let half = |b: &SharedBlock| -> u32 {
            if raw_mid.is_some() {
                (b.query.end - b.query.start) / 2
            } else {
                0
            }
        };
        let starts = |b: &SharedBlock| -> Vec<u32> {
            hits.iter()
                .zip(&b.targets)
                .map(|(h, t)| self.idx.tx_raw_pos(h.tx, t.start + half(b)).unwrap())
                .collect()
        };
        let mut prev = raw_mid.map(|_| starts(&blocks[0]));
        for pair in blocks.windows(2) {
            if let (Some(r), Some(lo)) = (raw_mid, prev.as_mut()) {
                let hi = starts(&pair[1]);
                let texts = hits
                    .iter()
                    .enumerate()
                    .map(|(i, h)| self.idx.tx_raw_slice(h.tx, lo[i], hi[i]).unwrap())
                    .collect();
                add(seg(
                    texts,
                    &r.codes[r.pos(pair[0].query.start + half(&pair[0]))
                        ..r.pos(pair[1].query.start + half(&pair[1]))],
                    align::Mode::Global,
                    false,
                    None,
                )
                .map(|(d, _, _)| d));
                *lo = hi;
                continue;
            }
            let (left, right) = (&pair[0], &pair[1]);
            let lo0 = left.targets[0].end;
            let len0 = right.targets[0].start - lo0;
            if hits.iter().enumerate().skip(1).all(|(i, h)| {
                let lo = left.targets[i].end;
                let len = right.targets[i].start - lo;
                len == len0 && self.idx.tx_codes_eq(hits[0].tx, lo0, h.tx, lo, len)
            }) {
                continue;
            }
            let texts = hits
                .iter()
                .enumerate()
                .map(|(i, h)| {
                    self.idx
                        .tx_range(h.tx, left.targets[i].end, right.targets[i].start, false)
                        .unwrap()
                })
                .collect();
            add(seg(
                texts,
                &codes[left.query.end as usize..right.query.start as usize],
                align::Mode::Global,
                false,
                None,
            )
            .map(|(d, _, _)| d));
        }
        // Left end: read[0, q0) against transcript windows ending at each candidate's anchor.
        let first = &blocks[0];
        let h0 = half(first);
        let q0 = first.query.start + h0;
        if q0 > 0 {
            let spans: Vec<(u32, u32, u32)> = hits
                .iter()
                .zip(&first.targets)
                .map(|(h, t)| (h.tx, (t.start + h0).saturating_sub(tb.window(q0)), t.start + h0))
                .collect();
            let short: Vec<bool> = spans.iter().map(|&(_, lo, _)| lo == 0).collect();
            add(match self.raw_ranges(tb, raw, &spans) {
                Some((r, texts)) => end(texts, &r.codes[..r.pos(q0)], true, &short),
                None => end(
                    spans
                        .iter()
                        .map(|&(t, lo, hi)| self.idx.tx_range(t, lo, hi, false).unwrap())
                        .collect(),
                    &codes[..q0 as usize],
                    true,
                    &short,
                ),
            });
        }
        // Right end: read[last block end, plen) against windows from each candidate's block end.
        // Blocks are exact in parsed space only: if the candidates' raw bases differ inside the
        // last block (homopolymer run lengths), the raw read end is aligned through it.
        let last = blocks.last().unwrap();
        let t0 = &last.targets[0];
        let through = raw.is_some()
            && tb.raw_ends
            && self
                .idx
                .tx_range(hits[0].tx, t0.start, t0.end, true)
                .is_some_and(|r0| {
                    hits.iter().zip(&last.targets).skip(1).any(|(h, t)| {
                        self.idx
                            .tx_range(h.tx, t.start, t.end, true)
                            .is_some_and(|r| !r.sequence_eq(&r0))
                    })
                });
        let hl = half(last);
        let ql = if hl > 0 {
            last.query.start + hl
        } else if through {
            last.query.start
        } else {
            last.query.end
        };
        let rest = plen - ql;
        let spans: Vec<(u32, u32, u32)> = hits
            .iter()
            .zip(&last.targets)
            .map(|(h, t)| {
                let tlen = self.idx.tx_plen()[h.tx as usize];
                let s = if hl > 0 {
                    t.start + hl
                } else if through {
                    t.start
                } else {
                    t.end
                };
                (h.tx, s, (s + tb.window(rest)).min(tlen))
            })
            .collect();
        // The read's poly(A) tail, searched for past the last block (raw positions from ql).
        let a_code = lcp_core::encode::CODE[b'A' as usize] & 3;
        // Whether each window's end (its length from the read end's start, in the tail's space)
        // is where the tail starts: within the tolerance of its start, or further but only by
        // annotated A's (an end that includes part of the tail) up to the run's end.
        let ends_at_tail = |texts: &[index::PackedRange<'_>], cut: usize, run_end: usize| -> Vec<bool> {
            let tol = TAIL_TOL + cut / 10;
            texts
                .iter()
                .map(|t| {
                    let n = t.len();
                    n.abs_diff(cut) <= tol
                        || (n > cut && n <= run_end + tol && t.count_other(cut, a_code) * 10 <= n - cut)
                })
                .collect()
        };
        let tail = match raw {
            _ if !tb.tail_trim && tb.tail_bonus <= 0 => None,
            Some(r) => {
                let (r0, r1) = (r.pos(ql), r.pos(last.query.end));
                tail_start(&r.codes[r1..], a_code, tb.tail_run)
                    .map(|(c, e)| (Some(r), r0, r1 - r0 + c, r1 - r0 + e))
            }
            // Without HPC the parsed read is the raw read (and `ql` the last block's end).
            None if !self.seeder.hpc => {
                tail_start(&codes[ql as usize..], a_code, tb.tail_run).map(|(c, e)| (None, ql as usize, c, e))
            }
            None => None,
        };
        let short: Vec<bool> = spans
            .iter()
            .map(|&(tx, _, hi)| hi == self.idx.tx_plen()[tx as usize])
            .collect();
        let ends = match self.raw_ranges(tb, raw, &spans) {
            Some((r, texts)) => {
                let rest = &r.codes[r.pos(ql)..];
                let cut = tail.map_or(rest.len(), |(_, _, c, _)| c);
                let run_end = tail.map_or(cut, |(_, _, _, e)| e);
                let mut at_end = tail.map(|_| ends_at_tail(&texts, cut, run_end));
                // A candidate with its own A-run where the read's tail starts explains the
                // tail as well as one that ends there. One that ends inside or just after
                // its run (an annotated tail) ends at the tail; one that continues gains the
                // bonus too, and what the read's bases after the run gain against its
                // continuation (if that reaches the bonus: the read went on along it).
                if let Some(at) = at_end.as_mut()
                    && tb.tail_templated
                {
                    let tol = TAIL_TOL + cut / 10;
                    let mut d = vec![0i64; texts.len()];
                    let mut cont = Vec::with_capacity(texts.len());
                    for (i, (t, &(tx, lo, _))) in texts.iter().zip(&spans).enumerate() {
                        let raw_lo = self.idx.tx_raw_pos(tx, lo).unwrap();
                        let n = t.len();
                        let run = (!t.is_empty())
                            .then(|| {
                                own_run(
                                    &t.decode(),
                                    a_code,
                                    tb.tail_run,
                                    cut.saturating_sub(tol),
                                    cut + tol,
                                )
                            })
                            .flatten();
                        let mut from = n;
                        if let Some((_, e)) = run {
                            if short[i] && n - e <= TAIL_TOL {
                                at[i] = true;
                            } else if !at[i] {
                                d[i] = -i64::from(tb.tail_bonus);
                                from = e;
                            }
                        }
                        cont.push(
                            self.idx
                                .tx_raw_slice(tx, raw_lo + from as u32, raw_lo + n as u32)
                                .unwrap(),
                        );
                    }
                    let after = &rest[run_end..];
                    if tb.tail_trim
                        && !after.is_empty()
                        && let Some((g, _, _)) = seg(cont, after, align::Mode::Prefix, false, None)
                    {
                        for (x, g) in d.iter_mut().zip(g) {
                            if -g >= i64::from(tb.tail_bonus) {
                                *x += g;
                            }
                        }
                    }
                    add(Some(d));
                }
                let read = if tb.tail_trim { &rest[..cut] } else { rest };
                add(end(texts, read, false, &short));
                at_end
            }
            None => {
                let texts: Vec<index::PackedRange<'_>> = spans
                    .iter()
                    .map(|&(t, lo, hi)| self.idx.tx_range(t, lo, hi, false).unwrap())
                    .collect();
                // The tail's start and run end in parsed space: the last parsed base starting at
                // or before each.
                let parsed = |x: usize, r: &RawRead| {
                    let (mut lo, mut hi) = (ql, plen);
                    while lo < hi {
                        let mid = (lo + hi).div_ceil(2);
                        if r.pos(mid) <= x { lo = mid } else { hi = mid - 1 }
                    }
                    (lo - ql) as usize
                };
                let cut = tail.map_or((plen - ql) as usize, |(r, r0, c, _)| {
                    r.map_or(c, |r| parsed(r0 + c, r))
                });
                let run_end = tail.map_or(cut, |(r, r0, _, e)| r.map_or(e, |r| parsed(r0 + e, r)));
                let at_end = tail.map(|_| ends_at_tail(&texts, cut, run_end));
                let read = &codes[ql as usize..];
                let read = if tb.tail_trim { &read[..cut] } else { read };
                add(end(texts, read, false, &short));
                at_end
            }
        };
        // Candidates whose 3' end is where the read's tail starts gain the bonus.
        if let Some(at_end) = ends
            && tb.tail_bonus > 0
        {
            for ((c, &(tx, _, hi)), &ok) in cost.iter_mut().zip(&spans).zip(&at_end) {
                if ok && hi == self.idx.tx_plen()[tx as usize] {
                    *c -= i64::from(tb.tail_bonus);
                }
            }
        }
        let unit = if tb.affine.is_some() { tb.lambda } else { tb.kappa };
        Some((cost, unit))
    }

    /// Transcript occurrences (transcript, first base) of the seed at index nodes `nodes`, until
    /// `f` returns true; whether it did.
    fn node_occurrences(&self, nodes: std::ops::Range<u32>, mut f: impl FnMut(u32, u32) -> bool) -> bool {
        let flags = self.idx.node_flags();
        for node in nodes {
            if flags[node as usize] & index::FLAG_MASKED != 0 {
                continue;
            }
            let off = self.idx.node_offset()[node as usize];
            let (ts, ps) = self.idx.occurrences(self.idx.node_unitig()[node as usize]);
            if ts.iter().zip(ps).any(|(&t, &p)| f(t, p + off)) {
                return true;
            }
        }
        false
    }

    /// The [`AntiSeed`]s of transcript `x`, sorted by position.
    fn anti_seeds(&self, x: u32) -> &[AntiSeed] {
        self.anti[x as usize].get_or_init(|| {
            let tlen = self.idx.tx_plen()[x as usize];
            let mut codes = self.idx.tx_codes(x, 0, tlen);
            codes.reverse();
            codes.iter_mut().for_each(|c| *c = 3 - *c);
            let mut seeds = Vec::new();
            self.seeder.seed_codes_into(&codes, &mut seeds, &mut Vec::new());
            let mut out: Vec<AntiSeed> = seeds
                .iter()
                .filter_map(|s| {
                    let r = self.idx.lookup(s.label);
                    let found = !r.is_empty()
                        && r.len() <= self.p.anchor.h_max
                        && self.node_occurrences(r.clone(), |t, _| t != x);
                    found.then_some(AntiSeed {
                        x: tlen - s.end,
                        len: s.end - s.start,
                        nodes: (r.start, r.end),
                    })
                })
                .collect();
            out.sort_unstable_by_key(|a| a.x);
            out.into_boxed_slice()
        })
    }

    /// The range of diagonals of `h` and the transcript interval a read of `plen` parsed
    /// bases on them covers.
    fn hit_span(&self, h: &TxHit, plen: u32) -> ((i64, i64), (u32, u32)) {
        let tlen = i64::from(self.idx.tx_plen()[h.tx as usize]);
        let (d1, d2) = (
            i64::from(h.tstart) - i64::from(h.qstart),
            i64::from(h.tend) - i64::from(h.qend),
        );
        let d = (d1.min(d2), d1.max(d2));
        (
            d,
            (
                d.0.clamp(0, tlen) as u32,
                (d.1 + i64::from(plen)).clamp(0, tlen) as u32,
            ),
        )
    }

    /// The transcripts of `found` (transcript -> weight and diagonal range): the `k` of
    /// largest weight first, in that order, then the others in no order.
    fn mates(found: &mut Mates, k: usize) -> Vec<(u32, u32, i64, i64)> {
        let mut out: Vec<(u32, u32, i64, i64)> =
            found.drain().map(|(t, (w, lo, hi))| (t, w, lo, hi)).collect();
        let key = |&(t, w, _, _): &(u32, u32, i64, i64)| (std::cmp::Reverse(w), t);
        if k > 0 && out.len() > k {
            out.select_nth_unstable_by_key(k - 1, key);
            out[..k].sort_unstable_by_key(key);
        } else {
            out.sort_unstable_by_key(key);
        }
        out
    }

    /// Transcripts sharing sequence with the region of hit `h` on the opposite strand, with
    /// the diagonals the read's other orientation would lie on: the `k` sharing most first,
    /// in that order, then the others in no order.
    fn antisense_mates(&self, h: &TxHit, plen: u32, k: usize) -> Vec<(u32, u32, i64, i64)> {
        let seeds = self.anti_seeds(h.tx);
        if seeds.is_empty() {
            return Vec::new();
        }
        let ((dlo, dhi), (lo, hi)) = self.hit_span(h, plen);
        let seeds = &seeds[seeds.partition_point(|a| a.x < lo)..seeds.partition_point(|a| a.x < hi)];
        // A seed at read position q of this orientation lies at plen - q - len of the other;
        // its position here + its first base in the mate + its length is constant along one
        // antisense match.
        MATES.with(|found| {
            let found = &mut *found.borrow_mut();
            for a in seeds {
                self.node_occurrences(a.nodes.0..a.nodes.1, |t, p| {
                    if t != h.tx {
                        let d = i64::from(a.x + p + a.len) - i64::from(plen);
                        let e = found.entry(t).or_insert((0, d - dhi, d - dlo));
                        e.0 += a.len;
                        e.1 = e.1.min(d - dhi);
                        e.2 = e.2.max(d - dlo);
                    }
                    false
                });
            }
            Self::mates(found, k)
        })
    }

    fn map_strand(&self, (seq, raw_read): (Vec<u8>, Option<RawRead>), salt: u64) -> Strand {
        let seeds = self.seeder.seed_parsed(&seq);
        let (mut cands, _) = anchor::select(&self.idx, &self.node_fam, &seeds, &self.p.anchor);
        if self.p.merge_collinear {
            for c in &mut cands {
                anchor::merge_collinear(&self.idx, c);
            }
        }
        let mut hits: Vec<TxHit> = cands
            .iter()
            .flat_map(|c| self.chainer.chain(&self.idx, c))
            .collect();
        sort_distinct_hits(&mut hits, salt);
        Strand {
            hits,
            seq,
            codes: Default::default(),
            raw: raw_read,
            seeds,
        }
    }

    /// [`Self::parse`] of `raw` and of its reverse complement. A read of A/C/G/T only (almost
    /// all) is compressed once: the reverse strand's compression and codes mirror the forward
    /// one's.
    fn parse_both(seeder: &BaseSeeder, raw: &[u8]) -> [(Vec<u8>, Option<RawRead>); 2] {
        let acgt = || raw.iter().all(|&b| lcp_core::encode::CODE[b as usize] < 4);
        match seeder.parsed(raw) {
            Some(h) if acgt() => {
                let rc = h.revcomp_acgt();
                let codes = codes(raw);
                let rcodes = codes.iter().rev().map(|&c| 3 - c).collect();
                let done = |mut h: lcp_core::hpc::Hpc, codes| {
                    let seq = std::mem::take(&mut h.seq);
                    (seq, Some(RawRead { codes, hpc: h }))
                };
                [done(h, codes), done(rc, rcodes)]
            }
            _ => [Self::parse(seeder, raw), Self::parse(seeder, &revcomp(raw))],
        }
    }

    /// Parsed sequence and, with HPC, the raw read for raw-space read ends.
    fn parse(seeder: &BaseSeeder, raw: &[u8]) -> (Vec<u8>, Option<RawRead>) {
        match seeder.parsed(raw) {
            Some(mut h) => {
                let seq = std::mem::take(&mut h.seq);
                let codes = codes(raw);
                (seq, Some(RawRead { codes, hpc: h }))
            }
            None => (raw.to_vec(), None),
        }
    }

    /// Extend `h`'s span over the read flanks outside its outermost anchors (`codes`: the
    /// parsed read in `h`'s orientation). The transcript flank is a window of the read flank's
    /// length plus slack, cut short at the transcript's end. The extension depends on it only
    /// through the bases it read (and, if it read them all, the flank's length), so isoforms
    /// whose flanks agree there share it; `seen` holds the read's extensions.
    fn extend_ends(&self, e: &Ends, codes: &[u8], h: &mut TxHit, seen: &mut EndCache) {
        let plen = codes.len() as u32;
        let tlen = self.idx.tx_plen()[h.tx as usize];
        let window = |x: u32| x + e.pad + x / e.pad_div;
        let mut ext = |left: bool, q: u32, read: &[u8], text: Vec<u8>| -> (u32, u32) {
            stats::add(&stats::SPAN_ENDS, 1);
            let done = seen.entry((left, h.reverse, q)).or_default();
            let hit = done
                .iter()
                .find(|(p, all, _)| text.starts_with(p) && (!*all || text.len() == p.len()));
            if let Some((_, _, r)) = hit {
                return *r;
            }
            let (m, w) = (read.len() as u64, window(read.len() as u32) as usize);
            stats::add(&stats::SPAN_END_CALLS, 1);
            stats::add(&stats::SPAN_END_CELLS, m * (2 * e.band as u64 + w as u64 - m + 1));
            let w0 = align::dp_work();
            let (_, qi, tj, used) =
                align::affine_extend_window(read, &text, w, &e.sc, e.band, e.end_bonus, e.xdrop);
            if stats::on() {
                let w1 = align::dp_work();
                stats::add(&stats::DP_SPAN_CELLS, w1.0 - w0.0);
                stats::add(&stats::DP_SPAN_ROWS, w1.1 - w0.1);
            }
            let r = (qi as u32, tj as u32);
            done.push((text[..used].to_vec(), used == text.len(), r));
            r
        };
        if h.qstart > 0 && h.tstart > 0 {
            let lo = h.tstart.saturating_sub(window(h.qstart));
            let mut text = self.idx.tx_codes(h.tx, lo, h.tstart);
            text.reverse();
            let read: Vec<u8> = codes[..h.qstart as usize].iter().rev().copied().collect();
            let (qi, tj) = ext(true, h.qstart, &read, text);
            h.qstart -= qi;
            h.tstart -= tj;
        }
        let rest = plen - h.qend;
        if rest > 0 && h.tend < tlen {
            let hi = (h.tend + window(rest)).min(tlen);
            let text = self.idx.tx_codes(h.tx, h.tend, hi);
            let (qi, tj) = ext(false, h.qend, &codes[h.qend as usize..], text);
            h.qend += qi;
            h.tend += tj;
        }
    }

    fn span_frac(h: Option<&TxHit>, plen: usize) -> f32 {
        h.map_or(0.0, |h| (h.qend - h.qstart) as f32 / plen.max(1) as f32)
    }

    pub fn map(&self, raw: &[u8]) -> ReadResult {
        stats::add(&stats::READS, 1);
        // Direct RNA reads spell T as U; seeds hash bases, so normalize first.
        let dna: Vec<u8>;
        let raw = if raw.iter().any(|&b| b == b'U' || b == b'u') {
            dna = raw
                .iter()
                .map(|&b| if b == b'U' || b == b'u' { b'T' } else { b })
                .collect();
            &dna[..]
        } else {
            raw
        };
        // A skipped orientation has no seeds, hence no hits; it keeps its raw read for
        // `ReadResult::hpc` (the forward runs, also when only the reverse is mapped).
        let skip = |p: (Vec<u8>, Option<RawRead>)| Strand {
            hits: Vec::new(),
            seq: p.0,
            codes: Default::default(),
            raw: p.1,
            seeds: Vec::new(),
        };
        // Equal scores are ordered by a hash of the transcript salted with the read, so ties
        // are broken the same way on every run without favouring the same transcript.
        let salt = xxhash_rust::xxh3::xxh3_64(raw);
        let [pf, pr] = match self.p.strands {
            // The reverse orientation would be skipped: not parsed either.
            Strands::Forward => [Self::parse(&self.seeder, raw), (Vec::new(), None)],
            _ => Self::parse_both(&self.seeder, raw),
        };
        let (fwd, rc) = match self.p.strands {
            Strands::Both => (self.map_strand(pf, salt), self.map_strand(pr, salt)),
            Strands::Forward => (self.map_strand(pf, salt), skip(pr)),
            Strands::Reverse => (skip(pf), self.map_strand(pr, salt)),
        };
        let best = |h: &[TxHit]| h.first().map_or(f32::NEG_INFINITY, |x| x.score);
        // Whether a best chain is more than a chance match (no `min_score`: `none` says).
        let passes = |h: Option<&TxHit>, rescued: bool, none: bool| {
            h.is_some_and(|t| self.p.min_score.map_or(none, |m| m.passes(t, raw.len(), rescued)))
        };
        let low = |h: &[TxHit], rescued: bool| !h.is_empty() && !passes(h.first(), rescued, true);
        let (mut reverse, mut st, mut other) = if best(&rc.hits) > best(&fwd.hits) {
            (true, rc, fwd)
        } else {
            (false, fwd, rc)
        };
        let mut rescued = false;
        // Hierarchical rescue when the base-level chain covers too little of the read.
        if let Some(r) = &self.rescuer
            && Self::span_frac(st.hits.first(), st.seq.len()) < self.p.pi_cov
        {
            let rescue = |s: &Strand| {
                let mut h = r.rescue(
                    &self.idx,
                    &self.node_fam,
                    &s.seq,
                    &s.seeds,
                    &self.p.anchor,
                    &self.p.chain,
                );
                sort_hits(&mut h, salt);
                h
            };
            let a = rescue(&st);
            // The other orientation is re-seeded only if its base-level best comes within
            // `other_margin` of this one's (None: always) or is, by `min_score`, more than
            // a chance match itself: this one's lead means little when neither chain
            // covers the read.
            let try_other = match r.p.other_margin {
                None => true,
                Some(m) => other
                    .hits
                    .first()
                    .is_some_and(|o| o.score >= best(&st.hits) - m || passes(Some(o), false, false)),
            };
            let b = if try_other { rescue(&other) } else { Vec::new() };
            let (rev, hits, alt, seq) = if best(&b) > best(&a) {
                (!reverse, b, a, &other.seq)
            } else {
                (reverse, a, b, &st.seq)
            };
            if Self::span_frac(hits.first(), seq.len()) > Self::span_frac(st.hits.first(), st.seq.len()) {
                if rev != reverse {
                    std::mem::swap(&mut st, &mut other);
                }
                st.hits = hits;
                other.hits = alt;
                reverse = rev;
                rescued = true;
            }
        }
        // The rescue level's table of every transcript, if the mapper has one.
        let global = self
            .rescuer
            .as_ref()
            .and_then(|r| r.p.global.filter(|_| r.has_global()).map(|g| (r, g)));
        // A read both steps leave weak (a short best chain, or a rescued one scoring near the
        // threshold) is seeded against every transcript: the candidate genes of the rescue
        // come from the read's base-level seeds, and a read with few of those may not name
        // its own gene.
        if let Some((r, g)) = global {
            let t = self.p.min_score.unwrap_or_default().at(raw.len());
            let weak = st.hits.first().is_none_or(|h| {
                Self::span_frac(Some(h), st.seq.len()) < self.p.pi_cov
                    || (rescued && h.score < t + g.weak_margin)
            });
            if weak {
                let seed = |s: &Strand| {
                    let mut h = r.global_hits(&self.idx, &s.seq, &self.p.chain, g.weak_cands);
                    sort_hits(&mut h, salt);
                    h
                };
                let (a, b) = (seed(&st), seed(&other));
                if rescued {
                    // Rescue-level chains on both sides: the higher score per transcript.
                    let join = |hits: &mut Vec<TxHit>, mut new: Vec<TxHit>| {
                        new.retain(|n| hits.iter().all(|h| h.tx != n.tx || h.score < n.score));
                        hits.retain(|h| new.iter().all(|n| n.tx != h.tx));
                        hits.extend(new);
                        sort_hits(hits, salt);
                    };
                    join(&mut st.hits, a);
                    join(&mut other.hits, b);
                    if best(&other.hits) > best(&st.hits) {
                        std::mem::swap(&mut st, &mut other);
                        reverse = !reverse;
                    }
                } else {
                    // As the rescue: its chains replace base-level ones if they span more.
                    let (rev, hits, alt, seq) = if best(&b) > best(&a) {
                        (!reverse, b, a, &other.seq)
                    } else {
                        (reverse, a, b, &st.seq)
                    };
                    if Self::span_frac(hits.first(), seq.len())
                        > Self::span_frac(st.hits.first(), st.seq.len())
                    {
                        if rev != reverse {
                            std::mem::swap(&mut st, &mut other);
                        }
                        st.hits = hits;
                        other.hits = alt;
                        reverse = rev;
                        rescued = true;
                    }
                }
            } else if rescued {
                // Transcripts the rescue's candidate genes left out join this orientation's
                // hits.
                let mut h = r.global_hits(&self.idx, &st.seq, &self.p.chain, g.extra_cands);
                h.retain(|n| st.hits.iter().all(|h| h.tx != n.tx));
                if !h.is_empty() {
                    st.hits.extend(h);
                    sort_hits(&mut st.hits, salt);
                }
            }
        }
        let plen = st.seq.len() as u32;
        // Transcripts the read's own seeds did not nominate, found through the best hit's
        // reference region instead (error-free, unlike the read), on the opposite strand: for
        // the read's other orientation, whose seeds are independent of this one's.
        let mut anti = false;
        // The best hits were chained again in their own windows at the rescue level.
        let mut targeted = false;
        if let (Some(r), Some(top)) = (&self.rescuer, st.hits.first().cloned()) {
            // Chains of `new` replace those of the same transcripts unless they score lower.
            let merge = |hits: &mut Vec<TxHit>, mut new: Vec<TxHit>| {
                new.retain(|n| hits.iter().all(|h| h.tx != n.tx || h.score < n.score));
                if new.is_empty() {
                    return false;
                }
                hits.retain(|h| new.iter().all(|n| n.tx != h.tx));
                hits.extend(new);
                sort_hits(hits, salt);
                true
            };
            // An orientation's hits within the tie-break margin, or `w` × the best score, of
            // its best (`top`), chained again at the rescue level in their own windows.
            let rechain = |hits: &[TxHit], seq: &[u8], top: &TxHit, w: f32| {
                let margin = self
                    .p
                    .tiebreak
                    .as_ref()
                    .map_or(0.0, |tb| tb.margin_abs.max(tb.margin_rel * top.score.abs()))
                    .max(w * top.score.abs());
                let txs: Vec<(u32, i64, i64)> = hits
                    .iter()
                    .take_while(|h| h.score >= top.score - margin)
                    .take(rescue::RECHAIN_MAX)
                    .map(|h| {
                        let ((lo, hi), _) = self.hit_span(h, plen);
                        (h.tx, lo, hi)
                    })
                    .collect();
                r.rescue_txs(&self.idx, seq, &txs, &self.p.chain)
            };
            // A read whose base-level chain anchors enough of it, with no rival on the other
            // orientation, is settled: nothing is nominated for it.
            let settled = !rescued
                && r.p.settled.is_some_and(|s| {
                    top.anchors >= s.anchors
                        && top.anchored as f32 >= s.anchored * plen as f32
                        && other.hits.first().is_none_or(|o| o.score < s.other * top.score)
                });
            let mates = if r.p.antisense > 0 && matches!(self.p.strands, Strands::Both) && !settled {
                self.antisense_mates(&top, plen, r.p.antisense)
            } else {
                Vec::new()
            };
            if let Some(w) = r.p.rechain
                && !rescued
                && !settled
            {
                // What is nominated below is chained at the rescue level, whose chains score
                // higher than base-level ones: so are this orientation's best hits, to compare
                // like with like.
                let new = rechain(&st.hits, &st.seq, &top, w);
                targeted = merge(&mut st.hits, new);
                // Transcripts the base-level candidates left out join the hits.
                if let Some((_, g)) = global {
                    let new = r.global_hits(&self.idx, &st.seq, &self.p.chain, g.extra_cands);
                    merge(&mut st.hits, new);
                }
            }
            if !mates.is_empty() {
                let txs: Vec<(u32, i64, i64)> = mates
                    .iter()
                    .take(r.p.antisense)
                    .map(|&(t, _, lo, hi)| (t, lo, hi))
                    .collect();
                stats::add(&stats::ANTISENSE_TX, txs.len() as u64);
                merge(
                    &mut other.hits,
                    r.rescue_txs(&self.idx, &other.seq, &txs, &self.p.chain),
                );
                // The other orientation's best hit shares sequence with this one's region:
                // the two are compared by alignment, beyond the tie-break's orientation margin.
                let floor = r.p.antisense_min * st.hits.first().map_or(0.0, |a| a.score);
                anti = other
                    .hits
                    .first()
                    .is_some_and(|o| o.score >= floor && mates.iter().any(|m| m.0 == o.tx));
            }
            // The other orientation, if it competes, gets what this one got: its best hits
            // chained at the rescue level, then the transcripts its candidates left out.
            if let Some((_, g)) = global
                && g.other
                && !settled
                && let (Some(a), Some(b)) = (st.hits.first(), other.hits.first().cloned())
                && let Some(tb) = &self.p.tiebreak
                && (anti || b.score >= a.score - tb.margin_abs.max(tb.orient_margin_rel * a.score.abs()))
            {
                if let Some(w) = r.p.rechain
                    && !rescued
                {
                    let new = rechain(&other.hits, &other.seq, &b, w);
                    merge(&mut other.hits, new);
                }
                let new = r.global_hits(&self.idx, &other.seq, &self.p.chain, g.extra_cands);
                merge(&mut other.hits, new);
            }
        }
        let rescued_any = rescued || targeted;
        let reject = |low: bool| self.p.min_score.is_some_and(|m| low && !m.keep);
        // Nothing below would raise the best score: skip the alignments of a read that is
        // rejected anyway.
        if reject(low(&st.hits, rescued_any)) {
            return ReadResult {
                len: raw.len() as u32,
                plen,
                reverse,
                hits: Vec::new(),
                low_coverage: true,
                rescued,
                low_score: true,
                hpc: if reverse { other.raw } else { st.raw }.map(|r| r.hpc),
            };
        }
        let mut both = false;
        if let Some(tb) = &self.p.tiebreak {
            let margin = |s: f32| tb.margin_abs.max(tb.orient_margin_rel * s.abs());
            both = matches!((st.hits.first(), other.hits.first()),
                    (Some(a), Some(b)) if b.score >= a.score - margin(a.score));
            let codes = || strand_codes(&st.codes, &st.seq);
            self.tiebreak(tb, codes, st.raw.as_ref(), &mut st.hits, salt);
            // The two best hits are one sequence on opposite strands: the read fits both
            // equally, however far apart seed luck leaves the two chains.
            let same = |a: &Strand, b: &Strand| match (a.hits.first(), b.hits.first()) {
                (Some(a), Some(b)) => self.same_antisense(tb, plen, a, b),
                _ => false,
            };
            let mut mirror = tb.orient_same && same(&st, &other);
            both |= mirror || anti;
            if both {
                // Near-tied orientations (antisense transcripts sharing sequence): re-rank each
                // orientation, then place the other's scores relative to this one by the
                // whole-read edit distances of the two best transcripts.
                stats::add(&stats::BOTH_ORIENTATIONS, 1);
                let top = other.hits[0].tx;
                let codes = || strand_codes(&other.codes, &other.seq);
                self.tiebreak(tb, codes, other.raw.as_ref(), &mut other.hits, salt);
                if tb.orient_same && other.hits[0].tx != top {
                    mirror = same(&st, &other);
                }
                let (qa, qb) = (
                    strand_codes(&st.codes, &st.seq),
                    strand_codes(&other.codes, &other.seq),
                );
                // Costs of the two best hits and the score charged per unit.
                let (ca, cb, unit) = match &tb.affine {
                    // Two banded alignments of one sequence differ only by their heuristics.
                    _ if mirror => {
                        stats::add(&stats::SAME_ANTISENSE, 1);
                        (0, 0, 0.0)
                    }
                    Some(sc) if tb.orient_clipped => (
                        -i64::from(self.clipped(tb, sc, qa, st.raw.as_ref(), &st.hits[0])),
                        -i64::from(self.clipped(tb, sc, qb, other.raw.as_ref(), &other.hits[0])),
                        tb.lambda,
                    ),
                    _ => {
                        // Both whole reads: trimming each to its own hit's anchors would score
                        // the two orientations over different read spans.
                        let da = self.whole_read_k(tb, qa, &st.hits[0], u32::MAX);
                        let db = self.whole_read_k(tb, qb, &other.hits[0], da.saturating_add(tb.cap));
                        (i64::from(da), i64::from(db), tb.kappa)
                    }
                };
                // This orientation's lead over the other, in score points: alignment, then the
                // poly(A) tail where one orientation shows it and the other doesn't.
                let mut lead = unit * (cb - ca) as f32;
                if tb.orient_tail > 0 {
                    // Without HPC the parsed read is the raw read.
                    let tailed = |s: &Strand| match &s.raw {
                        Some(r) => ends_in_tail(&r.codes, tb.tail_run),
                        None => !self.seeder.hpc && ends_in_tail(strand_codes(&s.codes, &s.seq), tb.tail_run),
                    };
                    lead += tb.orient_tail as f32
                        * (f32::from(u8::from(tailed(&st))) - f32::from(u8::from(tailed(&other))));
                }
                // A close call is decided by exact local alignments instead of the clipped
                // ones (`TieBreak::orient_exact`); the tail term stands.
                if let Some(sc) = &tb.affine
                    && !mirror
                    && tb.orient_clipped
                    && tb.orient_exact > 0.0
                    && lead.abs() <= tb.orient_exact
                {
                    let w = ORIENT_EXACT_BAND + plen / 20;
                    let tail = lead - unit * (cb - ca) as f32;
                    let xa = self.local_exact(tb, sc, qa, st.raw.as_ref(), &st.hits[0], w);
                    let xb = self.local_exact(tb, sc, qb, other.raw.as_ref(), &other.hits[0], w);
                    lead = unit * (xa - xb) as f32 + tail;
                }
                if lead < 0.0 {
                    std::mem::swap(&mut st, &mut other);
                    reverse = !reverse;
                }
                let target = st.hits[0].score - lead.abs();
                let shift = target - other.hits[0].score;
                for h in &mut other.hits {
                    h.score += shift;
                    h.reverse = !reverse;
                }
            }
        }
        for h in &mut st.hits {
            h.reverse = reverse;
        }
        let mut hits = st.hits;
        if both {
            hits.append(&mut other.hits);
            sort_hits(&mut hits, salt);
        }
        // A weak rescued best hit must also align (`MinScore::verify`).
        let mut unverified = false;
        if rescued
            && let Some(ms) = self.p.min_score
            && let Some(v) = ms.verify
            && let (Some(tb), Some(top)) = (&self.p.tiebreak, hits.first())
            && let Some(sc) = &tb.affine
        {
            let span = (top.qend - top.qstart) as f32;
            if span < self.p.pi_cov * plen as f32 || top.score < ms.at(raw.len()) + v.margin {
                let (c, r) = if top.reverse == reverse {
                    (strand_codes(&st.codes, &st.seq), st.raw.as_ref())
                } else {
                    (strand_codes(&other.codes, &other.seq), other.raw.as_ref())
                };
                unverified = (self.clipped(tb, sc, c, r, top) as f32) < v.base + v.per_base * span;
            }
        }
        let mut low_coverage = true;
        // Tie-broken orientations can leave the lower-scoring one on top.
        let low_score = low(&hits, rescued_any) || unverified;
        if reject(low_score) {
            hits.clear();
        }
        if let Some(top) = hits.first() {
            let floor = top.score - (self.p.delta_rep * top.score.abs()).max(self.p.rep_abs);
            low_coverage = ((top.qend - top.qstart) as f32) < self.p.pi_cov * plen as f32;
            hits.retain(|h| h.score >= floor);
        }
        if let Some(e) = &self.p.ends {
            let mut seen = rustc_hash::FxHashMap::default();
            let floor = hits.first().map_or(0.0, |t| t.score - e.within * t.score.abs());
            for h in hits.iter_mut().filter(|h| h.score >= floor) {
                // The other orientation's hits are reported only if both were tie-broken.
                let c = if h.reverse == reverse {
                    strand_codes(&st.codes, &st.seq)
                } else {
                    strand_codes(&other.codes, &other.seq)
                };
                self.extend_ends(e, c, h, &mut seen);
                h.extended = true;
            }
        }
        ReadResult {
            len: raw.len() as u32,
            plen,
            reverse,
            hits,
            low_coverage,
            rescued,
            low_score,
            // `st` is the orientation `reverse` names (swapped together), `other` the opposite.
            hpc: if reverse { other.raw } else { st.raw }.map(|r| r.hpc),
        }
    }
}

/// Diagonals searched on each side of a chain's by the exact comparison of near-tied read
/// orientations ([`TieBreak::orient_exact`]), plus a twentieth of the parsed read.
pub const ORIENT_EXACT_BAND: u32 = 100;
/// Bases between a transcript's 3' end and a read's tail start treated as the same end.
const TAIL_TOL: usize = 15;
/// Most read bases after a poly(A) tail (adapter and the like).
const TAIL_AFTER: usize = 80;
/// Raw read bases searched for a poly(A) tail when comparing read orientations.
const TAIL_WINDOW: usize = 200;

/// Whether a raw read (2-bit codes) ends in a poly(A) tail: a `run`-base window of A's within its
/// last `TAIL_WINDOW` bases, followed by at most `TAIL_AFTER` bases.
fn ends_in_tail(codes: &[u8], run: usize) -> bool {
    let a = lcp_core::encode::CODE[b'A' as usize] & 3;
    tail_start(&codes[codes.len().saturating_sub(TAIL_WINDOW)..], a, run).is_some()
}

/// The first A-run of `text` (raw codes) starting in `[lo, hi]`: a window of `run` bases with at
/// most one other base. Returns its first A and the end of the A-run from there (as
/// [`tail_start`]).
fn own_run(text: &[u8], a: u8, run: usize, lo: usize, hi: usize) -> Option<(usize, usize)> {
    let i = (lo..=hi.min(text.len().checked_sub(run)?))
        .find(|&i| text[i..i + run].iter().filter(|&&c| c != a).count() <= 1)?;
    let s = i + text[i..].iter().position(|&c| c == a)?;
    let (mut end, mut non) = (s, 0);
    for (k, &c) in text.iter().enumerate().skip(s) {
        non += usize::from(c != a);
        if c == a && non * 10 <= k + 1 - s {
            end = k + 1;
        }
    }
    Some((s, end))
}

/// The first poly(A) tail in `rest` (raw codes): a window of `run` bases with at most one other
/// base whose A-run is followed by at most `TAIL_AFTER` bases. Returns the window's first A and
/// the end of the A-run from there (the longest stretch ending in A with at most one other base
/// per ten).
fn tail_start(rest: &[u8], a: u8, run: usize) -> Option<(usize, usize)> {
    let others = |w: &[u8]| w.iter().filter(|&&c| c != a).count();
    let (mut i, mut other) = (0, others(rest.get(..run)?));
    loop {
        if other <= 1 {
            let s = i + rest[i..].iter().position(|&c| c == a).unwrap_or(0);
            let (mut end, mut non) = (s, 0);
            for (k, &c) in rest.iter().enumerate().skip(s) {
                non += usize::from(c != a);
                if c == a && non * 10 <= k + 1 - s {
                    end = k + 1;
                }
            }
            if rest.len() - end <= TAIL_AFTER {
                return Some((s, end));
            }
            // A tail is followed only by adapter: an A-run with more of the read after it is
            // transcript sequence (or internal priming), not the molecule's 3' end. The tail
            // may still come after it; the windows inside the run are part of it.
            i = end.max(i + 1);
            other = others(rest.get(i..i + run)?);
        } else {
            if i + run >= rest.len() {
                return None;
            }
            other = other + usize::from(rest[i + run] != a) - usize::from(rest[i] != a);
            i += 1;
        }
    }
}

/// Best score first; equal scores in the order of `tie_key(salt, tx)`, a per-read
/// pseudo-random permutation of the transcripts.
fn sort_hits(hits: &mut [TxHit], salt: u64) {
    hits.sort_unstable_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then_with(|| tie_key(salt, a.tx).cmp(&tie_key(salt, b.tx)))
            .then(a.tx.cmp(&b.tx))
    });
}

/// [`sort_hits`] for one orientation's chains, whose transcripts are distinct (a chainer gives
/// one hit per transcript of a candidate, and candidates are families, which share no
/// transcript): each hit's key is computed once and each hit moved once. Distinct transcripts
/// make the order total, so it is `sort_hits`' order.
fn sort_distinct_hits(hits: &mut [TxHit], salt: u64) {
    // `f32::total_cmp` order as an unsigned key, reversed (best score first).
    let desc = |s: f32| {
        let b = s.to_bits() as i32;
        !(((b ^ ((((b >> 31) as u32) >> 1) as i32)) as u32) ^ 0x8000_0000)
    };
    hits.sort_by_cached_key(|h| (desc(h.score), tie_key(salt, h.tx), h.tx));
    // Hits `sort_hits` can't order (same score and transcript) would now be adjacent.
    debug_assert!(
        hits.windows(2)
            .all(|w| (w[0].score.to_bits(), w[0].tx) != (w[1].score.to_bits(), w[1].tx))
    );
}

/// splitmix64 of the read salt and a transcript id.
fn tie_key(salt: u64, tx: u32) -> u64 {
    let mut z = salt ^ u64::from(tx).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod banded_tests {
    use super::*;
    use proptest::prelude::*;

    /// The sequential definition of [`Mapper::banded`]: one `distance_with` per distinct text,
    /// in order, at the threshold `k` = best so far + `cap`; a repeated text takes its first
    /// occurrence's distance re-capped at the current `k`.
    fn sequential<T: AsRef<[u8]>>(read: &[u8], texts: &[T], mode: align::Mode, cap: u32) -> Vec<u32> {
        let pat = align::Pattern::new(read);
        let mut seen: rustc_hash::FxHashMap<&[u8], u32> = rustc_hash::FxHashMap::default();
        let mut best = u32::MAX;
        texts
            .iter()
            .map(|t| {
                let t = t.as_ref();
                let k = best.saturating_add(cap);
                if let Some(&d) = seen.get(t) {
                    return if d <= k { d } else { k.saturating_add(1) };
                }
                let d = align::distance_with(&pat, t, mode, k).unwrap_or(k.saturating_add(1));
                best = best.min(d);
                seen.insert(t, d);
                d
            })
            .collect()
    }

    /// `s` with substitutions (kind 0), insertions (1) and deletions (2) at positions taken
    /// modulo its length.
    fn mutate(s: &[u8], edits: &[(usize, u8, u8)]) -> Vec<u8> {
        let mut v = s.to_vec();
        for &(pos, b, kind) in edits {
            if v.is_empty() {
                if kind == 1 {
                    v.push(b);
                }
                continue;
            }
            let i = pos % v.len();
            match kind {
                0 => v[i] = b,
                1 => v.insert(i, b),
                _ => {
                    v.remove(i);
                }
            }
        }
        v
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]

        /// `banded` (distinct texts in batches, re-capped in order) equals the sequential
        /// definition in every mode at every SIMD level, and so does `distinct` on the texts
        /// without their repeats: texts are mutated excerpts of the read's source (as isoform
        /// regions are), some empty or unrelated, with repeats, and caps range from 0 to near
        /// `u32::MAX` (where `k` saturates).
        #[test]
        fn banded_matches_sequential(
            base in prop::collection::vec(0u8..4, 0..300),
            read_edits in prop::collection::vec((0usize..300, 0u8..4, 0u8..3), 0..30),
            variants in prop::collection::vec(
                (
                    prop::collection::vec((0usize..300, 0u8..4, 0u8..3), 0..40),
                    0usize..40,
                    0usize..40,
                    any::<bool>(),
                ),
                1..24,
            ),
            repeats in prop::collection::vec((0usize..64, 0usize..64), 0..8),
            cap in prop_oneof![0u32..6, 6u32..40, (u32::MAX - 3)..=u32::MAX],
        ) {
            let read = mutate(&base, &read_edits);
            let mut texts: Vec<Vec<u8>> = variants
                .iter()
                .map(|(edits, a, b, unrelated)| {
                    if *unrelated {
                        edits.iter().map(|e| e.1).collect()
                    } else {
                        let t = mutate(&base, edits);
                        let lo = (*a).min(t.len());
                        let hi = t.len().saturating_sub(*b).max(lo);
                        t[lo..hi].to_vec()
                    }
                })
                .collect();
            for (from, to) in repeats {
                let t = texts[from % texts.len()].clone();
                texts.insert(to % (texts.len() + 1), t);
            }
            // The texts without repeats (first occurrences, in order), for `distinct`.
            let mut once = rustc_hash::FxHashSet::default();
            let uniq: Vec<&[u8]> = texts.iter().map(|t| t.as_slice()).filter(|t| once.insert(*t)).collect();
            for mode in [align::Mode::Global, align::Mode::Prefix, align::Mode::Infix] {
                let want = sequential(&read, &texts, mode, cap);
                let want_uniq = sequential(&read, &uniq, mode, cap);
                for level in lcp_core::simd::supported() {
                    let got = align::at_level(level, || {
                        Mapper::<OracleChainer>::banded(&read, &texts, mode, cap, false)
                    });
                    prop_assert_eq!(&got, &want, "{:?} at {}", mode, level);
                    let got = align::at_level(level, || {
                        Mapper::<OracleChainer>::banded(&read, &uniq, mode, cap, true)
                    });
                    prop_assert_eq!(&got, &want_uniq, "{:?} at {} (distinct)", mode, level);
                }
            }
        }
    }
}

#[cfg(test)]
mod segmented_tests {
    use super::*;

    fn hit(tx: u32, score: f32, chain: &[(u32, u32, u32)]) -> TxHit {
        chain::hit_from_path(tx, score, 0.0, chain.iter().map(|&(q, t, l)| (q, q + l, t, l)))
    }

    fn mapper(texts: &[&[u8]], hpc: bool) -> Mapper<OracleChainer> {
        let txs: Vec<_> = texts
            .iter()
            .enumerate()
            .map(|(i, seq)| index::TxInput {
                name: format!("t{i}"),
                gene: "g".into(),
                seq: seq.to_vec(),
            })
            .collect();
        let (data, _) = index::build(
            txs,
            &index::BuildParams {
                spec: "lcp:K=3".into(),
                hpc,
                f_max: 100,
                allow_collisions: false,
                collapse: index::Collapse::Gene,
                infer_genes: None,
            },
        )
        .unwrap();
        Mapper::new(
            Arc::new(Index::from_data(&data).unwrap()),
            OracleChainer {
                p: Default::default(),
            },
            MapParams::default(),
        )
        .unwrap()
    }

    #[test]
    fn complete_intervals_and_overlapping_blocks() {
        assert!(shared_blocks(&[hit(0, 1., &[(2, 2, 3)]), hit(1, 1., &[(2, 2, 4)])]).is_none());
        let same_start = [hit(0, 1., &[(2, 4, 3), (2, 8, 4)]), hit(1, 1., &[(2, 9, 4)])];
        assert_eq!(
            shared_blocks(&same_start).unwrap(),
            vec![SharedBlock {
                query: 2..6,
                targets: vec![8..12, 9..13],
            }]
        );
        let hits = [
            hit(0, 1., &[(2, 4, 5), (4, 6, 2), (6, 8, 5), (11, 15, 2)]),
            hit(1, 1., &[(2, 5, 5), (4, 7, 2), (6, 9, 5), (11, 16, 2)]),
        ];
        assert_eq!(
            shared_blocks(&hits).unwrap(),
            vec![
                SharedBlock {
                    query: 2..11,
                    targets: vec![4..13, 5..14]
                },
                SharedBlock {
                    query: 11..13,
                    targets: vec![15..17, 16..18]
                },
            ]
        );
        // An overlapping anchor on another diagonal, or one whose target runs backwards, is
        // dropped; the block before it stays.
        let one = |q: std::ops::Range<u32>, t: std::ops::Range<u32>| SharedBlock {
            query: q,
            targets: vec![t],
        };
        assert_eq!(
            shared_blocks(&[hit(0, 1., &[(2, 2, 5), (4, 5, 5)])]).unwrap(),
            vec![one(2..7, 2..7)]
        );
        assert_eq!(
            shared_blocks(&[hit(0, 1., &[(2, 2, 3), (6, 4, 3)])]).unwrap(),
            vec![one(2..5, 2..5)]
        );
    }

    #[test]
    fn ordered_intersection_skips_unshared_groups_and_handles_unsorted_lengths() {
        let hits = [
            hit(
                0,
                1.,
                &[(0, 0, 1), (2, 2, 4), (2, 2, 3), (2, 2, 3), (8, 8, 2), (12, 12, 2)],
            ),
            hit(1, 1., &[(1, 2, 1), (2, 3, 3), (2, 3, 4), (7, 8, 1), (12, 15, 2)]),
            hit(2, 1., &[(2, 4, 4), (2, 4, 3), (10, 12, 1), (12, 16, 2)]),
        ];
        assert_eq!(
            shared_blocks(&hits).unwrap(),
            vec![
                SharedBlock {
                    query: 2..6,
                    targets: vec![2..6, 3..7, 4..8]
                },
                SharedBlock {
                    query: 12..14,
                    targets: vec![12..14, 15..17, 16..18]
                },
            ]
        );
        let mut empty = hits[1].clone();
        empty.chain.clear();
        assert!(shared_blocks(&[hits[0].clone(), empty]).is_none());
    }

    #[test]
    fn paper_example_segments_and_per_transcript_scores() {
        let q = b"ACGTATGCACGATCTCAGTA";
        let texts: [&[u8]; 4] = [
            q,
            b"GACGTATGCACGATCTCAGT",
            b"ACGTATGCACATGATCTCAGT",
            b"ACGTACGTGCACGATCTCAGT",
        ];
        let m = mapper(&texts, false);
        let mut hits = vec![
            hit(0, 15., &[(7, 7, 3), (10, 10, 8)]),
            hit(1, 15., &[(7, 8, 3), (10, 11, 8)]),
            hit(2, 12., &[(7, 7, 3), (10, 12, 8)]),
            hit(3, 11., &[(7, 9, 3), (10, 12, 8)]),
        ];
        let b = shared_blocks(&hits).unwrap();
        assert_eq!(b[0].query.end..b[1].query.start, 10..10);
        let gaps: Vec<_> = hits
            .iter()
            .enumerate()
            .map(|(i, h)| m.idx.tx_codes(h.tx, b[0].targets[i].end, b[1].targets[i].start))
            .collect();
        assert_eq!(gaps, vec![vec![], vec![], codes(b"AT"), vec![]]);
        let sc = align::Affine::ONT;
        // Per-segment alignments of the example, without the flank rule.
        let tb = TieBreak {
            flank_min: 0,
            ..Default::default()
        };
        let mut totals = vec![0; 4];
        for (i, h) in hits.iter().enumerate() {
            let left = m.idx.tx_codes(h.tx, 0, b[0].targets[i].start);
            let left = align::affine_extend(
                &codes(&q[..7]).into_iter().rev().collect::<Vec<_>>(),
                &left.into_iter().rev().collect::<Vec<_>>(),
                &sc,
                tb.band,
                tb.end_bonus,
                tb.xdrop,
            );
            let interior = align::affine(&[], &gaps[i], align::Mode::Global, &sc, tb.band);
            let right_text = m.idx.tx_codes(h.tx, b[1].targets[i].end, m.idx.tx_plen()[i]);
            let right = align::affine_extend(
                &codes(&q[18..]),
                &right_text,
                &sc,
                tb.band,
                tb.end_bonus,
                tb.xdrop,
            );
            assert_eq!(
                (left, interior, right),
                [(19, 0, 9), (19, 0, 2), (19, -8, 2), (11, 0, 2)][i]
            );
            totals[i] = left + interior + right;
        }
        assert_eq!(totals, [28, 21, 13, 13]);
        assert_eq!(
            m.segmented(&tb, &codes(q), None, &hits).unwrap().0,
            [-28, -21, -13, -13]
        );
        let c = codes(q);
        assert_eq!(m.tiebreak(&tb, || &c, None, &mut hits, 0), 4);
        assert_eq!(
            hits.iter().map(|h| h.score).collect::<Vec<_>>(),
            [15., 8., 0., 0.]
        );
    }

    #[test]
    fn tail_start_finds_the_first_poly_a_window() {
        let a = codes(b"A")[0];
        assert_eq!(tail_start(&codes(b"CGTCAAAAAAAAAAAAAAAGT"), a, 15), Some((4, 19)));
        // One base in a window may be a sequencing error.
        assert_eq!(tail_start(&codes(b"CGAAAAAAACAAAAAAAGT"), a, 15), Some((2, 17)));
        assert_eq!(tail_start(&codes(b"CGTCAAAAAAAAAAGTCAAAA"), a, 15), None);
    }

    #[test]
    fn an_a_run_inside_the_read_does_not_hide_its_tail() {
        // A 3' UTR's A-run, more of the UTR, then the tail and an adapter.
        let utr = b"CGTCGATCGTACGATCGTCAGCTAGCTGACTGCATGCATG".repeat(3);
        let adapter: &[u8] = b"GTACTCTGCGTTGATACCACTGCTT";
        let read = [&b"A".repeat(16)[..], &utr, &b"A".repeat(30), adapter].concat();
        let a = codes(b"A")[0];
        assert_eq!(tail_start(&codes(&read), a, 15).map(|t| t.0), Some(136));
        let body = b"ACGTCGATCGTACGATCGTCAGCTAGCTGACTGCATGCAT";
        assert!(ends_in_tail(&codes(&[&body[..], &read].concat()), 15));
    }

    #[test]
    fn only_a_poly_a_near_the_read_end_marks_its_orientation() {
        let body = b"ACGTCGATCGTACGATCGTCAGCTAGCTGACTGCATGCAT".repeat(10);
        let tail = [b"A".repeat(20), b"GTACTCTGCGTTGATACCACTGCTT".to_vec()].concat();
        assert!(ends_in_tail(&codes(&[&body[..], &tail].concat()), 15));
        // A tail at the read's start is the other orientation's 3' end.
        assert!(!ends_in_tail(&codes(&[&tail[..], &body].concat()), 15));
        // An A-run with more than `TAIL_AFTER` read bases after it is transcript sequence.
        assert!(!ends_in_tail(
            &codes(&[&body[..], &tail, &body[..100]].concat()),
            15
        ));
    }

    #[test]
    fn poly_a_tail_is_trimmed_and_favours_the_isoform_ending_there() {
        // t1 continues past t0's 3' end; the read stops there and goes on into its tail.
        let body: &[u8] = b"ACGTCGATCGTACGATCGTCAGCTAGCTGACTGCATGCAT";
        let t1: Vec<u8> = [body, b"GCTCGTCGTCAGCTCGTCGATCGCT"].concat();
        let q: Vec<u8> = [body, b"AAAAAAAAAAAAAAAAAAAAAAAA", b"GTCGCTCGTGCA"].concat();
        // Without HPC there is no separate raw read: the tail is found in the parsed one.
        for hpc in [true, false] {
            let m = mapper(&[body, &t1], hpc);
            let (seq, raw) = Mapper::<OracleChainer>::parse(&m.seeder, &q);
            assert_eq!(raw.is_some(), hpc);
            let n = 30;
            let hits = [hit(0, 10., &[(0, 0, n)]), hit(1, 10., &[(0, 0, n)])];
            let cost = |tail_bonus| {
                let tb = TieBreak {
                    tail_bonus,
                    ..Default::default()
                };
                m.segmented(&tb, &codes(&seq), raw.as_ref(), &hits).unwrap().0
            };
            let (off, on) = (cost(0), cost(10));
            assert_eq!(off[0], off[1]);
            assert_eq!(on[0], off[0] - 10);
            assert_eq!(on[1], off[1]);
        }
    }

    #[test]
    fn chance_matches_in_a_longer_flank_do_not_separate_isoforms() {
        // t1 has 30 bases of 5' flank that t0 lacks; a read starting with a primer aligns a
        // few bases into it by chance, a read from t1 aligns all of it.
        let body: &[u8] = b"ACGTCGATCGTACGATCGTCAGCTAGCTGACTGCATGCAT";
        let flank: &[u8] = b"GTCAGTACGTAGCTGATCGATGCTAGTCGA";
        let t1: Vec<u8> = [flank, body].concat();
        let m = mapper(&[body, &t1], false);
        let primer: &[u8] = b"TTTCTGTTGGTGCTGATATTGCTGGG";
        let cost = |read: &[u8], off: u32, flank_min| {
            let hits = [hit(0, 10., &[(off, 0, 30)]), hit(1, 10., &[(off, 30, 30)])];
            let tb = TieBreak {
                flank_min,
                ..Default::default()
            };
            m.segmented(&tb, &codes(read), None, &hits).unwrap().0
        };
        let primed: Vec<u8> = [primer, body].concat();
        let c = cost(&primed, 26, 10);
        assert_eq!(c[0], c[1]);
        let from_t1: Vec<u8> = [flank, body].concat();
        let c = cost(&from_t1, 30, 10);
        assert!(c[1] < c[0]);
        // Ends two bases apart: the short exact extension still counts.
        let m = mapper(&[body, &[b"GT", body].concat()], false);
        let hits = [hit(0, 10., &[(2, 0, 30)]), hit(1, 10., &[(2, 2, 30)])];
        let tb = TieBreak {
            flank_min: 10,
            ..Default::default()
        };
        let read: Vec<u8> = [b"GT", body].concat();
        let c = m.segmented(&tb, &codes(&read), None, &hits).unwrap().0;
        assert!(c[1] < c[0]);
    }

    #[test]
    fn score_scaling_and_no_common_anchor_fallback() {
        let m = mapper(&[b"ACGTAC", b"ACAC"], false);
        let q = b"ACGTAC";
        let c = codes(q);
        let mut hits = [
            hit(0, 15., &[(0, 0, 2), (4, 4, 2)]),
            hit(1, 12., &[(0, 0, 2), (4, 2, 2)]),
        ];
        let tb = TieBreak {
            lambda: 0.5,
            ..Default::default()
        };
        m.tiebreak(&tb, || &c, None, &mut hits, 0);
        // Scores become the tie set's best chain score minus each deficit.
        assert_eq!((hits[0].score, hits[1].score), (15., 9.));
        let m = mapper(&[q, q], false);
        for grouped in [false, true] {
            let mut hits = [hit(0, 15., &[(0, 0, 2)]), hit(1, 12., &[(3, 3, 2)])];
            let tb = TieBreak {
                group_fallback: grouped,
                ..Default::default()
            };
            assert!(m.segmented(&tb, &c, None, &hits).is_none());
            m.tiebreak(&tb, || &c, None, &mut hits, 0);
            assert_eq!((hits[0].score, hits[1].score), (15., 15.));
        }
    }

    #[test]
    fn empty_and_deleted_interior_gaps() {
        let m = mapper(&[b"ACGTAC", b"ACAC"], false);
        let hits = [
            hit(0, 4., &[(0, 0, 2), (4, 4, 2)]),
            hit(1, 4., &[(0, 0, 2), (4, 2, 2)]),
        ];
        let tb = TieBreak::default();
        assert_eq!(
            m.segmented(&tb, &codes(b"ACGTAC"), None, &hits).unwrap().0,
            [-4, 8]
        );
        // Identical gaps contribute no relative cost, including empty query/target gaps.
        assert_eq!(
            m.segmented(&tb, &codes(b"ACGTAC"), None, &[hits[0].clone(), hits[0].clone()])
                .unwrap()
                .0,
            [0, 0]
        );
    }

    #[test]
    fn raw_hpc_end_boundaries() {
        let q = b"AACCCGGTTTAA";
        let raw = RawRead {
            codes: codes(q),
            hpc: lcp_core::hpc::Hpc::new(q),
        };
        let m = mapper(&[q, b"ACGTA"], true);
        let hits = [hit(0, 3., &[(1, 1, 2)]), hit(1, 3., &[(1, 1, 2)])];
        let blocks = shared_blocks(&hits).unwrap();
        assert_eq!(
            (raw.pos(blocks[0].query.start), raw.pos(blocks[0].query.end)),
            (2, 7)
        );
        let mut tb = TieBreak {
            flank_min: 0,
            ..Default::default()
        };
        let parsed = codes(&raw.hpc.seq);
        let costs = m.segmented(&tb, &parsed, Some(&raw), &hits).unwrap().0;
        // The block's raw bases differ (CCCGG vs CG), so the right end runs through it.
        // Raw left AA vs A: 9 vs 2; raw right CCCGGTTTAA vs CGTA: 25 vs 2.
        assert_eq!(costs, [-34, -4]);
        tb.raw_ends = false;
        assert_eq!(m.segmented(&tb, &parsed, Some(&raw), &hits).unwrap().0, [0, 0]);
    }

    proptest::proptest! {
        /// Sorting one orientation's hits (distinct transcripts) by cached keys gives
        /// `sort_hits`' order, with tied, signed-zero, infinite and NaN scores.
        #[test]
        fn distinct_hit_sort_is_sort_hits(
            scores in proptest::collection::vec(proptest::sample::select(vec![
                0.0, -0.0, 1.5, -3.0, 100.0, 7.25, f32::INFINITY, f32::NEG_INFINITY, f32::NAN,
                -f32::NAN,
            ]), 0..80),
            salt: u64,
            mul in 0u32..u32::MAX,
        ) {
            // Distinct transcripts: i -> i * odd is a bijection.
            let hits: Vec<TxHit> = scores
                .iter()
                .enumerate()
                .map(|(i, &s)| hit((i as u32).wrapping_mul(mul | 1), s, &[(0, 0, 1)]))
                .collect();
            let key = |h: &[TxHit]| h.iter().map(|h| (h.tx, h.score.to_bits())).collect::<Vec<_>>();
            let (mut a, mut b) = (hits.clone(), hits);
            sort_hits(&mut a, salt);
            sort_distinct_hits(&mut b, salt);
            proptest::prop_assert_eq!(key(&a), key(&b));
        }
    }
}

#[cfg(test)]
mod min_score_tests {
    use super::*;

    fn hit(score: f32, anchors: usize) -> TxHit {
        let path = (0..anchors as u32).map(|i| (10 * i, 10 * i + 5, 10 * i, 5));
        chain::hit_from_path(0, score, 0.0, path)
    }

    #[test]
    fn threshold_grows_with_the_log_of_read_length() {
        let m = MinScore::default();
        assert_eq!(m.at(1000), m.at_1kb);
        assert!((m.at(2718) - (m.at_1kb + m.per_ln_len)).abs() < 1e-2);
        assert!(m.at(100) < m.at(1000) && m.at(1000) < m.at(10_000));
        // A zero-length read has a finite threshold.
        assert!(m.at(0).is_finite());
    }

    #[test]
    fn the_default_judges_by_score_alone() {
        let m = MinScore::default();
        let t = m.at(1000);
        assert!(m.passes(&hit(t, 1), 1000, false) && m.passes(&hit(t, 1), 1000, true));
        assert!(!m.passes(&hit(t - 1.0, 1), 1000, false));
        assert!(!m.passes(&hit(t - 1.0, 5), 1000, false));
        assert!(m.passes(&hit(t, 5), 1000, false));
    }

    #[test]
    fn several_base_level_anchors_pass_at_any_score() {
        let m = MinScore {
            anchors: 2,
            ..Default::default()
        };
        let t = m.at(1000);
        assert!(!m.passes(&hit(t - 1.0, 1), 1000, false));
        // Two anchors pass under the threshold, but not when they are rescue seeds.
        assert_eq!(hit(1.0, 2).anchors, 2);
        assert!(m.passes(&hit(t - 1.0, 2), 1000, false));
        assert!(!m.passes(&hit(t - 1.0, 2), 1000, true));
    }
}
