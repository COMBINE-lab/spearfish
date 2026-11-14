//! Read mapping: seeding both strands, anchor lookup, candidate-gene selection, chaining, and
//! result extraction.

pub mod align;
pub mod anchor;
pub mod chain;
pub mod mates;
pub mod refine;
pub mod rescue;
pub mod stats;
pub mod tables;

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
    /// With `orient_tail`: most read bases after a poly(A) tail that starts where the best
    /// hit's transcript ends, for the orientation comparison. A cDNA read that kept its adapters
    /// carries the RT primer and the sequencing adapter after its tail, more than the
    /// [`TAIL_AFTER`] bases any A-run may have after it (`TAIL_AFTER` or less: that alone).
    #[serde(default = "tail_after")]
    pub orient_tail_after: u32,
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
    /// With `orient_same`: such a pair is compared by exact local alignments of the whole read
    /// (as [`TieBreak::orient_exact`] does) instead of being reported as tied. Under the chains
    /// the two transcripts are one sequence, which scores the same against both read
    /// orientations whatever the chains; the scores differ only by what the read holds beyond
    /// that sequence: bases of one transcript's own (a short antisense isoform's first bases)
    /// and how far each transcript reaches into the read's ends. Where the two are one
    /// sequence under the whole read as well, both alignments score the same and are skipped.
    /// The pair stands at the better of the two chain scores, whichever orientation leads.
    #[serde(default)]
    pub orient_same_exact: bool,
    /// With `orient_same_exact`: a read end earns `end_bonus` in those alignments only where
    /// both transcripts have sequence up to it (as `end_bonus_cut` within an orientation).
    /// Where a read ends at one transcript's end and the other goes on, the other could reach
    /// the read's end only through a junk base or two aligned into its continuation, which is
    /// no evidence for it.
    #[serde(default)]
    pub orient_same_end_cut: bool,
    /// Apply `orient_same_exact` and `orient_tail_after` to the read orientation through the
    /// quantification score only. The mapping score, the orientation, the reported hits and
    /// their ranks are then decided as without the two; each hit's quantification score is
    /// moved by what the two would have moved its mapping score. A quantifier weighs the read
    /// by that evidence, while every mapping decision stays as it was.
    #[serde(default)]
    pub orient_quant: bool,
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
    /// With `raw_mid`, a shared block longer than twice this many parsed bases is cut this
    /// far inside each of its edges instead of at its midpoint, and its bases between the two
    /// cuts are compared on their own (0: every block is cut at its midpoint). Candidates
    /// differ next to a block's edges, not deep inside it, and a short stretch that is
    /// identical in most of them is compared once for all.
    #[serde(default)]
    pub mid_inset: u32,
    /// The clipped alignment score (orientation comparison, anchor groups) follows the chain
    /// piece by piece and keeps the best run of pieces (off: one banded alignment of the
    /// whole chained span, as in 0.1.0, which a long noisy chain's indels drift out of).
    #[serde(default)]
    pub clip_pieces: bool,
    /// A candidate with its own A-run where the read's poly(A) tail starts explains the tail
    /// as well as one that ends there: it gains `tail_bonus` too, and what the read's bases
    /// after the run gain against its continuation.
    #[serde(default)]
    pub tail_templated: bool,
    /// Read ends: where the candidates' windows include one stopping at its transcript's end
    /// and longer ones, `end_bonus` counts only for reaching the read's end within the
    /// shortest such window. Past a transcript's end that candidate has no bases, so reaching
    /// the read's end there (through a few junk or adapter bases aligned into a longer
    /// candidate's flank) is no evidence for the longer one.
    #[serde(default)]
    pub end_bonus_cut: bool,
    /// Quantification score only: points added to the hits whose annotated 3' end the read's
    /// right end reaches, for a read without a poly(A) tail in either orientation (0: off).
    /// The mapping score, the best set and the orientation are unchanged; the adjustment is
    /// reported next to the score ([`ReadResult::end_fit`]) for a quantifier's weights. See
    /// `Mapper::end_fit`.
    #[serde(default)]
    pub end_fit: f32,
    /// With `end_fit`: the measurement error of where a read ends on a transcript, in parsed
    /// bases. A hit fits its 3' end if its right-end alignment stops within this many bases of
    /// it and reaches as far into the read as the other candidates' do, within as many;
    /// candidates whose annotated ends lie within it of each other get no relative adjustment.
    #[serde(default)]
    pub end_fit_tol: u32,
    /// With `end_fit`: measure where the read ends on each candidate from its chain instead of
    /// aligning the read's right flank: the read is taken to go on along the transcript from
    /// the last anchor, base for base, to the read's end or the transcript's, whichever comes
    /// first. No alignment is run for the end fit, which presets that do not extend spans would
    /// otherwise run for it alone. Candidates whose spans were extended are measured from
    /// their extension's alignment either way.
    #[serde(default)]
    pub end_fit_chain: bool,
}

/// Default [`TieBreak::end_fit`] in score points (2 points weigh a factor e at a score
/// denominator of 2).
pub const END_FIT: f32 = 5.0;
/// Default [`TieBreak::end_fit_tol`] in parsed bases.
pub const END_FIT_TOL: u32 = 4;

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
            orient_tail_after: 120,
            orient_margin_rel: 0.35,
            orient_same: true,
            orient_same_exact: true,
            orient_same_end_cut: true,
            orient_quant: true,
            orient_exact: 150.0,
            raw_mid: true,
            mid_inset: 8,
            clip_pieces: true,
            tail_templated: true,
            flank_min: 0,
            end_bonus_cut: true,
            end_fit: END_FIT,
            end_fit_tol: END_FIT_TOL,
            end_fit_chain: false,
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
    /// With both `strands`: map the two orientations in turn where the read's poly(A) tail
    /// names its strand (None: always both).
    pub tail_first: Option<TailFirst>,
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

/// Staged mapping of a read's two orientations. A read whose raw sequence ends in a poly(A)
/// tail in exactly one orientation (the tail itself, or a leading poly(T), its reverse
/// complement) is mapped in that orientation first, and in the other only if the first's best
/// chain spans less than `span` of the read. The tail names the molecule's strand, so the other
/// orientation's seeding and chaining, half the work of such a read, is rarely needed.
///
/// A read mapped in one orientation reports no hits in the other: where transcripts share
/// sequence antisense, the hits a near-tied opposite orientation would have added are the ones
/// the tail speaks against. Reads without a tail, or with one in both orientations, are mapped
/// as without staging, so a library without tails maps as before (bar the few reads whose
/// sequence happens to end in such a run and that also map well in that orientation).
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TailFirst {
    /// The tailed orientation stands alone if its best chain spans at least this fraction of
    /// the (parsed) read.
    pub span: f32,
    /// Shortest poly(A) tail that orders the orientations (a window as
    /// [`TieBreak::tail_run`]'s). Longer than the tie-break's 15: there a tail adds points to an
    /// alignment, here it decides what is aligned, and 20 leaves out more than half of the A-rich
    /// ends of tail-less reads while keeping all but 0.1% of the tails that 15 finds.
    pub run: usize,
}

impl Default for TailFirst {
    fn default() -> Self {
        TailFirst { span: 0.8, run: 20 }
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
type EndCache = rustc_hash::FxHashMap<(bool, bool, u32), Vec<(Flank, bool, (u32, u32))>>;

/// The transcript flank prefix an end extension read: decoded (left ends, read backwards), or
/// where it lies in a transcript (right ends: transcript, start, length), compared packed.
enum Flank {
    Codes(Vec<u8>),
    Tx(u32, u32, u32),
}

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
    /// The read as given (sense: direct RNA, reads oriented before mapping). A read of the
    /// other strand is not always left unmapped: it may have a weak hit in this orientation.
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
            tail_first: None,
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
    /// Adjustment of each hit's quantification score, in the order of `hits`: the end fit
    /// ([`TieBreak::end_fit`]) and the orientation evidence kept out of the mapping score
    /// ([`TieBreak::orient_quant`]); empty if no hit is adjusted. A hit's quantification
    /// score is its mapping score plus this.
    pub end_fit: Vec<f32>,
}

pub struct Mapper<C: Chainer> {
    pub idx: Arc<Index>,
    pub seeder: BaseSeeder,
    pub chainer: C,
    pub p: MapParams,
    rescuer: Option<Rescuer>,
    /// [`anchor::node_families`] of `idx`.
    node_fam: Vec<u32>,
    /// Every transcript's antisense mates.
    anti: mates::Mates,
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

/// Orders a read's two orientations by `lead`, the first's over the second's in score points:
/// the one it favours comes first (`reverse` follows it), and the other's hits are shifted so
/// that its best lies `|lead|` below the first's best. With `level` the first's hits are first
/// raised so that its best is at least that.
fn rank_orientations(st: &mut Strand, other: &mut Strand, reverse: &mut bool, lead: f32, level: Option<f32>) {
    if lead < 0.0 {
        std::mem::swap(st, other);
        *reverse = !*reverse;
    }
    if let Some(level) = level
        && st.hits[0].score < level
    {
        let up = level - st.hits[0].score;
        st.hits.iter_mut().for_each(|h| h.score += up);
    }
    let target = st.hits[0].score - lead.abs();
    let shift = target - other.hits[0].score;
    for h in &mut other.hits {
        h.score += shift;
        h.reverse = !*reverse;
    }
}

/// The best scores `rank_orientations` leaves the two orientations, in the order given, from
/// their best scores `tops` before it.
fn ranked_tops(tops: (f32, f32), lead: f32, level: Option<f32>) -> (f32, f32) {
    let up = |s: f32| level.map_or(s, |l| s.max(l));
    if lead < 0.0 {
        let b = up(tops.1);
        (b - lead.abs(), b)
    } else {
        let a = up(tops.0);
        (a, a - lead.abs())
    }
}

/// A strand's `codes`, made from its `seq` on first use (a function of the two fields rather
/// than a method, so the strand's hits can be borrowed mutably meanwhile).
fn strand_codes<'a>(codes: &'a std::cell::OnceCell<Vec<u8>>, seq: &[u8]) -> &'a [u8] {
    codes.get_or_init(|| lcp_core::encode::codes(seq))
}

impl<C: Chainer> Mapper<C> {
    /// Whether the mapper made its rescue tables itself, the index keeping none it could
    /// use (false also when it needs none).
    pub fn made_rescue_tables(&self) -> bool {
        self.rescuer.as_ref().is_some_and(|r| r.tables_made())
    }

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
            anti: mates::Mates::load(&idx, p.anchor.h_max),
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
    /// With `whole`, the transcripts are compared projected to the whole read instead of the
    /// chains' span: the read then scores the same against both, as far as its alignments stay
    /// within the projections.
    fn same_antisense(&self, tb: &TieBreak, plen: u32, a: &TxHit, b: &TxHit, whole: bool) -> bool {
        // A hit's transcript interval projected to the read span of both chains (or the whole
        // read), with `pad` more bases on each side.
        let region = |h: &TxHit, o: &TxHit, pad: u32| {
            let tlen = self.idx.tx_plen()[h.tx as usize];
            let (lo, hi) = match whole {
                true => (0, plen),
                false => (
                    h.qstart.min(plen.saturating_sub(o.qend)),
                    h.qend.max(plen.saturating_sub(o.qstart)),
                ),
            };
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

    /// Whether a read orientation ends in a poly(A) tail, for the orientation comparison: an
    /// A-run followed by at most [`TAIL_AFTER`] read bases, or by at most
    /// [`TieBreak::orient_tail_after`] where the A-run starts at the end of the best hit's
    /// transcript (projected through the hit's last anchor). A cDNA read that kept its
    /// adapters carries the RT primer and the sequencing adapter after its tail, more than
    /// `TAIL_AFTER`; an A-run that far from the read's end that the transcript goes on past is
    /// the transcript's own sequence.
    fn orient_tailed(&self, tb: &TieBreak, s: &Strand, after: u32) -> bool {
        // Without HPC the parsed read is the raw read.
        let (codes, raw) = match &s.raw {
            Some(r) => (&r.codes[..], Some(r)),
            None if !self.seeder.hpc => (strand_codes(&s.codes, &s.seq), None),
            None => return false,
        };
        let a = lcp_core::encode::CODE[b'A' as usize] & 3;
        let w0 = codes.len().saturating_sub(TAIL_WINDOW);
        if tail_start_within(&codes[w0..], a, tb.tail_run, TAIL_AFTER).is_some() {
            return true;
        }
        let after = after as usize;
        let (Some(h), true) = (s.hits.first(), after > TAIL_AFTER) else {
            return false;
        };
        let Some((c, _)) = tail_start_within(&codes[w0..], a, tb.tail_run, after) else {
            return false;
        };
        let (rq, rt, tlen) = match raw {
            Some(r) => match self.idx.tx_raw_pos(h.tx, h.tend) {
                Some(rt) => (
                    r.pos(h.qend),
                    rt as usize,
                    self.idx.tx_len()[h.tx as usize] as usize,
                ),
                None => return false,
            },
            None => (
                h.qend as usize,
                h.tend as usize,
                self.idx.tx_plen()[h.tx as usize] as usize,
            ),
        };
        let delta = (w0 + c) as i64 - rq as i64;
        let tol = (TAIL_TOL + delta.max(0) as usize / 10) as i64;
        rt as i64 + delta + tol >= tlen as i64
    }

    /// A mirrored pair's lead in alignment points (`TieBreak::orient_same_exact`): the exact
    /// local alignment of the whole read in orientation `a` against its best hit, less that of
    /// orientation `b` against its own. Each orientation is `(strand, parsed codes)`, and both
    /// best hits are one sequence under the chains, so the difference is what the read holds
    /// beyond it. `None` where the hits are one sequence under the whole read too: the two
    /// alignments then score the same.
    fn mirror_exact(
        &self,
        tb: &TieBreak,
        plen: u32,
        a: (&Strand, &[u8]),
        b: (&Strand, &[u8]),
    ) -> Option<i32> {
        let sc = tb.affine.as_ref()?;
        let (ha, hb) = (a.0.hits.first()?, b.0.hits.first()?);
        if self.same_antisense(tb, plen, ha, hb, true) {
            return None;
        }
        let w = ORIENT_EXACT_BAND + plen / 20;
        let (ba, bb) = self.mirror_end_bonus(tb, plen, ha, hb);
        let xa = self.local_exact(tb, sc, a.1, a.0.raw.as_ref(), ha, w, ba);
        let xb = self.local_exact(tb, sc, b.1, b.0.raw.as_ref(), hb, w, bb);
        Some(xa - xb)
    }

    /// The end bonuses (read start, read end) of a mirrored pair's two alignments, `a`'s and
    /// `b`'s. `a`'s read start is `b`'s read end and the other way round. With
    /// `TieBreak::orient_same_end_cut` a read end earns the bonus only where both transcripts,
    /// projected through their hits, reach it; otherwise every read end can.
    fn mirror_end_bonus(&self, tb: &TieBreak, plen: u32, a: &TxHit, b: &TxHit) -> ((i32, i32), (i32, i32)) {
        if !tb.orient_same_end_cut {
            let both = (tb.end_bonus, tb.end_bonus);
            return (both, both);
        }
        let p = i64::from(plen);
        let tlen = |h: &TxHit| i64::from(self.idx.tx_plen()[h.tx as usize]);
        // Transcript bases before the hit's read start, and after its read end, beyond the
        // read bases there.
        let before = |h: &TxHit| i64::from(h.tstart) - i64::from(h.qstart);
        let after = |h: &TxHit| (tlen(h) - i64::from(h.tend)) - (p - i64::from(h.qend));
        let at_a_start = before(a) >= 0 && after(b) >= 0;
        let at_a_end = after(a) >= 0 && before(b) >= 0;
        let bonus = |reached: bool| if reached { tb.end_bonus } else { 0 };
        (
            (bonus(at_a_start), bonus(at_a_end)),
            (bonus(at_a_end), bonus(at_a_start)),
        )
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
            if tb.clip_pieces && q >= pq + 2 * tb.band as u32 && t > pt && q < h.qend && t < h.tend {
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

    /// Local affine score of the whole read (ends free to clip, `bonus.0` for reaching the
    /// read's start and `bonus.1` for its end) against one hit's transcript, over every
    /// diagonal within `w` of the chain's: the chain only places the search, it forces nothing
    /// into the alignment.
    #[allow(clippy::too_many_arguments)]
    fn local_exact(
        &self,
        tb: &TieBreak,
        sc: &align::Affine,
        codes: &[u8],
        raw: Option<&RawRead>,
        h: &TxHit,
        w: u32,
        bonus: (i32, i32),
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
        align::local_band_ends(read, &text, sc, dlo, dhi, bonus)
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

    /// Ask the caches for the raw positions [`Self::raw_ranges`] looks up for `spans`: their
    /// starts (`lo`) and ends, or only their ends where the starts were asked for already.
    fn prefetch_spans(&self, tb: &TieBreak, raw: Option<&RawRead>, spans: &[(u32, u32, u32)], lo: bool) {
        if raw.is_none() || !tb.raw_ends {
            return;
        }
        let q: Vec<(u32, u32)> = if lo {
            spans.iter().flat_map(|&(t, lo, hi)| [(t, lo), (t, hi)]).collect()
        } else {
            spans.iter().map(|&(t, _, hi)| (t, hi)).collect()
        };
        self.idx.prefetch_raw(&q);
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
                        align::affine_extend_sorted_cut_bonus(
                            read,
                            &tx,
                            sc,
                            tb.band,
                            tb.end_bonus,
                            tb.xdrop,
                            diag,
                            window,
                            cut.unwrap_or(usize::MAX),
                            tb.flank_min > 0,
                            tb.end_bonus_cut,
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
            let Some(cut) =
                cut.filter(|&c| (tb.flank_min > 0 || tb.end_bonus_cut) && texts.iter().any(|t| t.len() > c))
            else {
                return seg(texts, read, align::Mode::Prefix, rev, None).map(|(d, _, _)| d);
            };
            // With `end_bonus_cut` alone: the end bonus counts only within the cut.
            if tb.flank_min <= 0 {
                return seg(texts, read, align::Mode::Prefix, rev, Some(cut)).map(|(d, _, _)| d);
            }
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
        // Where a block is cut, from its start: the segment before it ends at the first cut and
        // the one after it starts at the second. A long block is cut `MID_INSET` bases inside
        // each edge, and its bases between the cuts are compared on their own; a shorter one
        // at its midpoint.
        let inset = tb.mid_inset;
        let cuts = |b: &SharedBlock| -> (u32, u32) {
            let len = b.query.end - b.query.start;
            if raw_mid.is_none() {
                (0, 0)
            } else if inset == 0 || len <= 2 * inset {
                (len / 2, len / 2)
            } else {
                (inset, len - inset)
            }
        };
        // Each candidate's raw position `off` bases into a block: counted on from its raw
        // position `from.1` at `from.0` bases into the block `from.2` (the same block or the
        // one before), whose bases up to this one are read anyway.
        let starts = |b: &SharedBlock, off: u32, from: Option<(u32, &[u32], &SharedBlock)>| -> Vec<u32> {
            hits.iter()
                .zip(&b.targets)
                .enumerate()
                .map(|(i, (h, t))| {
                    let p = t.start + off;
                    match from.map(|(o, raw, pb)| (pb.targets[i].start + o, raw[i])) {
                        Some((p0, r0)) if p0 <= p => self.idx.tx_raw_after(h.tx, r0, p - p0).unwrap(),
                        _ => self.idx.tx_raw_pos(h.tx, p).unwrap(),
                    }
                })
                .collect()
        };
        // The candidates' costs for the read's raw bases between two cuts.
        let between = |r: &RawRead, q: (u32, u32), lo: &[u32], hi: &[u32]| {
            let texts = hits
                .iter()
                .enumerate()
                .map(|(i, h)| self.idx.tx_raw_slice(h.tx, lo[i], hi[i]).unwrap())
                .collect();
            seg(
                texts,
                &r.codes[r.pos(q.0)..r.pos(q.1)],
                align::Mode::Global,
                false,
                None,
            )
            .map(|(d, _, _)| d)
        };
        // A block's bases between its cuts: the candidates' raw positions at its second cut
        // and their costs.
        let inside = |r: &RawRead, b: &SharedBlock, at: Vec<u32>| {
            let (c0, c1) = cuts(b);
            if c0 == c1 {
                return (at, None);
            }
            let to = starts(b, c1, Some((c0, &at, b)));
            let d = between(r, (b.query.start + c0, b.query.start + c1), &at, &to);
            (to, d)
        };
        // The first cut's raw positions are looked up for every candidate at once: ask for their
        // lines first (later cuts count on from these, through bases read anyway).
        if raw_mid.is_some() {
            let c0 = cuts(&blocks[0]).0;
            let q: Vec<(u32, u32)> = hits
                .iter()
                .zip(&blocks[0].targets)
                .map(|(h, t)| (h.tx, t.start + c0))
                .collect();
            self.idx.prefetch_raw(&q);
        }
        let mut prev = raw_mid.map(|r| {
            let at = starts(&blocks[0], cuts(&blocks[0]).0, None);
            let (at, d) = inside(r, &blocks[0], at);
            add(d);
            at
        });
        for pair in blocks.windows(2) {
            if let (Some(r), Some(lo)) = (raw_mid, prev.as_mut()) {
                let (c0, c1) = (cuts(&pair[0]).1, cuts(&pair[1]).0);
                let hi = starts(&pair[1], c1, Some((c0, lo, &pair[0])));
                add(between(
                    r,
                    (pair[0].query.start + c0, pair[1].query.start + c1),
                    lo,
                    &hi,
                ));
                let (at, d) = inside(r, &pair[1], hi);
                add(d);
                *lo = at;
                continue;
            }
            let (left, right) = (&pair[0], &pair[1]);
            let lo0 = left.targets[0].end;
            let len0 = right.targets[0].start - lo0;
            // Every candidate's segment is compared with the first's: ask for all their lines
            // before the first comparison waits on one.
            let q: Vec<(u32, u32)> = hits
                .iter()
                .zip(&left.targets)
                .map(|(h, t)| (h.tx, t.end))
                .collect();
            self.idx.prefetch_codes(&q);
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
        let h0 = cuts(first).0;
        let q0 = first.query.start + h0;
        if q0 > 0 {
            let spans: Vec<(u32, u32, u32)> = hits
                .iter()
                .zip(&first.targets)
                .map(|(h, t)| (h.tx, (t.start + h0).saturating_sub(tb.window(q0)), t.start + h0))
                .collect();
            let short: Vec<bool> = spans.iter().map(|&(_, lo, _)| lo == 0).collect();
            self.prefetch_spans(tb, raw, &spans, true);
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
        // The last block's raw ranges (compared below) and the right-end windows starting there.
        if raw.is_some() && tb.raw_ends {
            let q: Vec<(u32, u32)> = hits
                .iter()
                .zip(&last.targets)
                .flat_map(|(h, t)| [(h.tx, t.start), (h.tx, t.end)])
                .collect();
            self.idx.prefetch_raw(&q);
        }
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
        let hl = cuts(last).1;
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
        self.prefetch_spans(tb, raw, &spans, false);
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

    /// Transcripts sharing sequence with the region of hit `h` on the opposite strand, with
    /// the diagonals the read's other orientation would lie on: those sharing at least `min`
    /// of the read's `plen` parsed bases, the `k` sharing most first, in that order, then the
    /// others in no order.
    fn antisense_mates(&self, h: &TxHit, plen: u32, k: usize, min: f32) -> Vec<(u32, u32, i64, i64)> {
        let segs = self.anti.of(&self.idx, &self.seeder, self.p.anchor.h_max, h.tx);
        if segs.is_empty() {
            return Vec::new();
        }
        let ((dlo, dhi), (lo, hi)) = self.hit_span(h, plen);
        let min = ((min * plen as f32) as u64).max(1);
        stats::add(&stats::MATE_CALLS, 1);
        // Per distinct stretch, the bases it shares within the region (a stretch partly
        // inside it counts in proportion), if it reaches into the region at all.
        let mut stretches = Vec::new();
        let first = segs.stretches(&mut stretches);
        let within: Vec<Option<u64>> = stretches
            .iter()
            .map(|&(x0, x1, w)| {
                let (a, b) = (x0.max(lo), x1.min(hi));
                (a < b).then(|| u64::from(w) * u64::from(b - a) / u64::from(x1 - x0))
            })
            .collect();
        let mut out: Vec<(u32, u32, i64, i64)> = Vec::new();
        // The mate at hand: its shared bases and its range of diagonals (a seed at read
        // position q of this orientation lies at plen - q - len of the other).
        let (mut t, mut w, mut d) = (u32::MAX, 0u64, (i64::MAX, i64::MIN));
        let mut n = 0u64;
        segs.for_each(first, |u, id, c| {
            n += 1;
            if u != t {
                if w >= min {
                    out.push((t, w as u32, d.0, d.1));
                }
                (t, w, d) = (u, 0, (i64::MAX, i64::MIN));
            }
            if let Some(&Some(x)) = within.get(id) {
                w += x;
                let c = i64::from(c) - i64::from(plen);
                d = (d.0.min(c - dhi), d.1.max(c - dlo));
            }
        });
        if w >= min {
            out.push((t, w as u32, d.0, d.1));
        }
        stats::add(&stats::MATE_SEEDS, n);
        stats::add(&stats::MATE_FOUND, out.len() as u64);
        let key = |&(t, w, _, _): &(u32, u32, i64, i64)| (std::cmp::Reverse(w), t);
        if k > 0 && out.len() > k {
            out.select_nth_unstable_by_key(k - 1, key);
            out[..k].sort_unstable_by_key(key);
        } else {
            out.sort_unstable_by_key(key);
        }
        out
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
        let window = |x: u32| x + e.pad + x / e.pad_div;
        if h.qstart > 0 && h.tstart > 0 {
            let lo = h.tstart.saturating_sub(window(h.qstart));
            let mut text = self.idx.tx_codes(h.tx, lo, h.tstart);
            text.reverse();
            let read: Vec<u8> = codes[..h.qstart as usize].iter().rev().copied().collect();
            let (qi, tj) = Self::extend_one(e, seen, (true, h.reverse, h.qstart), &read, text);
            h.qstart -= qi;
            h.tstart -= tj;
        }
        (h.qend, h.tend) = self.right_reach(e, codes, h, seen);
    }

    /// Whether [`Self::right_reach`] of `h` aligns (rather than returning `h`'s own end, or
    /// its span extension's).
    fn right_aligns(&self, codes_len: u32, h: &TxHit) -> bool {
        !h.extended && codes_len > h.qend && h.tend < self.idx.tx_plen()[h.tx as usize]
    }

    /// Where `h`'s span extension over the read's right flank ends: the read and transcript
    /// positions of its last aligned base (as [`Self::extend_ends`] extends the right end).
    fn right_reach(&self, e: &Ends, codes: &[u8], h: &TxHit, seen: &mut EndCache) -> (u32, u32) {
        let plen = codes.len() as u32;
        let tlen = self.idx.tx_plen()[h.tx as usize];
        let rest = plen - h.qend;
        if rest > 0 && h.tend < tlen {
            stats::add(&stats::SPAN_ENDS, 1);
            let hi = (h.tend + rest + e.pad + rest / e.pad_div).min(tlen);
            let len = hi - h.tend;
            // The flank is compared packed with those of earlier extensions from the same read
            // position, and decoded only if none read a prefix of it.
            let done = seen.entry((false, h.reverse, h.qend)).or_default();
            let hit = done.iter().find(|(p, all, _)| match *p {
                Flank::Tx(t, lo, used) => {
                    used <= len && (!*all || len == used) && self.idx.tx_codes_eq(h.tx, h.tend, t, lo, used)
                }
                Flank::Codes(_) => false,
            });
            if let Some((_, _, (qi, tj))) = hit {
                return (h.qend + qi, h.tend + tj);
            }
            let text = self.idx.tx_codes(h.tx, h.tend, hi);
            let (r, used) = Self::extend_dp(e, &codes[h.qend as usize..], &text);
            done.push((Flank::Tx(h.tx, h.tend, used as u32), used == text.len(), r));
            (h.qend + r.0, h.tend + r.1)
        } else {
            (h.qend, h.tend)
        }
    }

    /// One read-end extension of `read` into `text` (see [`Self::extend_ends`]), cached in
    /// `seen` under `key` (left side, reverse, read position): the read and text bases it covers.
    fn extend_one(
        e: &Ends,
        seen: &mut EndCache,
        key: (bool, bool, u32),
        read: &[u8],
        text: Vec<u8>,
    ) -> (u32, u32) {
        stats::add(&stats::SPAN_ENDS, 1);
        let done = seen.entry(key).or_default();
        let hit = done.iter().find(|(p, all, _)| match p {
            Flank::Codes(p) => text.starts_with(p) && (!*all || text.len() == p.len()),
            Flank::Tx(..) => false,
        });
        if let Some((_, _, r)) = hit {
            return *r;
        }
        let (r, used) = Self::extend_dp(e, read, &text);
        let all = used == text.len();
        let mut text = text;
        text.truncate(used);
        done.push((Flank::Codes(text), all, r));
        r
    }

    /// The X-drop extension of `read` into `text` (a window of the read's length plus slack,
    /// cut short at the transcript's end): the read and text bases it covers, and the text
    /// bases it read.
    fn extend_dp(e: &Ends, read: &[u8], text: &[u8]) -> ((u32, u32), usize) {
        let window = |x: u32| x + e.pad + x / e.pad_div;
        let (m, w) = (read.len() as u64, window(read.len() as u32) as usize);
        stats::add(&stats::SPAN_END_CALLS, 1);
        stats::add(&stats::SPAN_END_CELLS, m * (2 * e.band as u64 + w as u64 - m + 1));
        let w0 = align::dp_work();
        let (_, qi, tj, used) =
            align::affine_extend_window(read, text, w, &e.sc, e.band, e.end_bonus, e.xdrop);
        if stats::on() {
            let w1 = align::dp_work();
            stats::add(&stats::DP_SPAN_CELLS, w1.0 - w0.0);
            stats::add(&stats::DP_SPAN_ROWS, w1.1 - w0.1);
        }
        ((qi as u32, tj as u32), used)
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
        // With `tail_first`: the orientation ending in a poly(A) tail, if exactly one does.
        let tailed = match (&self.p.tail_first, self.p.strands) {
            (Some(t), Strands::Both) => {
                let tail = |p: &(Vec<u8>, Option<RawRead>)| match &p.1 {
                    Some(r) => ends_in_tail(&r.codes, t.run, TAIL_AFTER),
                    // Without HPC the parsed read is the raw read.
                    None if self.seeder.hpc => false,
                    None => ends_in_tail(
                        &codes(&p.0[p.0.len().saturating_sub(TAIL_WINDOW)..]),
                        t.run,
                        TAIL_AFTER,
                    ),
                };
                match (tail(&pf), tail(&pr)) {
                    (true, false) => Some((false, t.span)),
                    (false, true) => Some((true, t.span)),
                    _ => None,
                }
            }
            _ => None,
        };
        let (fwd, rc) = match self.p.strands {
            Strands::Both => match tailed {
                None => (self.map_strand(pf, salt), self.map_strand(pr, salt)),
                Some((reverse, span)) => {
                    stats::add(&stats::TAIL_FIRST, 1);
                    let (first, second) = if reverse { (pr, pf) } else { (pf, pr) };
                    let first = self.map_strand(first, salt);
                    let second = if Self::span_frac(first.hits.first(), first.seq.len()) >= span {
                        stats::add(&stats::TAIL_FIRST_ALONE, 1);
                        skip(second)
                    } else {
                        self.map_strand(second, salt)
                    };
                    if reverse { (second, first) } else { (first, second) }
                }
            },
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
                stats::add(&stats::WEAK, 1);
                // The candidates of both orientations: one whose shared seeds fall far short
                // of the best candidate's (of either orientation) is not chained, past the
                // first few of its orientation.
                let (mut va, mut vb) = (
                    r.global_votes(&st.seq, g.weak_cands),
                    r.global_votes(&other.seq, g.weak_cands),
                );
                let most = va.first().map_or(0, |c| c.0).max(vb.first().map_or(0, |c| c.0));
                for v in [&mut va, &mut vb] {
                    let keep = v
                        .iter()
                        .enumerate()
                        .take_while(|&(k, c)| {
                            k < g.weak_keep || f64::from(c.0) >= f64::from(g.weak_rel) * f64::from(most)
                        })
                        .count();
                    v.truncate(keep);
                }
                let seed = |s: &Strand| {
                    let v = if std::ptr::eq(s, &st) { &va } else { &vb };
                    let mut h = r.global_chain(&self.idx, &s.seq, &self.p.chain, v);
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
                stats::add(&stats::RECHAIN_TX, txs.len() as u64);
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
            stats::add(&stats::UNSETTLED, u64::from(!settled));
            let mates = if r.p.antisense > 0 && matches!(self.p.strands, Strands::Both) && !settled {
                self.antisense_mates(&top, plen, r.p.antisense, r.p.antisense_share)
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
            if !settled
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
                if let Some((_, g)) = global {
                    let new = r.global_hits(&self.idx, &other.seq, &self.p.chain, g.extra_cands);
                    merge(&mut other.hits, new);
                }
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
                end_fit: Vec::new(),
            };
        }
        let mut both = false;
        // With `TieBreak::orient_quant`: the orientation (`reverse`) that came first, and what
        // the orientation rules it leaves out of the mapping score add to the quantification
        // scores of that orientation's hits and of the other's.
        let mut orient_q = None;
        if let Some(tb) = &self.p.tiebreak {
            let margin = |s: f32| tb.margin_abs.max(tb.orient_margin_rel * s.abs());
            both = matches!((st.hits.first(), other.hits.first()),
                    (Some(a), Some(b)) if b.score >= a.score - margin(a.score));
            let codes = || strand_codes(&st.codes, &st.seq);
            self.tiebreak(tb, codes, st.raw.as_ref(), &mut st.hits, salt);
            // The two best hits are one sequence on opposite strands: the read fits both
            // equally, however far apart seed luck leaves the two chains.
            let same = |a: &Strand, b: &Strand| match (a.hits.first(), b.hits.first()) {
                (Some(a), Some(b)) => self.same_antisense(tb, plen, a, b, false),
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
                // poly(A) tail where one orientation shows it and the other doesn't. `lead`
                // takes a tail limit (`TieBreak::orient_tail_after`, or `TAIL_AFTER` as
                // before it) and whether a mirrored pair is compared (`orient_same_exact`).
                let aligned = unit * (cb - ca) as f32;
                let tail = |after: u32| {
                    if tb.orient_tail == 0 {
                        return 0.0;
                    }
                    let tailed = |s: &Strand| f32::from(u8::from(self.orient_tailed(tb, s, after)));
                    tb.orient_tail as f32 * (tailed(&st) - tailed(&other))
                };
                // One sequence on opposite strands: exact local alignments of the whole read
                // tell the pair apart by what lies beyond it (`TieBreak::orient_same_exact`).
                let same_x = (mirror && tb.orient_same_exact)
                    .then(|| self.mirror_exact(tb, plen, (&st, qa), (&other, qb)))
                    .flatten();
                // A close call is decided by exact local alignments instead of the clipped
                // ones (`TieBreak::orient_exact`); the tail term stands. Aligned once, for
                // either lead that asks.
                let mut exact = None;
                let mut lead = |after: u32, same: bool| {
                    let t = tail(after);
                    let mut lead = aligned + t;
                    if same && let Some(x) = same_x {
                        lead += tb.lambda * x as f32;
                    }
                    if let Some(sc) = &tb.affine
                        && !mirror
                        && tb.orient_clipped
                        && tb.orient_exact > 0.0
                        && lead.abs() <= tb.orient_exact
                    {
                        let x = *exact.get_or_insert_with(|| {
                            let w = ORIENT_EXACT_BAND + plen / 20;
                            let both = (tb.end_bonus, tb.end_bonus);
                            let xa = self.local_exact(tb, sc, qa, st.raw.as_ref(), &st.hits[0], w, both);
                            let xb =
                                self.local_exact(tb, sc, qb, other.raw.as_ref(), &other.hits[0], w, both);
                            xa - xb
                        });
                        lead = unit * x as f32 + t;
                    }
                    lead
                };
                let full = lead(tb.orient_tail_after, true);
                // One sequence scores the same whichever orientation's chain found it: the pair
                // stands at the better chain's score, not at the one seed luck left lower.
                let tops = (st.hits[0].score, other.hits[0].score);
                let level = (mirror && tb.orient_same_exact).then(|| tops.0.max(tops.1));
                if tb.orient_quant {
                    // The mapping decision as without the two rules; the quantification
                    // score takes the difference they make to each orientation's scores.
                    let plain = lead(TAIL_AFTER as u32, false);
                    let (a, b) = ranked_tops(tops, full, level);
                    let (a0, b0) = ranked_tops(tops, plain, None);
                    orient_q = Some((reverse, a - a0, b - b0));
                    rank_orientations(&mut st, &mut other, &mut reverse, plain, None);
                } else {
                    rank_orientations(&mut st, &mut other, &mut reverse, full, level);
                }
            }
        }
        for h in &mut st.hits {
            h.reverse = reverse;
        }
        let mut hits = std::mem::take(&mut st.hits);
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
        let mut seen = rustc_hash::FxHashMap::default();
        if let Some(e) = &self.p.ends {
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
        // The quantification score's end fit: after everything that decides the hits and
        // their scores, which it leaves as they are.
        let mut end_fit = match &self.p.tiebreak {
            Some(tb) if tb.end_fit > 0.0 && hits.len() >= 2 => {
                self.end_fit(tb, &hits, reverse, (&st, &other), &mut seen)
            }
            _ => Vec::new(),
        };
        if let Some((first, a, b)) = orient_q
            && (a != 0.0 || b != 0.0)
            && !hits.is_empty()
        {
            end_fit.resize(hits.len(), 0.0);
            for (q, h) in end_fit.iter_mut().zip(&hits) {
                *q += if h.reverse == first { a } else { b };
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
            end_fit,
        }
    }

    /// End-fit adjustments of the quantification scores of `hits` (reported, best first;
    /// `reverse`: the orientation of `st`, the first of `strands`), or none.
    ///
    /// Without a poly(A) tail a read carries no strand or 3'-end evidence in its mapping
    /// score, and isoforms that differ only in where they end tie there. Where the read ends on
    /// each candidate is measured from its right-end alignment (the span extension's, run here
    /// for candidates not extended, unless its result cannot change the outcome; with
    /// [`TieBreak::end_fit_chain`], the chain's gapless continuation to the read's end instead):
    /// `q`, the read position of its last aligned base, and `d`, the transcript bases left from
    /// there to the candidate's annotated 3' end. A
    /// candidate fits if `d` is at most `end_fit_tol` and its alignment reaches as far into
    /// the read (in its own orientation) as the farthest-reaching candidate's, within
    /// `end_fit_tol`: a shorter isoform that ends inside the read does not fit.
    ///
    /// Each candidate's annotated end is placed on the read, relative to that orientation's
    /// read end (`q + d` minus the farthest `q`); ends within `end_fit_tol` of a fitting
    /// candidate's are not told apart from it and gain as it does (transitively). The
    /// candidates are the hits within the tie-break margin of the best, of either orientation,
    /// so near-tied orientations are weighed by which one's transcript ends where the read
    /// does. Nothing is adjusted if the read has a tail in either orientation (the tail is
    /// in the mapping score), if no candidate fits, or if every candidate gains.
    fn end_fit(
        &self,
        tb: &TieBreak,
        hits: &[TxHit],
        reverse: bool,
        (st, other): (&Strand, &Strand),
        seen: &mut EndCache,
    ) -> Vec<f32> {
        let best = hits[0].score;
        let margin = tb.margin_abs.max(tb.margin_rel * best.abs());
        let n = hits
            .iter()
            .take(tb.max_cands)
            .take_while(|h| h.score >= best - margin)
            .count();
        if n < 2 {
            return Vec::new();
        }
        let tol = tb.end_fit_tol;
        let e = self.p.ends.clone().unwrap_or_else(|| Ends {
            sc: tb.affine.unwrap_or(align::Affine::ONT),
            ..Default::default()
        });
        let plen = st.seq.len() as u32;
        let tlen = |h: &TxHit| self.idx.tx_plen()[h.tx as usize];
        // Only a candidate whose 3' end lies within reach of its right-end extension can fit.
        let reachable = |h: &TxHit| {
            let rest = plen.saturating_sub(h.qend);
            tlen(h) - h.tend <= rest + e.pad + rest / e.pad_div + tol
        };
        if !hits[..n].iter().any(reachable) {
            return Vec::new();
        }
        // Without HPC the parsed read is the raw read; an orientation that was not mapped has
        // no sequence.
        let tailed = |s: &Strand| match &s.raw {
            Some(r) => ends_in_tail(&r.codes, tb.tail_run, TAIL_AFTER),
            None => {
                !self.seeder.hpc
                    && !s.seq.is_empty()
                    && ends_in_tail(strand_codes(&s.codes, &s.seq), tb.tail_run, TAIL_AFTER)
            }
        };
        if tailed(st) || tailed(other) {
            return Vec::new();
        }
        let side = |h: &TxHit| usize::from(h.reverse != reverse);
        let strand = |h: &TxHit| if side(h) == 0 { st } else { other };
        let mut reach_of = |h: &TxHit| {
            let s = strand(h);
            if h.extended {
                (h.qend, h.tend)
            } else if tb.end_fit_chain {
                let k = (s.seq.len() as u32).saturating_sub(h.qend).min(tlen(h) - h.tend);
                (h.qend + k, h.tend + k)
            } else {
                self.right_reach(&e, strand_codes(&s.codes, &s.seq), h, seen)
            }
        };
        // The reachable candidates first: if none of their alignments ends within the tolerance
        // of its 3' end, none fits, and the others need no alignment. (Candidates without an
        // alignment to run come along: they cost nothing; with `end_fit_chain`, none has one.)
        let mut reach: Vec<Option<(u32, u32)>> = hits[..n]
            .iter()
            .map(|h| {
                let slen = strand(h).seq.len() as u32;
                (reachable(h) || tb.end_fit_chain || !self.right_aligns(slen, h)).then(|| reach_of(h))
            })
            .collect();
        if !hits[..n]
            .iter()
            .zip(&reach)
            .any(|(h, r)| r.is_some_and(|(_, t)| tlen(h) - t <= tol))
        {
            return Vec::new();
        }
        let tol = i64::from(tol);
        // A candidate left is not reachable, so it cannot fit. Where its alignment can end is
        // bounded by the band (`Ends::band` cells either side of the diagonals the window's
        // slack spans): its read position `q` lies within `qend` + min(flank, window + band),
        // and the read position of its annotated end, `q + d`, at least `qend` + the
        // transcript bases left − min(window, band + slack). Such a candidate is aligned only
        // if that bound lets it raise the farthest reach of its orientation, or lie within the
        // tolerance of a gaining candidate's end; otherwise its alignment cannot change the
        // outcome. (None: aligned anyway, e.g. if its own orientation's flank is longer.)
        let bound: Vec<Option<(u32, i64)>> = hits[..n]
            .iter()
            .zip(&reach)
            .map(|(h, r)| {
                let rest = (strand(h).seq.len() as u32).checked_sub(h.qend)?;
                let w = rest + e.pad + rest / e.pad_div;
                let tl = tlen(h) - h.tend;
                (r.is_none() && i64::from(tl) > i64::from(w) + tol).then(|| {
                    let l = w.min(tl);
                    let band = u32::try_from(e.band).unwrap_or(u32::MAX);
                    let qub = h.qend + rest.min(l.saturating_add(band));
                    let slack = i64::from(w - rest);
                    let la = i64::from(h.qend) + i64::from(tl) - i64::from(l).min(i64::from(band) + slack);
                    (qub, la)
                })
            })
            .collect();
        for (h, (r, b)) in hits[..n].iter().zip(reach.iter_mut().zip(&bound)) {
            if r.is_none() && b.is_none() {
                *r = Some(reach_of(h));
            }
        }
        let gain = loop {
            let (mut far, mut known, mut qhi) = ([0u32; 2], [false; 2], [0u32; 2]);
            for ((h, r), b) in hits[..n].iter().zip(&reach).zip(&bound) {
                match (r, b) {
                    (Some((q, _)), _) => {
                        far[side(h)] = far[side(h)].max(*q);
                        known[side(h)] = true;
                    }
                    (None, Some((qub, _))) => qhi[side(h)] = qhi[side(h)].max(*qub),
                    (None, None) => unreachable!("aligned above"),
                }
            }
            // The farthest reach of an orientation with aligned candidates must be final.
            let mut grew = false;
            for ((h, r), b) in hits[..n].iter().zip(reach.iter_mut()).zip(&bound) {
                if r.is_none() && known[side(h)] && b.is_some_and(|(qub, _)| qub > far[side(h)]) {
                    *r = Some(reach_of(h));
                    grew = true;
                }
            }
            if grew {
                continue;
            }
            // Bases from the last aligned base to the annotated end, and that end's place
            // relative to the read's end, of the aligned candidates.
            let (d, at): (Vec<i64>, Vec<i64>) = hits[..n]
                .iter()
                .zip(&reach)
                .map(|(h, r)| {
                    r.map_or((i64::MAX, i64::MAX), |(q, t)| {
                        let d = i64::from(tlen(h)) - i64::from(t);
                        (d, i64::from(q) + d - i64::from(far[side(h)]))
                    })
                })
                .unzip();
            let mut gain: Vec<bool> = (0..n)
                .map(|i| {
                    reach[i].is_some_and(|(q, _)| {
                        d[i] <= tol && i64::from(q) + tol >= i64::from(far[side(&hits[i])])
                    })
                })
                .collect();
            if !gain.contains(&true) {
                return Vec::new();
            }
            loop {
                let mut grew = false;
                for j in 0..n {
                    if !gain[j]
                        && reach[j].is_some()
                        && (0..n).any(|i| gain[i] && (at[i] - at[j]).abs() <= tol)
                    {
                        gain[j] = true;
                        grew = true;
                    }
                }
                if !grew {
                    break;
                }
            }
            // A candidate not aligned gains only if its end can lie within the tolerance of a
            // gaining one's.
            let top = (0..n)
                .filter(|&i| gain[i])
                .map(|i| at[i])
                .max()
                .unwrap_or(i64::MIN);
            for ((h, r), b) in hits[..n].iter().zip(reach.iter_mut()).zip(&bound) {
                if let (None, Some((_, la))) = (&r, b) {
                    let s = side(h);
                    let lb = la - i64::from(if known[s] { far[s] } else { qhi[s] });
                    if lb <= top + tol {
                        *r = Some(reach_of(h));
                        grew = true;
                    }
                }
            }
            if !grew {
                break gain;
            }
        };
        if gain.iter().all(|&g| g) {
            return Vec::new();
        }
        stats::add(&stats::END_FIT, 1);
        let mut out = vec![0.0; hits.len()];
        for (o, &g) in out.iter_mut().zip(&gain) {
            if g {
                *o = tb.end_fit;
            }
        }
        out
    }
}

/// Diagonals searched on each side of a chain's by the exact comparison of near-tied read
/// orientations ([`TieBreak::orient_exact`]), plus a twentieth of the parsed read.
pub const ORIENT_EXACT_BAND: u32 = 100;
/// Bases between a transcript's 3' end and a read's tail start treated as the same end.
const TAIL_TOL: usize = 15;
/// Most read bases after a poly(A) tail (adapter and the like).
pub const TAIL_AFTER: usize = 80;

/// Serde default of [`TieBreak::orient_tail_after`]: options recorded without it map as they
/// did then.
fn tail_after() -> u32 {
    TAIL_AFTER as u32
}
/// Raw read bases searched for a poly(A) tail when comparing read orientations.
const TAIL_WINDOW: usize = 200;

/// Whether a raw read (2-bit codes) ends in a poly(A) tail: a `run`-base window of A's within its
/// last `TAIL_WINDOW` bases, followed by at most `after` bases.
fn ends_in_tail(codes: &[u8], run: usize, after: usize) -> bool {
    let a = lcp_core::encode::CODE[b'A' as usize] & 3;
    tail_start_within(&codes[codes.len().saturating_sub(TAIL_WINDOW)..], a, run, after).is_some()
}

/// The first A-run of `text` (raw codes) starting in `[lo, hi]`: a window of `run` bases with at
/// most one other base. Returns its first A and the end of the A-run from there (as
/// [`tail_start`]).
fn own_run(text: &[u8], a: u8, run: usize, lo: usize, hi: usize) -> Option<(usize, usize)> {
    let last = hi.min(text.len().checked_sub(run)?);
    if lo > last {
        return None;
    }
    // The other bases in the window at `i`, counted on as it slides instead of again at each
    // start: a read end's windows are searched for every candidate of a tie set.
    let mut other = text[lo..lo + run].iter().filter(|&&c| c != a).count();
    let mut i = lo;
    while other > 1 {
        if i == last {
            return None;
        }
        other = other + usize::from(text[i + run] != a) - usize::from(text[i] != a);
        i += 1;
    }
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
    tail_start_within(rest, a, run, TAIL_AFTER)
}

/// [`tail_start`] with at most `after` read bases after the tail.
fn tail_start_within(rest: &[u8], a: u8, run: usize, after: usize) -> Option<(usize, usize)> {
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
            if rest.len() - end <= after {
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
            end_bonus_cut: false,
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
    fn own_run_finds_the_first_window_as_a_scan_of_every_start_does() {
        // The window search as one scan per start, then the run's end as `own_run` finds it.
        fn scan(text: &[u8], a: u8, run: usize, lo: usize, hi: usize) -> Option<(usize, usize)> {
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
        let a = codes(b"A")[0];
        let check = |text: &[u8], run: usize, lo: usize, hi: usize| {
            assert_eq!(
                own_run(text, a, run, lo, hi),
                scan(text, a, run, lo, hi),
                "{text:?} run {run} in [{lo}, {hi}]"
            );
        };
        // A text shorter than the window, a run at either end, no A at all, an empty range.
        for t in [
            &b"AAAA"[..],
            b"AAAAAAAAAAAAAAACGT",
            b"CGTCGTAAAAAAAAAAAAAAA",
            b"CGTCGTCGTCGTCGTCGTCG",
            b"",
        ] {
            let t = codes(t);
            for run in [0, 1, 5, 15, 30] {
                for (lo, hi) in [(0, 0), (0, 40), (3, 2), (5, 10), (20, 100)] {
                    check(&t, run, lo, hi);
                }
            }
        }
        // Random A-rich texts (a fixed xorshift sequence).
        let mut x = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        for _ in 0..2000 {
            let n = (next() % 80) as usize;
            let rich = next() % 4;
            let t: Vec<u8> = (0..n)
                .map(|_| if next() % 4 <= rich { a } else { (next() % 4) as u8 })
                .collect();
            let run = (next() % 20) as usize;
            let lo = (next() % 60) as usize;
            let hi = lo + (next() % 40) as usize;
            check(&t, run, lo, hi);
        }
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
        assert!(ends_in_tail(&codes(&[&body[..], &read].concat()), 15, TAIL_AFTER));
    }

    #[test]
    fn only_a_poly_a_near_the_read_end_marks_its_orientation() {
        let body = b"ACGTCGATCGTACGATCGTCAGCTAGCTGACTGCATGCAT".repeat(10);
        let tail = [b"A".repeat(20), b"GTACTCTGCGTTGATACCACTGCTT".to_vec()].concat();
        assert!(ends_in_tail(&codes(&[&body[..], &tail].concat()), 15, TAIL_AFTER));
        // A tail at the read's start is the other orientation's 3' end.
        assert!(!ends_in_tail(
            &codes(&[&tail[..], &body].concat()),
            15,
            TAIL_AFTER
        ));
        // An A-run with more than `after` read bases after it is transcript sequence.
        assert!(!ends_in_tail(
            &codes(&[&body[..], &tail, &body[..100]].concat()),
            15,
            120
        ));
        // The orientation comparison accepts a longer adapter after the tail than the tie-break.
        let long = [&body[..], &tail, &body[..60]].concat();
        assert!(!ends_in_tail(&codes(&long), 15, TAIL_AFTER));
        assert!(ends_in_tail(&codes(&long), 15, 120));
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
    fn end_fit_favours_the_isoform_whose_end_the_tailless_read_reaches() {
        // A pseudo-random body; t0 ends with it, t1 continues 200 bases, t2 only 3.
        let mut x = 7u64;
        let mut base = || {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 62) as usize]
        };
        let body: Vec<u8> = (0..300).map(|_| base()).collect();
        let more: Vec<u8> = (0..200).map(|_| base()).collect();
        let t1: Vec<u8> = [&body[..], &more].concat();
        let t2: Vec<u8> = [&body[..], &more[..3]].concat();
        let m = mapper(&[&body, &t1, &t2], false);
        let tb = TieBreak::default();
        // Hits of a read of `body[100..end]` (plus `tail`), aligned to `end` on each candidate.
        let fit = |txs: &[u32], end: u32, tail: &[u8]| {
            let read: Vec<u8> = [&body[100..end as usize], tail].concat();
            let (seq, raw) = Mapper::<OracleChainer>::parse(&m.seeder, &read);
            let strand = |seq: Vec<u8>, raw| Strand {
                hits: Vec::new(),
                seq,
                codes: Default::default(),
                raw,
                seeds: Vec::new(),
            };
            let rc = revcomp(&read);
            let (st, other) = (strand(seq, raw), strand(rc, None));
            let hits: Vec<TxHit> = txs
                .iter()
                .map(|&t| {
                    let mut h = hit(t, 100.0, &[(0, 100, 40)]);
                    (h.qend, h.tend, h.extended) = (end - 100, end, true);
                    h
                })
                .collect();
            m.end_fit(&tb, &hits, false, (&st, &other), &mut Default::default())
        };
        // The read stops at t0's end: t0 gains, t1 (200 bases on) does not.
        assert_eq!(fit(&[0, 1], 300, b""), vec![tb.end_fit, 0.0]);
        assert_eq!(fit(&[1, 0], 300, b""), vec![0.0, tb.end_fit]);
        // A poly(A) tail is the mapping score's evidence: nothing here.
        assert!(fit(&[0, 1], 300, &[b'A'; 30]).is_empty());
        // The read stops inside both: neither fits.
        assert!(fit(&[0, 1], 250, b"").is_empty());
        // Ends 3 bases apart are within the measurement error: no relative adjustment, and t2
        // keeps t1 apart from t0 by no more than it is itself.
        assert!(fit(&[0, 2], 300, b"").is_empty());
        assert_eq!(fit(&[0, 2, 1], 300, b""), vec![tb.end_fit, tb.end_fit, 0.0]);
        // One candidate only: nothing to weigh.
        assert!(fit(&[0], 300, b"").is_empty());
    }

    #[test]
    fn end_fit_of_unextended_hits_by_alignment_and_by_chain() {
        // As above: t0 ends with the body, t1 continues 200 bases, t2 only 3. The hits end at
        // their last anchor, 40 bases before the read's end (which is t0's), so the end fit
        // places their ends itself: t1's lies far past the read's end (no alignment needed to
        // rule it out), t0's and t2's within reach.
        let mut x = 7u64;
        let mut base = || {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 62) as usize]
        };
        let body: Vec<u8> = (0..300).map(|_| base()).collect();
        let more: Vec<u8> = (0..200).map(|_| base()).collect();
        let t1: Vec<u8> = [&body[..], &more].concat();
        let t2: Vec<u8> = [&body[..], &more[..3]].concat();
        let m = mapper(&[&body, &t1, &t2], false);
        let fit = |txs: &[u32], end: u32, chain: bool| {
            let tb = TieBreak {
                end_fit_chain: chain,
                ..Default::default()
            };
            let read: Vec<u8> = body[100..end as usize].to_vec();
            let (seq, raw) = Mapper::<OracleChainer>::parse(&m.seeder, &read);
            let strand = |seq: Vec<u8>, raw| Strand {
                hits: Vec::new(),
                seq,
                codes: Default::default(),
                raw,
                seeds: Vec::new(),
            };
            let rc = revcomp(&read);
            let (st, other) = (strand(seq, raw), strand(rc, None));
            let hits: Vec<TxHit> = txs
                .iter()
                .map(|&t| {
                    let mut h = hit(t, 100.0, &[(0, 100, 40)]);
                    (h.qend, h.tend) = (end - 140, end - 40);
                    h
                })
                .collect();
            m.end_fit(&tb, &hits, false, (&st, &other), &mut Default::default())
        };
        for chain in [false, true] {
            let g = TieBreak::default().end_fit;
            assert_eq!(fit(&[0, 1], 300, chain), vec![g, 0.0]);
            assert_eq!(fit(&[1, 0, 2], 300, chain), vec![0.0, g, g]);
            // The read stops inside every candidate: none fits.
            assert!(fit(&[0, 1], 250, chain).is_empty());
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
            end_bonus_cut: false,
            ..Default::default()
        };
        let read: Vec<u8> = [b"GT", body].concat();
        let c = m.segmented(&tb, &codes(&read), None, &hits).unwrap().0;
        assert!(c[1] < c[0]);
    }

    #[test]
    fn end_bonus_counts_only_where_every_candidate_has_sequence() {
        // t0 ends where t1 goes on. A read from t0 with one unaligned base after its end
        // reaches its end in t1 only, through a mismatch; that is no evidence for t1.
        let body: &[u8] = b"ACGTCGATCGTACGATCGTCAGCTAGCTGACTGCATGCAT";
        let t1: Vec<u8> = [body, b"GTCAGTACGT"].concat();
        let m = mapper(&[body, &t1], false);
        let hits = [hit(0, 10., &[(0, 0, 30)]), hit(1, 10., &[(0, 0, 30)])];
        let cost = |read: &[u8], end_bonus_cut| {
            let tb = TieBreak {
                end_bonus_cut,
                ..Default::default()
            };
            m.segmented(&tb, &codes(read), None, &hits).unwrap().0
        };
        let junk: Vec<u8> = [body, b"C"].concat();
        let (off, on) = (cost(&junk, false), cost(&junk, true));
        // Without the rule t1 wins by the end bonus less a mismatch.
        let tb = TieBreak::default();
        let mis = i64::from(tb.affine.as_ref().unwrap().mis);
        assert_eq!(off[0] - off[1], i64::from(tb.end_bonus) - mis);
        assert_eq!(on[0], on[1]);
        // A read that goes on along t1 still prefers it, by its matched bases.
        let longer: Vec<u8> = [body, b"GTC"].concat();
        let (off, on) = (cost(&longer, false), cost(&longer, true));
        assert!(on[1] < on[0]);
        assert_eq!(on[0] - on[1], off[0] - off[1] - i64::from(tb.end_bonus));
    }

    /// A deterministic sequence without long repeats (xorshift).
    fn random_seq(n: usize, mut x: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x >> 33) as usize & 3]
            })
            .collect()
    }

    /// One read orientation with the given hits (the parsed read is the raw read, no HPC).
    fn strand(seq: &[u8], hits: Vec<TxHit>) -> Strand {
        Strand {
            hits,
            seq: seq.to_vec(),
            codes: Default::default(),
            raw: None,
            seeds: Vec::new(),
        }
    }

    /// The mirrored-pair comparison of read orientation `a` (codes `qa`) against `b`.
    fn mirror(m: &Mapper<OracleChainer>, tb: &TieBreak, a: &Strand, b: &Strand) -> Option<i32> {
        let (qa, qb) = (codes(&a.seq), codes(&b.seq));
        m.mirror_exact(tb, qa.len() as u32, (a, &qa), (b, &qb))
    }

    #[test]
    fn a_mirrored_pair_is_told_apart_by_what_lies_beyond_the_shared_sequence() {
        // t1 is 8 bases of its own, then the reverse complement of t0[100..250): a short
        // antisense isoform. A read of t1 fits both orientations equally under the chains.
        let t0 = random_seq(300, 7);
        let own: &[u8] = b"GTACCAAA";
        let t1 = [own, &revcomp(&t0[100..250])].concat();
        let m = mapper(&[&t0, &t1], false);
        let tb = TieBreak::default();
        let fwd = strand(&t1, vec![hit(1, 50., &[(8, 8, 150)])]);
        let rev = strand(&revcomp(&t1), vec![hit(0, 50., &[(0, 100, 150)])]);
        // Exact alignments of the whole read: t1's own bases decide, for t1.
        let lead = mirror(&m, &tb, &fwd, &rev).unwrap();
        let mat = tb.affine.as_ref().unwrap().mat;
        assert!(lead > 0 && lead <= 8 * mat + tb.end_bonus, "{lead}");
        // The comparison does not depend on which orientation is asked first.
        assert_eq!(mirror(&m, &tb, &rev, &fwd), Some(-lead));
        // Without the own bases the read is one sequence with both transcripts under its whole
        // length: both alignments score the same, and they are skipped.
        let read = revcomp(&t0[100..250]);
        let fwd = strand(&read, vec![hit(1, 50., &[(0, 8, 150)])]);
        let rev = strand(&revcomp(&read), vec![hit(0, 50., &[(0, 100, 150)])]);
        assert_eq!(mirror(&m, &tb, &fwd, &rev), None);
        let both = (tb.end_bonus, tb.end_bonus);
        let sc = tb.affine.as_ref().unwrap();
        let x = |s: &Strand| m.local_exact(&tb, sc, &codes(&s.seq), None, &s.hits[0], 100, both);
        assert_eq!(x(&fwd), x(&rev));
    }

    #[test]
    fn a_mirrored_pair_earns_the_end_bonus_only_where_both_transcripts_reach_the_read_end() {
        // t1 is t0[100..250) reverse complemented, wholly inside t0. A read of t1 with one junk
        // base after t1's 3' end reaches that end only in t0's orientation, through the junk
        // base aligned as a mismatch into t0's continuation.
        let mut t0 = random_seq(300, 11);
        t0[99] = b'A';
        let t1 = revcomp(&t0[100..250]);
        let m = mapper(&[&t0, &t1], false);
        let read = [&t1[..], b"C"].concat();
        let fwd = strand(&read, vec![hit(1, 50., &[(0, 0, 150)])]);
        let rev = strand(&revcomp(&read), vec![hit(0, 50., &[(1, 100, 150)])]);
        let lead = |orient_same_end_cut| {
            let tb = TieBreak {
                orient_same_end_cut,
                ..Default::default()
            };
            mirror(&m, &tb, &fwd, &rev).unwrap()
        };
        // Without the rule t0 wins by the end bonus less a mismatch.
        let tb = TieBreak::default();
        let mis = tb.affine.as_ref().unwrap().mis;
        assert_eq!(lead(false), -(tb.end_bonus - mis));
        // With it neither alignment earns a bonus at that read end: a tie.
        assert_eq!(lead(true), 0);
        // Where both transcripts reach both read ends, every read end earns it.
        let (a, b) = (&fwd.hits[0], &rev.hits[0]);
        let e = tb.end_bonus;
        assert_eq!(m.mirror_end_bonus(&tb, 151, a, b), ((e, 0), (0, e)));
        let (a, b) = (hit(1, 50., &[(0, 0, 150)]), hit(0, 50., &[(0, 100, 150)]));
        assert_eq!(m.mirror_end_bonus(&tb, 150, &a, &b), ((e, e), (e, e)));
    }

    #[test]
    fn a_long_adapter_after_the_tail_counts_only_at_the_transcript_end() {
        // A cDNA read that kept its adapters: the tail, then 100 bases of primer and adapter.
        let body = random_seq(400, 13);
        let adapter = b"GTACTCTGCGTTGATCCCACTGCTT".repeat(4);
        let read = [&body[..], &b"A".repeat(20), &adapter].concat();
        let ends_at_tail = body.clone();
        let goes_on = [&body[..], &random_seq(300, 17)].concat();
        let m = mapper(&[&ends_at_tail, &goes_on], false);
        let tailed = |tx: u32, orient_tail_after: u32, read: &[u8]| {
            let tb = TieBreak::default();
            m.orient_tailed(
                &tb,
                &strand(read, vec![hit(tx, 50., &[(0, 0, 400)])]),
                orient_tail_after,
            )
        };
        // The A-run starts at the hit transcript's 3' end: a tail, with up to 120 bases after.
        assert!(tailed(0, 120, &read));
        assert!(!tailed(0, TAIL_AFTER as u32, &read));
        // Where the transcript goes on past it, an A-run that far from the read's end is not.
        assert!(!tailed(1, 120, &read));
        // With at most `TAIL_AFTER` bases after it, any A-run marks the tail, as before.
        let short = [&body[..], &b"A".repeat(20), &adapter[..60]].concat();
        assert!(tailed(1, 120, &short) && tailed(1, TAIL_AFTER as u32, &short));
    }

    #[test]
    fn a_mirrored_pair_stands_at_the_better_chain_score() {
        let hits = |tx: u32, scores: &[f32]| -> Vec<TxHit> {
            scores.iter().map(|&s| hit(tx, s, &[(0, 0, 10)])).collect()
        };
        let rank = |lead: f32, level: Option<f32>| {
            let mut st = strand(b"ACGT", hits(0, &[50., 48.]));
            let mut other = strand(b"ACGT", hits(1, &[40., 35.]));
            let mut reverse = false;
            rank_orientations(&mut st, &mut other, &mut reverse, lead, level);
            let sc = |s: &Strand| s.hits.iter().map(|h| (h.tx, h.score)).collect::<Vec<_>>();
            // The other's hits are marked with the other orientation.
            assert!(other.hits.iter().all(|h| h.reverse == !reverse));
            (sc(&st), sc(&other), reverse)
        };
        // The alignments favour the orientation seed luck left 10 points lower, by 6: it comes
        // first, at its own chain score or, for a mirrored pair, at the better one.
        let (st, other, reverse) = rank(-6.0, None);
        assert!(reverse);
        assert_eq!(st, [(1, 40.), (1, 35.)]);
        assert_eq!(other, [(0, 34.), (0, 32.)]);
        let (st, other, _) = rank(-6.0, Some(50.));
        assert_eq!(st, [(1, 50.), (1, 45.)]);
        assert_eq!(other, [(0, 44.), (0, 42.)]);
        // A lead for the first orientation keeps the order; the level is then its own.
        let (st, other, reverse) = rank(3.0, Some(50.));
        assert!(!reverse);
        assert_eq!(st, [(0, 50.), (0, 48.)]);
        assert_eq!(other, [(1, 47.), (1, 42.)]);
    }

    #[test]
    fn the_quantification_score_takes_the_orientation_scores_the_rules_would_give() {
        // `ranked_tops` gives each orientation the best score `rank_orientations` leaves it
        // (in the order given), so the quantification score can carry what the orientation
        // rules would move without moving the mapping scores.
        for lead in [-31.0, -6.0, 0.0, 3.0, 30.0] {
            for level in [None, Some(40.), Some(50.), Some(60.)] {
                let mut st = strand(b"ACGT", vec![hit(0, 50., &[(0, 0, 10)])]);
                let mut other = strand(b"ACGT", vec![hit(1, 40., &[(0, 0, 10)])]);
                let mut reverse = false;
                rank_orientations(&mut st, &mut other, &mut reverse, lead, level);
                // `st` and `other` swap with the order; tx 0 is the first given.
                let top = |tx: u32| [&st, &other].iter().find(|s| s.hits[0].tx == tx).unwrap().hits[0].score;
                assert_eq!(
                    ranked_tops((50., 40.), lead, level),
                    (top(0), top(1)),
                    "{lead} {level:?}"
                );
            }
        }
        // A mirrored pair the rules part by 16 points (own bases), against a tie: the leading
        // orientation's quantification scores stay, the other's fall by the lead.
        let (a, b) = ranked_tops((50., 46.), 16.0, Some(50.));
        let (a0, b0) = ranked_tops((50., 46.), 0.0, None);
        assert_eq!((a - a0, b - b0), (0.0, -16.0));
        // The rules reverse the order (a tail after a long adapter, 30 points, against the
        // other orientation's 6-point alignment lead).
        let (a, b) = ranked_tops((50., 48.), -24.0, None);
        let (a0, b0) = ranked_tops((50., 48.), 6.0, None);
        assert_eq!((a - a0, b - b0), (-26.0, 4.0));
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
            end_bonus_cut: false,
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
