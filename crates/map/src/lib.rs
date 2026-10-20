//! Read mapping: seeding both strands, anchor lookup, candidate-gene selection, chaining, and
//! result extraction (design doc, "Read processing" and "Chaining").

pub mod align;
pub mod anchor;
pub mod chain;
pub mod refine;
pub mod rescue;
mod segmented;
pub mod stats;

use std::sync::Arc;

use anyhow::{Context, Result};
use index::Index;
use lcp_core::encode::revcomp;
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
    /// Whole-read distances cover only the read between the tied hits' outermost anchors.
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
    /// followed by at most `TAIL_AFTER` read bases.
    pub tail_run: usize,
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
            group_fallback: false,
            whole_trim: true,
            orient_clipped: true,
            raw_ends: true,
            tail_bonus: 10,
            tail_trim: true,
            tail_run: 15,
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
            pi_cov: 0.5,
            strands: Strands::Both,
            merge_collinear: true,
            ends: Some(Ends::default()),
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
    pub anchors: u32,
    pub candidates: u32,
    /// The best chain spans less than `pi_cov` of the (parsed) read.
    pub low_coverage: bool,
    /// Transcripts re-ranked by the tie-break (0 if it didn't run).
    pub tiebreak_size: u32,
    /// The hits come from the hierarchical rescue.
    pub rescued: bool,
}

pub struct Mapper<C: Chainer> {
    pub idx: Arc<Index>,
    pub seeder: BaseSeeder,
    pub chainer: C,
    pub p: MapParams,
    rescuer: Option<Rescuer>,
    /// [`anchor::node_families`] of `idx`.
    node_fam: Vec<u32>,
    #[cfg(test)]
    reference_tiebreak: bool,
}

/// One orientation of a read: hits (sorted), the parsed sequence and its base seeds.
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
/// translations agree in every chain. Otherwise no valid anchored segmentation exists.
/// Adjacent intervals remain separate: a target insertion may lie between them.
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

struct Strand {
    hits: Vec<TxHit>,
    seq: Vec<u8>,
    raw: Option<RawRead>,
    seeds: Vec<seeds::Seed>,
    anchors: u32,
    candidates: u32,
}

impl<C: Chainer> Mapper<C> {
    /// A mapper seeding exactly as the index was built.
    pub fn new(idx: Arc<Index>, chainer: C, p: MapParams) -> Result<Self> {
        idx.check_compatible()?;
        let spec = idx.meta("spec").context("index lacks its seed spec")?;
        let hpc = idx.meta("hpc").context("index lacks its HPC flag")? == "true";
        let seeder = BaseSeeder::new(spec, hpc)?;
        let rescuer = p.rescue.clone().map(|r| Rescuer::new(&idx, r)).transpose()?;
        let node_fam = anchor::node_families(&idx);
        Ok(Mapper {
            #[cfg(test)]
            reference_tiebreak: false,
            node_fam,
            idx,
            seeder,
            chainer,
            p,
            rescuer,
        })
    }

    /// Re-rank the near-tie set of `hits` (sorted, best first) by edit distance between the
    /// parsed read `seq` and each transcript. Returns the tie-set size. `salt` orders equal
    /// scores (see [`sort_hits`]).
    /// `prune` is false when orientation reconciliation can still shift the reporting floor.
    #[allow(clippy::too_many_arguments)]
    fn tiebreak(
        &self,
        tb: &TieBreak,
        seq: &[u8],
        raw: Option<&RawRead>,
        hits: &mut [TxHit],
        salt: u64,
        prune: bool,
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
        let codes: Vec<u8> = seq
            .iter()
            .map(|&b| lcp_core::encode::CODE[b as usize] & 3)
            .collect();
        // Costs (edits, or negated affine scores) and the score charged per unit.
        let seg = if tb.segmented {
            if let Some(mut plan) = self.affine_plan(tb, &codes, raw, &hits[..n]) {
                let filter = prune.then_some(segmented::Retention {
                    baseline: best,
                    delta: self.p.delta_rep,
                    lambda: tb.lambda,
                });
                Some((plan.evaluate(tb, &hits[..n], filter), tb.lambda))
            } else {
                self.segmented_reference(tb, &codes, raw, &hits[..n])
                    .map(|(cost, unit)| (cost.into_iter().map(Some).collect(), unit))
            }
        } else {
            None
        };
        // Alignment deficit relative to the tie set's best alignment, per hit.
        let penalty: Vec<f32> = match seg {
            Some((cost, unit)) => {
                let cmin = *cost.iter().flatten().min().unwrap();
                // Missing costs are certified below the final reporting floor. Keep their
                // slots (and the original tie-set metadata) until the normal retain step.
                cost.iter()
                    .map(|&c| c.map_or(f32::INFINITY, |c| unit * (c - cmin) as f32))
                    .collect()
            }
            None => {
                stats::add(&stats::TIEBREAK_WHOLE_READ, 1);
                if tb.segmented && tb.group_fallback {
                    self.grouped(tb, &codes, raw, &hits[..n])
                } else {
                    let d = self.whole_read(tb, &codes, &hits[..n]);
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
                // Members share a complete interval; inconsistent placements may still fall back.
                if let Some((cost, unit)) = self.segmented(tb, codes, raw, &sub) {
                    let cmin = *cost.iter().min().unwrap();
                    for (k, &i) in members.iter().enumerate() {
                        penalty[i] = unit * (cost[k] - cmin) as f32;
                    }
                    rep = members[cost.iter().position(|&c| c == cmin).unwrap()];
                }
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
    /// later ones are capped at (best so far + `cap`) + 1.
    fn banded(read: &[u8], texts: &[Vec<u8>], mode: align::Mode, cap: u32) -> Vec<u32> {
        // Identical texts (isoforms agreeing over the region) are aligned once. The band only
        // tightens (k = best so far + cap), so a repeat's result is the first occurrence's
        // re-capped at the current k: exact if d <= k, else k + 1, as a fresh call returns.
        let pat = align::Pattern::new(read);
        if mode == align::Mode::Infix {
            return Self::banded_infix(&pat, read.len(), texts, cap);
        }
        let mut seen: rustc_hash::FxHashMap<&[u8], u32> = rustc_hash::FxHashMap::default();
        let mut best = u32::MAX;
        texts
            .iter()
            .map(|t| {
                let k = best.saturating_add(cap);
                if let Some(&d) = seen.get(t.as_slice()) {
                    return if d <= k { d } else { k.saturating_add(1) };
                }
                stats::add(&stats::MYERS_CALLS, 1);
                stats::add(
                    &stats::MYERS_CELLS,
                    (read.len() as u64).div_ceil(64) * t.len() as u64,
                );
                let d = align::distance_with(&pat, t, mode, k).unwrap_or(k.saturating_add(1));
                best = best.min(d);
                seen.insert(t.as_slice(), d);
                d
            })
            .collect()
    }

    /// [`Self::banded`] in infix mode, aligning distinct texts in batches of eight
    /// ([`align::infix_distances`]) with the threshold in force at each batch's first text:
    /// later texts' thresholds are no larger (the best only falls), so a batch result is exact
    /// wherever the sequential pass's would be, and re-capping it at each text's own threshold,
    /// in order, gives the sequential results exactly.
    fn banded_infix(pat: &align::Pattern, m: usize, texts: &[Vec<u8>], cap: u32) -> Vec<u32> {
        // Each text's distinct-text index, and the distinct texts in first-occurrence order.
        let mut first: rustc_hash::FxHashMap<&[u8], u32> = rustc_hash::FxHashMap::default();
        let mut uniq: Vec<&[u8]> = Vec::new();
        let mut at: Vec<usize> = Vec::new();
        let id: Vec<u32> = texts
            .iter()
            .enumerate()
            .map(|(i, t)| {
                *first.entry(t.as_slice()).or_insert_with(|| {
                    uniq.push(t.as_slice());
                    at.push(i);
                    uniq.len() as u32 - 1
                })
            })
            .collect();
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
                raw.extend(align::infix_distances(pat, batch, k));
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
        )
    }

    /// A hit's transcript region projected to the ends of the read span `[lo, hi)`, padded.
    fn hit_region(&self, tb: &TieBreak, lo: u32, hi: u32, h: &TxHit) -> Vec<u8> {
        let tlen = self.idx.tx_plen()[h.tx as usize];
        let (left, right) = (h.qstart - lo, hi - h.qend);
        let slack = |ext: u32| tb.pad + ext / tb.pad_div;
        let lo = h.tstart.saturating_sub(left + slack(left));
        let hi = (h.tend + right + slack(right)).min(tlen);
        self.idx.tx_codes(h.tx, lo, hi)
    }

    /// Affine score of the read against one hit with both read ends clipped: the chained span
    /// aligned globally, each read end extended from it into the padded transcript flank.
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
        let slack = |ext: u32| tb.pad + ext / tb.pad_div;
        let (qs, qe) = (h.qstart as usize, h.qend as usize);
        let w0 = align::dp_work();
        let mid = align::affine(
            &codes[qs..qe],
            &self.idx.tx_codes(h.tx, h.tstart, h.tend),
            align::Mode::Global,
            sc,
            tb.band,
        );
        let w1 = align::dp_work();
        let left = if qs > 0 {
            let span = [(
                h.tx,
                h.tstart.saturating_sub(h.qstart + slack(h.qstart)),
                h.tstart,
            )];
            let (read, text) = match self.raw_texts(tb, raw, &span) {
                Some((r, mut t)) => (&r.codes[..r.pos(h.qstart)], t.pop().unwrap()),
                None => (&codes[..qs], self.idx.tx_codes(h.tx, span[0].1, h.tstart)),
            };
            let read: Vec<u8> = read.iter().rev().copied().collect();
            let text: Vec<u8> = text.into_iter().rev().collect();
            align::affine_extend(&read, &text, sc, tb.band, tb.end_bonus, tb.xdrop)
        } else {
            tb.end_bonus
        };
        let rest = plen - h.qend;
        let right = if rest > 0 {
            let span = [(h.tx, h.tend, (h.tend + rest + slack(rest)).min(tlen))];
            let (read, text) = match self.raw_texts(tb, raw, &span) {
                Some((r, mut t)) => (&r.codes[r.pos(h.qend)..], t.pop().unwrap()),
                None => (&codes[qe..], self.idx.tx_codes(h.tx, h.tend, span[0].2)),
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
    /// excluded from every DP slice. `None` if no complete query interval is shared or if
    /// the placements cannot form consistent, ordered shared blocks.
    fn segmented(
        &self,
        tb: &TieBreak,
        codes: &[u8],
        raw: Option<&RawRead>,
        hits: &[TxHit],
    ) -> Option<(Vec<i64>, f32)> {
        if let Some(mut plan) = self.affine_plan(tb, codes, raw, hits) {
            let costs = plan
                .evaluate(tb, hits, None)
                .into_iter()
                .map(Option::unwrap)
                .collect();
            Some((costs, tb.lambda))
        } else {
            self.segmented_reference(tb, codes, raw, hits)
        }
    }

    // Original evaluator: differential oracle and fallback for unsupported scoring paths.
    fn segmented_reference(
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
        let seg = |texts: Vec<Vec<u8>>,
                   read: &[u8],
                   mode: align::Mode,
                   rev: bool,
                   cut: Option<usize>|
         -> Option<(Vec<i64>, Vec<i64>, Vec<usize>)> {
            if texts.windows(2).all(|w| w[0] == w[1]) {
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
            let (read, texts): (Vec<u8>, Vec<Vec<u8>>) = if rev {
                (
                    read.iter().rev().copied().collect(),
                    texts.into_iter().map(|t| t.into_iter().rev().collect()).collect(),
                )
            } else {
                (read.to_vec(), texts)
            };
            match &tb.affine {
                Some(sc) => {
                    // Tied isoforms often share a segment's sequence: score each distinct text
                    // once. With `affine_within`, a capped unit-cost (Myers) pass comes first
                    // and only texts within that many edits of the best get an affine score;
                    // the rest are charged one mismatch (mat + mis) per extra edit.
                    let mut order: Vec<usize> = (0..texts.len()).collect();
                    order.sort_unstable_by(|&a, &b| texts[a].cmp(&texts[b]));
                    let mut uniq: Vec<usize> = Vec::with_capacity(order.len());
                    let mut rep = vec![0usize; texts.len()];
                    for &i in &order {
                        if uniq.last().is_none_or(|&u| texts[u] != texts[i]) {
                            uniq.push(i);
                        }
                        rep[i] = uniq.len() - 1;
                    }
                    let unit = match tb.affine_within {
                        Some(w) if uniq.len() > 1 => {
                            let ut: Vec<Vec<u8>> = uniq.iter().map(|&u| texts[u].clone()).collect();
                            Some((Self::banded(&read, &ut, mode, w), w))
                        }
                        _ => None,
                    };
                    let dmin = unit.as_ref().map_or(0, |(d, _)| *d.iter().min().unwrap());
                    if stats::on() {
                        // Text bases, and those outside the common prefix with the previous
                        // distinct text in sorted order (what prefix sharing would compute).
                        let mut prev: Option<&Vec<u8>> = None;
                        for &u in &uniq {
                            let t = &texts[u];
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
                        let (m, n) = (read.len() as u64, texts[uniq[k]].len() as u64);
                        stats::add(&stats::AFFINE_CALLS, 1);
                        let drift = if diag && n >= m { 0 } else { m.abs_diff(n) };
                        let width = (tb.band as u64)
                            .saturating_mul(2)
                            .saturating_add(drift)
                            .saturating_add(1);
                        let cells = m.saturating_mul(width.min(n + 1));
                        stats::add(&stats::AFFINE_CELLS, cells);
                        if mode == align::Mode::Prefix {
                            stats::add(&stats::AFFINE_END_CELLS, cells);
                        }
                    }
                    let tx: Vec<&[u8]> = todo.iter().map(|&k| texts[uniq[k]].as_slice()).collect();
                    let w0 = align::dp_work();
                    let xs = if mode == align::Mode::Prefix && tb.clip_ends {
                        // The band is laid out for the window the read end would get against a
                        // transcript that doesn't end inside it.
                        let window = read.len() + tb.pad as usize + read.len() / tb.pad_div as usize;
                        align::affine_extend_sorted_cut(
                            &read,
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
                        align::affine_sorted(&read, &tx, mode, sc, tb.band)
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
                    for (c, x) in delta.iter_mut().zip(Self::banded(&read, &texts, mode, tb.cap)) {
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
        let end = |texts: Vec<Vec<u8>>, read: &[u8], rev: bool, short: &[bool]| -> Option<Vec<i64>> {
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
        // Strict interior gaps [previous anchor end, next anchor start).
        for pair in blocks.windows(2) {
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
                        .tx_codes(h.tx, left.targets[i].end, right.targets[i].start)
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
        let q0 = first.query.start;
        if q0 > 0 {
            let spans: Vec<(u32, u32, u32)> = hits
                .iter()
                .zip(&first.targets)
                .map(|(h, t)| {
                    (
                        h.tx,
                        t.start.saturating_sub(q0 + tb.pad + q0 / tb.pad_div),
                        t.start,
                    )
                })
                .collect();
            let short: Vec<bool> = spans.iter().map(|&(_, lo, _)| lo == 0).collect();
            add(match self.raw_texts(tb, raw, &spans) {
                Some((r, texts)) => end(texts, &r.codes[..r.pos(q0)], true, &short),
                None => end(
                    spans
                        .iter()
                        .map(|&(t, lo, hi)| self.idx.tx_codes(t, lo, hi))
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
            && hits.iter().zip(&last.targets).skip(1).any(|(h, t)| {
                self.idx
                    .tx_raw_eq(hits[0].tx, t0.start, t0.end, h.tx, t.start, t.end)
                    == Some(false)
            });
        let ql = if through { last.query.start } else { last.query.end };
        let rest = plen - ql;
        let spans: Vec<(u32, u32, u32)> = hits
            .iter()
            .zip(&last.targets)
            .map(|(h, t)| {
                let tlen = self.idx.tx_plen()[h.tx as usize];
                let s = if through { t.start } else { t.end };
                (h.tx, s, (s + rest + tb.pad + rest / tb.pad_div).min(tlen))
            })
            .collect();
        // The read's poly(A) tail, searched for past the last block (raw positions from ql).
        let a_code = lcp_core::encode::CODE[b'A' as usize] & 3;
        // Whether each window's end (its length from the read end's start, in the tail's space)
        // is where the tail starts: within the tolerance of its start, or further but only by
        // annotated A's (an end that includes part of the tail) up to the run's end.
        let ends_at_tail = |texts: &[Vec<u8>], cut: usize, run_end: usize| -> Vec<bool> {
            let tol = TAIL_TOL + cut / 10;
            texts
                .iter()
                .map(|t| {
                    let n = t.len();
                    n.abs_diff(cut) <= tol
                        || (n > cut
                            && n <= run_end + tol
                            && t[cut..].iter().filter(|&&c| c != a_code).count() * 10 <= n - cut)
                })
                .collect()
        };
        let tail = match raw {
            Some(r) if tb.tail_trim || tb.tail_bonus > 0 => {
                let (r0, r1) = (r.pos(ql), r.pos(last.query.end));
                tail_start(&r.codes[r1..], a_code, tb.tail_run)
                    .map(|(c, e)| (r, r0, r1 - r0 + c, r1 - r0 + e))
            }
            _ => None,
        };
        let short: Vec<bool> = spans
            .iter()
            .map(|&(tx, _, hi)| hi == self.idx.tx_plen()[tx as usize])
            .collect();
        let ends = match self.raw_texts(tb, raw, &spans) {
            Some((r, texts)) => {
                let rest = &r.codes[r.pos(ql)..];
                let cut = tail.map_or(rest.len(), |(_, _, c, _)| c);
                let run_end = tail.map_or(cut, |(_, _, _, e)| e);
                let at_end = tail.map(|_| ends_at_tail(&texts, cut, run_end));
                let read = if tb.tail_trim { &rest[..cut] } else { rest };
                add(end(texts, read, false, &short));
                at_end
            }
            None => {
                let texts: Vec<Vec<u8>> = spans
                    .iter()
                    .map(|&(t, lo, hi)| self.idx.tx_codes(t, lo, hi))
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
                let cut = tail.map_or((plen - ql) as usize, |(r, r0, c, _)| parsed(r0 + c, r));
                let run_end = tail.map_or(cut, |(r, r0, _, e)| parsed(r0 + e, r));
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

    fn map_strand(&self, (seq, raw_read): (Vec<u8>, Option<RawRead>), salt: u64) -> Strand {
        let seeds = self.seeder.seed_parsed(&seq);
        let (mut cands, n_anchors) = anchor::select(&self.idx, &self.node_fam, &seeds, &self.p.anchor);
        if self.p.merge_collinear {
            for c in &mut cands {
                anchor::merge_collinear(&self.idx, c);
            }
        }
        let mut hits: Vec<TxHit> = cands
            .iter()
            .flat_map(|c| self.chainer.chain(&self.idx, c))
            .collect();
        sort_hits(&mut hits, salt);
        Strand {
            hits,
            seq,
            raw: raw_read,
            seeds,
            anchors: n_anchors,
            candidates: cands.len() as u32,
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
                let codes: Vec<u8> = raw
                    .iter()
                    .map(|&b| lcp_core::encode::CODE[b as usize] & 3)
                    .collect();
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
                let codes = raw
                    .iter()
                    .map(|&b| lcp_core::encode::CODE[b as usize] & 3)
                    .collect();
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
        let skip = |p: (Vec<u8>, Option<RawRead>)| Strand {
            hits: Vec::new(),
            seq: p.0,
            raw: None,
            seeds: Vec::new(),
            anchors: 0,
            candidates: 0,
        };
        // Equal scores are ordered by a hash of the transcript salted with the read, so ties
        // are broken the same way on every run without favouring the same transcript.
        let salt = xxhash_rust::xxh3::xxh3_64(raw);
        let [pf, pr] = Self::parse_both(&self.seeder, raw);
        let (fwd, rc) = match self.p.strands {
            Strands::Both => (self.map_strand(pf, salt), self.map_strand(pr, salt)),
            Strands::Forward => (self.map_strand(pf, salt), skip(pr)),
            Strands::Reverse => (skip(pf), self.map_strand(pr, salt)),
        };
        let best = |h: &[TxHit]| h.first().map_or(f32::NEG_INFINITY, |x| x.score);
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
            let mut a = r.rescue(&self.idx, &st.seq, &st.seeds, &self.p.anchor, &self.p.chain);
            // The other orientation is re-seeded only if its base-level best comes within
            // `other_margin` of this one's (None: always).
            let try_other = match r.p.other_margin {
                None => true,
                Some(m) => other.hits.first().is_some_and(|o| o.score >= best(&st.hits) - m),
            };
            let mut b = if try_other {
                r.rescue(&self.idx, &other.seq, &other.seeds, &self.p.anchor, &self.p.chain)
            } else {
                Vec::new()
            };
            sort_hits(&mut a, salt);
            sort_hits(&mut b, salt);
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
        let plen = st.seq.len() as u32;
        let mut tiebreak_size = 0;
        let mut both = false;
        if let Some(tb) = &self.p.tiebreak {
            let margin = |s: f32| tb.margin_abs.max(tb.margin_rel * s.abs());
            both = matches!((st.hits.first(), other.hits.first()),
                    (Some(a), Some(b)) if b.score >= a.score - margin(a.score));
            tiebreak_size = self.tiebreak(tb, &st.seq, st.raw.as_ref(), &mut st.hits, salt, !both);
            if both {
                // Near-tied orientations (antisense transcripts sharing sequence): re-rank each
                // orientation, then place the other's scores relative to this one by the
                // whole-read edit distances of the two best transcripts.
                stats::add(&stats::BOTH_ORIENTATIONS, 1);
                tiebreak_size +=
                    self.tiebreak(tb, &other.seq, other.raw.as_ref(), &mut other.hits, salt, false);
                let codes = |s: &[u8]| -> Vec<u8> {
                    s.iter()
                        .map(|&c| lcp_core::encode::CODE[c as usize] & 3)
                        .collect()
                };
                // Costs of the two best hits and the score charged per unit.
                let (ca, cb, unit) = match &tb.affine {
                    Some(sc) if tb.orient_clipped => (
                        -i64::from(self.clipped(tb, sc, &codes(&st.seq), st.raw.as_ref(), &st.hits[0])),
                        -i64::from(self.clipped(
                            tb,
                            sc,
                            &codes(&other.seq),
                            other.raw.as_ref(),
                            &other.hits[0],
                        )),
                        tb.lambda,
                    ),
                    _ => {
                        let da = self.whole_read(tb, &codes(&st.seq), &st.hits[..1])[0];
                        let db = self.whole_read_k(
                            tb,
                            &codes(&other.seq),
                            &other.hits[0],
                            da.saturating_add(tb.cap),
                        );
                        (i64::from(da), i64::from(db), tb.kappa)
                    }
                };
                if cb < ca {
                    std::mem::swap(&mut st, &mut other);
                    reverse = !reverse;
                }
                let target = st.hits[0].score - unit * ca.abs_diff(cb) as f32;
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
        let mut low_coverage = true;
        if let Some(top) = hits.first() {
            let floor = top.score - self.p.delta_rep * top.score.abs();
            low_coverage = ((top.qend - top.qstart) as f32) < self.p.pi_cov * plen as f32;
            hits.retain(|h| h.score >= floor);
        }
        if let Some(e) = &self.p.ends {
            let codes = |s: &[u8]| -> Vec<u8> {
                s.iter()
                    .map(|&c| lcp_core::encode::CODE[c as usize] & 3)
                    .collect()
            };
            let this = codes(&st.seq);
            let that = if both { codes(&other.seq) } else { Vec::new() };
            let mut seen = rustc_hash::FxHashMap::default();
            let floor = hits.first().map_or(0.0, |t| t.score - e.within * t.score.abs());
            for h in hits.iter_mut().filter(|h| h.score >= floor) {
                let c = if h.reverse == reverse { &this } else { &that };
                self.extend_ends(e, c, h, &mut seen);
                h.extended = true;
            }
        }
        ReadResult {
            len: raw.len() as u32,
            plen,
            reverse,
            hits,
            anchors: st.anchors,
            candidates: st.candidates,
            low_coverage,
            tiebreak_size,
            rescued,
        }
    }
}

/// Bases between a transcript's 3' end and a read's tail start treated as the same end.
const TAIL_TOL: usize = 15;
/// Most read bases after a poly(A) tail (adapter and the like).
const TAIL_AFTER: usize = 80;

/// The first poly(A) window of `run` bases in `rest` (raw codes): its first A, and the end of the A-run from
/// there (the longest stretch ending in A with at most one other base per ten).
fn tail_start(rest: &[u8], a: u8, run: usize) -> Option<(usize, usize)> {
    if rest.len() < run {
        return None;
    }
    let mut other = rest[..run].iter().filter(|&&c| c != a).count();
    for i in 0..=rest.len() - run {
        if i > 0 {
            other = other + usize::from(rest[i + run - 1] != a) - usize::from(rest[i - 1] != a);
        }
        if other <= 1 {
            let s = i + rest[i..].iter().position(|&c| c == a).unwrap_or(0);
            let (mut end, mut non) = (s, 0);
            for (k, &c) in rest.iter().enumerate().skip(s) {
                non += usize::from(c != a);
                if c == a && non * 10 <= k + 1 - s {
                    end = k + 1;
                }
            }
            // A tail is followed only by adapter: an A-run with more of the read after it is
            // transcript sequence (or internal priming), not the molecule's 3' end.
            if rest.len() - end > TAIL_AFTER {
                return None;
            }
            return Some((s, end));
        }
    }
    None
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

/// splitmix64 of the read salt and a transcript id.
fn tie_key(salt: u64, tx: u32) -> u64 {
    let mut z = salt ^ u64::from(tx).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod segmented_tests {
    use super::*;

    pub(super) fn hit(tx: u32, score: f32, chain: &[(u32, u32, u32)]) -> TxHit {
        chain::hit_from_path(tx, score, 0.0, chain.iter().map(|&(q, t, l)| (q, q + l, t, l)))
    }

    pub(super) fn mapper(texts: &[&[u8]], hpc: bool) -> Mapper<OracleChainer> {
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
            &txs,
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
            Arc::new(Index::from_bytes(&data.to_bytes()).unwrap()),
            OracleChainer {
                p: Default::default(),
            },
            MapParams::default(),
        )
        .unwrap()
    }

    pub(super) fn codes(s: &[u8]) -> Vec<u8> {
        s.iter()
            .map(|&b| lcp_core::encode::CODE[b as usize] & 3)
            .collect()
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
        assert_eq!(m.tiebreak(&tb, q, None, &mut hits, 0, false), 4);
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
    fn poly_a_tail_is_trimmed_and_favours_the_isoform_ending_there() {
        // t1 continues past t0's 3' end; the read stops there and goes on into its tail.
        let body: &[u8] = b"ACGTCGATCGTACGATCGTCAGCTAGCTGACTGCATGCAT";
        let t1: Vec<u8> = [body, b"GCTCGTCGTCAGCTCGTCGATCGCT"].concat();
        let m = mapper(&[body, &t1], true);
        let q: Vec<u8> = [body, b"AAAAAAAAAAAAAAAAAAAAAAAA", b"GTCGCTCGTGCA"].concat();
        let (seq, raw) = Mapper::<OracleChainer>::parse(&m.seeder, &q);
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
        let mut hits = [
            hit(0, 15., &[(0, 0, 2), (4, 4, 2)]),
            hit(1, 12., &[(0, 0, 2), (4, 2, 2)]),
        ];
        let tb = TieBreak {
            lambda: 0.5,
            ..Default::default()
        };
        m.tiebreak(&tb, q, None, &mut hits, 0, false);
        // Scores become the tie set's best chain score minus each deficit.
        assert_eq!((hits[0].score, hits[1].score), (15., 9.));
        let m = mapper(&[q, q], false);
        for grouped in [false, true] {
            let mut hits = [hit(0, 15., &[(0, 0, 2)]), hit(1, 12., &[(3, 3, 2)])];
            let tb = TieBreak {
                group_fallback: grouped,
                ..Default::default()
            };
            assert!(m.segmented(&tb, &codes(q), None, &hits).is_none());
            m.tiebreak(&tb, q, None, &mut hits, 0, false);
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
}
