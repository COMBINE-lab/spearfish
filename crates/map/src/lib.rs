//! Read mapping: seeding both strands, anchor lookup, candidate-gene selection, chaining, and
//! result extraction (design doc, "Read processing" and "Chaining").

pub mod align;
pub mod anchor;
pub mod chain;
pub mod refine;
pub mod rescue;
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
#[derive(Clone, Debug)]
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
    /// whole read across group representatives) instead of aligning every hit's whole read.
    pub group_fallback: bool,
    /// With `clip_ends`: stop extending a read end once a row falls this far below the best
    /// score so far (X-drop; `i32::MAX`: never).
    pub xdrop: i32,
    /// Near-tied read orientations: compare their best hits by clipped affine alignment
    /// scores (read ends past the transcript are clipped) instead of whole-read edit distance,
    /// which charges a short transcript for the read's adapters.
    pub orient_clipped: bool,
    /// HPC indexes with raw transcripts: align read ends in raw (uncompressed) space, where
    /// homopolymers at the read ends (e.g. a template-switch oligo's Gs) can't merge into the
    /// transcript's first or last run.
    pub raw_ends: bool,
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
            group_fallback: false,
            orient_clipped: true,
            raw_ends: true,
        }
    }
}

#[derive(Clone, Debug)]
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
#[derive(Clone, Debug)]
pub struct Ends {
    pub sc: align::Affine,
    pub band: usize,
    pub end_bonus: i32,
    pub xdrop: i32,
    /// Transcript flank searched beyond the read flank's length: `pad + flank / pad_div`.
    pub pad: u32,
    pub pad_div: u32,
}

/// End extensions of one read: (left side, reverse, read position, transcript flank) to the
/// read and transcript bases they cover.
type EndCache = rustc_hash::FxHashMap<(bool, bool, u32, Vec<u8>), (u32, u32)>;

impl Default for Ends {
    fn default() -> Self {
        Ends {
            sc: align::Affine::ONT,
            band: 16,
            end_bonus: 5,
            xdrop: 50,
            pad: 20,
            pad_div: 10,
        }
    }
}

/// Which orientations of a read to map against the transcripts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strands {
    Both,
    /// The read as given (sense: direct RNA).
    Forward,
    /// Only its reverse complement.
    Reverse,
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
        let spec = idx.meta("spec").context("index lacks its seed spec")?;
        let hpc = idx.meta("hpc").context("index lacks its HPC flag")? == "true";
        let seeder = BaseSeeder::new(spec, hpc)?;
        let rescuer = p.rescue.clone().map(|r| Rescuer::new(&idx, r)).transpose()?;
        Ok(Mapper {
            idx,
            seeder,
            chainer,
            p,
            rescuer,
        })
    }

    /// Re-rank the near-tie set of `hits` (sorted, best first) by edit distance between the
    /// parsed read `seq` and each transcript. Returns the tie-set size.
    fn tiebreak(&self, tb: &TieBreak, seq: &[u8], raw: Option<&RawRead>, hits: &mut [TxHit]) -> u32 {
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
            self.segmented(tb, &codes, raw, &hits[..n])
        } else {
            None
        };
        // Score below the tie set's best, per hit.
        let penalty: Vec<f32> = match seg {
            Some((cost, unit)) => {
                let cmin = *cost.iter().min().unwrap();
                cost.iter().map(|&c| unit * (c - cmin) as f32).collect()
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
        hits.sort_by(|a, b| b.score.total_cmp(&a.score).then(a.tx.cmp(&b.tx)));
        n as u32
    }

    /// Penalties for a tie set whose hits share no common anchor (typically paralogs): the
    /// hits are grouped greedily so each group shares an anchor, scored by [`Self::segmented`]
    /// within groups, and the groups are placed relative to each other by the whole-read edit
    /// distance of each group's best hit (`kappa` per edit), as the two read orientations are.
    /// Aligns one whole-read region per group instead of one per hit.
    fn grouped(&self, tb: &TieBreak, codes: &[u8], raw: Option<&RawRead>, hits: &[TxHit]) -> Vec<f32> {
        let mut groups: Vec<(rustc_hash::FxHashSet<u32>, Vec<usize>)> = Vec::new();
        for (i, h) in hits.iter().enumerate() {
            let qs: rustc_hash::FxHashSet<u32> = h.chain.iter().map(|a| a.0).collect();
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
                // Members share an anchor by construction, so this is always Some.
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
            let rep_hits: Vec<TxHit> = reps.iter().map(|&i| hits[i].clone()).collect();
            let d = self.whole_read(tb, codes, &rep_hits);
            let dmin = *d.iter().min().unwrap();
            for ((_, members), &dg) in groups.iter().zip(&d) {
                let off = tb.kappa * (dg - dmin) as f32;
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
                let d = align::distance(read, t, mode, k).unwrap_or(k.saturating_add(1));
                best = best.min(d);
                seen.insert(t.as_slice(), d);
                d
            })
            .collect()
    }

    /// Whole read against each candidate's covered transcript region (projected to the read
    /// ends, padded), infix mode.
    fn whole_read(&self, tb: &TieBreak, codes: &[u8], hits: &[TxHit]) -> Vec<u32> {
        let plen = codes.len() as u32;
        let texts: Vec<Vec<u8>> = hits.iter().map(|h| self.hit_region(tb, plen, h)).collect();
        Self::banded(codes, &texts, align::Mode::Infix, tb.cap)
    }

    /// A hit's transcript region projected to the read ends (read length `plen`), padded.
    fn hit_region(&self, tb: &TieBreak, plen: u32, h: &TxHit) -> Vec<u8> {
        let tlen = self.idx.tx_plen()[h.tx as usize];
        let (left, right) = (h.qstart, plen - h.qend);
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
        let mid = align::affine(
            &codes[qs..qe],
            &self.idx.tx_codes(h.tx, h.tstart, h.tend),
            align::Mode::Global,
            sc,
            tb.band,
        );
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
        mid + left + right
    }

    /// Whole-read infix distance against one hit's region, capped: `k + 1` if above `k`.
    fn whole_read_k(&self, tb: &TieBreak, codes: &[u8], h: &TxHit, k: u32) -> u32 {
        let text = self.hit_region(tb, codes.len() as u32, h);
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
    /// aligned against windows anchored at the outermost shared anchors. `None` if the
    /// candidates share no anchor.
    fn segmented(
        &self,
        tb: &TieBreak,
        codes: &[u8],
        raw: Option<&RawRead>,
        hits: &[TxHit],
    ) -> Option<(Vec<i64>, f32)> {
        // Query starts present in every chain, with each candidate's transcript start.
        let mut shared: Vec<u32> = hits[0].chain.iter().map(|a| a.0).collect();
        for h in &hits[1..] {
            let qs: rustc_hash::FxHashSet<u32> = h.chain.iter().map(|a| a.0).collect();
            shared.retain(|q| qs.contains(q));
        }
        if shared.is_empty() {
            return None;
        }
        let tpos: Vec<Vec<u32>> = hits
            .iter()
            .map(|h| {
                let mut it = h.chain.iter().peekable();
                shared
                    .iter()
                    .map(|&q| {
                        while it.peek().is_some_and(|a| a.0 < q) {
                            it.next();
                        }
                        it.peek().unwrap().1
                    })
                    .collect()
            })
            .collect();
        let plen = codes.len() as u32;
        let mut cost = vec![0i64; hits.len()];
        let mut add = |texts: Vec<Vec<u8>>, read: &[u8], mode: align::Mode, rev: bool| {
            if texts.windows(2).all(|w| w[0] == w[1]) {
                return; // same cost for every candidate
            }
            stats::add(&stats::TIEBREAK_SEGMENTS, 1);
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
                    let mut best: Option<(u32, i64)> = None; // (unit distance, affine score)
                    for (k, &u) in uniq.iter().enumerate() {
                        if let Some((d, w)) = &unit
                            && d[k] > dmin + w
                        {
                            continue;
                        }
                        let (m, n) = (read.len() as u64, texts[u].len() as u64);
                        stats::add(&stats::AFFINE_CALLS, 1);
                        let cells = m * (2 * tb.band as u64 + m.abs_diff(n) + 1).min(n + 1);
                        stats::add(&stats::AFFINE_CELLS, cells);
                        if mode == align::Mode::Prefix {
                            stats::add(&stats::AFFINE_END_CELLS, cells);
                        }
                        let x = i64::from(if mode == align::Mode::Prefix && tb.clip_ends {
                            align::affine_extend(&read, &texts[u], sc, tb.band, tb.end_bonus, tb.xdrop)
                        } else {
                            align::affine(&read, &texts[u], mode, sc, tb.band)
                        });
                        score[k] = Some(x);
                        let du = unit.as_ref().map_or(0, |(d, _)| d[k]);
                        if best.is_none_or(|(bd, bx)| (du, -x) < (bd, -bx)) {
                            best = Some((du, x));
                        }
                    }
                    let per_edit = i64::from(sc.mat + sc.mis);
                    for (i, c) in cost.iter_mut().enumerate() {
                        let k = rep[i];
                        let x = score[k].unwrap_or_else(|| {
                            let (bd, bx) = best.unwrap();
                            let d = unit.as_ref().unwrap().0[k];
                            bx - per_edit * i64::from(d - bd)
                        });
                        *c -= x;
                    }
                }
                None => {
                    for (c, x) in cost.iter_mut().zip(Self::banded(&read, &texts, mode, tb.cap)) {
                        *c += i64::from(x);
                    }
                }
            }
        };
        // Interior segments [shared_k, shared_k+1).
        for k in 0..shared.len() - 1 {
            let texts = hits
                .iter()
                .zip(&tpos)
                .map(|(h, t)| self.idx.tx_codes(h.tx, t[k], t[k + 1].max(t[k])))
                .collect();
            add(
                texts,
                &codes[shared[k] as usize..shared[k + 1] as usize],
                align::Mode::Global,
                false,
            );
        }
        // Left end: read[0, q0) against transcript windows ending at each candidate's anchor.
        let q0 = shared[0];
        if q0 > 0 {
            let spans: Vec<(u32, u32, u32)> = hits
                .iter()
                .zip(&tpos)
                .map(|(h, t)| (h.tx, t[0].saturating_sub(q0 + tb.pad + q0 / tb.pad_div), t[0]))
                .collect();
            match self.raw_texts(tb, raw, &spans) {
                Some((r, texts)) => add(texts, &r.codes[..r.pos(q0)], align::Mode::Prefix, true),
                None => add(
                    spans
                        .iter()
                        .map(|&(t, lo, hi)| self.idx.tx_codes(t, lo, hi))
                        .collect(),
                    &codes[..q0 as usize],
                    align::Mode::Prefix,
                    true,
                ),
            }
        }
        // Right end: read[q_last, plen) against windows starting at the last shared anchor.
        let ql = *shared.last().unwrap();
        let rest = plen - ql;
        let spans: Vec<(u32, u32, u32)> = hits
            .iter()
            .zip(&tpos)
            .map(|(h, t)| {
                let tlen = self.idx.tx_plen()[h.tx as usize];
                let s = *t.last().unwrap();
                (h.tx, s, (s + rest + tb.pad + rest / tb.pad_div).min(tlen))
            })
            .collect();
        match self.raw_texts(tb, raw, &spans) {
            Some((r, texts)) => add(texts, &r.codes[r.pos(ql)..], align::Mode::Prefix, false),
            None => add(
                spans
                    .iter()
                    .map(|&(t, lo, hi)| self.idx.tx_codes(t, lo, hi))
                    .collect(),
                &codes[ql as usize..],
                align::Mode::Prefix,
                false,
            ),
        }
        let unit = if tb.affine.is_some() { tb.lambda } else { tb.kappa };
        Some((cost, unit))
    }

    fn map_strand(&self, raw: &[u8]) -> Strand {
        let (seq, raw_read) = Self::parse(&self.seeder, raw);
        let seeds = self.seeder.seed_parsed(&seq);
        let anchors = anchor::anchors(&self.idx, &seeds, &self.p.anchor);
        let n_anchors = anchors.len() as u32;
        let mut cands = anchor::candidates(&self.idx, anchors, &self.p.anchor);
        if self.p.merge_collinear {
            for c in &mut cands {
                anchor::merge_collinear(&self.idx, c);
            }
        }
        let mut hits: Vec<TxHit> = cands
            .iter()
            .flat_map(|c| self.chainer.chain(&self.idx, c))
            .collect();
        sort_hits(&mut hits);
        Strand {
            hits,
            seq,
            raw: raw_read,
            seeds,
            anchors: n_anchors,
            candidates: cands.len() as u32,
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
    /// parsed read in `h`'s orientation). Isoforms sharing a flank's sequence share its
    /// extension, so `seen` caches extensions by (side, orientation, read position, text).
    fn extend_ends(&self, e: &Ends, codes: &[u8], h: &mut TxHit, seen: &mut EndCache) {
        let plen = codes.len() as u32;
        let tlen = self.idx.tx_plen()[h.tx as usize];
        let slack = |x: u32| e.pad + x / e.pad_div;
        let mut ext = |left: bool, q: u32, read: Vec<u8>, text: Vec<u8>| -> (u32, u32) {
            stats::add(&stats::SPAN_ENDS, 1);
            *seen.entry((left, h.reverse, q, text)).or_insert_with_key(|k| {
                let (m, n) = (read.len() as u64, k.3.len() as u64);
                stats::add(&stats::SPAN_END_CALLS, 1);
                stats::add(
                    &stats::SPAN_END_CELLS,
                    m * (2 * e.band as u64 + m.abs_diff(n) + 1).min(n + 1),
                );
                let (_, qi, tj) = align::affine_extend_end(&read, &k.3, &e.sc, e.band, e.end_bonus, e.xdrop);
                (qi as u32, tj as u32)
            })
        };
        if h.qstart > 0 && h.tstart > 0 {
            let lo = h.tstart.saturating_sub(h.qstart + slack(h.qstart));
            let mut text = self.idx.tx_codes(h.tx, lo, h.tstart);
            text.reverse();
            let read = codes[..h.qstart as usize].iter().rev().copied().collect();
            let (qi, tj) = ext(true, h.qstart, read, text);
            h.qstart -= qi;
            h.tstart -= tj;
        }
        let rest = plen - h.qend;
        if rest > 0 && h.tend < tlen {
            let hi = (h.tend + rest + slack(rest)).min(tlen);
            let text = self.idx.tx_codes(h.tx, h.tend, hi);
            let (qi, tj) = ext(false, h.qend, codes[h.qend as usize..].to_vec(), text);
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
        let skip = |s: &Self, r: &[u8]| Strand {
            hits: Vec::new(),
            seq: Self::parse(&s.seeder, r).0,
            raw: None,
            seeds: Vec::new(),
            anchors: 0,
            candidates: 0,
        };
        let (fwd, rc) = match self.p.strands {
            Strands::Both => (self.map_strand(raw), self.map_strand(&revcomp(raw))),
            Strands::Forward => (self.map_strand(raw), skip(self, &revcomp(raw))),
            Strands::Reverse => (skip(self, raw), self.map_strand(&revcomp(raw))),
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
            sort_hits(&mut a);
            sort_hits(&mut b);
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
            tiebreak_size = self.tiebreak(tb, &st.seq, st.raw.as_ref(), &mut st.hits);
            if both {
                // Near-tied orientations (antisense transcripts sharing sequence): re-rank each
                // orientation, then place the other's scores relative to this one by the
                // whole-read edit distances of the two best transcripts.
                stats::add(&stats::BOTH_ORIENTATIONS, 1);
                tiebreak_size += self.tiebreak(tb, &other.seq, other.raw.as_ref(), &mut other.hits);
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
            sort_hits(&mut hits);
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
            for h in &mut hits {
                let c = if h.reverse == reverse { &this } else { &that };
                self.extend_ends(e, c, h, &mut seen);
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

fn sort_hits(hits: &mut [TxHit]) {
    hits.sort_unstable_by(|a, b| b.score.total_cmp(&a.score).then(a.tx.cmp(&b.tx)));
}
