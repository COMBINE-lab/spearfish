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
    /// With `clip_ends`: stop extending a read end once a row falls this far below the best
    /// score so far (X-drop; `i32::MAX`: never).
    pub xdrop: i32,
}

impl Default for TieBreak {
    fn default() -> Self {
        TieBreak {
            margin_abs: 30.0,
            margin_rel: 0.05,
            kappa: 6.0,
            pad: 50,
            pad_div: 6,
            max_cands: 256,
            cap: 16,
            segmented: true,
            affine: Some(align::Affine::default()),
            lambda: 1.0,
            band: 16,
            affine_within: None,
            clip_ends: true,
            end_bonus: 5,
            xdrop: 100,
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
struct Strand {
    hits: Vec<TxHit>,
    seq: Vec<u8>,
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
    fn tiebreak(&self, tb: &TieBreak, seq: &[u8], hits: &mut [TxHit]) -> u32 {
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
            self.segmented(tb, &codes, &hits[..n])
        } else {
            None
        };
        let (cost, unit) = seg.unwrap_or_else(|| {
            stats::add(&stats::TIEBREAK_WHOLE_READ, 1);
            let d = self.whole_read(tb, &codes, &hits[..n]);
            (d.into_iter().map(i64::from).collect(), tb.kappa)
        });
        let cmin = *cost.iter().min().unwrap();
        for (h, &c) in hits[..n].iter_mut().zip(&cost) {
            h.score = best - unit * (c - cmin) as f32;
        }
        hits.sort_by(|a, b| b.score.total_cmp(&a.score).then(a.tx.cmp(&b.tx)));
        n as u32
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

    /// Distances through anchors shared by every candidate: segments between consecutive
    /// shared anchors are aligned globally, but only where the candidates' transcript
    /// sequences differ (identical segments cost every candidate the same); the read ends are
    /// aligned against windows anchored at the outermost shared anchors. `None` if the
    /// candidates share no anchor.
    fn segmented(&self, tb: &TieBreak, codes: &[u8], hits: &[TxHit]) -> Option<(Vec<i64>, f32)> {
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
            let texts = hits
                .iter()
                .zip(&tpos)
                .map(|(h, t)| {
                    self.idx
                        .tx_codes(h.tx, t[0].saturating_sub(q0 + tb.pad + q0 / tb.pad_div), t[0])
                })
                .collect();
            add(texts, &codes[..q0 as usize], align::Mode::Prefix, true);
        }
        // Right end: read[q_last, plen) against windows starting at the last shared anchor.
        let ql = *shared.last().unwrap();
        let rest = plen - ql;
        let texts = hits
            .iter()
            .zip(&tpos)
            .map(|(h, t)| {
                let tlen = self.idx.tx_plen()[h.tx as usize];
                let s = *t.last().unwrap();
                self.idx
                    .tx_codes(h.tx, s, (s + rest + tb.pad + rest / tb.pad_div).min(tlen))
            })
            .collect();
        add(texts, &codes[ql as usize..], align::Mode::Prefix, false);
        let unit = if tb.affine.is_some() { tb.lambda } else { tb.kappa };
        Some((cost, unit))
    }

    fn map_strand(&self, raw: &[u8]) -> Strand {
        let seq = self.seeder.parsed(raw).map_or_else(|| raw.to_vec(), |h| h.seq);
        let seeds = self.seeder.seed_parsed(&seq);
        let anchors = anchor::anchors(&self.idx, &seeds, &self.p.anchor);
        let n_anchors = anchors.len() as u32;
        let cands = anchor::candidates(&self.idx, anchors, &self.p.anchor);
        let mut hits: Vec<TxHit> = cands
            .iter()
            .flat_map(|c| self.chainer.chain(&self.idx, c))
            .collect();
        sort_hits(&mut hits);
        Strand {
            hits,
            seq,
            seeds,
            anchors: n_anchors,
            candidates: cands.len() as u32,
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
            seq: s.seeder.parsed(r).map_or_else(|| r.to_vec(), |h| h.seq),
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
            tiebreak_size = self.tiebreak(tb, &st.seq, &mut st.hits);
            if both {
                // Near-tied orientations (antisense transcripts sharing sequence): re-rank each
                // orientation, then place the other's scores relative to this one by the
                // whole-read edit distances of the two best transcripts.
                stats::add(&stats::BOTH_ORIENTATIONS, 1);
                tiebreak_size += self.tiebreak(tb, &other.seq, &mut other.hits);
                let codes = |s: &[u8]| -> Vec<u8> {
                    s.iter()
                        .map(|&c| lcp_core::encode::CODE[c as usize] & 3)
                        .collect()
                };
                let da = self.whole_read(tb, &codes(&st.seq), &st.hits[..1])[0];
                let db = self.whole_read_k(tb, &codes(&other.seq), &other.hits[0], da.saturating_add(tb.cap));
                if db < da {
                    std::mem::swap(&mut st, &mut other);
                    reverse = !reverse;
                }
                let target = st.hits[0].score - tb.kappa * da.abs_diff(db) as f32;
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
