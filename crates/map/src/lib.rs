//! Read mapping: seeding both strands, anchor lookup, candidate-gene selection, chaining, and
//! result extraction (design doc, "Read processing" and "Chaining").

pub mod align;
pub mod anchor;
pub mod chain;

use std::sync::Arc;

use anyhow::{Context, Result};
use index::Index;
use lcp_core::encode::revcomp;
use seeds::BaseSeeder;

pub use anchor::{Anchor, AnchorParams, Candidate};
pub use chain::{ChainParams, Chainer, OracleChainer, TxHit};

/// Targeted base-level tie-break among near-best transcripts.
#[derive(Clone, Debug)]
pub struct TieBreak {
    /// Transcripts within max(margin_abs, margin_rel × best) of the best chain score form
    /// the tie set.
    pub margin_abs: f32,
    pub margin_rel: f32,
    /// Score charged per extra edit relative to the tie set's best alignment.
    pub kappa: f32,
    /// Transcript bases added beyond the read's projected ends.
    pub pad: u32,
    /// Align at most this many transcripts per read.
    pub max_cands: usize,
}

impl Default for TieBreak {
    fn default() -> Self {
        TieBreak {
            margin_abs: 30.0,
            margin_rel: 0.05,
            kappa: 6.0,
            pad: 50,
            max_cands: 16,
        }
    }
}

#[derive(Clone, Debug)]
pub struct MapParams {
    pub anchor: AnchorParams,
    pub chain: ChainParams,
    /// Re-rank near-ties by edit distance to the read (None: chain scores only).
    pub tiebreak: Option<TieBreak>,
    /// Report transcripts scoring within this fraction of the best.
    pub delta_rep: f32,
    /// Flag reads whose best chain anchors less than this fraction of the read.
    pub pi_cov: f32,
}

impl Default for MapParams {
    fn default() -> Self {
        MapParams {
            anchor: AnchorParams::default(),
            chain: ChainParams::default(),
            tiebreak: Some(TieBreak::default()),
            delta_rep: 0.3,
            pi_cov: 0.6,
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
    /// Reported transcripts, best first.
    pub hits: Vec<TxHit>,
    pub anchors: u32,
    pub candidates: u32,
    /// The best chain anchors less than `pi_cov` of the read.
    pub low_coverage: bool,
    /// Transcripts re-ranked by the tie-break (0 if it didn't run).
    pub tiebreak_size: u32,
}

pub struct Mapper<C: Chainer> {
    pub idx: Arc<Index>,
    pub seeder: BaseSeeder,
    pub chainer: C,
    pub p: MapParams,
}

impl<C: Chainer> Mapper<C> {
    /// A mapper seeding exactly as the index was built.
    pub fn new(idx: Arc<Index>, chainer: C, p: MapParams) -> Result<Self> {
        let spec = idx.meta("spec").context("index lacks its seed spec")?;
        let hpc = idx.meta("hpc").context("index lacks its HPC flag")? == "true";
        let seeder = BaseSeeder::new(spec, hpc)?;
        Ok(Mapper {
            idx,
            seeder,
            chainer,
            p,
        })
    }

    /// Re-rank the near-tie set of `hits` (sorted, best first) by edit distance between the
    /// parsed read `seq` and each transcript's covered region. Returns the tie-set size.
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
        let codes: Vec<u8> = seq
            .iter()
            .map(|&b| lcp_core::encode::CODE[b as usize] & 3)
            .collect();
        let plen = seq.len() as u32;
        let dist: Vec<u32> = hits[..n]
            .iter()
            .map(|h| {
                let tlen = self.idx.tx_plen()[h.tx as usize];
                // Project the read's ends onto the transcript, with slack for indels.
                let (left, right) = (h.qstart, plen - h.qend);
                let slack = |ext: u32| tb.pad + ext / 6;
                let lo = h.tstart.saturating_sub(left + slack(left));
                let hi = (h.tend + right + slack(right)).min(tlen);
                align::infix_distance(&codes, &self.idx.tx_codes(h.tx, lo, hi))
            })
            .collect();
        let dmin = *dist.iter().min().unwrap();
        for (h, &d) in hits[..n].iter_mut().zip(&dist) {
            h.score = best - tb.kappa * (d - dmin) as f32;
        }
        hits.sort_by(|a, b| b.score.total_cmp(&a.score).then(a.tx.cmp(&b.tx)));
        n as u32
    }

    fn map_strand(&self, raw: &[u8]) -> (Vec<TxHit>, u32, u32, u32) {
        let h = self.seeder.parsed(raw);
        let seq: &[u8] = h.as_ref().map_or(raw, |h| &h.seq);
        let seeds = self.seeder.seed_parsed(seq);
        let anchors = anchor::anchors(&self.idx, &seeds, &self.p.anchor);
        let n_anchors = anchors.len() as u32;
        let cands = anchor::candidates(&self.idx, anchors, &self.p.anchor);
        let mut hits: Vec<TxHit> = cands
            .iter()
            .flat_map(|c| self.chainer.chain(&self.idx, c))
            .collect();
        hits.sort_unstable_by(|a, b| b.score.total_cmp(&a.score).then(a.tx.cmp(&b.tx)));
        (hits, n_anchors, cands.len() as u32, seq.len() as u32)
    }

    pub fn map(&self, raw: &[u8]) -> ReadResult {
        let fwd = self.map_strand(raw);
        let rc = self.map_strand(&revcomp(raw));
        let best = |h: &[TxHit]| h.first().map_or(f32::NEG_INFINITY, |x| x.score);
        let (reverse, (mut hits, anchors, candidates, plen)) = if best(&rc.0) > best(&fwd.0) {
            (true, rc)
        } else {
            (false, fwd)
        };
        let mut tiebreak_size = 0;
        if let Some(tb) = &self.p.tiebreak {
            let strand = if reverse { revcomp(raw) } else { raw.to_vec() };
            let seq = self.seeder.parsed(&strand).map_or(strand, |h| h.seq);
            tiebreak_size = self.tiebreak(tb, &seq, &mut hits);
        }
        let mut low_coverage = true;
        if let Some(top) = hits.first() {
            let floor = top.score - self.p.delta_rep * top.score.abs();
            low_coverage = (top.anchored as f32) < self.p.pi_cov * plen as f32;
            hits.retain(|h| h.score >= floor);
        }
        ReadResult {
            len: raw.len() as u32,
            plen,
            reverse,
            hits,
            anchors,
            candidates,
            low_coverage,
            tiebreak_size,
        }
    }
}
