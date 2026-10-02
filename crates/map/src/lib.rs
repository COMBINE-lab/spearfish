//! Read mapping: seeding both strands, anchor lookup, candidate-gene selection, chaining, and
//! result extraction (design doc, "Read processing" and "Chaining").

pub mod anchor;
pub mod chain;

use std::sync::Arc;

use anyhow::{Context, Result};
use index::Index;
use lcp_core::encode::revcomp;
use seeds::BaseSeeder;

pub use anchor::{Anchor, AnchorParams, Candidate};
pub use chain::{ChainParams, Chainer, OracleChainer, TxHit};

#[derive(Clone, Debug)]
pub struct MapParams {
    pub anchor: AnchorParams,
    pub chain: ChainParams,
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
        }
    }
}
