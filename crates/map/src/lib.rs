//! Read mapping: seeding both strands, anchor lookup, candidate-gene selection, chaining, and
//! result extraction (design doc, "Read processing" and "Chaining").

pub mod align;
pub mod anchor;
pub mod chain;
pub mod rescue;

use std::sync::Arc;

use anyhow::{Context, Result};
use index::Index;
use lcp_core::encode::revcomp;
use seeds::BaseSeeder;

pub use anchor::{Anchor, AnchorParams, Candidate};
pub use chain::{ChainParams, Chainer, OracleChainer, TxHit};
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
    /// Transcript bases added beyond the read's projected ends.
    pub pad: u32,
    /// Align at most this many transcripts per read.
    pub max_cands: usize,
    /// Band: distances beyond the best candidate's + cap are not computed exactly.
    pub cap: u32,
    /// Align only where tied transcripts differ (false: whole read per candidate).
    pub segmented: bool,
}

impl Default for TieBreak {
    fn default() -> Self {
        TieBreak {
            margin_abs: 30.0,
            margin_rel: 0.05,
            kappa: 6.0,
            pad: 50,
            max_cands: 16,
            cap: 16,
            segmented: true,
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
        let codes: Vec<u8> = seq
            .iter()
            .map(|&b| lcp_core::encode::CODE[b as usize] & 3)
            .collect();
        let seg = if tb.segmented {
            self.segmented(tb, &codes, &hits[..n])
        } else {
            None
        };
        let dist = seg.unwrap_or_else(|| self.whole_read(tb, &codes, &hits[..n]));
        let dmin = *dist.iter().min().unwrap();
        for (h, &d) in hits[..n].iter_mut().zip(&dist) {
            h.score = best - tb.kappa * (d - dmin) as f32;
        }
        hits.sort_by(|a, b| b.score.total_cmp(&a.score).then(a.tx.cmp(&b.tx)));
        n as u32
    }

    /// Distances of `read` against several texts with a shared band: the first is exact,
    /// later ones are capped at (best so far + `cap`) + 1.
    fn banded(read: &[u8], texts: &[Vec<u8>], mode: align::Mode, cap: u32) -> Vec<u32> {
        let mut best = u32::MAX;
        texts
            .iter()
            .map(|t| {
                let k = best.saturating_add(cap);
                let d = align::distance(read, t, mode, k).unwrap_or(k.saturating_add(1));
                best = best.min(d);
                d
            })
            .collect()
    }

    /// Whole read against each candidate's covered transcript region (projected to the read
    /// ends, padded), infix mode.
    fn whole_read(&self, tb: &TieBreak, codes: &[u8], hits: &[TxHit]) -> Vec<u32> {
        let plen = codes.len() as u32;
        let texts: Vec<Vec<u8>> = hits
            .iter()
            .map(|h| {
                let tlen = self.idx.tx_plen()[h.tx as usize];
                let (left, right) = (h.qstart, plen - h.qend);
                let slack = |ext: u32| tb.pad + ext / 6;
                let lo = h.tstart.saturating_sub(left + slack(left));
                let hi = (h.tend + right + slack(right)).min(tlen);
                self.idx.tx_codes(h.tx, lo, hi)
            })
            .collect();
        Self::banded(codes, &texts, align::Mode::Infix, tb.cap)
    }

    /// Distances through anchors shared by every candidate: segments between consecutive
    /// shared anchors are aligned globally, but only where the candidates' transcript
    /// sequences differ (identical segments cost every candidate the same); the read ends are
    /// aligned against windows anchored at the outermost shared anchors. `None` if the
    /// candidates share no anchor.
    fn segmented(&self, tb: &TieBreak, codes: &[u8], hits: &[TxHit]) -> Option<Vec<u32>> {
        // Query starts present in every chain, with each candidate's transcript start.
        let mut shared: Vec<u32> = hits[0].chain.iter().map(|a| a.0).collect();
        for h in &hits[1..] {
            let qs: std::collections::HashSet<u32> = h.chain.iter().map(|a| a.0).collect();
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
        let mut dist = vec![0u32; hits.len()];
        let mut add = |texts: Vec<Vec<u8>>, read: &[u8], mode: align::Mode, rev: bool| {
            if texts.windows(2).all(|w| w[0] == w[1]) {
                return; // same cost for every candidate
            }
            let (read, texts): (Vec<u8>, Vec<Vec<u8>>) = if rev {
                (
                    read.iter().rev().copied().collect(),
                    texts.into_iter().map(|t| t.into_iter().rev().collect()).collect(),
                )
            } else {
                (read.to_vec(), texts)
            };
            for (d, x) in dist.iter_mut().zip(Self::banded(&read, &texts, mode, tb.cap)) {
                *d += x;
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
                        .tx_codes(h.tx, t[0].saturating_sub(q0 + tb.pad + q0 / 6), t[0])
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
                    .tx_codes(h.tx, s, (s + rest + tb.pad + rest / 6).min(tlen))
            })
            .collect();
        add(texts, &codes[ql as usize..], align::Mode::Prefix, false);
        Some(dist)
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
        let fwd = self.map_strand(raw);
        let rc = self.map_strand(&revcomp(raw));
        let best = |h: &[TxHit]| h.first().map_or(f32::NEG_INFINITY, |x| x.score);
        let (reverse, mut st, other) = if best(&rc.hits) > best(&fwd.hits) {
            (true, rc, fwd)
        } else {
            (false, fwd, rc)
        };
        let (mut reverse, mut rescued) = (reverse, false);
        // Hierarchical rescue when the base-level chain covers too little of the read.
        if let Some(r) = &self.rescuer
            && Self::span_frac(st.hits.first(), st.seq.len()) < self.p.pi_cov
        {
            let mut a = r.rescue(&self.idx, &st.seq, &st.seeds, &self.p.anchor, &self.p.chain);
            let mut b = r.rescue(&self.idx, &other.seq, &other.seeds, &self.p.anchor, &self.p.chain);
            sort_hits(&mut a);
            sort_hits(&mut b);
            let (rev, hits, seq) = if best(&b) > best(&a) {
                (!reverse, b, &other.seq)
            } else {
                (reverse, a, &st.seq)
            };
            if Self::span_frac(hits.first(), seq.len()) > Self::span_frac(st.hits.first(), st.seq.len()) {
                if rev != reverse {
                    st = Strand { hits, ..other };
                } else {
                    st.hits = hits;
                }
                reverse = rev;
                rescued = true;
            }
        }
        let plen = st.seq.len() as u32;
        let mut hits = st.hits;
        let mut tiebreak_size = 0;
        if let Some(tb) = &self.p.tiebreak {
            tiebreak_size = self.tiebreak(tb, &st.seq, &mut hits);
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
