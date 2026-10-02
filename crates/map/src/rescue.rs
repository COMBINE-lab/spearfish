//! Hierarchical rescue for reads whose base-level chain covers too little of the read
//! (typically low-accuracy reads, whose base cores rarely survive). The read and each
//! candidate transcript are parsed down to a lower LCP level, whose shorter cores survive
//! errors far better but are too unspecific for genome-wide lookup; restricted to candidate
//! transcripts, they are informative. Candidates are the transcripts with the most base-level
//! anchors in each gene with any anchor on the strand; transcript sequences come from the
//! index (parsed, 2-bit), so no extra index is needed.

use index::Index;
use seeds::{BaseSeeder, Seed};

use crate::anchor::{self, AnchorParams};
use crate::chain::{ChainParams, Pt, TxHit, best_chain};

#[derive(Clone, Debug)]
pub struct RescueParams {
    /// Seed scheme for rescue seeds (its highest level), applied to parsed sequences.
    pub spec: String,
    /// Candidate genes per strand (best base-level score first).
    pub max_genes: usize,
    /// Skip transcripts longer than this (parsed bases).
    pub max_tx_len: u32,
    /// Per gene, re-seed at most this many transcripts: those with the most base anchors.
    pub max_tx_per_gene: usize,
    /// Re-seed only the window of each transcript the read can occupy (projected from its
    /// base anchors, padded by `pad` + 1/`pad_div` of the read), and compute each distinct
    /// window sequence once per read (isoforms usually agree around the read). `None`:
    /// whole transcripts.
    pub window: Option<(u32, u32)>,
}

impl Default for RescueParams {
    fn default() -> Self {
        RescueParams {
            spec: "lcp:K=2".into(),
            max_genes: 16,
            max_tx_len: 50_000,
            max_tx_per_gene: 128,
            window: Some((100, 5)),
        }
    }
}

pub struct Rescuer {
    pub p: RescueParams,
    seeder: BaseSeeder,
}

/// Per-thread buffers, reused across reads.
#[derive(Default)]
struct Scratch {
    read: Vec<Seed>,
    tx: Vec<Seed>,
    tmp: Vec<(u64, u32, u32)>,
    tseq: Vec<u8>,
    pts: Vec<Pt>,
    /// Transcript -> (base anchors, min and max anchor diagonal ts - qs).
    per_tx: rustc_hash::FxHashMap<u32, (u32, i64, i64)>,
    txs: Vec<(u32, u32, i64, i64)>,
    /// Window sequence (hash, length) -> index into `results`, per read orientation.
    seen: rustc_hash::FxHashMap<(u64, u32), u32>,
    /// Per distinct window: its best chain, or None.
    results: Vec<Option<WindowChain>>,
}

/// A rescue chain in one transcript window: (score, max discrepancy, anchored bases, chain
/// as (qs, qe, ts relative to the window, len)).
type WindowChain = (f32, f32, u32, Vec<(u32, u32, u32, u32)>);

thread_local! {
    static RESCUE: std::cell::RefCell<Scratch> = std::cell::RefCell::new(Scratch::default());
}

impl Rescuer {
    pub fn new(idx: &Index, p: RescueParams) -> anyhow::Result<Self> {
        let _ = idx;
        // Sequences are already parsed (HPC if the index used it): no HPC here.
        let seeder = BaseSeeder::new(&p.spec, false)?;
        Ok(Rescuer { p, seeder })
    }

    /// Rescue hits for one orientation of a read: `seq` is the parsed read, `base` its
    /// base-level seeds. Hits are unsorted.
    pub fn rescue(
        &self,
        idx: &Index,
        seq: &[u8],
        base: &[Seed],
        ap: &AnchorParams,
        cp: &ChainParams,
    ) -> Vec<TxHit> {
        RESCUE.with(|sc| self.rescue_with(idx, seq, base, ap, cp, &mut sc.borrow_mut()))
    }

    /// Best chain of the read's rescue seeds (`sc.read`, sorted by label) against one
    /// transcript sequence: (score, max discrepancy, anchored bases, chain as (qs, qe, ts,
    /// len)), or None.
    fn best_in(&self, tseq: &[u8], cp: &ChainParams, sc: &mut Scratch) -> Option<WindowChain> {
        self.seeder.seed_parsed_into(tseq, &mut sc.tx, &mut sc.tmp);
        sc.tx.sort_unstable_by_key(|s| (s.label, s.start));
        // Merge join on labels.
        sc.pts.clear();
        let (r, x) = (&sc.read, &sc.tx);
        let (mut i, mut j) = (0, 0);
        while i < r.len() && j < x.len() {
            match r[i].label.cmp(&x[j].label) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    let l = r[i].label;
                    let j0 = j;
                    while j < x.len() && x[j].label == l {
                        j += 1;
                    }
                    while i < r.len() && r[i].label == l {
                        let len = r[i].end - r[i].start;
                        for xs in &x[j0..j] {
                            sc.pts.push(Pt {
                                qs: r[i].start,
                                qe: r[i].end,
                                ts: xs.start,
                                len,
                                w: len as f32,
                            });
                        }
                        i += 1;
                    }
                }
            }
        }
        sc.pts.sort_unstable_by_key(|p| (p.qs, p.ts));
        let (score, chain, maxd) = best_chain(cp, &sc.pts)?;
        let pts = &sc.pts;
        let mut anchored = 0;
        let mut to = 0;
        for &k in &chain {
            let from = pts[k].qs.max(to);
            if pts[k].qe > from {
                anchored += pts[k].qe - from;
            }
            to = to.max(pts[k].qe);
        }
        Some((
            score,
            maxd,
            anchored,
            chain
                .iter()
                .map(|&k| (pts[k].qs, pts[k].qe, pts[k].ts, pts[k].len))
                .collect(),
        ))
    }

    #[allow(clippy::too_many_arguments)]
    fn rescue_with(
        &self,
        idx: &Index,
        seq: &[u8],
        base: &[Seed],
        ap: &AnchorParams,
        cp: &ChainParams,
        sc: &mut Scratch,
    ) -> Vec<TxHit> {
        let anchors = anchor::anchors(idx, base, ap);
        let all = AnchorParams {
            f_keep: self.p.max_genes,
            f_max: self.p.max_genes,
            alpha_fam: 0.0,
            ..ap.clone()
        };
        let cands = anchor::candidates(idx, anchors, &all);
        if cands.is_empty() {
            return Vec::new();
        }
        crate::stats::add(&crate::stats::RESCUE_STRANDS, 1);
        crate::stats::add(&crate::stats::RESCUE_GENES, cands.len() as u64);
        // Read seeds sorted by label, for merge joins with each transcript's seeds.
        self.seeder.seed_parsed_into(seq, &mut sc.read, &mut sc.tmp);
        sc.read.sort_unstable_by_key(|s| (s.label, s.start));
        let mut hits = Vec::new();
        let plen = seq.len() as i64;
        sc.seen.clear();
        sc.results.clear();
        for c in &cands {
            // Transcripts of this gene that received base anchors, most anchors first, with
            // the range of their anchors' diagonals.
            sc.per_tx.clear();
            for a in &c.anchors {
                let off = idx.node_offset()[a.node as usize] as i64;
                let (ts, ps) = idx.occurrences(idx.node_unitig()[a.node as usize]);
                for (&t, &p) in ts.iter().zip(ps) {
                    let d = p as i64 + off - a.qs as i64;
                    let e = sc.per_tx.entry(t).or_insert((0, d, d));
                    e.0 += 1;
                    e.1 = e.1.min(d);
                    e.2 = e.2.max(d);
                }
            }
            sc.txs.clear();
            sc.txs
                .extend(sc.per_tx.iter().map(|(&t, &(n, lo, hi))| (t, n, lo, hi)));
            sc.txs
                .sort_unstable_by_key(|&(t, n, _, _)| (std::cmp::Reverse(n), t));
            sc.txs.truncate(self.p.max_tx_per_gene);
            for i in 0..sc.txs.len() {
                let (t, _, dlo, dhi) = sc.txs[i];
                let tlen = idx.tx_plen()[t as usize];
                if tlen > self.p.max_tx_len {
                    continue;
                }
                let (lo, hi) = match self.p.window {
                    Some((pad, div)) => {
                        let slack = pad as i64 + plen / div.max(1) as i64;
                        (
                            (dlo - slack).clamp(0, tlen as i64) as u32,
                            (dhi + plen + slack).clamp(0, tlen as i64) as u32,
                        )
                    }
                    None => (0, tlen),
                };
                sc.tseq.clear();
                idx.tx_codes_into(t, lo, hi, &mut sc.tseq, true);
                let key = (xxhash_rust::xxh3::xxh3_64(&sc.tseq), sc.tseq.len() as u32);
                let r = match sc.seen.get(&key) {
                    Some(&r) => r,
                    None => {
                        crate::stats::add(&crate::stats::RESCUE_TX, 1);
                        crate::stats::add(&crate::stats::RESCUE_TX_BASES, (hi - lo) as u64);
                        let tseq = std::mem::take(&mut sc.tseq);
                        let res = self.best_in(&tseq, cp, sc);
                        sc.tseq = tseq;
                        let r = sc.results.len() as u32;
                        sc.results.push(res);
                        sc.seen.insert(key, r);
                        r
                    }
                };
                let Some((score, maxd, anchored, chain)) = &sc.results[r as usize] else {
                    continue;
                };
                let (first, last) = (chain[0], chain[chain.len() - 1]);
                hits.push(TxHit {
                    tx: t,
                    score: *score,
                    qstart: first.0,
                    qend: last.1,
                    tstart: first.2 + lo,
                    tend: last.2 + lo + last.3,
                    anchors: chain.len() as u32,
                    anchored: *anchored,
                    max_discrepancy: *maxd,
                    chain: chain.iter().map(|&(qs, _, ts, len)| (qs, ts + lo, len)).collect(),
                    reverse: false,
                });
            }
        }
        hits
    }
}
