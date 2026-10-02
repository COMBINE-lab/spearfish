//! Hierarchical rescue for reads whose base-level chain covers too little of the read
//! (typically low-accuracy reads, whose base cores rarely survive). The read and each
//! candidate transcript are parsed down to a lower LCP level, whose shorter cores survive
//! errors far better but are too unspecific for genome-wide lookup; restricted to candidate
//! transcripts, they are informative. Candidates are the transcripts with the most base-level
//! anchors in each gene with any anchor on the strand; transcript sequences come from the
//! index (parsed, 2-bit), so no extra index is needed.

use std::collections::HashMap;

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
}

impl Default for RescueParams {
    fn default() -> Self {
        RescueParams {
            spec: "lcp:K=2".into(),
            max_genes: 16,
            max_tx_len: 50_000,
            max_tx_per_gene: 128,
        }
    }
}

pub struct Rescuer {
    pub p: RescueParams,
    seeder: BaseSeeder,
}

const ACGT: [u8; 4] = *b"ACGT";

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
        let anchors = anchor::anchors(idx, base, ap);
        let all = AnchorParams {
            f_keep: self.p.max_genes,
            alpha_fam: 0.0,
            ..ap.clone()
        };
        let cands = anchor::candidates(idx, anchors, &all);
        if cands.is_empty() {
            return Vec::new();
        }
        let read_seeds = self.seeder.seed_parsed(seq);
        let mut hits = Vec::new();
        let mut tx_index: HashMap<u64, Vec<u32>> = HashMap::new();
        let mut pts: Vec<Pt> = Vec::new();
        let mut per_tx: HashMap<u32, u32> = HashMap::new();
        for c in &cands {
            // Transcripts of this gene that received base anchors, most anchors first.
            per_tx.clear();
            for a in &c.anchors {
                for &t in idx.occurrences(idx.node_unitig()[a.node as usize]).0 {
                    *per_tx.entry(t).or_default() += 1;
                }
            }
            let mut txs: Vec<(u32, u32)> = per_tx.iter().map(|(&t, &n)| (t, n)).collect();
            txs.sort_unstable_by_key(|&(t, n)| (std::cmp::Reverse(n), t));
            txs.truncate(self.p.max_tx_per_gene);
            for &(t, _) in &txs {
                let tlen = idx.tx_plen()[t as usize];
                if tlen > self.p.max_tx_len {
                    continue;
                }
                let tseq: Vec<u8> = idx
                    .tx_codes(t, 0, tlen)
                    .iter()
                    .map(|&c| ACGT[c as usize])
                    .collect();
                tx_index.clear();
                for s in self.seeder.seed_parsed(&tseq) {
                    tx_index.entry(s.label).or_default().push(s.start);
                }
                pts.clear();
                for s in &read_seeds {
                    if let Some(ps) = tx_index.get(&s.label) {
                        let len = s.end - s.start;
                        for &ts in ps {
                            pts.push(Pt {
                                qs: s.start,
                                qe: s.end,
                                ts,
                                len,
                                w: len as f32,
                            });
                        }
                    }
                }
                pts.sort_unstable_by_key(|p| (p.qs, p.ts));
                let Some((score, chain, maxd)) = best_chain(cp, &pts) else {
                    continue;
                };
                let (first, last) = (pts[chain[0]], pts[*chain.last().unwrap()]);
                let mut anchored = 0;
                let mut to = 0;
                for &k in &chain {
                    let from = pts[k].qs.max(to);
                    if pts[k].qe > from {
                        anchored += pts[k].qe - from;
                    }
                    to = to.max(pts[k].qe);
                }
                hits.push(TxHit {
                    tx: t,
                    score,
                    qstart: first.qs,
                    qend: last.qe,
                    tstart: first.ts,
                    tend: last.ts + last.len,
                    anchors: chain.len() as u32,
                    anchored,
                    max_discrepancy: maxd,
                    chain: chain
                        .iter()
                        .map(|&k| (pts[k].qs, pts[k].ts, pts[k].len))
                        .collect(),
                });
            }
        }
        hits
    }
}
