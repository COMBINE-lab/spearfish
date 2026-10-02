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

/// Per-thread buffers, reused across reads.
#[derive(Default)]
struct Scratch {
    read: Vec<Seed>,
    tx: Vec<Seed>,
    tmp: Vec<(u64, u32, u32)>,
    tseq: Vec<u8>,
    pts: Vec<Pt>,
    per_tx: rustc_hash::FxHashMap<u32, u32>,
    txs: Vec<(u32, u32)>,
}

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
            alpha_fam: 0.0,
            ..ap.clone()
        };
        let cands = anchor::candidates(idx, anchors, &all);
        if cands.is_empty() {
            return Vec::new();
        }
        // Read seeds sorted by label, for merge joins with each transcript's seeds.
        self.seeder.seed_parsed_into(seq, &mut sc.read, &mut sc.tmp);
        sc.read.sort_unstable_by_key(|s| (s.label, s.start));
        let mut hits = Vec::new();
        for c in &cands {
            // Transcripts of this gene that received base anchors, most anchors first.
            sc.per_tx.clear();
            for a in &c.anchors {
                for &t in idx.occurrences(idx.node_unitig()[a.node as usize]).0 {
                    *sc.per_tx.entry(t).or_default() += 1;
                }
            }
            sc.txs.clear();
            sc.txs.extend(sc.per_tx.iter().map(|(&t, &n)| (t, n)));
            sc.txs.sort_unstable_by_key(|&(t, n)| (std::cmp::Reverse(n), t));
            sc.txs.truncate(self.p.max_tx_per_gene);
            for i in 0..sc.txs.len() {
                let t = sc.txs[i].0;
                let tlen = idx.tx_plen()[t as usize];
                if tlen > self.p.max_tx_len {
                    continue;
                }
                sc.tseq.clear();
                idx.tx_codes_into(t, 0, tlen, &mut sc.tseq, true);
                self.seeder.seed_parsed_into(&sc.tseq, &mut sc.tx, &mut sc.tmp);
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
                let Some((score, chain, maxd)) = best_chain(cp, &sc.pts) else {
                    continue;
                };
                let pts = &sc.pts;
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
