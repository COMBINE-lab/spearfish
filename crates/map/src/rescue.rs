//! Hierarchical rescue for reads whose base-level chain covers too little of the read
//! (typically low-accuracy reads, whose base cores rarely survive). The read and each
//! candidate transcript are parsed down to a lower LCP level, whose shorter cores survive
//! errors far better but are too unspecific for genome-wide lookup; restricted to candidate
//! transcripts, they are informative. Candidates are the strand's best genes by base-level
//! score (those scoring at least `alpha` times the best gene, at most
//! [`RescueParams::max_genes_for`] of them) and, in each, the `max_tx_per_gene` transcripts
//! with the most base anchors, skipping those longer than `max_tx_len`. With `window` set,
//! only the part of each transcript the read can occupy is re-seeded, and identical windows
//! are seeded once per read. Transcript sequences come from the index (parsed, 2-bit), so no
//! extra index is needed.

use index::Index;
use seeds::{BaseSeeder, Seed};

use crate::anchor::{self, AnchorParams};
use crate::chain::{ChainParams, Pt, TxHit, best_chains};

#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct RescueParams {
    /// Seed scheme for rescue seeds (its highest level), applied to parsed sequences.
    pub spec: String,
    /// Candidate genes per strand (best base-level score first).
    pub max_genes: usize,
    /// Reads shorter than this (parsed bases) get proportionally more candidate genes, up to
    /// `GENES_LEN_CAP` × `max_genes`: re-seeding a gene costs work in proportion to the read's
    /// length, and a short read's few base anchors rank its true gene among many genes that
    /// share one common label with it (0: `max_genes` for every read).
    pub genes_len: u32,
    /// Skip transcripts longer than this (parsed bases).
    pub max_tx_len: u32,
    /// Per gene, re-seed at most this many transcripts: those with the most base anchors.
    pub max_tx_per_gene: usize,
    /// Re-seed only the window of each transcript the read can occupy (projected from its
    /// base anchors, padded by `pad` + 1/`pad_div` of the read), and compute each distinct
    /// window sequence once per read (isoforms usually agree around the read). `None`:
    /// whole transcripts.
    pub window: Option<(u32, u32)>,
    /// Re-seed the other read orientation only if its base-level best score is within this
    /// margin of the chosen orientation's (None: always).
    pub other_margin: Option<f32>,
    /// Candidate genes must score at least this fraction of the best gene's base-level
    /// score (0: any gene with an anchor).
    pub alpha: f32,
}

impl Default for RescueParams {
    fn default() -> Self {
        RescueParams {
            spec: "lcp:K=2".into(),
            max_genes: 16,
            genes_len: 350,
            max_tx_len: 50_000,
            max_tx_per_gene: 128,
            window: Some((10, 50)),
            other_margin: Some(30.0),
            alpha: 0.5,
        }
    }
}

/// Most candidate genes a short read gets, as a multiple of `max_genes` (see `genes_len`).
pub const GENES_LEN_CAP: usize = 4;

impl RescueParams {
    /// Candidate genes for a read of `plen` parsed bases.
    pub fn max_genes_for(&self, plen: usize) -> usize {
        if self.genes_len == 0 {
            return self.max_genes;
        }
        (self.max_genes * self.genes_len as usize / plen.max(1))
            .clamp(self.max_genes, GENES_LEN_CAP * self.max_genes)
    }
}

pub struct Rescuer {
    pub p: RescueParams,
    seeder: BaseSeeder,
}

/// Per-thread buffers, reused across reads.
#[derive(Default)]
struct Scratch {
    /// The read's rescue seeds, in start order (starts are distinct).
    read: Vec<Seed>,
    /// Read seed label -> its first seed in `read`; `read_next` links each seed to the next one
    /// with its label (`u32::MAX`: none).
    read_at: rustc_hash::FxHashMap<u64, u32>,
    read_next: Vec<u32>,
    tx: Vec<Seed>,
    tmp: Vec<(u64, u32, u32)>,
    tseq: Vec<u8>,
    pts: Vec<Pt>,
    /// Transcript -> (base anchors, min and max anchor diagonal ts - qs).
    per_tx: rustc_hash::FxHashMap<u32, (u32, i64, i64)>,
    txs: Vec<(u32, u32, i64, i64)>,
    /// Window sequence (hash, length) -> distinct window index, per read orientation.
    seen: rustc_hash::FxHashMap<(u64, u32), u32>,
    /// Points of all distinct windows, and each window's range in it.
    pool: Vec<Pt>,
    wins: Vec<(u32, u32)>,
    /// Per (transcript, window start, distinct window), in visiting order.
    reqs: Vec<(u32, u32, u32)>,
}

/// A rescue chain in one transcript window: (score, max discrepancy, anchored bases, chain
/// as (qs, qe, ts relative to the window, len)).
type WindowChain = (f32, f32, u32, Vec<(u32, u32, u32, u32)>);

/// Rescue points chained per batch (320 kB of `Pt`s): a read whose windows hold more is
/// chained in several batches rather than leaving a large pool on its thread.
const POOL_CAP: usize = 1 << 14;

/// Best chain of each window `wins` (ranges of `pool`), as [`WindowChain`]s.
fn window_chains(cp: &ChainParams, pool: &[Pt], wins: &[(u32, u32)]) -> Vec<Option<WindowChain>> {
    let lists: Vec<&[Pt]> = wins
        .iter()
        .map(|&(a, n)| &pool[a as usize..(a + n) as usize])
        .collect();
    best_chains(cp, &lists)
        .into_iter()
        .zip(&lists)
        .map(|(res, pts)| {
            let (score, chain, maxd) = res?;
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
        })
        .collect()
}

thread_local! {
    static RESCUE: std::cell::RefCell<Scratch> = std::cell::RefCell::new(Scratch::default());
}

impl Rescuer {
    pub fn new(p: RescueParams) -> anyhow::Result<Self> {
        // Sequences are already parsed (HPC if the index used it): no HPC here.
        let seeder = BaseSeeder::new(&p.spec, false)?;
        Ok(Rescuer { p, seeder })
    }

    /// Rescue hits for one orientation of a read: `seq` is the parsed read, `base` its
    /// base-level seeds, `node_fam` the index's [`anchor::node_families`]. Hits are unsorted.
    pub fn rescue(
        &self,
        idx: &Index,
        node_fam: &[u32],
        seq: &[u8],
        base: &[Seed],
        ap: &AnchorParams,
        cp: &ChainParams,
    ) -> Vec<TxHit> {
        RESCUE.with(|sc| self.rescue_with(idx, node_fam, seq, base, ap, cp, &mut sc.borrow_mut()))
    }

    /// The read's rescue seeds into `sc.read`, with the label table `sc.read_at`/`read_next`
    /// that [`Self::points_in`] probes with each transcript window's seeds.
    fn seed_read(&self, seq: &[u8], sc: &mut Scratch) {
        self.seeder.seed_parsed_into(seq, &mut sc.read, &mut sc.tmp);
        debug_assert!(sc.read.windows(2).all(|w| w[0].start < w[1].start));
        sc.read_at.clear();
        sc.read_next.clear();
        sc.read_next.resize(sc.read.len(), u32::MAX);
        // Back to front, so that each label's chain runs in start order.
        for (i, s) in sc.read.iter().enumerate().rev() {
            if let Some(next) = sc.read_at.insert(s.label, i as u32) {
                sc.read_next[i] = next;
            }
        }
    }

    /// Points of the read's rescue seeds (`sc.read`, see [`Self::seed_read`]) against one
    /// transcript window sequence, sorted by (qs, ts), appended to `sc.pool`. Seeds have
    /// distinct starts on either side, so (qs, ts) is unique and fixes the order.
    fn points_in(&self, tseq: &[u8], sc: &mut Scratch) {
        self.seeder.seed_parsed_into(tseq, &mut sc.tx, &mut sc.tmp);
        // Hash join on labels: each window seed probes the read's label table.
        sc.pts.clear();
        for x in &sc.tx {
            let mut i = sc.read_at.get(&x.label).copied().unwrap_or(u32::MAX);
            while i != u32::MAX {
                let r = sc.read[i as usize];
                let len = r.end - r.start;
                sc.pts.push(Pt {
                    qs: r.start,
                    qe: r.end,
                    ts: x.start,
                    len,
                    w: len as f32,
                });
                i = sc.read_next[i as usize];
            }
        }
        sc.pts.sort_unstable_by_key(|p| (p.qs, p.ts));
        sc.pool.extend_from_slice(&sc.pts);
    }

    #[allow(clippy::too_many_arguments)]
    fn rescue_with(
        &self,
        idx: &Index,
        node_fam: &[u32],
        seq: &[u8],
        base: &[Seed],
        ap: &AnchorParams,
        cp: &ChainParams,
        sc: &mut Scratch,
    ) -> Vec<TxHit> {
        let genes = self.p.max_genes_for(seq.len());
        let all = AnchorParams {
            f_keep: genes,
            f_max: genes,
            alpha_fam: self.p.alpha,
            ..ap.clone()
        };
        let (cands, _) = anchor::select(idx, node_fam, base, &all);
        if cands.is_empty() {
            return Vec::new();
        }
        crate::stats::add(&crate::stats::RESCUE_STRANDS, 1);
        crate::stats::add(&crate::stats::RESCUE_GENES, cands.len() as u64);
        self.seed_read(seq, sc);
        let mut hits = Vec::new();
        let plen = seq.len() as i64;
        sc.seen.clear();
        sc.pool.clear();
        sc.wins.clear();
        sc.reqs.clear();
        // Window chains, computed in batches of at most POOL_CAP points so a read with many
        // or long windows doesn't leave a large pool behind on its thread.
        let mut results: Vec<Option<WindowChain>> = Vec::new();
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
                        let start = sc.pool.len() as u32;
                        self.points_in(&tseq, sc);
                        sc.tseq = tseq;
                        let r = (results.len() + sc.wins.len()) as u32;
                        sc.wins.push((start, sc.pool.len() as u32 - start));
                        sc.seen.insert(key, r);
                        if sc.pool.len() >= POOL_CAP {
                            results.extend(window_chains(cp, &sc.pool, &sc.wins));
                            sc.pool.clear();
                            sc.wins.clear();
                        }
                        r
                    }
                };
                sc.reqs.push((t, lo, r));
            }
        }
        // Best chain per distinct window, isoform windows sharing their DP prefixes.
        results.extend(window_chains(cp, &sc.pool, &sc.wins));
        if sc.pool.capacity() > POOL_CAP * 2 {
            sc.pool = Vec::with_capacity(POOL_CAP);
        }
        for &(t, lo, r) in &sc.reqs {
            let Some((score, maxd, anchored, chain)) = &results[r as usize] else {
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
                extended: false,
            });
        }
        hits
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn short_reads_get_more_candidate_genes_up_to_the_cap() {
        let p = RescueParams {
            genes_len: 500,
            ..RescueParams::default()
        };
        assert_eq!(p.max_genes_for(2000), p.max_genes);
        assert_eq!(p.max_genes_for(500), p.max_genes);
        assert_eq!(p.max_genes_for(250), 2 * p.max_genes);
        assert_eq!(p.max_genes_for(10), GENES_LEN_CAP * p.max_genes);
        let off = RescueParams { genes_len: 0, ..p };
        assert_eq!(off.max_genes_for(10), off.max_genes);
    }

    /// `points_in` equals a nested-loop join of the read's and the window's seeds on labels.
    #[test]
    fn points_in_matches_nested_loop_join() {
        let r = Rescuer::new(RescueParams::default()).unwrap();
        let mut x = 0x2545f4914f6cdd1du64;
        let mut rand = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    b"ACGT"[(x >> 62) as usize]
                })
                .collect()
        };
        let mut sc = Scratch::default();
        for _ in 0..50 {
            // A window repeating a segment, and a read sharing an excerpt of it and the repeat
            // amid unrelated sequence: labels match several times on both sides.
            let mut tseq = rand(2000);
            let rep = tseq[200..400].to_vec();
            tseq.extend_from_slice(&rep);
            let mut read = rand(300);
            read.extend_from_slice(&tseq[100..1300]);
            read.extend_from_slice(&rep);
            read.extend(rand(300));
            r.seed_read(&read, &mut sc);
            // Points of earlier windows stay in place.
            sc.pool.clear();
            sc.pool.push(Pt {
                qs: 7,
                qe: 9,
                ts: 0,
                len: 2,
                w: 2.0,
            });
            r.points_in(&tseq, &mut sc);
            let got: Vec<_> = sc
                .pool
                .iter()
                .map(|p| (p.qs, p.ts, p.qe, p.len, p.w.to_bits()))
                .collect();
            let (rs, ws) = (r.seeder.seed_parsed(&read), r.seeder.seed_parsed(&tseq));
            let mut want: Vec<_> = rs
                .iter()
                .flat_map(|a| {
                    let len = a.end - a.start;
                    ws.iter()
                        .filter(move |b| b.label == a.label)
                        .map(move |b| (a.start, b.start, a.end, len, (len as f32).to_bits()))
                })
                .collect();
            want.sort_unstable();
            // Points are unique in (qs, ts), and some seeds match several times on each side.
            assert!(want.windows(2).all(|w| (w[0].0, w[0].1) < (w[1].0, w[1].1)));
            let distinct = |f: fn(&(u32, u32, u32, u32, u32)) -> u32| {
                want.iter().map(f).collect::<std::collections::HashSet<_>>().len()
            };
            assert!(distinct(|p| p.0) < want.len() && distinct(|p| p.1) < want.len());
            assert_eq!(got[0], (7, 0, 9, 2, 2f32.to_bits()));
            assert_eq!(got[1..], want[..]);
        }
    }
}
