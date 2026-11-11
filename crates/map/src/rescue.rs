//! Hierarchical rescue for reads whose base-level chain covers too little of the read
//! (typically low-accuracy reads, whose base cores rarely survive). The read and each
//! candidate transcript are parsed down to a lower LCP level, whose shorter cores survive
//! errors far better but are too unspecific for genome-wide lookup; restricted to candidate
//! transcripts, they are informative. Candidates are the strand's best genes by base-level
//! score (those scoring at least `alpha` times the best gene, at most
//! [`RescueParams::max_genes_for`] of them) and, in each, the `max_tx_per_gene` transcripts
//! with the most base anchors. With `window` set, only the part of each transcript the read
//! can occupy is re-seeded (parts longer than `max_tx_len` are skipped); identical windows
//! (or whole transcripts) are seeded once per read orientation. Transcript sequences come from the
//! index (parsed, 2-bit), so no extra index is needed.
//!
//! The same seeds serve reads the base level did place, when the placement is in doubt
//! (every read but the [`Settled`] ones):
//! - the read's best hits are chained again at the rescue level
//!   ([`RescueParams::rechain`]), which ranks similar transcripts on many short seeds
//!   instead of the few long ones that survived the read's errors;
//! - the read is seeded against the whole transcriptome ([`RescueParams::global`]), by a
//!   table of runs of consecutive rescue seeds, for candidates no base-level seed named;
//! - the other read orientation is chained in transcripts that share sequence,
//!   reverse-complemented, with the best hit ([`RescueParams::antisense`]).
//!
//! For these every transcript's rescue seeds are made once and kept
//! ([`RescueParams::store`], [`Rescuer::build_tables`]).

use std::sync::Arc;

use index::Index;
use seeds::{BaseSeeder, Seed};

use crate::anchor::{self, AnchorParams};
use crate::chain::{ChainParams, Pt, TxHit, best_chains};
use crate::tables::{self, Global, LOOKUPS, Store, Window, mix, tuple};

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
    /// Skip transcript windows longer than this (parsed bases).
    pub max_tx_len: u32,
    /// `max_tx_len` bounds the window re-seeded (off: the whole transcript's length, whatever
    /// part of it the window is, as in 0.1.0: a read from a transcript longer than the cap
    /// then has no rescue there).
    #[serde(default)]
    pub cap_window: bool,
    /// Per gene, re-seed at most this many transcripts: those with the most base anchors.
    pub max_tx_per_gene: usize,
    /// Re-seed only the window of each transcript the read can occupy (projected from its
    /// base anchors, padded by `pad` + 1/`pad_div` of the read), and compute each distinct
    /// window sequence once per read orientation (isoforms usually agree around the read).
    /// `None`: whole transcripts.
    pub window: Option<(u32, u32)>,
    /// Re-seed the other read orientation only if its base-level best score is within this
    /// margin of the chosen orientation's or passes `MapParams::min_score` by itself
    /// (None: always).
    pub other_margin: Option<f32>,
    /// Candidate genes must score at least this fraction of the best gene's base-level
    /// score (0: any gene with an anchor).
    pub alpha: f32,
    /// Also chain the other read orientation in up to this many transcripts that share
    /// base-level seeds, reverse-complemented, with the best hit's own region (antisense
    /// transcripts sharing sequence; 0: off).
    #[serde(default)]
    pub antisense: usize,
    /// The other orientation's best hit, when it is one of those transcripts, is compared
    /// with this one's by alignment if it scores at least this fraction of it, whatever the
    /// tie-break's orientation margin (0: always).
    #[serde(default)]
    pub antisense_min: f32,
    /// Those transcripts share at least this fraction of the read's parsed bases with the
    /// best hit's region (0: any shared seed).
    #[serde(default)]
    pub antisense_share: f32,
    /// Take every transcript's rescue seeds from a table (read from the index, or made once
    /// when the mapper is built if the index has none) instead of parsing each transcript
    /// window again for every read that needs it. The table keeps stretches of seeds that
    /// transcripts share once; without it only the occasional low-coverage read pays for
    /// rescue seeds, with it every read may use them cheaply.
    #[serde(default)]
    pub store: bool,
    /// With `store`, a window's stored seeds become chaining points only within this many
    /// window slacks (`window`'s pad + read share) of the diagonals the request names (0:
    /// every seed of the window that matches a read seed). A request names the diagonals of
    /// a chain or of shared table keys, and matches far from them are repeats that cost
    /// chaining work without joining the chain.
    #[serde(default)]
    pub band: f32,
    /// Chain a read's best base-level hits again at the rescue level, in their own windows:
    /// those scoring within this fraction of the best (or the tie-break margin, if wider), at
    /// most [`RECHAIN_MAX`] of them (None: only reads the hierarchical rescue takes get
    /// rescue chains). Base seeds are long, so a read with errors keeps few of them, and
    /// which of several similar transcripts they favor is partly luck; the rescue level's
    /// many short seeds rank the same transcripts on far more evidence.
    #[serde(default)]
    pub rechain: Option<f32>,
    /// Reads whose base-level chains are taken as they stand: nothing below (`rechain`,
    /// `antisense`, `global`) is tried for them (None: for no read).
    #[serde(default)]
    pub settled: Option<Settled>,
    /// Seed reads against every transcript at the rescue level (None: only within candidates
    /// the base-level seeds name).
    #[serde(default)]
    pub global: Option<GlobalParams>,
}

/// A read is settled when its best base-level chain covers the read well and has no rival
/// on the other orientation. Most accurate reads are, and skip the work that only doubtful
/// reads repay.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Settled {
    /// The best chain's anchors cover at least this fraction of the parsed read,
    pub anchored: f32,
    /// it has at least this many anchors,
    pub anchors: u32,
    /// and the other orientation's best chain scores under this fraction of it.
    pub other: f32,
}

impl Default for Settled {
    fn default() -> Self {
        Settled {
            anchored: 0.3,
            anchors: 3,
            other: 0.65,
        }
    }
}

/// Most hits of one read orientation chained again ([`RescueParams::rechain`]).
pub const RECHAIN_MAX: usize = 256;

/// Rescue seeding against the whole transcriptome, by a table of every transcript's runs of
/// consecutive rescue seeds (made once when the mapper is built). One rescue seed is too
/// unspecific to look up; a short run of them in order is not. The table names transcripts
/// the read's base-level seeds never reached: paralogs and pseudogenes of the candidates,
/// and the true origin of a read whose base seeds all fell to errors.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct GlobalParams {
    /// Consecutive rescue seeds per table key.
    pub tuple: usize,
    /// Keys occurring more often than this in the transcriptome are not looked up.
    pub max_occ: u32,
    /// A transcript is a candidate if it shares at least this many keys with the read on
    /// nearby diagonals.
    pub min_shared: u32,
    /// A weak read is seeded against the whole transcriptome in both orientations: one
    /// whose best chain spans less than `MapParams::pi_cov` of it or, when it comes from the
    /// rescue, scores within this margin of the `MapParams::min_score` threshold.
    pub weak_margin: f32,
    /// Transcripts chained per orientation of a weak read: those sharing the most keys.
    pub weak_cands: usize,
    /// Of those, past the first `weak_keep` of an orientation, only the ones sharing at
    /// least `weak_rel` times as many keys as the best candidate of either orientation
    /// (0: all `weak_cands`). The true transcript is rarely far down the list, and chaining
    /// the tail of it is most of a weak read's cost.
    #[serde(default)]
    pub weak_rel: f32,
    #[serde(default)]
    pub weak_keep: usize,
    /// Transcripts chained for a read the rescue took that is not weak, next to the
    /// transcripts of the rescue's candidate genes (0: none). Few are needed: these reads
    /// have a good hit, and more candidates mostly add equally good alternates of large
    /// families.
    pub extra_cands: usize,
}

impl Default for GlobalParams {
    fn default() -> Self {
        GlobalParams {
            tuple: 3,
            max_occ: 5000,
            min_shared: 2,
            weak_margin: 10.0,
            weak_cands: 64,
            weak_rel: 0.25,
            weak_keep: 8,
            extra_cands: 8,
        }
    }
}

impl Default for RescueParams {
    fn default() -> Self {
        RescueParams {
            spec: "lcp:K=2".into(),
            max_genes: 16,
            genes_len: 350,
            max_tx_len: 50_000,
            cap_window: true,
            max_tx_per_gene: 128,
            window: Some((10, 50)),
            other_margin: Some(30.0),
            alpha: 0.5,
            antisense: 16,
            antisense_min: 0.5,
            antisense_share: 0.05,
            store: true,
            band: 4.0,
            rechain: Some(0.5),
            settled: Some(Settled::default()),
            global: Some(GlobalParams::default()),
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
    /// The table of [`RescueParams::global`].
    global: Option<Global>,
    /// The seeds of [`RescueParams::store`].
    store: Option<Store>,
    /// Whether the tables had to be made, the index keeping none of this rescue scheme
    /// ([`Self::build_tables`]).
    tables_made: bool,
    /// Distinct per rescuer: a thread's [`Session`]s are those of one rescuer.
    id: u64,
}

/// Read seed label -> its first seed in the read: an open-addressing table sized for the
/// read (a transcript window probes it once per seed).
#[derive(Default)]
struct LabelTable {
    /// (label, first seed), `u32::MAX` for an empty slot; a power of two of them.
    slots: Vec<(u64, u32)>,
    shift: u32,
}

impl LabelTable {
    /// Empty the table, for `n` labels.
    fn reset(&mut self, n: usize) {
        let cap = (2 * n).next_power_of_two().max(16);
        self.slots.clear();
        self.slots.resize(cap, (0, u32::MAX));
        self.shift = 64 - cap.trailing_zeros();
    }

    fn slot(&self, label: u64) -> usize {
        (label.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> self.shift) as usize
    }

    /// The first seed of `label` (`u32::MAX`: none).
    #[inline]
    fn get(&self, label: u64) -> u32 {
        let mask = self.slots.len() - 1;
        let mut i = self.slot(label);
        loop {
            let (l, first) = self.slots[i];
            if first == u32::MAX || l == label {
                return first;
            }
            i = (i + 1) & mask;
        }
    }

    /// Make `first` the first seed of `label`; returns the one it replaces (`u32::MAX`: none).
    fn insert(&mut self, label: u64, first: u32) -> u32 {
        let mask = self.slots.len() - 1;
        let mut i = self.slot(label);
        loop {
            let (l, was) = self.slots[i];
            if was == u32::MAX || l == label {
                self.slots[i] = (label, first);
                return was;
            }
            i = (i + 1) & mask;
        }
    }
}

/// A candidate of [`Rescuer::global_votes`]: (shared keys within one diagonal window,
/// transcript, lowest and highest diagonal).
type Vote = (u32, u32, i64, i64);

/// One read orientation's rescue state: its seeds and the windows chained so far, which
/// later calls for the same orientation reuse.
#[derive(Default)]
struct Session {
    /// The orientation's parsed sequence (hash, length) and the rescuer's id; `None`: unused.
    key: Option<(u64, usize, u64)>,
    /// The read's rescue seeds, in start order (starts are distinct).
    read: Vec<Seed>,
    /// Read seed label -> its first seed in `read`; `read_next` links each seed to the next one
    /// with its label (`u32::MAX`: none).
    read_at: LabelTable,
    read_next: Vec<u32>,
    /// Distinct window (hash of its points, their number) -> its index in `results`.
    seen: rustc_hash::FxHashMap<(u64, u32), u32>,
    /// With stored seeds: a window's seeds ([`Window::key`], their chunks) -> its index in
    /// `results`, and where its points' transcript positions count from, relative to its
    /// first seed. Isoforms mostly repeat a window seed for seed, and then need no join.
    seen_seeds: rustc_hash::FxHashMap<(u64, u32), (u32, i64)>,
    /// Transcript -> the first and last of its windows in `tx_wins`.
    by_tx: rustc_hash::FxHashMap<u32, (u32, u32)>,
    /// The windows requested (start, end, the transcript position its chain's target
    /// coordinates are relative to, index in `results`), each with the next one requested of
    /// its transcript (`u32::MAX`: none).
    tx_wins: Vec<((u32, u32, u32, u32), u32)>,
    /// The best chain of each distinct window (those of the pending batch follow).
    results: Vec<Option<WindowChain>>,
    /// The candidates of [`Rescuer::global_votes`] and how many were asked for: fewer than
    /// that are all there are, and the best of them are a prefix.
    voted: Option<(usize, Vec<Vote>)>,
}

/// Per-thread buffers, reused across reads.
#[derive(Default)]
struct Scratch {
    /// The current read orientation, and the one before it (usually the read's other one).
    s: Session,
    alt: Session,
    /// A window's rescue seeds as (label, start).
    tx: Vec<(u64, u32)>,
    tmp: Vec<(u64, u32, u32)>,
    /// A transcript window, as 2-bit codes.
    tseq: Vec<u8>,
    /// One window's points as (seed index in `read` << 32 | window start).
    keys: Vec<u64>,
    /// Transcript -> (base anchors, min and max anchor diagonal ts - qs).
    per_tx: rustc_hash::FxHashMap<u32, (u32, i64, i64)>,
    txs: Vec<(u32, u32, i64, i64)>,
    /// Points of the pending batch of distinct windows, and each window's range in it.
    pool: Vec<Pt>,
    wins: Vec<(u32, u32)>,
    /// Per (transcript, window start, distinct window), in visiting order.
    reqs: Vec<(u32, u32, u32)>,
    /// The postings of a read's table keys ([`Rescuer::votes`]), and how many fell to each
    /// value of a hash of their transcripts.
    votes: Vec<u64>,
    counts: Vec<u8>,
}

/// Most postings (eight bytes each) a thread keeps room for between reads.
const VOTES_KEEP: usize = 1 << 16;

/// Added to a posting's diagonal in [`Rescuer::votes`], which keeps it positive.
const OFF: i64 = 1 << 30;

/// [`Rescuer::votes`] takes first the hash values with the most postings: enough of them to
/// hold this many times the candidates asked for.
const VOTES_FIRST: usize = 2;

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
        static NEXT_ID: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        Ok(Rescuer {
            p,
            seeder,
            global: None,
            store: None,
            tables_made: false,
            id: NEXT_ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
        })
    }

    /// Whether the tables were made for this mapper rather than read from the index.
    pub fn tables_made(&self) -> bool {
        self.tables_made
    }

    /// Whether reads can be seeded against the whole transcriptome ([`RescueParams::global`],
    /// once [`Self::build_tables`] made its table).
    pub fn has_global(&self) -> bool {
        self.global.is_some()
    }

    /// The tables of [`RescueParams::global`] and [`RescueParams::store`] (nothing to do
    /// without either): those kept in the index, if they are of this rescue scheme, or else
    /// every transcript parsed at the rescue level, once.
    pub fn build_tables(&mut self, idx: &Arc<Index>) {
        let (tuple_k, store) = (self.p.global.map(|g| g.tuple.max(1)), self.p.store);
        if tuple_k.is_none() && !store {
            return;
        }
        let (s, g) = match tables::load(idx, &self.p.spec, tuple_k) {
            Some(kept) => kept,
            None => {
                self.tables_made = true;
                let nt = std::thread::available_parallelism()
                    .map_or(8, |x| x.get())
                    .min(64);
                let codes =
                    |t: u32, out: &mut Vec<u8>| idx.tx_codes_into(t, 0, idx.tx_plen()[t as usize], out);
                tables::build(&self.seeder, idx.n_transcripts(), &codes, tuple_k, nt)
            }
        };
        self.global = g;
        if store {
            self.store = Some(s);
        }
    }

    /// Rescue hits of one read orientation (`seq`, parsed) in the `n` transcripts sharing the
    /// most table keys with it on nearby diagonals ([`RescueParams::global`]; none without
    /// the table). Hits are unsorted.
    pub fn global_hits(&self, idx: &Index, seq: &[u8], cp: &ChainParams, n: usize) -> Vec<TxHit> {
        let cand = self.global_votes(seq, n);
        self.global_chain(idx, seq, cp, &cand)
    }

    /// The candidates of [`Self::global_hits`]: (shared keys within one diagonal window,
    /// transcript, lowest and highest diagonal), most keys first.
    pub fn global_votes(&self, seq: &[u8], n: usize) -> Vec<(u32, u32, i64, i64)> {
        let (Some(g), Some(table)) = (&self.p.global, &self.global) else {
            return Vec::new();
        };
        if n == 0 {
            return Vec::new();
        }
        RESCUE.with(|sc| {
            let sc = &mut *sc.borrow_mut();
            self.begin(seq, sc);
            if let Some((asked, cand)) = &sc.s.voted
                && (n <= *asked || cand.len() < *asked)
            {
                return cand[..n.min(cand.len())].to_vec();
            }
            let cand = Self::votes(g, table, seq.len(), n, sc);
            sc.s.voted = Some((n, cand.clone()));
            cand
        })
    }

    /// Rescue hits of one read orientation (`seq`, parsed) in the candidates `cand` of
    /// [`Self::global_votes`]. Hits are unsorted.
    pub fn global_chain(
        &self,
        idx: &Index,
        seq: &[u8],
        cp: &ChainParams,
        cand: &[(u32, u32, i64, i64)],
    ) -> Vec<TxHit> {
        if cand.is_empty() {
            return Vec::new();
        }
        RESCUE.with(|sc| {
            let sc = &mut *sc.borrow_mut();
            self.begin(seq, sc);
            crate::stats::add(&crate::stats::RESCUE_STRANDS, 1);
            for (k, &(_, t, lo, hi)) in cand.iter().enumerate() {
                let next = cand.get(k + 1).map(|&(_, t, lo, hi)| (t, lo, hi));
                self.ahead(idx, seq.len(), next.as_ref(), cand.get(k + 2).map(|x| x.1));
                self.request(idx, cp, seq.len(), (t, lo, hi), sc);
            }
            self.finish(cp, sc)
        })
    }

    /// The `n` transcripts sharing the most table keys with a read of `plen` parsed bases
    /// (`sc.s.read` holds its seeds): (keys within one diagonal window, transcript, lowest and
    /// highest diagonal), most keys first.
    fn votes(
        g: &GlobalParams,
        table: &Global,
        plen: usize,
        n: usize,
        sc: &mut Scratch,
    ) -> Vec<(u32, u32, i64, i64)> {
        let max_occ = g.max_occ as usize;
        let (read, tk) = (&sc.s.read, g.tuple.max(1));
        // Every posting of the read's keys, as transcript << 32 | diagonal.
        let keys = &mut sc.votes;
        keys.clear();
        let mut batch = [0u64; LOOKUPS];
        for from in (0..(read.len() + 1).saturating_sub(tk)).step_by(LOOKUPS) {
            let m = LOOKUPS.min(read.len() + 1 - tk - from);
            for (k, w) in batch.iter_mut().zip(read[from..].windows(tk)) {
                *k = tuple(w.iter().map(|x| mix(x.label)));
            }
            table.get_batch(&batch[..m], |i, posts| {
                if posts.len() <= max_occ {
                    let d = OFF - i64::from(read[from + i].start);
                    let at = keys.len();
                    keys.resize(at + posts.len(), 0);
                    let mut out = keys[at..].iter_mut();
                    posts.for_each(|t, p| {
                        if let Some(o) = out.next() {
                            *o = (u64::from(t) << 32) | (i64::from(p) + d) as u64;
                        }
                    });
                    let unread = out.len();
                    keys.truncate(keys.len() - unread);
                }
            });
        }
        crate::stats::add(&crate::stats::VOTE_CALLS, 1);
        crate::stats::add(
            &crate::stats::VOTE_KEYS,
            (read.len() + 1).saturating_sub(tk) as u64,
        );
        crate::stats::add(&crate::stats::VOTE_POSTS, keys.len() as u64);
        let w = 32 + plen as u64 / 8;
        let by = |&(c, t, _, _): &(u32, u32, i64, i64)| (std::cmp::Reverse(c), t);
        let mut cand: Vec<(u32, u32, i64, i64)> = Vec::new();
        // Most postings are the only one of their transcript, which then cannot be a
        // candidate: count postings by a hash of the transcript, and sort only those that
        // have, or may have, enough company.
        if g.min_shared > 1 && !keys.is_empty() {
            let bits = (4 * keys.len())
                .next_power_of_two()
                .trailing_zeros()
                .clamp(10, 20);
            let slot = |k: u64| ((k >> 32).wrapping_mul(0x9e37_79b9_7f4a_7c15) >> (64 - bits)) as usize;
            let count = &mut sc.counts;
            count.clear();
            count.resize(1 << bits, 0u8);
            // reached[c]: the hash values with at least `c` postings (255: or more).
            let mut reached = [0u32; 256];
            for &k in keys.iter() {
                let c = &mut count[slot(k)];
                let up = c.saturating_add(1);
                reached[usize::from(up)] += u32::from(up != *c);
                *c = up;
            }
            let need = g.min_shared.min(255) as u8;
            // A transcript has at most as many postings as its hash value, and a candidate no
            // more shared keys than postings. So the transcripts of the fullest hash values go
            // first: once they make `n` candidates, those of a hash value with fewer postings
            // than the least of them has shared keys cannot be among the best `n`, and their
            // postings are never sorted. (An equal count can be, by its transcript's rank.)
            let want = n.saturating_mul(VOTES_FIRST);
            let first = (need..=255)
                .rev()
                .find(|&c| reached[usize::from(c)] as usize >= want)
                .unwrap_or(need);
            let mut top = 0;
            for i in 0..keys.len() {
                let k = keys[i];
                keys.swap(i, top);
                top += usize::from(count[slot(k)] >= first);
            }
            Self::windows(&mut keys[..top], w, g.min_shared, &mut cand);
            let mut floor = need;
            if cand.len() >= n {
                cand.select_nth_unstable_by_key(n - 1, by);
                cand.truncate(n);
                floor = floor.max(cand[n - 1].0.min(255) as u8);
            }
            let mut kept = top;
            if floor < first {
                for i in top..keys.len() {
                    let k = keys[i];
                    keys[kept] = k;
                    kept += usize::from(count[slot(k)] >= floor);
                }
                Self::windows(&mut keys[top..kept], w, g.min_shared, &mut cand);
            }
            crate::stats::add(&crate::stats::VOTE_KEPT, kept as u64);
        } else {
            crate::stats::add(&crate::stats::VOTE_KEPT, keys.len() as u64);
            Self::windows(keys, w, g.min_shared, &mut cand);
        }
        crate::stats::add(&crate::stats::VOTE_CANDS, cand.len() as u64);
        if keys.capacity() > VOTES_KEEP {
            *keys = Vec::new();
            sc.counts = Vec::new();
        }
        // Transcripts are distinct, so the order is total: the best `n`, then those in order.
        if cand.len() > n {
            cand.select_nth_unstable_by_key(n - 1, by);
            cand.truncate(n);
        }
        cand.sort_unstable_by_key(by);
        cand
    }

    /// Adds to `cand` each transcript of `keys` (transcript << 32 | diagonal, which this
    /// sorts) with at least `min_shared` of them within `w` diagonals: (the most within `w`,
    /// transcript, lowest and highest diagonal of those).
    fn windows(keys: &mut [u64], w: u64, min_shared: u32, cand: &mut Vec<(u32, u32, i64, i64)>) {
        keys.sort_unstable();
        let mut i = 0;
        while i < keys.len() {
            let t = keys[i] >> 32;
            let mut j = i;
            while j < keys.len() && keys[j] >> 32 == t {
                j += 1;
            }
            let (mut bestc, mut lo, mut hi) = (0u32, 0, 0);
            let mut a = i;
            for b in i..j {
                while keys[b] - keys[a] > w {
                    a += 1;
                }
                if (b - a + 1) as u32 > bestc {
                    bestc = (b - a + 1) as u32;
                    lo = (keys[a] & 0xffff_ffff) as i64 - OFF;
                    hi = (keys[b] & 0xffff_ffff) as i64 - OFF;
                }
            }
            if bestc >= min_shared {
                cand.push((bestc, t as u32, lo, hi));
            }
            i = j;
        }
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

    /// A rescue seed label as the read's label table keys it: as the store keeps the
    /// transcripts' when there is one.
    fn label_key(&self, label: u64) -> u64 {
        if self.store.is_some() {
            u64::from(mix(label))
        } else {
            label
        }
    }

    /// The read's rescue seeds into `sc.s.read`, with the label table `read_at`/`read_next`
    /// that [`Self::points_in`] probes with each transcript window's seeds.
    fn seed_read(&self, seq: &[u8], sc: &mut Scratch) {
        let s = &mut sc.s;
        self.seeder.seed_parsed_into(seq, &mut s.read, &mut sc.tmp);
        debug_assert!(s.read.windows(2).all(|w| w[0].start < w[1].start));
        s.read_at.reset(s.read.len());
        s.read_next.clear();
        s.read_next.resize(s.read.len(), u32::MAX);
        // Back to front, so that each label's chain runs in start order.
        for (i, x) in s.read.iter().enumerate().rev() {
            s.read_next[i] = s.read_at.insert(self.label_key(x.label), i as u32);
        }
    }

    /// Keys of the points of the read's rescue seeds (`sc.s.read`, see [`Self::seed_read`])
    /// against one transcript window `tseq` (2-bit codes), into `sc.keys`: (read seed index,
    /// window start), sorted. Read seeds are in start order and seeds have distinct starts on
    /// either side, so this is the points' (qs, ts) order.
    fn point_keys(&self, tseq: &[u8], sc: &mut Scratch) {
        // The join needs only the window seeds' labels and starts.
        self.seeder.seed_starts_codes_into(tseq, &mut sc.tx);
        // Hash join on labels: each window seed probes the read's label table.
        sc.keys.clear();
        for &(label, start) in &sc.tx {
            let mut i = sc.s.read_at.get(label);
            while i != u32::MAX {
                sc.keys.push((u64::from(i) << 32) | u64::from(start));
                i = sc.s.read_next[i as usize];
            }
        }
        sc.keys.sort_unstable();
    }

    /// [`Self::point_keys`] of a transcript window starting at `lo`, by its stored seeds `win`.
    fn stored_point_keys(win: &Window, lo: u32, band: (i64, i64), sc: &mut Scratch) {
        sc.keys.clear();
        let (s, keys) = (&sc.s, &mut sc.keys);
        win.for_each(|label, start| {
            let mut i = s.read_at.get(u64::from(label));
            while i != u32::MAX {
                let d = i64::from(start) - i64::from(s.read[i as usize].start);
                if band.0 <= d && d <= band.1 {
                    keys.push((u64::from(i) << 32) | u64::from(start - lo));
                }
                i = s.read_next[i as usize];
            }
        });
        sc.keys.sort_unstable();
    }

    /// The points of `sc.keys` appended to `sc.pool`.
    fn push_points(sc: &mut Scratch) {
        sc.pool.extend(sc.keys.iter().map(|&k| {
            let r = sc.s.read[(k >> 32) as usize];
            let len = r.end - r.start;
            Pt {
                qs: r.start,
                qe: r.end,
                ts: k as u32,
                len,
                w: len as f32,
            }
        }));
    }

    /// Points of the read's rescue seeds against one transcript window `tseq` (2-bit codes),
    /// sorted by (qs, ts), appended to `sc.pool`.
    #[cfg(test)]
    fn points_in(&self, tseq: &[u8], sc: &mut Scratch) {
        self.point_keys(tseq, sc);
        Self::push_points(sc);
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
        self.begin(seq, sc);
        // Window chains are computed in batches of at most POOL_CAP points so a read with many
        // or long windows doesn't leave a large pool behind on its thread.
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
                let next = sc.txs.get(i + 1).map(|&(t, _, lo, hi)| (t, lo, hi));
                self.ahead(idx, seq.len(), next.as_ref(), sc.txs.get(i + 2).map(|x| x.0));
                self.request(idx, cp, seq.len(), (t, dlo, dhi), sc);
            }
        }
        self.finish(cp, sc)
    }

    /// Rescue hits of one orientation of a read (`seq`, parsed) in the transcripts `txs`, each
    /// with the range of diagonals (transcript position - read position) the read may lie on.
    /// Hits are unsorted.
    pub fn rescue_txs(
        &self,
        idx: &Index,
        seq: &[u8],
        txs: &[(u32, i64, i64)],
        cp: &ChainParams,
    ) -> Vec<TxHit> {
        if txs.is_empty() {
            return Vec::new();
        }
        RESCUE.with(|sc| {
            let sc = &mut *sc.borrow_mut();
            self.begin(seq, sc);
            for (k, &tx) in txs.iter().enumerate() {
                self.ahead(idx, seq.len(), txs.get(k + 1), txs.get(k + 2).map(|x| x.0));
                self.request(idx, cp, seq.len(), tx, sc);
            }
            self.finish(cp, sc)
        })
    }

    /// Start a call for a read orientation: its rescue seeds and the windows chained for it
    /// so far, if it is one of the last two orientations seen (a read's two, usually).
    fn begin(&self, seq: &[u8], sc: &mut Scratch) {
        let key = Some((xxhash_rust::xxh3::xxh3_64(seq), seq.len(), self.id));
        if sc.s.key != key {
            std::mem::swap(&mut sc.s, &mut sc.alt);
        }
        if sc.s.key != key {
            self.seed_read(seq, sc);
            let s = &mut sc.s;
            s.key = key;
            s.seen.clear();
            s.seen_seeds.clear();
            s.by_tx.clear();
            s.tx_wins.clear();
            s.results.clear();
            s.voted = None;
        }
        debug_assert!(sc.pool.is_empty() && sc.wins.is_empty());
        sc.reqs.clear();
    }

    /// The transcript window `(length of t, lo, hi)` that a read of `plen` parsed bases on
    /// the diagonals `dlo..=dhi` of transcript `t` can occupy.
    fn window_of(&self, idx: &Index, plen: usize, (t, dlo, dhi): (u32, i64, i64)) -> (u32, u32, u32) {
        let plen = plen as i64;
        let tlen = idx.tx_plen()[t as usize];
        match self.p.window {
            Some((pad, div)) => {
                let slack = pad as i64 + plen / div.max(1) as i64;
                (
                    tlen,
                    (dlo - slack).clamp(0, tlen as i64) as u32,
                    (dhi + plen + slack).clamp(0, tlen as i64) as u32,
                )
            }
            None => (tlen, 0, tlen),
        }
    }

    /// Ask the caches for the stored seeds of the windows requested next (`next`) and after
    /// it (`after`): a request's first reads are of memory no earlier one touched.
    fn ahead(&self, idx: &Index, plen: usize, next: Option<&(u32, i64, i64)>, after: Option<u32>) {
        let Some(store) = &self.store else {
            return;
        };
        if let Some(t) = after {
            store.prefetch_tx(t);
        }
        if let Some(&w) = next {
            let (tlen, lo, _) = self.window_of(idx, plen, w);
            store.prefetch_window(w.0, tlen, lo);
        }
    }

    /// Queue the window of transcript `t` that a read of `plen` parsed bases on the diagonals
    /// `dlo..=dhi` can occupy, chaining a full batch of windows.
    fn request(&self, idx: &Index, cp: &ChainParams, plen: usize, w: (u32, i64, i64), sc: &mut Scratch) {
        let t = w.0;
        let (tlen, lo, hi) = self.window_of(idx, plen, w);
        crate::stats::add(&crate::stats::REQUESTS, 1);
        // The cap bounds the bases re-seeded: a window of a longer transcript is within it.
        if (if self.p.cap_window { hi - lo } else { tlen }) > self.p.max_tx_len {
            return;
        }
        // A window within one already chained for this transcript has no better chain.
        let mut k = sc.s.by_tx.get(&t).map_or(u32::MAX, |x| x.0);
        while k != u32::MAX {
            let (w, next) = sc.s.tx_wins[k as usize];
            if w.0 <= lo && hi <= w.1 {
                crate::stats::add(&crate::stats::REQ_INSIDE, 1);
                sc.reqs.push((t, w.2, w.3));
                return;
            }
            k = next;
        }
        // The window's key, and the transcript position its points are relative to.
        let mut seeds = None;
        let (key, base) = match &self.store {
            Some(store) => {
                let win = store.window(t, tlen, lo, hi);
                // The same seeds as an earlier window's make the same points.
                let s0 = win.first_start().unwrap_or(lo);
                // Only the stored seeds a read on the requested diagonals can chain through
                // become points: `band` window slacks on either side of them.
                let bm = f64::from(self.p.band);
                let band = if bm > 0.0 {
                    let slack = self.p.window.map_or(0.0, |(pad, div)| {
                        f64::from(pad) + (plen / div.max(1) as usize) as f64
                    });
                    let b = (bm * slack) as i64;
                    (w.1 - b, w.2 + b)
                } else {
                    (i64::MIN, i64::MAX)
                };
                let bh = if bm > 0.0 {
                    ((band.0 - i64::from(s0)) as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15)
                        ^ ((band.1 - i64::from(s0)) as u64).wrapping_mul(0xc2b2_ae3d_27d4_eb4f)
                } else {
                    0
                };
                let skey = (win.key(s0) ^ bh, win.chunks() as u32);
                if let Some(&(r, d)) = sc.s.seen_seeds.get(&skey) {
                    crate::stats::add(&crate::stats::REQ_SAME_SEEDS, 1);
                    let base = (i64::from(s0) + d).max(0) as u32;
                    Self::note(t, (lo, hi, base, r), sc);
                    return;
                }
                seeds = Some((skey, s0));
                Self::stored_point_keys(&win, lo, band, sc);
                // Equal points up to a shift along the transcript have equal chains: the
                // points are kept relative to the first one's transcript position.
                let first = sc.keys.iter().map(|&k| k as u32).min().unwrap_or(0);
                let mut h = 0x9e37_79b9_7f4a_7c15u64;
                for k in &mut sc.keys {
                    *k -= u64::from(first);
                    h = (h.rotate_left(23) ^ *k).wrapping_mul(0xff51_afd7_ed55_8ccd);
                }
                ((h, sc.keys.len() as u32), lo + first)
            }
            None => {
                sc.tseq.clear();
                idx.tx_codes_into(t, lo, hi, &mut sc.tseq);
                ((xxhash_rust::xxh3::xxh3_64(&sc.tseq), sc.tseq.len() as u32), lo)
            }
        };
        let r = match sc.s.seen.get(&key) {
            Some(&r) => {
                crate::stats::add(&crate::stats::REQ_SAME_POINTS, 1);
                r
            }
            None => {
                crate::stats::add(&crate::stats::RESCUE_TX, 1);
                crate::stats::add(&crate::stats::RESCUE_TX_BASES, (hi - lo) as u64);
                if self.store.is_none() {
                    let tseq = std::mem::take(&mut sc.tseq);
                    self.point_keys(&tseq, sc);
                    sc.tseq = tseq;
                }
                let start = sc.pool.len() as u32;
                Self::push_points(sc);
                let r = (sc.s.results.len() + sc.wins.len()) as u32;
                sc.wins.push((start, sc.pool.len() as u32 - start));
                sc.s.seen.insert(key, r);
                if sc.pool.len() >= POOL_CAP {
                    Self::flush(cp, sc);
                }
                r
            }
        };
        if let Some((skey, s0)) = seeds {
            sc.s.seen_seeds.insert(skey, (r, i64::from(base) - i64::from(s0)));
        }
        Self::note(t, (lo, hi, base, r), sc);
    }

    /// Record the window `w` of transcript `t` (see [`Session::tx_wins`]) as requested.
    fn note(t: u32, w: (u32, u32, u32, u32), sc: &mut Scratch) {
        let s = &mut sc.s;
        let k = s.tx_wins.len() as u32;
        s.tx_wins.push((w, u32::MAX));
        match s.by_tx.entry(t) {
            std::collections::hash_map::Entry::Occupied(mut e) => {
                let last = &mut e.get_mut().1;
                s.tx_wins[*last as usize].1 = k;
                *last = k;
            }
            std::collections::hash_map::Entry::Vacant(e) => {
                e.insert((k, k));
            }
        }
        sc.reqs.push((t, w.2, w.3));
    }

    /// Chain the pending batch of windows.
    fn flush(cp: &ChainParams, sc: &mut Scratch) {
        sc.s.results.extend(window_chains(cp, &sc.pool, &sc.wins));
        sc.pool.clear();
        sc.wins.clear();
    }

    /// Hits of the windows queued since [`Self::begin`]: the best chain per distinct window,
    /// isoform windows sharing their DP prefixes.
    fn finish(&self, cp: &ChainParams, sc: &mut Scratch) -> Vec<TxHit> {
        Self::flush(cp, sc);
        if sc.pool.capacity() > POOL_CAP * 2 {
            sc.pool = Vec::with_capacity(POOL_CAP);
        }
        let mut hits = Vec::new();
        for &(t, lo, r) in &sc.reqs {
            let Some((score, maxd, anchored, chain)) = &sc.s.results[r as usize] else {
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
            let codes: Vec<u8> = tseq.iter().map(|&b| lcp_core::encode::CODE[b as usize]).collect();
            r.points_in(&codes, &mut sc);
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
