//! Hierarchical rescue for reads whose base-level chain covers too little of the read
//! (typically low-accuracy reads, whose base cores rarely survive). The read and each
//! candidate transcript are parsed down to a lower LCP level, whose shorter cores survive
//! errors far better but are too unspecific for genome-wide lookup; restricted to candidate
//! transcripts, they are informative. Candidates are the strand's best genes by base-level
//! score (those scoring at least `alpha` times the best gene, at most
//! [`RescueParams::max_genes_for`] of them) and, in each, the `max_tx_per_gene` transcripts
//! with the most base anchors, skipping those longer than `max_tx_len` (windows, with `window`). With `window` set,
//! only the part of each transcript the read can occupy is re-seeded; identical windows (or
//! whole transcripts) are seeded once per read orientation. Transcript sequences come from the
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
    /// Skip transcripts (with `window`, transcript windows) longer than this (parsed bases).
    pub max_tx_len: u32,
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
    /// Keep every transcript's rescue seeds in memory, made once when the mapper is built,
    /// instead of parsing each transcript window again for every read that needs it. Costs
    /// eight bytes per rescue seed of the transcriptome; without it only the occasional
    /// low-coverage read pays for rescue seeds, with it every read may use them cheaply.
    #[serde(default)]
    pub store: bool,
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
            other: 0.4,
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
    /// Transcripts chained for any other read that is not [`Settled`], next to its
    /// base-level candidates (0: none). Few are needed: these reads have a good hit, and
    /// more candidates mostly add equally good alternates of large families.
    pub extra_cands: usize,
    /// The other read orientation gets the same when its best chain competes with this
    /// one's (within the tie-break's orientation margin).
    pub other: bool,
}

impl Default for GlobalParams {
    fn default() -> Self {
        GlobalParams {
            tuple: 3,
            max_occ: 5000,
            min_shared: 2,
            weak_margin: 10.0,
            weak_cands: 64,
            extra_cands: 8,
            other: true,
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
            max_tx_per_gene: 128,
            window: Some((10, 50)),
            other_margin: Some(30.0),
            alpha: 0.5,
            antisense: 16,
            store: true,
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
}

/// Every transcript's rescue seeds as (label hash, start), in start order per transcript.
/// Transcripts are kept [`PIECE_TX`] to a piece, as the builder's threads made them.
pub struct Store {
    pieces: Vec<StorePiece>,
}

/// Transcripts per piece of the tables' build (a power of two).
const PIECE_TX: usize = 256;

#[derive(Default)]
struct StorePiece {
    /// Where each transcript's seeds start in `seeds` (one more entry closes the last).
    off: Vec<u32>,
    seeds: Vec<(u32, u32)>,
}

impl Store {
    /// The seeds of transcript `t` that start in `[lo, hi)`.
    fn window(&self, t: u32, lo: u32, hi: u32) -> &[(u32, u32)] {
        let (p, i) = (&self.pieces[t as usize / PIECE_TX], t as usize % PIECE_TX);
        let s = &p.seeds[p.off[i] as usize..p.off[i + 1] as usize];
        &s[s.partition_point(|x| x.1 < lo)..s.partition_point(|x| x.1 < hi)]
    }
}

/// A rescue seed label as the store keeps it.
fn mix(label: u64) -> u32 {
    (label.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 32) as u32
}

/// The table key of a run of consecutive rescue seeds, from their label hashes (the store
/// has no more of a label).
fn tuple(it: impl Iterator<Item = u32>) -> u64 {
    it.fold(0x9e37_79b9_7f4a_7c15, |h, l| {
        (h ^ u64::from(l))
            .wrapping_mul(0xff51_afd7_ed55_8ccd)
            .rotate_left(29)
    })
}

/// Every transcript's runs of consecutive rescue seeds: tuple key -> (transcript, start of
/// the run's first seed). Keys are split by their top eight bits.
pub struct Global {
    buckets: Vec<Bucket>,
}

/// The keys of one value of the top eight bits.
#[derive(Default)]
struct Bucket {
    /// The distinct keys in order, and where each one's postings start (one more entry
    /// closes the last).
    keys: Vec<u64>,
    first: Vec<u32>,
    /// Where the keys of each value of the next `DIR_BITS` bits start in `keys`.
    dir: Vec<u32>,
    posts: Vec<(u32, u32)>,
}

const BUCKETS: usize = 256;
const DIR_BITS: u32 = 12;
/// Buckets sorted at a time while the table is built: the unsorted postings of that many
/// are all the build holds beyond its input and the table.
const GROUP: usize = 32;

/// A run of rescue seeds while the table is built: (key, transcript, start).
type Run = (u64, u32, u32);

impl Bucket {
    /// The bucket of `runs`, all of its keys.
    fn new(runs: &mut [Run]) -> Self {
        assert!(runs.len() < u32::MAX as usize, "too many rescue seed postings");
        runs.sort_unstable();
        let (mut keys, mut first) = (Vec::new(), Vec::new());
        let mut dir = vec![0u32; (1 << DIR_BITS) + 1];
        let mut at = 0u32;
        for run in runs.chunk_by(|a, b| a.0 == b.0) {
            keys.push(run[0].0);
            first.push(at);
            at += run.len() as u32;
            dir[Self::dir_of(run[0].0) + 1] += 1;
        }
        first.push(at);
        for d in 0..1 << DIR_BITS {
            dir[d + 1] += dir[d];
        }
        keys.shrink_to_fit();
        first.shrink_to_fit();
        Bucket {
            keys,
            first,
            dir,
            posts: runs.iter().map(|&(_, t, s)| (t, s)).collect(),
        }
    }

    fn dir_of(key: u64) -> usize {
        (key >> (56 - DIR_BITS)) as usize & ((1 << DIR_BITS) - 1)
    }
}

impl Global {
    /// The postings of `key`.
    fn get(&self, key: u64) -> &[(u32, u32)] {
        let b = &self.buckets[(key >> 56) as usize];
        let d = Bucket::dir_of(key);
        let (lo, hi) = (b.dir[d] as usize, b.dir[d + 1] as usize);
        match b.keys[lo..hi].binary_search(&key) {
            Ok(i) => &b.posts[b.first[lo + i] as usize..b.first[lo + i + 1] as usize],
            Err(_) => &[],
        }
    }

    /// The table of the runs of `tk` seeds in the store's `pieces` (`hist`: each piece's runs
    /// per bucket), on `nt` threads. Buckets are made [`GROUP`] at a time: the seeds are
    /// read once per group, and only one group's postings are ever held unsorted.
    fn build(pieces: &[StorePiece], hist: &[Vec<u32>], tk: usize, nt: usize) -> Self {
        let size = |b: usize| hist.iter().map(|h| h[b] as usize).sum::<usize>();
        let most = (0..BUCKETS)
            .step_by(GROUP)
            .map(|b0| (b0..b0 + GROUP).map(size).sum())
            .max()
            .unwrap_or(0);
        let mut runs: Vec<Run> = vec![(0, 0, 0); most];
        let mut buckets: Vec<(usize, Bucket)> = Vec::with_capacity(BUCKETS);
        for b0 in (0..BUCKETS).step_by(GROUP) {
            // Each bucket's runs, cut into the parts the pieces fill.
            let mut of_bucket: Vec<&mut [Run]> = Vec::with_capacity(GROUP);
            let mut rest = runs.as_mut_slice();
            for b in b0..b0 + GROUP {
                let (a, r) = rest.split_at_mut(size(b));
                of_bucket.push(a);
                rest = r;
            }
            let mut parts: Vec<Vec<&mut [Run]>> = hist.iter().map(|_| Vec::with_capacity(GROUP)).collect();
            for (b, v) in of_bucket.iter_mut().enumerate() {
                let mut rest = &mut **v;
                for (h, part) in hist.iter().zip(&mut parts) {
                    let (a, r) = rest.split_at_mut(h[b0 + b] as usize);
                    part.push(a);
                    rest = r;
                }
            }
            let work = std::sync::Mutex::new(parts.into_iter().enumerate());
            std::thread::scope(|sc| {
                for _ in 0..nt {
                    sc.spawn(|| {
                        loop {
                            let Some((p, mut part)) = work.lock().unwrap().next() else {
                                break;
                            };
                            let sp = &pieces[p];
                            let mut fill = [0usize; GROUP];
                            for (i, o) in sp.off.windows(2).enumerate() {
                                let t = (p * PIECE_TX + i) as u32;
                                for w in sp.seeds[o[0] as usize..o[1] as usize].windows(tk) {
                                    let key = tuple(w.iter().map(|x| x.0));
                                    let b = ((key >> 56) as usize).wrapping_sub(b0);
                                    if b < GROUP {
                                        part[b][fill[b]] = (key, t, w[0].1);
                                        fill[b] += 1;
                                    }
                                }
                            }
                        }
                    });
                }
            });
            let work = std::sync::Mutex::new(of_bucket.into_iter().enumerate());
            let made: Vec<(usize, Bucket)> = std::thread::scope(|sc| {
                let hs: Vec<_> = (0..nt.min(GROUP))
                    .map(|_| {
                        sc.spawn(|| {
                            let mut out = Vec::new();
                            loop {
                                let Some((b, runs)) = work.lock().unwrap().next() else {
                                    break;
                                };
                                out.push((b0 + b, Bucket::new(runs)));
                            }
                            out
                        })
                    })
                    .collect();
                hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
            });
            buckets.extend(made);
        }
        buckets.sort_unstable_by_key(|x| x.0);
        Global {
            buckets: buckets.into_iter().map(|x| x.1).collect(),
        }
    }
}

/// One read orientation's rescue state: its seeds and the windows chained so far, which
/// later calls for the same orientation reuse.
#[derive(Default)]
struct Session {
    /// The orientation's parsed sequence (hash, length); `None`: unused.
    key: Option<(u64, usize)>,
    /// The read's rescue seeds, in start order (starts are distinct).
    read: Vec<Seed>,
    /// Read seed label -> its first seed in `read`; `read_next` links each seed to the next one
    /// with its label (`u32::MAX`: none).
    read_at: rustc_hash::FxHashMap<u64, u32>,
    read_next: Vec<u32>,
    /// Distinct window (hash, length) -> its index in `results`.
    seen: rustc_hash::FxHashMap<(u64, u32), u32>,
    /// Transcript -> the windows requested (start, end, the transcript position its chain's
    /// target coordinates are relative to, index in `results`).
    by_tx: rustc_hash::FxHashMap<u32, Vec<(u32, u32, u32, u32)>>,
    /// The best chain of each distinct window (those of the pending batch follow).
    results: Vec<Option<WindowChain>>,
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
        Ok(Rescuer {
            p,
            seeder,
            global: None,
            store: None,
        })
    }

    /// Whether reads can be seeded against the whole transcriptome ([`RescueParams::global`],
    /// once [`Self::build_tables`] made its table).
    pub fn has_global(&self) -> bool {
        self.global.is_some()
    }

    /// Parse every transcript at the rescue level, once, for the table of
    /// [`RescueParams::global`] and the seeds of [`RescueParams::store`] (nothing to do
    /// without either).
    pub fn build_tables(&mut self, idx: &Index) {
        let (tuple_k, store) = (self.p.global.map(|g| g.tuple.max(1)), self.p.store);
        if tuple_k.is_none() && !store {
            return;
        }
        let n = idx.n_transcripts();
        let nt = std::thread::available_parallelism()
            .map_or(8, |x| x.get())
            .min(64);
        let seeder = &self.seeder;
        // Pieces of PIECE_TX transcripts, taken by whichever thread is free. A thread grows
        // its seeds in one buffer and keeps each piece as an exact copy, so that no piece
        // holds more memory than its seeds.
        let next = std::sync::atomic::AtomicUsize::new(0);
        let n_pieces = n.div_ceil(PIECE_TX);
        let mut made: Vec<(usize, StorePiece, Vec<u32>)> = std::thread::scope(|sc| {
            let hs: Vec<_> = (0..nt)
                .map(|_| {
                    sc.spawn(|| {
                        let mut out = Vec::new();
                        let (mut codes, mut sd) = (Vec::new(), Vec::new());
                        let (mut off, mut seeds) = (Vec::new(), Vec::new());
                        loop {
                            let p = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            if p >= n_pieces {
                                break;
                            }
                            off.clear();
                            seeds.clear();
                            off.push(0);
                            let mut hist = vec![0u32; BUCKETS];
                            for t in p * PIECE_TX..((p + 1) * PIECE_TX).min(n) {
                                codes.clear();
                                idx.tx_codes_into(t as u32, 0, idx.tx_plen()[t], &mut codes);
                                seeder.seed_starts_codes_into(&codes, &mut sd);
                                if !sd.is_sorted_by_key(|x| x.1) {
                                    sd.sort_unstable_by_key(|x| x.1);
                                }
                                let from = seeds.len();
                                seeds.extend(sd.iter().map(|&(l, s)| (mix(l), s)));
                                if let Some(tk) = tuple_k {
                                    for w in seeds[from..].windows(tk) {
                                        hist[(tuple(w.iter().map(|x| x.0)) >> 56) as usize] += 1;
                                    }
                                }
                                off.push(seeds.len() as u32);
                            }
                            let piece = StorePiece {
                                off: off.clone(),
                                seeds: seeds.clone(),
                            };
                            out.push((p, piece, hist));
                        }
                        out
                    })
                })
                .collect();
            hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
        });
        made.sort_unstable_by_key(|x| x.0);
        let (pieces, hist): (Vec<StorePiece>, Vec<Vec<u32>>) =
            made.into_iter().map(|(_, sp, h)| (sp, h)).unzip();
        if let Some(tk) = tuple_k {
            self.global = Some(Global::build(&pieces, &hist, tk, nt.min(32)));
        }
        if store {
            self.store = Some(Store { pieces });
        }
    }

    /// Rescue hits of one read orientation (`seq`, parsed) in the `n` transcripts sharing the
    /// most table keys with it on nearby diagonals ([`RescueParams::global`]; none without
    /// the table). Hits are unsorted.
    pub fn global_hits(&self, idx: &Index, seq: &[u8], cp: &ChainParams, n: usize) -> Vec<TxHit> {
        let (Some(g), Some(table)) = (&self.p.global, &self.global) else {
            return Vec::new();
        };
        if n == 0 {
            return Vec::new();
        }
        RESCUE.with(|sc| {
            let sc = &mut *sc.borrow_mut();
            self.begin(seq, sc);
            let mut cand = Self::votes(g, table, seq.len(), sc);
            cand.truncate(n);
            crate::stats::add(&crate::stats::RESCUE_STRANDS, 1);
            for &(_, t, lo, hi) in &cand {
                self.request(idx, cp, seq.len(), (t, lo, hi), sc);
            }
            self.finish(cp, sc)
        })
    }

    /// Transcripts sharing table keys with a read of `plen` parsed bases (`sc.s.read` holds
    /// its seeds): (keys within one diagonal window, transcript, lowest and highest diagonal),
    /// most keys first.
    fn votes(g: &GlobalParams, table: &Global, plen: usize, sc: &mut Scratch) -> Vec<(u32, u32, i64, i64)> {
        {
            const OFF: i64 = 1 << 30;
            let mut keys: Vec<u64> = Vec::new();
            for w in sc.s.read.windows(g.tuple.max(1)) {
                let posts = table.get(tuple(w.iter().map(|x| mix(x.label))));
                if posts.len() <= g.max_occ as usize {
                    let start = w[0].start;
                    keys.extend(
                        posts
                            .iter()
                            .map(|&(t, p)| (u64::from(t) << 32) | (p as i64 - start as i64 + OFF) as u64),
                    );
                }
            }
            keys.sort_unstable();
            let w = 32 + plen as u64 / 8;
            // Per transcript: the most seeds within `w` diagonals.
            let mut cand: Vec<(u32, u32, i64, i64)> = Vec::new();
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
                if bestc >= g.min_shared {
                    cand.push((bestc, t as u32, lo, hi));
                }
                i = j;
            }
            cand.sort_unstable_by_key(|&(c, t, _, _)| (std::cmp::Reverse(c), t));
            cand
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
        s.read_at.clear();
        s.read_next.clear();
        s.read_next.resize(s.read.len(), u32::MAX);
        // Back to front, so that each label's chain runs in start order.
        for (i, x) in s.read.iter().enumerate().rev() {
            if let Some(next) = s.read_at.insert(self.label_key(x.label), i as u32) {
                s.read_next[i] = next;
            }
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
            let mut i = sc.s.read_at.get(&label).copied().unwrap_or(u32::MAX);
            while i != u32::MAX {
                sc.keys.push((u64::from(i) << 32) | u64::from(start));
                i = sc.s.read_next[i as usize];
            }
        }
        sc.keys.sort_unstable();
    }

    /// [`Self::point_keys`] of the window `[lo, hi)` of transcript `t`, by its stored seeds.
    fn stored_point_keys(&self, store: &Store, t: u32, lo: u32, hi: u32, sc: &mut Scratch) {
        sc.keys.clear();
        for &(label, start) in store.window(t, lo, hi) {
            let mut i = sc.s.read_at.get(&u64::from(label)).copied().unwrap_or(u32::MAX);
            while i != u32::MAX {
                sc.keys.push((u64::from(i) << 32) | u64::from(start - lo));
                i = sc.s.read_next[i as usize];
            }
        }
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
            for &tx in txs {
                self.request(idx, cp, seq.len(), tx, sc);
            }
            self.finish(cp, sc)
        })
    }

    /// Start a call for a read orientation: its rescue seeds and the windows chained for it
    /// so far, if it is one of the last two orientations seen (a read's two, usually).
    fn begin(&self, seq: &[u8], sc: &mut Scratch) {
        let key = Some((xxhash_rust::xxh3::xxh3_64(seq), seq.len()));
        if sc.s.key != key {
            std::mem::swap(&mut sc.s, &mut sc.alt);
        }
        if sc.s.key != key {
            self.seed_read(seq, sc);
            let s = &mut sc.s;
            s.key = key;
            s.seen.clear();
            s.by_tx.clear();
            s.results.clear();
        }
        debug_assert!(sc.pool.is_empty() && sc.wins.is_empty());
        sc.reqs.clear();
    }

    /// Queue the window of transcript `t` that a read of `plen` parsed bases on the diagonals
    /// `dlo..=dhi` can occupy, chaining a full batch of windows.
    fn request(
        &self,
        idx: &Index,
        cp: &ChainParams,
        plen: usize,
        (t, dlo, dhi): (u32, i64, i64),
        sc: &mut Scratch,
    ) {
        let plen = plen as i64;
        let tlen = idx.tx_plen()[t as usize];
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
        // The cap bounds the bases re-seeded: a window of a longer transcript is within it.
        if hi - lo > self.p.max_tx_len {
            return;
        }
        // A window within one already chained for this transcript has no better chain.
        if let Some(&(_, _, base, r)) =
            sc.s.by_tx
                .get(&t)
                .and_then(|w| w.iter().find(|w| w.0 <= lo && hi <= w.1))
        {
            sc.reqs.push((t, base, r));
            return;
        }
        // The window's key, and the transcript position its points are relative to.
        let (key, base) = match &self.store {
            Some(store) => {
                self.stored_point_keys(store, t, lo, hi, sc);
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
            Some(&r) => r,
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
        sc.s.by_tx.entry(t).or_default().push((lo, hi, base, r));
        sc.reqs.push((t, base, r));
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
