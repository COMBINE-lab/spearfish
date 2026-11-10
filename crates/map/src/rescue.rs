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
    /// The other orientation's best hit, when it is one of those transcripts, is compared
    /// with this one's by alignment if it scores at least this fraction of it, whatever the
    /// tie-break's orientation margin (0: always).
    #[serde(default)]
    pub antisense_min: f32,
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
            other: 0.5,
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
            antisense_min: 0.5,
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
    /// The seeds of transcript `t` (of `tlen` parsed bases) that start in `[lo, hi)`.
    fn window(&self, t: u32, tlen: u32, lo: u32, hi: u32) -> &[(u32, u32)] {
        let (p, i) = (&self.pieces[t as usize / PIECE_TX], t as usize % PIECE_TX);
        let s = &p.seeds[p.off[i] as usize..p.off[i + 1] as usize];
        let a = first_from(s, lo, Self::guess(s.len(), tlen, lo));
        let b = first_from(s, hi, Self::guess(s.len(), tlen, hi).max(a));
        &s[a..b]
    }

    /// About where, among a transcript's `n` seeds, those starting at `pos` are: seeds are
    /// spread nearly evenly along a transcript.
    fn guess(n: usize, tlen: u32, pos: u32) -> usize {
        (u64::from(pos) * n as u64 / u64::from(tlen.max(1))) as usize
    }

    /// Ask the caches for what [`Self::prefetch_window`] and [`Self::window`] of `t` read
    /// first.
    fn prefetch_tx(&self, t: u32) {
        let (p, i) = (&self.pieces[t as usize / PIECE_TX], t as usize % PIECE_TX);
        index::prefetch(p.off.as_ptr().wrapping_add(i));
    }

    /// Ask the caches for the seeds [`Self::window`] searches from (reads the transcript's
    /// offsets: call it after [`Self::prefetch_tx`] has had time to land).
    fn prefetch_window(&self, t: u32, tlen: u32, lo: u32) {
        let (p, i) = (&self.pieces[t as usize / PIECE_TX], t as usize % PIECE_TX);
        let (a, b) = (p.off[i] as usize, p.off[i + 1] as usize);
        let at = p.seeds.as_ptr().wrapping_add(a + Self::guess(b - a, tlen, lo));
        index::prefetch(at);
        index::prefetch(at.wrapping_add(8));
    }
}

/// The first of the seeds `s` (in start order) that starts at `x` or later, searched outward
/// from index `g`: a few probes next to each other when the guess is close, where a binary
/// search of a long transcript's seeds misses the caches at every step.
fn first_from(s: &[(u32, u32)], x: u32, g: usize) -> usize {
    let n = s.len();
    let g = g.min(n);
    // The answer is in `lo..=hi`.
    let (mut lo, mut hi) = (0, n);
    let mut step = 1;
    if g < n && s[g].1 < x {
        lo = g + 1;
        while lo + step - 1 < n {
            let k = lo + step - 1;
            if s[k].1 < x {
                lo = k + 1;
                step *= 2;
            } else {
                hi = k;
                break;
            }
        }
    } else {
        hi = g;
        while step <= hi {
            let k = hi - step;
            if s[k].1 < x {
                lo = k + 1;
                break;
            }
            hi = k;
            step *= 2;
        }
    }
    lo + s[lo..hi].partition_point(|y| y.1 < x)
}

/// A hash of the seeds `w` of a transcript window that windows with the same seeds elsewhere
/// along a transcript share: starts count from the first seed's. Four interleaved chains, so
/// that the multiplications of one seed do not wait for those of the last.
fn window_hash(w: &[(u32, u32)]) -> u64 {
    let s0 = w.first().map_or(0, |x| x.1);
    let mut h = [
        0x9e37_79b9_7f4a_7c15u64,
        0xc2b2_ae3d_27d4_eb4f,
        0x1656_67b1_9e37_79f9,
        0x27d4_eb2f_1656_67c5,
    ];
    let step = |h: u64, x: &(u32, u32)| {
        (h ^ (u64::from(x.0) << 32 | u64::from(x.1 - s0)))
            .wrapping_mul(0xff51_afd7_ed55_8ccd)
            .rotate_left(29)
    };
    let (quads, rest) = w.as_chunks::<4>();
    for q in quads {
        for (h, x) in h.iter_mut().zip(q) {
            *h = step(*h, x);
        }
    }
    for (h, x) in h.iter_mut().zip(rest) {
        *h = step(*h, x);
    }
    h.iter().fold(w.len() as u64, |a, &b| {
        (a.rotate_left(23) ^ b).wrapping_mul(0x9e37_79b9_7f4a_7c15)
    })
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
const DIR_BITS: u32 = 16;
/// Buckets sorted at a time while the table is built: the unsorted postings of that many
/// are all the build holds beyond its input and the table.
const GROUP: usize = 32;

/// A run of rescue seeds while the table is built: (key, transcript, start).
type Run = (u64, u32, u32);

impl Bucket {
    /// The bucket of `runs`, all of its keys. A directory cell holds a handful of runs: they
    /// are dealt out by cell first, in place and in two steps (each to few enough places for
    /// the caches), which leaves little for the comparison sort.
    fn new(runs: &mut [Run]) -> Self {
        assert!(runs.len() < u32::MAX as usize, "too many rescue seed postings");
        let low = DIR_BITS / 2;
        deal(runs, 1 << (DIR_BITS - low), |k| Self::dir_of(k) >> low);
        for part in runs.chunk_by_mut(|a, b| Self::dir_of(a.0) >> low == Self::dir_of(b.0) >> low) {
            deal(part, 1 << low, |k| Self::dir_of(k) & ((1 << low) - 1));
            for cell in part.chunk_by_mut(|a, b| Self::dir_of(a.0) == Self::dir_of(b.0)) {
                cell.sort_unstable();
            }
        }
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

/// Reorder `runs` by `digit` of their keys (below `n`), in place: each run is swapped to the
/// next free place of its digit until the one that belongs here turns up.
fn deal(runs: &mut [Run], n: usize, digit: impl Fn(u64) -> usize) {
    let mut end = vec![0usize; n];
    for r in runs.iter() {
        end[digit(r.0)] += 1;
    }
    let mut at = Vec::with_capacity(n);
    let mut sum = 0;
    for e in &mut end {
        at.push(sum);
        sum += *e;
        *e = sum;
    }
    for d in 0..n {
        while at[d] < end[d] {
            let mut r = runs[at[d]];
            loop {
                let k = digit(r.0);
                if k == d {
                    break;
                }
                std::mem::swap(&mut r, &mut runs[at[k]]);
                at[k] += 1;
            }
            runs[at[d]] = r;
            at[d] += 1;
        }
    }
}

/// Table keys looked up side by side ([`Global::get_batch`]).
const LOOKUPS: usize = 32;

impl Global {
    /// The postings of each of `keys` (at most [`LOOKUPS`]) to `f`, with the key's index. A
    /// lookup reads four places in turn, each of them rarely in the caches; taking every
    /// step for all the keys before the next lets the caches fetch them side by side.
    fn get_batch<'a>(&'a self, keys: &[u64], mut f: impl FnMut(usize, &'a [(u32, u32)])) {
        debug_assert!(keys.len() <= LOOKUPS);
        let bucket = |k: u64| &self.buckets[(k >> 56) as usize];
        for &k in keys {
            index::prefetch(bucket(k).dir.as_ptr().wrapping_add(Bucket::dir_of(k)));
        }
        // Each key's range of `keys`, then of `posts` (empty: no such key).
        let mut at = [(0u32, 0u32); LOOKUPS];
        for (&k, at) in keys.iter().zip(&mut at) {
            let (b, d) = (bucket(k), Bucket::dir_of(k));
            *at = (b.dir[d], b.dir[d + 1]);
            index::prefetch(b.keys.as_ptr().wrapping_add(at.0 as usize));
            index::prefetch(b.keys.as_ptr().wrapping_add((at.1 as usize).saturating_sub(1)));
        }
        for (&k, at) in keys.iter().zip(&mut at) {
            let b = bucket(k);
            *at = match b.keys[at.0 as usize..at.1 as usize].binary_search(&k) {
                Ok(i) => {
                    let i = at.0 as usize + i;
                    index::prefetch(b.first.as_ptr().wrapping_add(i));
                    index::prefetch(b.first.as_ptr().wrapping_add(i + 1));
                    (i as u32, 1)
                }
                Err(_) => (0, 0),
            };
        }
        for (&k, at) in keys.iter().zip(&mut at) {
            if at.1 != 0 {
                let b = bucket(k);
                *at = (b.first[at.0 as usize], b.first[at.0 as usize + 1]);
                // The first lines of the postings; a longer list is read in order.
                for i in (at.0..at.1).step_by(8).take(4) {
                    index::prefetch(b.posts.as_ptr().wrapping_add(i as usize));
                }
            }
        }
        for (i, (&k, at)) in keys.iter().zip(&at).enumerate() {
            f(i, &bucket(k).posts[at.0 as usize..at.1 as usize]);
        }
    }

    /// The table of the runs of `tk` seeds in the store's `pieces` (`hist`: each piece's runs
    /// per bucket; `tops`: the bucket of each of its runs, in order), on `nt` threads. Buckets
    /// are made [`GROUP`] at a time: the seeds are read once per group, and only one group's
    /// postings are ever held unsorted.
    fn build(pieces: &[StorePiece], hist: &[Vec<u32>], tops: Vec<Vec<u8>>, tk: usize, nt: usize) -> Self {
        let mut tops = Some(tops);
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
            let of_piece = tops.as_deref().unwrap();
            std::thread::scope(|sc| {
                for _ in 0..nt {
                    sc.spawn(|| {
                        loop {
                            let Some((p, mut part)) = work.lock().unwrap().next() else {
                                break;
                            };
                            let sp = &pieces[p];
                            let mut top = of_piece[p].iter();
                            let mut fill = [0usize; GROUP];
                            for (i, o) in sp.off.windows(2).enumerate() {
                                let t = (p * PIECE_TX + i) as u32;
                                for w in sp.seeds[o[0] as usize..o[1] as usize].windows(tk) {
                                    let b = usize::from(*top.next().unwrap()).wrapping_sub(b0);
                                    if b < GROUP {
                                        let key = tuple(w.iter().map(|x| x.0));
                                        part[b][fill[b]] = (key, t, w[0].1);
                                        fill[b] += 1;
                                    }
                                }
                            }
                        }
                    });
                }
            });
            if b0 + GROUP >= BUCKETS {
                // Read for the last time: gone before the last buckets take their memory.
                tops = None;
            }
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
    read_at: LabelTable,
    read_next: Vec<u32>,
    /// Distinct window (hash of its points, their number) -> its index in `results`.
    seen: rustc_hash::FxHashMap<(u64, u32), u32>,
    /// With stored seeds: a window's seeds ([`window_hash`], their number) -> its index in
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
        let mut made: Vec<(usize, StorePiece, Vec<u32>, Vec<u8>)> = std::thread::scope(|sc| {
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
                            let mut top = Vec::new();
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
                                        let b = (tuple(w.iter().map(|x| x.0)) >> 56) as u8;
                                        hist[usize::from(b)] += 1;
                                        top.push(b);
                                    }
                                }
                                off.push(seeds.len() as u32);
                            }
                            let piece = StorePiece {
                                off: off.clone(),
                                seeds: seeds.clone(),
                            };
                            top.shrink_to_fit();
                            out.push((p, piece, hist, top));
                        }
                        out
                    })
                })
                .collect();
            hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
        });
        made.sort_unstable_by_key(|x| x.0);
        let (mut pieces, mut hist, mut tops) = (Vec::new(), Vec::new(), Vec::new());
        for (_, sp, h, top) in made {
            pieces.push(sp);
            hist.push(h);
            tops.push(top);
        }
        if let Some(tk) = tuple_k {
            self.global = Some(Global::build(&pieces, &hist, tops, tk, nt.min(32)));
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
            let cand = Self::votes(g, table, seq.len(), n, sc);
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
        const OFF: i64 = 1 << 30;
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
                    let start = read[from + i].start;
                    keys.extend(
                        posts
                            .iter()
                            .map(|&(t, p)| (u64::from(t) << 32) | (p as i64 - start as i64 + OFF) as u64),
                    );
                }
            });
        }
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
            for &k in keys.iter() {
                let c = &mut count[slot(k)];
                *c = c.saturating_add(1);
            }
            let need = g.min_shared.min(255) as u8;
            let mut kept = 0;
            for i in 0..keys.len() {
                let k = keys[i];
                keys[kept] = k;
                kept += usize::from(count[slot(k)] >= need);
            }
            keys.truncate(kept);
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
        if keys.capacity() > VOTES_KEEP {
            *keys = Vec::new();
            sc.counts = Vec::new();
        }
        // Transcripts are distinct, so the order is total: the best `n`, then those in order.
        let by = |&(c, t, _, _): &(u32, u32, i64, i64)| (std::cmp::Reverse(c), t);
        if cand.len() > n {
            cand.select_nth_unstable_by_key(n - 1, by);
            cand.truncate(n);
        }
        cand.sort_unstable_by_key(by);
        cand
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
    fn stored_point_keys(win: &[(u32, u32)], lo: u32, sc: &mut Scratch) {
        sc.keys.clear();
        for &(label, start) in win {
            let mut i = sc.s.read_at.get(u64::from(label));
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
        let key = Some((xxhash_rust::xxh3::xxh3_64(seq), seq.len()));
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
        // The cap bounds the bases re-seeded: a window of a longer transcript is within it.
        if hi - lo > self.p.max_tx_len {
            return;
        }
        // A window within one already chained for this transcript has no better chain.
        let mut k = sc.s.by_tx.get(&t).map_or(u32::MAX, |x| x.0);
        while k != u32::MAX {
            let (w, next) = sc.s.tx_wins[k as usize];
            if w.0 <= lo && hi <= w.1 {
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
                let (skey, s0) = (
                    (window_hash(win), win.len() as u32),
                    win.first().map_or(lo, |x| x.1),
                );
                if let Some(&(r, d)) = sc.s.seen_seeds.get(&skey) {
                    let base = (i64::from(s0) + d).max(0) as u32;
                    Self::note(t, (lo, hi, base, r), sc);
                    return;
                }
                seeds = Some((skey, s0));
                Self::stored_point_keys(win, lo, sc);
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
