//! The rescue tables: every transcript's rescue seeds ([`Store`]) and the table of their runs
//! ([`Global`]). Both are a few flat arrays, made once ([`build`]) and either kept in the
//! index file, where they are read in place ([`load`]), or made when the mapper is built.
//!
//! A transcriptome repeats itself: the isoforms of a gene share most of their sequence, and
//! so most of their rescue seeds. The store keeps the seeds in chunks cut where the seeds
//! themselves say so (at a seed whose label hash has its low bits clear), so that the same
//! stretch of seeds is cut the same way in every transcript holding it, and keeps each
//! distinct chunk once. The table's posting lists, sorted by transcript, are kept as byte
//! varints of the differences.

use std::ops::Deref;
use std::sync::Arc;

use index::Index;
use seeds::BaseSeeder;

/// An array the tables own, or one inside the index's backing memory (which the table
/// holding it keeps alive).
pub(crate) struct Arr<T> {
    ptr: *const T,
    len: usize,
    _own: Vec<T>,
}

// SAFETY: an `Arr` is a read-only view of memory that outlives it and never moves: its own
// vector's buffer, or the index its table holds.
unsafe impl<T: Send + Sync> Send for Arr<T> {}
unsafe impl<T: Send + Sync> Sync for Arr<T> {}

impl<T> Arr<T> {
    fn own(v: Vec<T>) -> Self {
        Arr {
            ptr: v.as_ptr(),
            len: v.len(),
            _own: v,
        }
    }

    /// A view of `s`, without its lifetime: the caller keeps the memory alive and in place
    /// for as long as the `Arr` lives.
    pub(crate) fn of(s: &[T]) -> Self {
        Arr {
            ptr: s.as_ptr(),
            len: s.len(),
            _own: Vec::new(),
        }
    }
}

impl<T> Deref for Arr<T> {
    type Target = [T];
    #[inline]
    fn deref(&self) -> &[T] {
        // SAFETY: `ptr` and `len` are those of a slice that lives as long as `self` (see
        // `own` and `of`).
        unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
    }
}

/// A rescue seed label as the tables keep it.
pub(crate) fn mix(label: u64) -> u32 {
    (label.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 32) as u32
}

/// The table key of a run of consecutive rescue seeds, from their label hashes (the store
/// has no more of a label).
pub(crate) fn tuple(it: impl Iterator<Item = u32>) -> u64 {
    it.fold(0x9e37_79b9_7f4a_7c15, |h, l| {
        (h ^ u64::from(l))
            .wrapping_mul(0xff51_afd7_ed55_8ccd)
            .rotate_left(29)
    })
}

/// A seed whose label hash has these bits clear starts a chunk: one seed in sixteen.
const CHUNK_MASK: u32 = 15;
/// Most seeds of a chunk (a chunk's length is kept in a byte).
const CHUNK_CAP: usize = 64;

/// Every transcript's rescue seeds as (label hash, start), in start order per transcript.
///
/// A transcript is a list of chunks, each named by where its seeds are kept and where its
/// first seed starts in the transcript; a chunk's seeds are their label hashes and their
/// starts counted from the first one's. Equal chunks, wherever they occur, are kept once.
pub struct Store {
    /// Where each transcript's chunks start in `occ` (one more entry closes the last).
    tx_off: Arr<u32>,
    /// Each transcript's chunks, in order: (where the chunk's seeds start in `hash` and
    /// `rel`, the transcript position of its first seed).
    occ: Arr<[u32; 2]>,
    /// The number of seeds of each of them.
    occ_len: Arr<u8>,
    /// The seeds of the distinct chunks.
    hash: Arr<u32>,
    rel: Arr<u16>,
    _idx: Option<Arc<Index>>,
}

/// The seeds of a transcript that start in some range: a run of the transcript's chunks, the
/// first and last of them in part.
#[derive(Clone, Copy)]
pub(crate) struct Window<'a> {
    st: &'a Store,
    /// The chunks holding the window's seeds are `occ[a..b]`;
    a: usize,
    b: usize,
    /// the window skips the first `skip` seeds of the first one, and takes those before
    /// `last` of the last one.
    skip: usize,
    last: usize,
}

impl Store {
    /// The seeds of transcript `t` (of `tlen` parsed bases) that start in `[lo, hi)`.
    pub(crate) fn window(&self, t: u32, tlen: u32, lo: u32, hi: u32) -> Window<'_> {
        let (o0, o1) = (
            self.tx_off[t as usize] as usize,
            self.tx_off[t as usize + 1] as usize,
        );
        let occ = &self.occ[o0..o1];
        // The first seed starting at `x` or later, as (chunk, seed of the chunk), searched
        // from chunk `g`: it is in the last chunk that starts before `x`, or opens the next.
        let at = |x: u32, g: usize| -> (usize, usize) {
            let j = first_from(occ, x, g);
            if j > 0 {
                let [a, pos] = occ[j - 1];
                let n = usize::from(self.occ_len[o0 + j - 1]);
                let k = self.rel[a as usize..a as usize + n].partition_point(|&r| pos + u32::from(r) < x);
                if k < n {
                    return (j - 1, k);
                }
            }
            (j, 0)
        };
        let a = at(lo, Self::guess(occ.len(), tlen, lo));
        let b = at(hi, Self::guess(occ.len(), tlen, hi).max(a.0));
        if b <= a {
            return Window {
                st: self,
                a: o0 + a.0,
                b: o0 + a.0,
                skip: 0,
                last: 0,
            };
        }
        let (b, last) = if b.1 > 0 {
            (b.0 + 1, b.1)
        } else {
            (b.0, usize::from(self.occ_len[o0 + b.0 - 1]))
        };
        Window {
            st: self,
            a: o0 + a.0,
            b: o0 + b,
            skip: a.1,
            last,
        }
    }

    /// About where, among a transcript's `n` chunks, the one holding `pos` is: seeds are
    /// spread nearly evenly along a transcript.
    fn guess(n: usize, tlen: u32, pos: u32) -> usize {
        (u64::from(pos) * n as u64 / u64::from(tlen.max(1))) as usize
    }

    /// Ask the caches for what [`Self::prefetch_window`] and [`Self::window`] of `t` read
    /// first.
    pub(crate) fn prefetch_tx(&self, t: u32) {
        index::prefetch(self.tx_off.as_ptr().wrapping_add(t as usize));
    }

    /// Ask the caches for the chunks [`Self::window`] searches from (reads the transcript's
    /// offsets: call it after [`Self::prefetch_tx`] has had time to land).
    pub(crate) fn prefetch_window(&self, t: u32, tlen: u32, lo: u32) {
        let (a, b) = (
            self.tx_off[t as usize] as usize,
            self.tx_off[t as usize + 1] as usize,
        );
        let g = a + Self::guess(b - a, tlen, lo);
        let at = self.occ.as_ptr().wrapping_add(g);
        index::prefetch(at);
        index::prefetch(at.wrapping_add(8));
        index::prefetch(self.occ_len.as_ptr().wrapping_add(g));
    }

    /// The seeds of transcript `t` as (label hash, start), in start order.
    pub(crate) fn seeds_into(&self, t: u32, out: &mut Vec<(u32, u32)>) {
        out.clear();
        for i in self.tx_off[t as usize] as usize..self.tx_off[t as usize + 1] as usize {
            let [at, pos] = self.occ[i];
            let (at, n) = (at as usize, usize::from(self.occ_len[i]));
            out.extend(
                self.hash[at..at + n]
                    .iter()
                    .zip(&self.rel[at..at + n])
                    .map(|(&h, &r)| (h, pos + u32::from(r))),
            );
        }
    }

    pub fn n_transcripts(&self) -> usize {
        self.tx_off.len().saturating_sub(1)
    }

    /// Bytes of the store's arrays.
    pub fn bytes(&self) -> usize {
        4 * self.tx_off.len() + 9 * self.occ.len() + 6 * self.hash.len()
    }
}

/// The first of the chunks `s` (in order) that starts at `x` or later, searched outward from
/// index `g`: a few probes next to each other when the guess is close, where a binary search
/// of a long transcript's chunks misses the caches at every step.
fn first_from(s: &[[u32; 2]], x: u32, g: usize) -> usize {
    let n = s.len();
    let g = g.min(n);
    // The answer is in `lo..=hi`.
    let (mut lo, mut hi) = (0, n);
    let mut step = 1;
    if g < n && s[g][1] < x {
        lo = g + 1;
        while lo + step - 1 < n {
            let k = lo + step - 1;
            if s[k][1] < x {
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
            if s[k][1] < x {
                lo = k + 1;
                break;
            }
            hi = k;
            step *= 2;
        }
    }
    lo + s[lo..hi].partition_point(|y| y[1] < x)
}

/// Chunks whose seeds [`Window::for_each`] asks the caches for ahead of the one it reads.
const AHEAD: usize = 3;

impl Window<'_> {
    pub(crate) fn is_empty(&self) -> bool {
        self.a == self.b
    }

    /// The seeds of chunk `i` that are the window's.
    #[inline]
    fn part(&self, i: usize) -> (usize, usize) {
        (
            if i == self.a { self.skip } else { 0 },
            if i + 1 == self.b {
                self.last
            } else {
                usize::from(self.st.occ_len[i])
            },
        )
    }

    /// Where the window's first seed starts (None: the window is empty).
    pub(crate) fn first_start(&self) -> Option<u32> {
        (!self.is_empty()).then(|| {
            let [at, pos] = self.st.occ[self.a];
            pos + u32::from(self.st.rel[at as usize + self.skip])
        })
    }

    /// The number of chunks the window's seeds are in.
    pub(crate) fn chunks(&self) -> usize {
        self.b - self.a
    }

    /// A hash of the window's seeds, their starts counted from `s0` (the first one's), that
    /// windows with the same seeds elsewhere along a transcript mostly share. A whole chunk
    /// counts by where its seeds are kept, which names them: the seeds of the window's two
    /// ends are all it reads. (Equal seeds cut into chunks differently have different
    /// hashes: equal hashes mean equal seeds, not the reverse.)
    pub(crate) fn key(&self, s0: u32) -> u64 {
        let st = self.st;
        let mut h = 0x9e37_79b9_7f4a_7c15u64;
        for i in self.a..self.b {
            let [at, pos] = st.occ[i];
            let whole = (i != self.a || self.skip == 0)
                && (i + 1 != self.b || self.last == usize::from(st.occ_len[i]));
            if whole {
                h = (h ^ (u64::from(at) << 32 | u64::from(pos - s0)))
                    .wrapping_mul(0xff51_afd7_ed55_8ccd)
                    .rotate_left(29);
            } else {
                let (lo, hi) = self.part(i);
                let at = at as usize;
                h = h.wrapping_add(0x1656_67b1_9e37_79f9);
                for (&l, &r) in st.hash[at + lo..at + hi].iter().zip(&st.rel[at + lo..at + hi]) {
                    h = (h ^ (u64::from(l) << 32 | u64::from(pos + u32::from(r) - s0)))
                        .wrapping_mul(0xc2b2_ae3d_27d4_eb4f)
                        .rotate_left(31);
                }
            }
        }
        h
    }

    /// Ask the caches for the seeds of chunk `i`, if it is one of the window's.
    #[inline]
    fn prefetch(&self, i: usize) {
        if i < self.b {
            let at = self.st.occ[i][0] as usize;
            let (h, r) = (
                self.st.hash.as_ptr().wrapping_add(at),
                self.st.rel.as_ptr().wrapping_add(at),
            );
            index::prefetch(h);
            index::prefetch(h.wrapping_add(16));
            index::prefetch(r);
            index::prefetch(r.wrapping_add(32));
        }
    }

    /// The window's seeds as (label hash, start), in order.
    #[inline]
    pub(crate) fn for_each(&self, mut f: impl FnMut(u32, u32)) {
        let st = self.st;
        for i in self.a..self.a + AHEAD {
            self.prefetch(i);
        }
        for i in self.a..self.b {
            self.prefetch(i + AHEAD);
            let [at, pos] = st.occ[i];
            let (lo, hi) = self.part(i);
            let at = at as usize;
            for (&l, &r) in st.hash[at + lo..at + hi].iter().zip(&st.rel[at + lo..at + hi]) {
                f(l, pos + u32::from(r));
            }
        }
    }
}

/// Every transcript's runs of consecutive rescue seeds: tuple key -> (transcript, start of
/// the run's first seed).
pub struct Global {
    /// Where the keys of each value of their top bits start in `keys` (one more entry
    /// closes the last); a power of two of them, about one per key.
    dir: Arr<u32>,
    /// What to shift a key by for its top bits.
    shift: u32,
    /// The distinct keys, in order.
    keys: Arr<u64>,
    /// Where each key's postings start in `posts`, counted from where those of its top eight
    /// bits start (`base`; one more entry closes the last).
    first: Arr<u32>,
    base: Arr<u64>,
    /// Per key: the number of its postings, then each as (transcript, start), in that order;
    /// the transcript is kept as its difference from the one before (from zero, for the
    /// first). All are byte varints (seven bits a byte, low bits first, the top bit set on
    /// all but the last). [`PAD`] zero bytes close the array.
    posts: Arr<u8>,
    _idx: Option<Arc<Index>>,
}

/// Zero bytes after the last posting: a decoder may load a word at any posting.
const PAD: usize = 8;

const BUCKETS: usize = 256;
/// Bits (below the top eight) the runs of a bucket are dealt out by before they are sorted.
const DEAL_BITS: u32 = 16;
/// Buckets sorted at a time while the table is built: the unsorted postings of that many
/// are all the build holds beyond its input and the table.
const GROUP: usize = 32;

/// A run of rescue seeds while the table is built: (key, transcript, start).
type Run = (u64, u32, u32);

/// Append `x` as a byte varint.
pub(crate) fn put(out: &mut Vec<u8>, mut x: u32) {
    while x >= 0x80 {
        out.push(x as u8 | 0x80);
        x >>= 7;
    }
    out.push(x as u8);
}

/// The byte varint at `b[p..]` and where the next one starts (past the end of `b`, if it
/// runs out).
#[inline]
pub(crate) fn varint(b: &[u8], mut p: usize) -> (u32, usize) {
    let (mut x, mut sh) = (0u32, 0);
    loop {
        let Some(&c) = b.get(p) else {
            return (x, b.len() + 1);
        };
        p += 1;
        x |= u32::from(c & 0x7f) << sh;
        if c < 0x80 || sh >= 28 {
            return (x, p);
        }
        sh += 7;
    }
}

/// The postings of one key.
#[derive(Clone, Copy)]
pub(crate) struct Posts<'a> {
    n: usize,
    /// From the first posting on.
    bytes: &'a [u8],
}

impl Posts<'_> {
    pub(crate) fn len(&self) -> usize {
        self.n
    }

    /// Each posting as (transcript, start), in order.
    #[inline]
    pub(crate) fn for_each(&self, mut f: impl FnMut(u32, u32)) {
        let b = self.bytes;
        let (mut p, mut t) = (0usize, 0u32);
        for _ in 0..self.n {
            let Some(w) = b.get(p..p + 4) else {
                break;
            };
            let w = u32::from_le_bytes([w[0], w[1], w[2], w[3]]);
            let s;
            if w & 0x0080_8080 == 0x0000_8000 {
                // The usual posting: a transcript close to the last one in one byte, and a
                // start in two.
                t = t.wrapping_add(w & 0x7f);
                s = (w >> 8) & 0x7f | (w >> 9) & 0x3f80;
                p += 3;
            } else {
                let (dt, q) = varint(b, p);
                let (x, q) = varint(b, q);
                t = t.wrapping_add(dt);
                s = x;
                p = q;
            }
            f(t, s);
        }
    }
}

/// One value of the keys' top eight bits while the table is built: its keys in order, where
/// each one's postings start, and the postings.
struct Part {
    keys: Vec<u64>,
    first: Vec<u32>,
    posts: Vec<u8>,
}

impl Part {
    /// The part of `runs`, all of its keys. A cell of the dealing holds a handful of runs:
    /// they are dealt out by cell first, in place and in two steps (each to few enough
    /// places for the caches), which leaves little for the comparison sort.
    fn new(runs: &mut [Run]) -> Self {
        let cell = |k: u64| (k >> (56 - DEAL_BITS)) as usize & ((1 << DEAL_BITS) - 1);
        let low = DEAL_BITS / 2;
        deal(runs, 1 << (DEAL_BITS - low), |k| cell(k) >> low);
        for part in runs.chunk_by_mut(|a, b| cell(a.0) >> low == cell(b.0) >> low) {
            deal(part, 1 << low, |k| cell(k) & ((1 << low) - 1));
            for c in part.chunk_by_mut(|a, b| cell(a.0) == cell(b.0)) {
                c.sort_unstable();
            }
        }
        let (mut keys, mut first) = (Vec::new(), Vec::new());
        let mut posts = Vec::with_capacity(runs.len() * 7 / 2);
        for run in runs.chunk_by(|a, b| a.0 == b.0) {
            assert!(
                posts.len() < u32::MAX as usize && run.len() < u32::MAX as usize,
                "too many rescue seed postings"
            );
            keys.push(run[0].0);
            first.push(posts.len() as u32);
            put(&mut posts, run.len() as u32);
            let mut last = 0;
            for &(_, t, s) in run {
                put(&mut posts, t - last);
                put(&mut posts, s);
                last = t;
            }
        }
        keys.shrink_to_fit();
        first.shrink_to_fit();
        posts.shrink_to_fit();
        Part { keys, first, posts }
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
pub(crate) const LOOKUPS: usize = 32;

/// Bits of the key directory for `n` keys: about one key per entry.
fn dir_bits(n: usize) -> u32 {
    (usize::BITS - n.max(2).leading_zeros() - 1).clamp(8, 30)
}

impl Global {
    /// The postings of each of `keys` (at most [`LOOKUPS`]) to `f`, with the key's index. A
    /// lookup reads four places in turn, each of them rarely in the caches; taking every
    /// step for all the keys before the next lets the caches fetch them side by side.
    pub(crate) fn get_batch<'a>(&'a self, keys: &[u64], mut f: impl FnMut(usize, Posts<'a>)) {
        debug_assert!(keys.len() <= LOOKUPS);
        let cell = |k: u64| (k >> self.shift) as usize;
        for &k in keys {
            index::prefetch(self.dir.as_ptr().wrapping_add(cell(k)));
        }
        // Each key's range of `keys`, then its place in them and whether it is there, then
        // where its postings start.
        let mut at = [(0usize, 0usize); LOOKUPS];
        for (&k, at) in keys.iter().zip(&mut at) {
            let d = cell(k);
            *at = (self.dir[d] as usize, self.dir[d + 1] as usize);
            index::prefetch(self.keys.as_ptr().wrapping_add(at.0));
            index::prefetch(self.keys.as_ptr().wrapping_add(at.1.saturating_sub(1)));
        }
        for (&k, at) in keys.iter().zip(&mut at) {
            *at = match self.keys.get(at.0..at.1).unwrap_or(&[]).binary_search(&k) {
                Ok(i) => {
                    index::prefetch(self.first.as_ptr().wrapping_add(at.0 + i));
                    (at.0 + i, 1)
                }
                Err(_) => (0, 0),
            };
        }
        for (&k, at) in keys.iter().zip(&mut at) {
            if at.1 != 0 {
                at.0 = (self.base[(k >> 56) as usize] + u64::from(self.first[at.0])) as usize;
                // The first lines of the postings; a longer list is read in order.
                let p = self.posts.as_ptr().wrapping_add(at.0);
                index::prefetch(p);
                index::prefetch(p.wrapping_add(64));
            }
        }
        for (i, at) in at[..keys.len()].iter().enumerate() {
            let none = Posts { n: 0, bytes: &[] };
            f(
                i,
                if at.1 != 0 {
                    let (n, p) = varint(&self.posts, at.0);
                    self.posts
                        .get(p..)
                        .map_or(none, |bytes| Posts { n: n as usize, bytes })
                } else {
                    none
                },
            );
        }
    }

    /// The table of the runs of `tk` seeds of `store`'s transcripts, taken `PIECE_TX` at a
    /// time (`hist`: the runs per bucket of each such piece; `tops`: the bucket of each of
    /// its runs, in order), on `nt` threads. Buckets are made [`GROUP`] at a time: the seeds
    /// are read once per group, and only one group's postings are ever held unsorted.
    fn build(store: &Store, hist: &[Vec<u32>], tops: Vec<Vec<u8>>, tk: usize, nt: usize) -> Self {
        let n = store.n_transcripts();
        let mut tops = Some(tops);
        let size = |b: usize| hist.iter().map(|h| h[b] as usize).sum::<usize>();
        let most = (0..BUCKETS)
            .step_by(GROUP)
            .map(|b0| (b0..b0 + GROUP).map(size).sum())
            .max()
            .unwrap_or(0);
        let mut runs: Vec<Run> = vec![(0, 0, 0); most];
        let mut parts: Vec<(usize, Part)> = Vec::with_capacity(BUCKETS);
        for b0 in (0..BUCKETS).step_by(GROUP) {
            // Each bucket's runs, cut into the parts the pieces fill.
            let mut of_bucket: Vec<&mut [Run]> = Vec::with_capacity(GROUP);
            let mut rest = runs.as_mut_slice();
            for b in b0..b0 + GROUP {
                let (a, r) = rest.split_at_mut(size(b));
                of_bucket.push(a);
                rest = r;
            }
            let mut fills: Vec<Vec<&mut [Run]>> = hist.iter().map(|_| Vec::with_capacity(GROUP)).collect();
            for (b, v) in of_bucket.iter_mut().enumerate() {
                let mut rest = &mut **v;
                for (h, fill) in hist.iter().zip(&mut fills) {
                    let (a, r) = rest.split_at_mut(h[b0 + b] as usize);
                    fill.push(a);
                    rest = r;
                }
            }
            let work = std::sync::Mutex::new(fills.into_iter().enumerate());
            let of_piece = tops.as_deref().unwrap();
            std::thread::scope(|sc| {
                for _ in 0..nt {
                    sc.spawn(|| {
                        let mut seeds = Vec::new();
                        loop {
                            let Some((p, mut fill)) = work.lock().unwrap().next() else {
                                break;
                            };
                            let mut top = of_piece[p].iter();
                            let mut filled = [0usize; GROUP];
                            for t in p * PIECE_TX..((p + 1) * PIECE_TX).min(n) {
                                store.seeds_into(t as u32, &mut seeds);
                                for w in seeds.windows(tk) {
                                    let b = usize::from(*top.next().unwrap()).wrapping_sub(b0);
                                    if b < GROUP {
                                        let key = tuple(w.iter().map(|x| x.0));
                                        fill[b][filled[b]] = (key, t as u32, w[0].1);
                                        filled[b] += 1;
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
            let made: Vec<(usize, Part)> = std::thread::scope(|sc| {
                let hs: Vec<_> = (0..nt.min(GROUP))
                    .map(|_| {
                        sc.spawn(|| {
                            let mut out = Vec::new();
                            loop {
                                let Some((b, runs)) = work.lock().unwrap().next() else {
                                    break;
                                };
                                out.push((b0 + b, Part::new(runs)));
                            }
                            out
                        })
                    })
                    .collect();
                hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
            });
            parts.extend(made);
        }
        drop(runs);
        parts.sort_unstable_by_key(|x| x.0);
        let nk: usize = parts.iter().map(|x| x.1.keys.len()).sum();
        let nb: usize = parts.iter().map(|x| x.1.posts.len()).sum();
        assert!(nk < u32::MAX as usize, "too many rescue seed runs");
        let (mut keys, mut first) = (Vec::with_capacity(nk), Vec::with_capacity(nk));
        let mut posts = Vec::with_capacity(nb + PAD);
        let mut base = Vec::with_capacity(BUCKETS + 1);
        for (_, part) in parts {
            base.push(posts.len() as u64);
            keys.extend_from_slice(&part.keys);
            first.extend_from_slice(&part.first);
            posts.extend_from_slice(&part.posts);
        }
        base.push(posts.len() as u64);
        posts.resize(posts.len() + PAD, 0);
        let bits = dir_bits(nk);
        let mut dir = vec![0u32; (1 << bits) + 1];
        for &k in &keys {
            dir[(k >> (64 - bits)) as usize + 1] += 1;
        }
        for d in 0..1 << bits {
            dir[d + 1] += dir[d];
        }
        Global {
            dir: Arr::own(dir),
            shift: 64 - bits,
            keys: Arr::own(keys),
            first: Arr::own(first),
            base: Arr::own(base),
            posts: Arr::own(posts),
            _idx: None,
        }
    }

    /// Bytes of the table's arrays.
    pub fn bytes(&self) -> usize {
        4 * self.dir.len() + 12 * self.keys.len() + 8 * self.base.len() + self.posts.len()
    }
}

/// Transcripts per piece of the tables' build (a power of two).
const PIECE_TX: usize = 256;

/// What a builder thread makes of one piece: its distinct chunks, in order of first
/// appearance, and its transcripts as lists of them.
#[derive(Default)]
struct Piece {
    /// The distinct chunks: each one's seeds (`len` of them, in `hash` and `rel` one after
    /// the other) and a hash of them.
    len: Vec<u8>,
    key: Vec<u64>,
    hash: Vec<u32>,
    rel: Vec<u16>,
    /// The chunks of each transcript, and each as (distinct chunk, transcript position).
    tx_n: Vec<u32>,
    occ: Vec<(u32, u32)>,
    /// The piece's runs per bucket of the table, and the bucket of each run in order.
    hist: Vec<u32>,
    top: Vec<u8>,
}

/// A hash of a chunk's seeds, starts counted from the first one's.
fn chunk_key(s: &[(u32, u32)]) -> u64 {
    let s0 = s[0].1;
    s.iter().fold(s.len() as u64 ^ 0xc2b2_ae3d_27d4_eb4f, |h, x| {
        (h ^ (u64::from(x.0) << 16 | u64::from(x.1 - s0)))
            .wrapping_mul(0xff51_afd7_ed55_8ccd)
            .rotate_left(29)
    })
}

/// Where a chunk's hash is looked up next, when its place is taken by another chunk.
fn next_key(k: u64) -> u64 {
    k.wrapping_mul(0x9e37_79b9_7f4a_7c15).wrapping_add(1)
}

/// The rescue tables of `n` transcripts, the seeds being `seeder`'s: the store, and with
/// `tuple` the table of the runs of that many seeds. `codes(t, out)` appends the parsed
/// transcript `t` to `out` as 2-bit codes. Made on `nt` threads; the result does not depend
/// on their number.
pub fn build(
    seeder: &BaseSeeder,
    n: usize,
    codes: &(impl Fn(u32, &mut Vec<u8>) + Sync),
    tuple_k: Option<usize>,
    nt: usize,
) -> (Store, Option<Global>) {
    let nt = nt.max(1);
    // Pieces of PIECE_TX transcripts, taken by whichever thread is free.
    let next = std::sync::atomic::AtomicUsize::new(0);
    let n_pieces = n.div_ceil(PIECE_TX);
    let mut made: Vec<(usize, Piece)> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..nt)
            .map(|_| {
                sc.spawn(|| {
                    let mut out = Vec::new();
                    let (mut tseq, mut sd) = (Vec::new(), Vec::new());
                    let mut seeds: Vec<(u32, u32)> = Vec::new();
                    // Chunk hash -> (distinct chunk, where its seeds start in the piece).
                    let mut ids: rustc_hash::FxHashMap<u64, (u32, u32)> = Default::default();
                    loop {
                        let p = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if p >= n_pieces {
                            break;
                        }
                        let mut pc = Piece {
                            hist: vec![0u32; BUCKETS],
                            ..Piece::default()
                        };
                        ids.clear();
                        for t in p * PIECE_TX..((p + 1) * PIECE_TX).min(n) {
                            tseq.clear();
                            codes(t as u32, &mut tseq);
                            seeder.seed_starts_codes_into(&tseq, &mut sd);
                            if !sd.is_sorted_by_key(|x| x.1) {
                                sd.sort_unstable_by_key(|x| x.1);
                            }
                            seeds.clear();
                            seeds.extend(sd.iter().map(|&(l, s)| (mix(l), s)));
                            if let Some(tk) = tuple_k {
                                for w in seeds.windows(tk) {
                                    let b = (tuple(w.iter().map(|x| x.0)) >> 56) as u8;
                                    pc.hist[usize::from(b)] += 1;
                                    pc.top.push(b);
                                }
                            }
                            let before = pc.occ.len();
                            let mut a = 0;
                            while a < seeds.len() {
                                let mut b = a + 1;
                                while b < seeds.len()
                                    && b - a < CHUNK_CAP
                                    && seeds[b].0 & CHUNK_MASK != 0
                                    && seeds[b].1 - seeds[a].1 <= u32::from(u16::MAX)
                                {
                                    b += 1;
                                }
                                let (s, s0) = (&seeds[a..b], seeds[a].1);
                                let same = |pc: &Piece, at: usize| {
                                    s.iter()
                                        .zip(pc.hash[at..].iter().zip(&pc.rel[at..]))
                                        .all(|(x, (&h, &r))| x.0 == h && x.1 - s0 == u32::from(r))
                                };
                                let mut k = chunk_key(s);
                                let key = k;
                                let id = loop {
                                    match ids.get(&k) {
                                        Some(&(id, at))
                                            if usize::from(pc.len[id as usize]) == s.len()
                                                && same(&pc, at as usize) =>
                                        {
                                            break id;
                                        }
                                        Some(_) => k = next_key(k),
                                        None => {
                                            let id = pc.len.len() as u32;
                                            ids.insert(k, (id, pc.hash.len() as u32));
                                            pc.len.push(s.len() as u8);
                                            pc.key.push(key);
                                            pc.hash.extend(s.iter().map(|x| x.0));
                                            pc.rel.extend(s.iter().map(|x| (x.1 - s0) as u16));
                                            break id;
                                        }
                                    }
                                };
                                pc.occ.push((id, s0));
                                a = b;
                            }
                            pc.tx_n.push((pc.occ.len() - before) as u32);
                        }
                        pc.hash.shrink_to_fit();
                        pc.rel.shrink_to_fit();
                        pc.occ.shrink_to_fit();
                        pc.top.shrink_to_fit();
                        out.push((p, pc));
                    }
                    out
                })
            })
            .collect();
        hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
    });
    made.sort_unstable_by_key(|x| x.0);
    // The pieces' chunks into one store, in order: a chunk is kept where it first appears.
    let n_occ: usize = made.iter().map(|x| x.1.occ.len()).sum();
    assert!(n_occ < u32::MAX as usize, "too many rescue seeds");
    let mut tx_off = Vec::with_capacity(n + 1);
    let mut occ: Vec<[u32; 2]> = Vec::with_capacity(n_occ);
    let mut occ_len: Vec<u8> = Vec::with_capacity(n_occ);
    let (mut hash, mut rel): (Vec<u32>, Vec<u16>) = (Vec::new(), Vec::new());
    // Chunk hash -> (where its seeds start, their number).
    let mut kept: rustc_hash::FxHashMap<u64, (u32, u8)> = Default::default();
    let (mut hist, mut tops) = (Vec::with_capacity(made.len()), Vec::with_capacity(made.len()));
    let mut at_of: Vec<u32> = Vec::new();
    for (_, pc) in made {
        at_of.clear();
        let mut o = 0;
        for (&len, &key) in pc.len.iter().zip(&pc.key) {
            let n = usize::from(len);
            let (h, r) = (&pc.hash[o..o + n], &pc.rel[o..o + n]);
            o += n;
            let mut k = key;
            at_of.push(loop {
                match kept.get(&k) {
                    Some(&(at, l))
                        if l == len && hash[at as usize..][..n] == *h && rel[at as usize..][..n] == *r =>
                    {
                        break at;
                    }
                    Some(_) => k = next_key(k),
                    None => {
                        assert!(hash.len() < u32::MAX as usize, "too many rescue seeds");
                        let at = hash.len() as u32;
                        kept.insert(k, (at, len));
                        hash.extend_from_slice(h);
                        rel.extend_from_slice(r);
                        break at;
                    }
                }
            });
        }
        let mut o = 0;
        for &m in &pc.tx_n {
            tx_off.push(occ.len() as u32);
            for &(id, pos) in &pc.occ[o..o + m as usize] {
                occ.push([at_of[id as usize], pos]);
                occ_len.push(pc.len[id as usize]);
            }
            o += m as usize;
        }
        hist.push(pc.hist);
        tops.push(pc.top);
    }
    tx_off.push(occ.len() as u32);
    drop(kept);
    hash.shrink_to_fit();
    rel.shrink_to_fit();
    let store = Store {
        tx_off: Arr::own(tx_off),
        occ: Arr::own(occ),
        occ_len: Arr::own(occ_len),
        hash: Arr::own(hash),
        rel: Arr::own(rel),
        _idx: None,
    };
    let global = tuple_k.map(|tk| Global::build(&store, &hist, tops, tk, nt.min(32)));
    (store, global)
}

/// The layout of the tables in an index: a reader takes only tables of its own.
const FORMAT: &str = "1";

/// The names of the tables' sections in an index.
mod section {
    pub const META: &str = "rescue_meta";
    pub const TX_OFF: &str = "rescue_tx_off";
    pub const OCC: &str = "rescue_occ";
    pub const OCC_LEN: &str = "rescue_occ_len";
    pub const HASH: &str = "rescue_hash";
    pub const REL: &str = "rescue_rel";
    pub const DIR: &str = "rescue_dir";
    pub const KEYS: &str = "rescue_keys";
    pub const FIRST: &str = "rescue_first";
    pub const BASE: &str = "rescue_base";
    pub const POSTS: &str = "rescue_posts";
}

/// Make the rescue tables of `data`'s transcripts, for the rescue scheme `spec` and runs of
/// `tuple_k` seeds, and keep them in it as sections of the index: a mapper of that scheme
/// then reads them in place ([`load`]) instead of making them each time it is built.
pub fn build_into(data: &mut index::IndexData, spec: &str, tuple_k: usize, nt: usize) -> anyhow::Result<()> {
    // Sequences are already parsed (HPC if the index used it): no HPC here.
    let seeder = BaseSeeder::new(spec, false)?;
    let tuple_k = tuple_k.max(1);
    let (store, global) = {
        let d = &*data;
        let codes = |t: u32, out: &mut Vec<u8>| d.tx_codes_into(t, out);
        build(&seeder, d.tx_plen.len(), &codes, Some(tuple_k), nt)
    };
    let meta = format!("format={FORMAT}\nspec={spec}\ntuple={tuple_k}\n");
    data.add_section(section::META, meta.as_bytes());
    data.add_section(section::TX_OFF, &store.tx_off);
    data.add_section(section::OCC, &store.occ);
    data.add_section(section::OCC_LEN, &store.occ_len);
    data.add_section(section::HASH, &store.hash);
    data.add_section(section::REL, &store.rel);
    drop(store);
    if let Some(g) = global {
        data.add_section(section::DIR, &g.dir);
        data.add_section(section::KEYS, &g.keys);
        data.add_section(section::FIRST, &g.first);
        data.add_section(section::BASE, &g.base);
        data.add_section(section::POSTS, &g.posts);
    }
    Ok(())
}

/// The tables kept in `idx` ([`build_into`]), read in place, if they are those of the rescue
/// scheme `spec` (and, with `tuple_k`, of runs of that many seeds; without it only the store
/// is taken). `None` if the index has none, or other ones, or their sizes do not agree with
/// each other: the caller then makes its own. Values are not checked: a lookup in damaged
/// tables may panic or find nonsense, like one in the rest of a damaged index.
pub fn load(idx: &Arc<Index>, spec: &str, tuple_k: Option<usize>) -> Option<(Store, Option<Global>)> {
    let meta = std::str::from_utf8(idx.section::<u8>(section::META)?).ok()?;
    let field = |k: &str| {
        meta.lines()
            .filter_map(|l| l.split_once('='))
            .find(|x| x.0 == k)
            .map(|x| x.1)
    };
    if field("format")? != FORMAT || field("spec")? != spec {
        return None;
    }
    let store = Store {
        tx_off: Arr::of(idx.section(section::TX_OFF)?),
        occ: Arr::of(idx.section(section::OCC)?),
        occ_len: Arr::of(idx.section(section::OCC_LEN)?),
        hash: Arr::of(idx.section(section::HASH)?),
        rel: Arr::of(idx.section(section::REL)?),
        _idx: Some(idx.clone()),
    };
    let n_occ = store.occ.len();
    let whole = store.tx_off.len() == idx.n_transcripts() + 1
        && store.tx_off[0] == 0
        && store.tx_off.windows(2).all(|w| w[0] <= w[1])
        && store.tx_off.last().is_some_and(|&x| x as usize == n_occ)
        && store.occ_len.len() == n_occ
        && store.rel.len() == store.hash.len();
    if !whole {
        return None;
    }
    let Some(tuple_k) = tuple_k else {
        return Some((store, None));
    };
    if field("tuple")?.parse() != Ok(tuple_k.max(1)) {
        return None;
    }
    let dir: &[u32] = idx.section(section::DIR)?;
    let cells = dir.len().checked_sub(1)?;
    let global = Global {
        dir: Arr::of(dir),
        shift: 64 - cells.max(1).trailing_zeros(),
        keys: Arr::of(idx.section(section::KEYS)?),
        first: Arr::of(idx.section(section::FIRST)?),
        base: Arr::of(idx.section(section::BASE)?),
        posts: Arr::of(idx.section(section::POSTS)?),
        _idx: Some(idx.clone()),
    };
    let n_keys = global.keys.len();
    let whole = cells.is_power_of_two()
        && (8..=30).contains(&cells.trailing_zeros())
        && global.dir[cells] as usize == n_keys
        && global.first.len() == n_keys
        && global.base.len() == BUCKETS + 1
        && global.base[BUCKETS] as usize + PAD == global.posts.len();
    whole.then_some((store, Some(global)))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tables kept in an index are the tables a mapper would make of it, and are taken only
    /// by a mapper of their rescue scheme.
    #[test]
    fn tables_kept_in_an_index_are_read_back() {
        let txs: Vec<index::TxInput> = transcripts()
            .iter()
            .enumerate()
            .filter(|(_, s)| s.len() > 100)
            .take(200)
            .map(|(i, s)| index::TxInput {
                name: format!("t{i}"),
                gene: format!("g{}", i % 4),
                seq: s.iter().map(|&c| b"ACGT"[c as usize]).collect(),
            })
            .collect();
        let params = index::BuildParams {
            spec: "lcp:K=3".into(),
            hpc: true,
            f_max: 20,
            allow_collisions: false,
            collapse: index::Collapse::Gene,
            infer_genes: None,
        };
        let (mut data, _) = index::build(txs, &params).unwrap();
        let plain = Arc::new(Index::from_data(&data).unwrap());
        assert!(load(&plain, "lcp:K=2", Some(3)).is_none());
        build_into(&mut data, "lcp:K=2", 3, 3).unwrap();
        let idx = Arc::new(Index::from_data(&data).unwrap());
        assert_eq!(idx.n_transcripts(), plain.n_transcripts());

        let seeder = BaseSeeder::new("lcp:K=2", false).unwrap();
        let codes = |t: u32, out: &mut Vec<u8>| idx.tx_codes_into(t, 0, idx.tx_plen()[t as usize], out);
        let (store, global) = build(&seeder, idx.n_transcripts(), &codes, Some(3), 2);
        let global = global.unwrap();
        let (s, g) = load(&idx, "lcp:K=2", Some(3)).unwrap();
        let g = g.unwrap();
        assert!(!store.hash.is_empty() && !global.keys.is_empty());
        assert_eq!(*s.tx_off, *store.tx_off);
        assert_eq!(*s.occ, *store.occ);
        assert_eq!(*s.occ_len, *store.occ_len);
        assert_eq!(*s.hash, *store.hash);
        assert_eq!(*s.rel, *store.rel);
        assert_eq!((*g.dir).to_vec(), (*global.dir).to_vec());
        assert_eq!(g.shift, global.shift);
        assert_eq!(*g.keys, *global.keys);
        assert_eq!(*g.first, *global.first);
        assert_eq!(*g.base, *global.base);
        assert_eq!(*g.posts, *global.posts);
        // The tables outlive the caller's hold on the index.
        drop(idx);
        assert_eq!(*s.hash, *store.hash);
        assert_eq!(*g.posts, *global.posts);

        // Another scheme or run length makes its own; the store alone serves any run length.
        let idx = Arc::new(Index::from_data(&data).unwrap());
        assert!(load(&idx, "lcp:K=1", Some(3)).is_none());
        assert!(load(&idx, "lcp:K=2", Some(2)).is_none());
        assert!(load(&idx, "lcp:K=2", None).is_some_and(|x| x.1.is_none()));
    }

    /// Transcripts that share stretches of sequence, as isoforms do, as 2-bit codes.
    fn transcripts() -> Vec<Vec<u8>> {
        let mut x = 0x2545_f491_4f6c_dd1du64;
        let mut rand = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    (x >> 62) as u8
                })
                .collect()
        };
        let exons: Vec<Vec<u8>> = (0..40).map(|i| rand(150 + 97 * (i % 7))).collect();
        let mut txs = Vec::new();
        for t in 0..700usize {
            let mut s = Vec::new();
            // A gene is ten exons; its isoforms skip some of them.
            let g = t % 4;
            for (e, exon) in exons[g * 10..(g + 1) * 10].iter().enumerate() {
                if (t / 4 + e) % 5 != 0 {
                    s.extend_from_slice(exon);
                }
            }
            if t % 9 == 0 {
                s.extend(rand(300));
            }
            txs.push(s);
        }
        // A long transcript of its own, an empty one and a short one.
        txs.push(rand(200_000));
        txs.push(Vec::new());
        txs.push(rand(30));
        txs
    }

    fn seeds_of(seeder: &BaseSeeder, codes: &[u8]) -> Vec<(u32, u32)> {
        let mut sd = Vec::new();
        seeder.seed_starts_codes_into(codes, &mut sd);
        sd.sort_unstable_by_key(|x| x.1);
        sd.iter().map(|&(l, s)| (mix(l), s)).collect()
    }

    #[test]
    fn store_and_table_hold_every_seed_and_run() {
        let seeder = BaseSeeder::new("lcp:K=2", false).unwrap();
        let txs = transcripts();
        let codes = |t: u32, out: &mut Vec<u8>| out.extend_from_slice(&txs[t as usize]);
        let (store, global) = build(&seeder, txs.len(), &codes, Some(3), 4);
        let global = global.unwrap();
        // The same tables whatever the number of threads.
        let (one, one_global) = build(&seeder, txs.len(), &codes, Some(3), 1);
        let one_global = one_global.unwrap();
        assert_eq!(*store.tx_off, *one.tx_off);
        assert_eq!(*store.occ, *one.occ);
        assert_eq!(*store.occ_len, *one.occ_len);
        assert_eq!(*store.hash, *one.hash);
        assert_eq!(*store.rel, *one.rel);
        assert_eq!(*global.keys, *one_global.keys);
        assert_eq!(*global.posts, *one_global.posts);

        let mut runs: std::collections::BTreeMap<u64, Vec<(u32, u32)>> = Default::default();
        let (mut got, mut total) = (Vec::new(), 0);
        let mut x = 0x9e37_79b9_7f4a_7c15u64;
        for (t, tx) in txs.iter().enumerate() {
            let want = seeds_of(&seeder, tx);
            store.seeds_into(t as u32, &mut got);
            assert_eq!(got, want);
            total += want.len();
            for w in want.windows(3) {
                runs.entry(tuple(w.iter().map(|s| s.0)))
                    .or_default()
                    .push((t as u32, w[0].1));
            }
            // Windows: those of a few ranges, and the edges.
            let tlen = tx.len() as u32;
            let mut ranges = vec![(0, tlen), (0, 0), (tlen, tlen), (tlen / 2, tlen / 2)];
            for _ in 0..6 {
                x = (x ^ (x >> 29)).wrapping_mul(0xff51_afd7_ed55_8ccd);
                let (a, b) = ((x >> 40) as u32 % (tlen + 1), (x >> 8) as u32 % (tlen + 1));
                ranges.push((a.min(b), a.max(b)));
                if let Some(s) = want.get((x % 97) as usize) {
                    ranges.push((s.1, s.1 + 1));
                    ranges.push((s.1 + 1, tlen.max(s.1 + 1)));
                }
            }
            for (lo, hi) in ranges {
                let w = store.window(t as u32, tlen, lo, hi);
                let mut seen = Vec::new();
                w.for_each(|l, s| seen.push((l, s)));
                let want: Vec<_> = want.iter().copied().filter(|s| lo <= s.1 && s.1 < hi).collect();
                assert_eq!(seen, want, "transcript {t} window {lo}..{hi}");
                assert_eq!(w.first_start(), want.first().map(|s| s.1));
                assert_eq!(w.is_empty(), want.is_empty());
            }
        }
        // Isoforms share their chunks.
        assert!(
            store.hash.len() * 4 < total,
            "{} of {total} seeds kept",
            store.hash.len()
        );

        // Every key's postings, and none for other keys.
        let keys: Vec<u64> = runs.keys().copied().collect();
        assert_eq!(*global.keys, keys[..]);
        let mut probes = keys.clone();
        probes.extend(keys.iter().map(|k| k ^ 1).filter(|k| !runs.contains_key(k)));
        for batch in probes.chunks(LOOKUPS) {
            global.get_batch(batch, |i, posts| {
                let mut got = Vec::new();
                posts.for_each(|t, s| got.push((t, s)));
                let mut want = runs.get(&batch[i]).cloned().unwrap_or_default();
                want.sort_unstable();
                assert_eq!(posts.len(), want.len());
                assert_eq!(got, want);
            });
        }
    }

    /// Equal seeds have equal window hashes wherever the chunks agree, and different seeds
    /// different ones.
    #[test]
    fn window_keys_tell_windows_apart() {
        let seeder = BaseSeeder::new("lcp:K=2", false).unwrap();
        let txs = transcripts();
        let codes = |t: u32, out: &mut Vec<u8>| out.extend_from_slice(&txs[t as usize]);
        let (store, _) = build(&seeder, txs.len(), &codes, None, 2);
        let mut by_key: std::collections::HashMap<u64, Vec<(u32, u32)>> = Default::default();
        let (mut windows, mut shared) = (0, 0);
        for (t, tx) in txs.iter().enumerate().take(700) {
            let tlen = tx.len() as u32;
            for lo in (0..tlen).step_by(211) {
                let w = store.window(t as u32, tlen, lo, lo + 400);
                let Some(s0) = w.first_start() else {
                    continue;
                };
                let mut seeds = Vec::new();
                w.for_each(|l, s| seeds.push((l, s - s0)));
                windows += 1;
                match by_key.entry(w.key(s0)) {
                    std::collections::hash_map::Entry::Occupied(e) => {
                        assert_eq!(*e.get(), seeds);
                        shared += 1;
                    }
                    std::collections::hash_map::Entry::Vacant(e) => {
                        e.insert(seeds);
                    }
                }
            }
        }
        assert!(shared * 2 > windows, "{shared} of {windows} windows share a key");
    }

    #[test]
    fn varints_round_trip() {
        let xs = [
            0u32,
            1,
            127,
            128,
            16383,
            16384,
            1 << 21,
            (1 << 28) - 1,
            1 << 28,
            u32::MAX,
        ];
        let mut b = Vec::new();
        for &x in &xs {
            put(&mut b, x);
        }
        let mut p = 0;
        for &x in &xs {
            let (y, q) = varint(&b, p);
            assert_eq!(y, x);
            p = q;
        }
        assert_eq!(p, b.len());
        assert_eq!(varint(&b, b.len()).1, b.len() + 1);
    }
}
