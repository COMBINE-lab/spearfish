//! The spearfish index: a colored, compacted graph over base-level seeds.
//!
//! Nodes are distinct seed labels; unitigs are maximal non-branching chains of seeds with
//! fixed relative offsets; each transcript is a walk over whole unitigs, stored as
//! occurrences `(transcript, position of the unitig's first seed)`. A seed in unitig `u` at
//! offset `δ` therefore sits at `pos + δ` in every transcript `t` with an occurrence
//! `(t, pos)` of `u`, which is all the chainer needs. Coordinates are in the parsed
//! (homopolymer-compressed, if HPC is on) transcript. See `build` for construction.

mod build;
pub mod store;

use std::path::Path;

use anyhow::Result;

pub use build::{BuildParams, BuildReport, Collapse, InferGenes, TxInput, build};
use store::{Store, Writer};

/// `unitig_family` value for masked or unused unitigs.
pub const NO_FAMILY: u32 = u32::MAX;
/// `node_flags` bit: the node's unitig is masked (color touches too many genes).
pub const FLAG_MASKED: u8 = 1;
/// Parsed bases between raw-position checkpoints (`tx_raw_ck`).
pub const RAW_CK: u32 = 64;

/// Base `i` of 2-bit packed `words`.
#[inline]
#[cfg_attr(not(test), allow(dead_code))]
fn base_at(words: &[u64], i: u32) -> u8 {
    ((words[(i / 32) as usize] >> (2 * (i % 32))) & 3) as u8
}

/// Bases `[lo, hi)` of 2-bit packed `words`, through `map`, appended to `out` a word at a time.
fn decode_into(words: &[u64], lo: u32, hi: u32, map: [u8; 4], out: &mut Vec<u8>) {
    let start = out.len();
    out.resize(start + (hi - lo) as usize, 0);
    let dst = &mut out[start..];
    let (mut i, mut d) = (lo, 0usize);
    while i < hi {
        let w = words[(i / 32) as usize] >> (2 * (i % 32));
        let take = ((32 - i % 32).min(hi - i)) as usize;
        for (k, x) in dst[d..d + take].iter_mut().enumerate() {
            *x = map[((w >> (2 * k)) & 3) as usize];
        }
        i += take as u32;
        d += take;
    }
}

/// The `q`-th position after `r` (q >= 1) where the base differs from the one before it, in
/// 2-bit packed `words` of `rlen` bases, or `rlen` if there are fewer: where `q` runs after the
/// run starting at `r` end. Changes are found a word at a time (XOR with the shifted word).
fn nth_run_start(words: &[u64], rlen: u32, r: u32, mut q: u32) -> u32 {
    let mut j0 = r + 1;
    while j0 < rlen {
        let wi = (j0 / 32) as usize;
        let w = words[wi];
        // The base before the word's first (no change at position 0: it is never asked for).
        let prev = if wi == 0 { w & 3 } else { words[wi - 1] >> 62 };
        let x = w ^ ((w << 2) | prev);
        // One bit per base (the low bit of its field): set where the base changes.
        let mut m = (x | (x >> 1)) & 0x5555_5555_5555_5555;
        m &= !0u64 << (2 * (j0 % 32));
        let end = rlen - wi as u32 * 32;
        if end < 32 {
            m &= (1u64 << (2 * end)) - 1;
        }
        let c = m.count_ones();
        if c >= q {
            for _ in 1..q {
                m &= m - 1;
            }
            return wi as u32 * 32 + m.trailing_zeros() / 2;
        }
        q -= c;
        j0 = (wi as u32 + 1) * 32;
    }
    rlen
}

/// Owned index contents, as produced by `build`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexData {
    pub meta: Vec<(String, String)>,
    pub labels: Vec<u64>,
    pub node_unitig: Vec<u32>,
    pub node_offset: Vec<u32>,
    pub node_len: Vec<u32>,
    pub node_flags: Vec<u8>,
    pub unitig_off: Vec<u64>,
    pub unitig_nodes: Vec<u32>,
    pub unitig_span: Vec<u32>,
    pub unitig_color: Vec<u32>,
    pub unitig_family: Vec<u32>,
    pub occ_off: Vec<u64>,
    pub occ_tx: Vec<u32>,
    pub occ_pos: Vec<u32>,
    pub color_off: Vec<u64>,
    pub color_tx: Vec<u32>,
    pub tx_gene: Vec<u32>,
    pub tx_family: Vec<u32>,
    pub tx_len: Vec<u32>,
    pub tx_plen: Vec<u32>,
    /// Parsed transcript sequences, 2-bit packed (32 bases per word, base i at bits
    /// 2(i mod 32)); each transcript starts on a word boundary: (word offsets, words).
    pub tx_seq: (Vec<u64>, Vec<u64>),
    /// HPC indexes: raw transcript sequences, packed as `tx_seq` (empty for a transcript whose
    /// non-ACGT bases make its packed runs differ from its parsed runs).
    pub tx_raw: (Vec<u64>, Vec<u64>),
    /// HPC indexes: raw position of every `RAW_CK`-th parsed base, per transcript
    /// (offsets, positions; empty where `tx_raw` is).
    pub tx_raw_ck: (Vec<u64>, Vec<u32>),
    pub tx_names: Vec<String>,
    pub gene_names: Vec<String>,
}

/// 2-bit pack an ASCII sequence (non-ACGT as A), 32 bases per word, base i at bits
/// 2(i mod 32).
pub(crate) fn pack2(seq: &[u8]) -> Vec<u64> {
    let mut words = vec![0u64; seq.len().div_ceil(32)];
    for (i, &b) in seq.iter().enumerate() {
        let c = lcp_core::encode::CODE[b as usize];
        let c = if c > 3 { 0 } else { c as u64 };
        words[i / 32] |= c << (2 * (i % 32));
    }
    words
}

fn pool(strings: &[String]) -> (Vec<u64>, Vec<u8>) {
    let mut off = vec![0u64];
    let mut bytes = Vec::new();
    for s in strings {
        bytes.extend_from_slice(s.as_bytes());
        off.push(bytes.len() as u64);
    }
    (off, bytes)
}

impl IndexData {
    fn writer(&self) -> Writer {
        let mut w = Writer::default();
        let meta: String = self.meta.iter().map(|(k, v)| format!("{k}={v}\n")).collect();
        w.add("meta", meta.as_bytes());
        w.add("labels", &self.labels);
        w.add("label_dir", &label_dir(&self.labels));
        w.add("node_unitig", &self.node_unitig);
        w.add("node_offset", &self.node_offset);
        w.add("node_len", &self.node_len);
        w.add("node_flags", &self.node_flags);
        w.add("unitig_off", &self.unitig_off);
        w.add("unitig_nodes", &self.unitig_nodes);
        w.add("unitig_span", &self.unitig_span);
        w.add("unitig_color", &self.unitig_color);
        w.add("unitig_family", &self.unitig_family);
        w.add("occ_off", &self.occ_off);
        w.add("occ_tx", &self.occ_tx);
        w.add("occ_pos", &self.occ_pos);
        w.add("color_off", &self.color_off);
        w.add("color_tx", &self.color_tx);
        w.add("tx_gene", &self.tx_gene);
        w.add("tx_family", &self.tx_family);
        w.add("tx_len", &self.tx_len);
        w.add("tx_plen", &self.tx_plen);
        w.add("tx_seq_off", &self.tx_seq.0);
        w.add("tx_seq", &self.tx_seq.1);
        if !self.tx_raw.0.is_empty() {
            w.add("tx_raw_off", &self.tx_raw.0);
            w.add("tx_raw", &self.tx_raw.1);
            w.add("tx_raw_ck_off", &self.tx_raw_ck.0);
            w.add("tx_raw_ck", &self.tx_raw_ck.1);
        }
        let (o, b) = pool(&self.tx_names);
        w.add("tx_name_off", &o);
        w.add("tx_names", &b);
        let (o, b) = pool(&self.gene_names);
        w.add("gene_name_off", &o);
        w.add("gene_names", &b);
        w
    }

    /// Write to `path`; returns bytes written.
    pub fn write(&self, path: &Path) -> Result<u64> {
        self.writer().write(path)
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.writer().to_bytes()
    }
}

/// A cache hint for the line holding `p` (no effect off x86-64; never faults).
#[inline(always)]
fn prefetch<T>(p: *const T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: prefetching is a hint: it reads nothing and cannot fault, whatever the address.
    #[allow(unused_unsafe)]
    unsafe {
        core::arch::x86_64::_mm_prefetch(p as *const i8, core::arch::x86_64::_MM_HINT_T0);
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = p;
}

/// Bits of the label directory for `n` nodes: about one node per bucket.
pub fn label_dir_bits(n: usize) -> u32 {
    (usize::BITS - n.max(2).leading_zeros() - 1).clamp(1, 30)
}

/// For each of the 2^bits buckets of the label's top bits, the first node whose label falls
/// in it (plus a final entry, the node count). Labels are uniform hashes, so a lookup reads
/// one directory entry and scans about one label.
pub fn label_dir(labels: &[u64]) -> Vec<u32> {
    let bits = label_dir_bits(labels.len());
    let mut dir = vec![0u32; (1usize << bits) + 1];
    let mut i = 0usize;
    for (b, d) in dir.iter_mut().enumerate() {
        while i < labels.len() && ((labels[i] >> (64 - bits)) as usize) < b {
            i += 1;
        }
        *d = i as u32;
    }
    dir
}

/// A section borrowed from the index's own backing memory (mmap or heap buffer), which never
/// moves or changes while the `Index` lives.
struct Raw<T> {
    ptr: *const T,
    len: usize,
}

impl<T> Clone for Raw<T> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<T> Copy for Raw<T> {}

// SAFETY: the pointed-to data is immutable and owned by the Index holding the Raw.
unsafe impl<T: Sync> Send for Raw<T> {}
unsafe impl<T: Sync> Sync for Raw<T> {}

impl<T> Raw<T> {
    fn of(s: &[T]) -> Self {
        Raw {
            ptr: s.as_ptr(),
            len: s.len(),
        }
    }

    /// SAFETY: the caller ties the returned lifetime to the owner of the backing memory.
    #[inline(always)]
    unsafe fn slice<'a>(self) -> &'a [T] {
        unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
    }
}

macro_rules! sections {
    ($($name:ident : $t:ty),* $(,)?) => {
        struct Sections {
            $($name: Raw<$t>,)*
        }

        impl Sections {
            fn load(store: &Store) -> Result<Self> {
                Ok(Sections { $($name: Raw::of(store.get::<$t>(stringify!($name))?),)* })
            }
        }

        impl Index {
            $(
                #[inline(always)]
                pub fn $name(&self) -> &[$t] {
                    // SAFETY: borrowed from self.store, which outlives the returned slice.
                    unsafe { self.s.$name.slice() }
                }
            )*
        }
    };
}

sections!(
    labels: u64,
    label_dir: u32,
    node_unitig: u32,
    node_offset: u32,
    node_len: u32,
    node_flags: u8,
    unitig_off: u64,
    unitig_nodes: u32,
    unitig_span: u32,
    unitig_color: u32,
    unitig_family: u32,
    occ_off: u64,
    occ_tx: u32,
    occ_pos: u32,
    color_off: u64,
    color_tx: u32,
    tx_gene: u32,
    tx_family: u32,
    tx_len: u32,
    tx_plen: u32,
    tx_seq_off: u64,
    tx_seq: u64,
    tx_name_off: u64,
    tx_names: u8,
    gene_name_off: u64,
    gene_names: u8,
);

/// Optional sections: raw transcript sequences (HPC indexes).
struct RawSections {
    off: Raw<u64>,
    words: Raw<u64>,
    ck_off: Raw<u64>,
    ck: Raw<u32>,
}

/// A loaded (memory-mapped) index.
pub struct Index {
    store: Store,
    s: Sections,
    raw: Option<RawSections>,
    meta: Vec<(String, String)>,
    dir_shift: u32,
}

impl Index {
    pub fn open(path: &Path) -> Result<Self> {
        Self::from_store(Store::open(path)?)
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        Self::from_store(Store::from_bytes(bytes)?)
    }

    fn from_store(store: Store) -> Result<Self> {
        let text = std::str::from_utf8(store.get::<u8>("meta")?)?;
        let meta = text
            .lines()
            .filter_map(|l| l.split_once('='))
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect();
        let s = Sections::load(&store)?;
        let raw = (|| -> Result<RawSections> {
            Ok(RawSections {
                off: Raw::of(store.get::<u64>("tx_raw_off")?),
                words: Raw::of(store.get::<u64>("tx_raw")?),
                ck_off: Raw::of(store.get::<u64>("tx_raw_ck_off")?),
                ck: Raw::of(store.get::<u32>("tx_raw_ck")?),
            })
        })()
        .ok();
        let bits = label_dir_bits(s.labels.len);
        anyhow::ensure!(
            s.label_dir.len == (1usize << bits) + 1,
            "label directory has the wrong size"
        );
        Ok(Index {
            store,
            s,
            raw,
            meta,
            dir_shift: 64 - bits,
        })
    }

    pub fn meta(&self, key: &str) -> Option<&str> {
        self.meta.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
    }

    pub fn n_nodes(&self) -> usize {
        self.labels().len()
    }

    pub fn n_unitigs(&self) -> usize {
        self.unitig_span().len()
    }

    pub fn n_transcripts(&self) -> usize {
        self.tx_gene().len()
    }

    /// Ask the caches for the directory entry `lookup(label)` reads first. Lookups are random
    /// accesses into large arrays; issuing these a few seeds ahead overlaps their misses.
    #[inline]
    pub fn prefetch_dir(&self, label: u64) {
        let b = (label >> self.dir_shift) as usize;
        prefetch(self.label_dir().as_ptr().wrapping_add(b));
    }

    /// Ask the caches for the labels and node flags `lookup(label)` scans (reads the directory
    /// entry, so call it after [`Self::prefetch_dir`] has had time to land).
    #[inline]
    pub fn prefetch_labels(&self, label: u64) {
        let b = (label >> self.dir_shift) as usize;
        let lo = self.label_dir()[b] as usize;
        prefetch(self.labels().as_ptr().wrapping_add(lo));
        prefetch(self.node_flags().as_ptr().wrapping_add(lo));
    }

    /// Ask the caches for node `node`'s unitig entry (read by family grouping and chaining).
    #[inline]
    pub fn prefetch_node(&self, node: u32) {
        prefetch(self.node_unitig().as_ptr().wrapping_add(node as usize));
    }

    /// Ask the caches for the family entry of `node`'s unitig (reads the unitig entry, so call
    /// it once [`Self::prefetch_node`] has had time to land).
    #[inline]
    pub fn prefetch_family(&self, node: u32) {
        let u = self.node_unitig()[node as usize] as usize;
        prefetch(self.unitig_family().as_ptr().wrapping_add(u));
    }

    /// Nodes carrying a seed label (empty if the index lacks it). With global collapsing at
    /// most one; with per-gene collapsing one per gene containing the label.
    #[inline]
    pub fn lookup(&self, label: u64) -> std::ops::Range<u32> {
        let b = (label >> self.dir_shift) as usize;
        let dir = self.label_dir();
        let (lo, hi) = (dir[b] as usize, dir[b + 1] as usize);
        let l = &self.labels()[lo..hi];
        let start = l.partition_point(|&x| x < label);
        let end = start + l[start..].partition_point(|&x| x == label);
        (lo + start) as u32..(lo + end) as u32
    }

    /// 2-bit codes (A=0, C=1, G=2, T=3) of the parsed transcript `t` over `[lo, hi)`.
    pub fn tx_codes(&self, t: u32, lo: u32, hi: u32) -> Vec<u8> {
        let mut v = Vec::new();
        self.tx_codes_into(t, lo, hi, &mut v, false);
        v
    }

    /// Append the parsed transcript `t` over `[lo, hi)` to `out`, as 2-bit codes or (if
    /// `ascii`) as ACGT.
    pub fn tx_codes_into(&self, t: u32, lo: u32, hi: u32, out: &mut Vec<u8>, ascii: bool) {
        let words = &self.tx_seq()[self.tx_seq_off()[t as usize] as usize..];
        let map: [u8; 4] = if ascii { *b"ACGT" } else { [0, 1, 2, 3] };
        decode_into(words, lo, hi, map, out);
    }

    /// Raw 2-bit codes of transcript `t` over the raw span of its parsed bases `[lo, hi)`
    /// (whole homopolymer runs), or None if the index has no raw sequence for `t`.
    pub fn tx_raw_codes(&self, t: u32, lo: u32, hi: u32) -> Option<Vec<u8>> {
        let (a, b) = (self.tx_raw_pos(t, lo)?, self.tx_raw_pos(t, hi)?);
        let (words, _) = self.raw_words(t)?;
        let mut v = Vec::new();
        decode_into(words, a, b, [0, 1, 2, 3], &mut v);
        Some(v)
    }

    /// Raw position of transcript `t` where parsed base `p` starts (the raw length at the
    /// parsed length), if the index has its raw sequence.
    pub fn tx_raw_pos(&self, t: u32, p: u32) -> Option<u32> {
        let (words, ck) = self.raw_words(t)?;
        let rlen = self.tx_len()[t as usize];
        // From the checkpoint, skip whole runs.
        let r = ck[(p / RAW_CK) as usize];
        Some(match p % RAW_CK {
            0 => r,
            q => nth_run_start(words, rlen, r, q),
        })
    }

    /// Transcript `t`'s packed raw bases and run checkpoints, if present.
    fn raw_words(&self, t: u32) -> Option<(&[u64], &[u32])> {
        let r = self.raw.as_ref()?;
        // SAFETY: borrowed from self.store, which outlives these slices.
        let (off, words, ck_off, ck) =
            unsafe { (r.off.slice(), r.words.slice(), r.ck_off.slice(), r.ck.slice()) };
        let t = t as usize;
        let ck = &ck[ck_off[t] as usize..ck_off[t + 1] as usize];
        if ck.is_empty() {
            return None;
        }
        Some((&words[off[t] as usize..off[t + 1] as usize], ck))
    }

    /// Nodes of unitig `u`, in order.
    pub fn nodes_of(&self, u: u32) -> &[u32] {
        let off = self.unitig_off();
        &self.unitig_nodes()[off[u as usize] as usize..off[u as usize + 1] as usize]
    }

    /// Occurrences of unitig `u`: transcripts (sorted) and the position of its first seed.
    pub fn occurrences(&self, u: u32) -> (&[u32], &[u32]) {
        let off = self.occ_off();
        let r = off[u as usize] as usize..off[u as usize + 1] as usize;
        (&self.occ_tx()[r.clone()], &self.occ_pos()[r])
    }

    /// Transcripts of color class `c` (sorted).
    pub fn color(&self, c: u32) -> &[u32] {
        let off = self.color_off();
        &self.color_tx()[off[c as usize] as usize..off[c as usize + 1] as usize]
    }

    fn pooled<'a>(off: &[u64], bytes: &'a [u8], i: usize) -> &'a str {
        std::str::from_utf8(&bytes[off[i] as usize..off[i + 1] as usize]).unwrap()
    }

    pub fn tx_name(&self, t: u32) -> &str {
        Self::pooled(self.tx_name_off(), self.tx_names(), t as usize)
    }

    pub fn gene_name(&self, g: u32) -> &str {
        Self::pooled(self.gene_name_off(), self.gene_names(), g as usize)
    }

    /// Section sizes in bytes, by name.
    pub fn sizes(&self) -> Vec<(String, usize)> {
        self.store.sizes()
    }

    /// Every transcript's seeds, reconstructed from its unitig walk as (node, position),
    /// sorted by position. Used to check the index against fresh seeding.
    pub fn transcript_seeds(&self) -> Vec<Vec<(u32, u32)>> {
        let mut out = vec![Vec::new(); self.n_transcripts()];
        for u in 0..self.n_unitigs() as u32 {
            let (ts, ps) = self.occurrences(u);
            for (&t, &p) in ts.iter().zip(ps) {
                out[t as usize].extend(
                    self.nodes_of(u)
                        .iter()
                        .map(|&v| (v, p + self.node_offset()[v as usize])),
                );
            }
        }
        for v in &mut out {
            v.sort_unstable_by_key(|&(n, p)| (p, n));
        }
        out
    }
}

#[cfg(test)]
mod packed_tests {
    use super::*;
    use proptest::prelude::*;

    fn pack(bases: &[u8]) -> Vec<u64> {
        let mut w = vec![0u64; bases.len().div_ceil(32).max(1)];
        for (i, &b) in bases.iter().enumerate() {
            w[i / 32] |= (b as u64) << (2 * (i % 32));
        }
        w
    }

    proptest! {
        /// Word-at-a-time decoding and run skipping equal the base-at-a-time loops.
        #[test]
        fn packed_ops_match_scalar(
            bases in prop::collection::vec(prop_oneof![0u8..4, Just(1u8)], 1..300),
            a in 0usize..300, b in 0usize..300, q in 1u32..64,
        ) {
            let w = pack(&bases);
            let n = bases.len() as u32;
            let (lo, hi) = ((a as u32).min(n), (b as u32).min(n));
            let (lo, hi) = (lo.min(hi), lo.max(hi));
            let mut got = Vec::new();
            decode_into(&w, lo, hi, *b"ACGT", &mut got);
            let want: Vec<u8> = (lo..hi).map(|i| b"ACGT"[base_at(&w, i) as usize]).collect();
            prop_assert_eq!(got, want);
            // A run start: position 0 or a change.
            let starts: Vec<u32> =
                (0..n).filter(|&i| i == 0 || bases[i as usize] != bases[i as usize - 1]).collect();
            let r = starts[a % starts.len()];
            let mut want = r;
            for _ in 0..q {
                if want >= n { break; }
                let c = base_at(&w, want);
                want += 1;
                while want < n && base_at(&w, want) == c {
                    want += 1;
                }
            }
            prop_assert_eq!(nth_run_start(&w, n, r, q), want);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Directory lookup equals a plain binary search, for present labels (with duplicates, as
    /// per-gene nodes give), absent labels, and the extremes.
    #[test]
    fn label_dir_lookup_matches_binary_search() {
        let mut x = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        for n in [1usize, 2, 3, 100, 5000] {
            let mut labels: Vec<u64> = (0..n).map(|_| next()).collect();
            // Duplicates and extremes.
            let dup: Vec<u64> = labels.iter().step_by(7).copied().collect();
            labels.extend(dup);
            labels.extend([0, u64::MAX, u64::MAX]);
            labels.sort_unstable();
            let dir = label_dir(&labels);
            let shift = 64 - label_dir_bits(labels.len());
            let lookup = |label: u64| {
                let b = (label >> shift) as usize;
                let (lo, hi) = (dir[b] as usize, dir[b + 1] as usize);
                let l = &labels[lo..hi];
                let s = l.partition_point(|&v| v < label);
                lo + s..lo + s + l[s..].partition_point(|&v| v == label)
            };
            let reference = |label: u64| {
                let s = labels.partition_point(|&v| v < label);
                s..s + labels[s..].partition_point(|&v| v == label)
            };
            for &q in labels.iter().chain(&[1, 12345, u64::MAX - 1]) {
                assert_eq!(lookup(q), reference(q), "n {n} label {q:x}");
            }
            for _ in 0..1000 {
                let q = next();
                assert_eq!(lookup(q), reference(q));
            }
        }
    }
}
