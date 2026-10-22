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

/// Bases `[lo, hi)` of 2-bit packed `words`, appended to `out` as codes a word at a time.
fn decode_into(words: &[u64], lo: u32, hi: u32, out: &mut Vec<u8>) {
    let start = out.len();
    out.resize(start + (hi - lo) as usize, 0);
    let dst = &mut out[start..];
    let (mut i, mut d) = (lo, 0usize);
    while i < hi {
        let w = words[(i / 32) as usize] >> (2 * (i % 32));
        let take = ((32 - i % 32).min(hi - i)) as usize;
        for (k, x) in dst[d..d + take].iter_mut().enumerate() {
            *x = ((w >> (2 * k)) & 3) as u8;
        }
        i += take as u32;
        d += take;
    }
}

/// 32 bases of 2-bit packed `words` starting at base `i` (bases past the end of `words` read
/// as 0).
#[inline]
fn word_at(words: &[u64], i: u32) -> u64 {
    let (q, s) = ((i / 32) as usize, 2 * (i % 32));
    let lo = words[q] >> s;
    if s == 0 {
        lo
    } else {
        lo | words.get(q + 1).map_or(0, |w| w << (64 - s))
    }
}

/// Whether bases `[alo, alo + len)` of packed `a` equal bases `[blo, blo + len)` of packed `b`,
/// compared 32 at a time without decoding.
fn packed_eq(a: &[u64], alo: u32, b: &[u64], blo: u32, len: u32) -> bool {
    let mut k = 0;
    while k + 32 <= len {
        if word_at(a, alo + k) != word_at(b, blo + k) {
            return false;
        }
        k += 32;
    }
    let r = len - k;
    r == 0 || (word_at(a, alo + k) ^ word_at(b, blo + k)) & ((1u64 << (2 * r)) - 1) == 0
}

/// A borrowed interval in an existing packed transcript, in parsed or raw space.
/// Constructed by [`Index::tx_range`]; raw intervals cover complete HPC runs.
#[derive(Clone, Copy, Debug)]
pub struct PackedRange<'a> {
    words: &'a [u64],
    lo: u32,
    hi: u32,
}

impl PackedRange<'_> {
    pub fn len(&self) -> usize {
        (self.hi - self.lo) as usize
    }

    pub fn is_empty(&self) -> bool {
        self.lo == self.hi
    }

    /// Exact sequence equality, including unaligned starts and partial final words.
    pub fn sequence_eq(&self, other: &Self) -> bool {
        self.len() == other.len() && packed_eq(self.words, self.lo, other.words, other.lo, self.hi - self.lo)
    }

    pub fn decode(&self) -> Vec<u8> {
        let mut out = Vec::new();
        decode_into(self.words, self.lo, self.hi, &mut out);
        out
    }

    /// Count non-`code` bases in a suffix without allocating a decoded sequence.
    pub fn count_other(&self, from: usize, code: u8) -> usize {
        assert!(from <= self.len());
        (self.lo + from as u32..self.hi)
            .filter(|&p| (self.words[(p / 32) as usize] >> (2 * (p % 32))) & 3 != u64::from(code))
            .count()
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
    /// What built the index, as JSON (written as section `manifest`; the caller fills it in).
    pub manifest: Option<String>,
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

/// Concatenate `parts` with offsets: part `i` is `data[off[i]..off[i + 1]]`.
pub(crate) fn concat<'a, T: Copy + 'a>(parts: impl Iterator<Item = &'a [T]> + Clone) -> (Vec<u64>, Vec<T>) {
    let (n, total) = parts.clone().fold((0, 0), |(n, t), p| (n + 1, t + p.len()));
    let mut off = Vec::with_capacity(n + 1);
    let mut data = Vec::with_capacity(total);
    off.push(0);
    for p in parts {
        data.extend_from_slice(p);
        off.push(data.len() as u64);
    }
    (off, data)
}

impl IndexData {
    /// The `meta` section's text.
    pub fn meta_text(&self) -> String {
        self.meta.iter().map(|(k, v)| format!("{k}={v}\n")).collect()
    }

    fn writer(&self) -> Writer<'_> {
        let mut w = Writer::default();
        w.add_owned("meta", self.meta_text().into_bytes());
        w.add("labels", &self.labels);
        w.add_owned("label_dir", label_dir(&self.labels));
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
        let (o, b) = concat(self.tx_names.iter().map(|s| s.as_bytes()));
        w.add_owned("tx_name_off", o);
        w.add_owned("tx_names", b);
        let (o, b) = concat(self.gene_names.iter().map(|s| s.as_bytes()));
        w.add_owned("gene_name_off", o);
        w.add_owned("gene_names", b);
        if let Some(m) = &self.manifest {
            w.add("manifest", m.as_bytes());
        }
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

/// A cache hint for the line holding `p` (no effect off x86-64 and aarch64; never faults).
#[inline(always)]
pub fn prefetch<T>(p: *const T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: prefetching is a hint: it reads nothing and cannot fault, whatever the address.
    #[allow(unused_unsafe)]
    unsafe {
        core::arch::x86_64::_mm_prefetch(p as *const i8, core::arch::x86_64::_MM_HINT_T0);
    }
    #[cfg(target_arch = "aarch64")]
    // SAFETY: as on x86-64: `prfm` is a hint that reads nothing and cannot fault.
    unsafe {
        core::arch::asm!("prfm pldl1keep, [{0}]", in(reg) p, options(nostack, readonly, preserves_flags));
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
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

    /// An index built in memory, laid out once (`from_bytes(&data.to_bytes())` copies twice).
    pub fn from_data(data: &IndexData) -> Result<Self> {
        let (words, len) = data.writer().to_words();
        Self::from_store(Store::from_words(words, len)?)
    }

    fn from_store(store: Store) -> Result<Self> {
        let text = std::str::from_utf8(store.get::<u8>("meta")?)?;
        let meta = text
            .lines()
            .filter_map(|l| l.split_once('='))
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect();
        let s = Sections::load(&store)?;
        let raw = if store.has("tx_raw_off") {
            Some(RawSections {
                off: Raw::of(store.get::<u64>("tx_raw_off")?),
                words: Raw::of(store.get::<u64>("tx_raw")?),
                ck_off: Raw::of(store.get::<u64>("tx_raw_ck_off")?),
                ck: Raw::of(store.get::<u32>("tx_raw_ck")?),
            })
        } else {
            None
        };
        let bits = label_dir_bits(s.labels.len);
        anyhow::ensure!(
            s.label_dir.len == (1usize << bits) + 1,
            "label directory has the wrong size"
        );
        Self::check_sizes(&s, raw.as_ref())?;
        Ok(Index {
            store,
            s,
            raw,
            meta,
            dir_shift: 64 - bits,
        })
    }

    /// Check that the sections' sizes agree with each other, so that a truncated or mismatched
    /// index fails to load instead of panicking while mapping. Values are not checked.
    fn check_sizes(s: &Sections, raw: Option<&RawSections>) -> Result<()> {
        // The last entry of an offsets section (None if it is empty).
        // SAFETY: borrowed from the store, which outlives this call.
        let end = |r: Raw<u64>| unsafe { r.slice() }.last().map(|&x| x as usize);
        let (nn, nu, nt) = (s.labels.len, s.unitig_span.len, s.tx_gene.len);
        let mut checks = vec![
            ("node_unitig", s.node_unitig.len == nn),
            ("node_offset", s.node_offset.len == nn),
            ("node_len", s.node_len.len == nn),
            ("node_flags", s.node_flags.len == nn),
            ("unitig_off", s.unitig_off.len == nu + 1),
            ("unitig_nodes", Some(s.unitig_nodes.len) == end(s.unitig_off)),
            ("unitig_color", s.unitig_color.len == nu),
            ("unitig_family", s.unitig_family.len == nu),
            ("occ_off", s.occ_off.len == nu + 1),
            ("occ_tx", Some(s.occ_tx.len) == end(s.occ_off)),
            ("occ_pos", Some(s.occ_pos.len) == end(s.occ_off)),
            ("color_tx", Some(s.color_tx.len) == end(s.color_off)),
            ("tx_family", s.tx_family.len == nt),
            ("tx_len", s.tx_len.len == nt),
            ("tx_plen", s.tx_plen.len == nt),
            ("tx_seq_off", s.tx_seq_off.len == nt + 1),
            ("tx_seq", Some(s.tx_seq.len) == end(s.tx_seq_off)),
            ("tx_name_off", s.tx_name_off.len == nt + 1),
            ("tx_names", Some(s.tx_names.len) == end(s.tx_name_off)),
            ("gene_names", Some(s.gene_names.len) == end(s.gene_name_off)),
        ];
        if let Some(r) = raw {
            checks.extend([
                ("tx_raw_off", r.off.len == nt + 1),
                ("tx_raw", Some(r.words.len) == end(r.off)),
                ("tx_raw_ck_off", r.ck_off.len == nt + 1),
                ("tx_raw_ck", Some(r.ck.len) == end(r.ck_off)),
            ]);
        }
        for (name, ok) in checks {
            anyhow::ensure!(ok, "index section {name} has the wrong size; rebuild the index");
        }
        Ok(())
    }

    pub fn meta(&self, key: &str) -> Option<&str> {
        self.meta.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
    }

    /// The `meta` section's text (`key=value` lines).
    pub fn meta_text(&self) -> &str {
        std::str::from_utf8(self.store.get::<u8>("meta").unwrap_or(&[])).unwrap_or("")
    }

    /// The manifest (JSON) recorded when the index was built, if any.
    pub fn manifest(&self) -> Option<&str> {
        let b = self.store.get::<u8>("manifest").ok()?;
        std::str::from_utf8(b).ok()
    }

    /// The index holds raw (uncompressed) transcript sequence.
    pub fn has_raw(&self) -> bool {
        self.raw.is_some()
    }

    /// Check that this build parses sequence as the index was built: an index built with
    /// another LCP parser version or hash seed loads, but would give wrong mappings.
    pub fn check_compatible(&self) -> Result<()> {
        let spec = self.meta("spec").unwrap_or("");
        if !spec.starts_with("lcp") {
            return Ok(());
        }
        if let Some(v) = self.meta("lcp_parser_version") {
            let ours = lcp_core::LCP_PARSER_VERSION.to_string();
            anyhow::ensure!(
                v == ours,
                "index built with LCP parser version {v}, this build parses with version {ours}; rebuild the index"
            );
        }
        if let Some(seed) = self.meta("lcp_seed") {
            let ours = format!("{:#x}", lcp_core::LcpConfig::default().seed);
            anyhow::ensure!(
                seed == ours,
                "index built with LCP seed {seed}, this build uses {ours}; rebuild the index"
            );
        }
        Ok(())
    }

    /// Every section's name, size in bytes and xxh3-64 checksum.
    pub fn section_checksums(&self) -> Vec<(String, usize, u64)> {
        self.store
            .sizes()
            .into_iter()
            .map(|(name, len)| {
                let h = xxhash_rust::xxh3::xxh3_64(self.store.get::<u8>(&name).unwrap_or(&[]));
                (name, len, h)
            })
            .collect()
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
        self.tx_codes_into(t, lo, hi, &mut v);
        v
    }

    /// Borrow a packed reference interval without copying its bases. Coordinates are
    /// parsed coordinates even when `raw` is true; conversion uses the run checkpoints.
    /// Returns None for a raw interval if that transcript has no stored raw sequence.
    pub fn tx_range(&self, t: u32, lo: u32, hi: u32, raw: bool) -> Option<PackedRange<'_>> {
        if raw {
            Some(PackedRange {
                words: self.raw_words(t)?.0,
                lo: self.tx_raw_pos(t, lo)?,
                hi: self.tx_raw_pos(t, hi)?,
            })
        } else {
            Some(PackedRange {
                words: &self.tx_seq()[self.tx_seq_off()[t as usize] as usize..],
                lo,
                hi,
            })
        }
    }

    /// Append the parsed transcript `t` over `[lo, hi)` to `out`, as 2-bit codes.
    pub fn tx_codes_into(&self, t: u32, lo: u32, hi: u32, out: &mut Vec<u8>) {
        let words = &self.tx_seq()[self.tx_seq_off()[t as usize] as usize..];
        decode_into(words, lo, hi, out);
    }

    /// Whether the parsed transcripts `t` from `lo` and `u` from `ulo` agree on `len` bases,
    /// compared on the packed sequence (no decoding).
    pub fn tx_codes_eq(&self, t: u32, lo: u32, u: u32, ulo: u32, len: u32) -> bool {
        let (seq, off) = (self.tx_seq(), self.tx_seq_off());
        packed_eq(
            &seq[off[t as usize] as usize..],
            lo,
            &seq[off[u as usize] as usize..],
            ulo,
            len,
        )
    }

    /// Raw 2-bit codes of transcript `t` over the raw span of its parsed bases `[lo, hi)`
    /// (whole homopolymer runs), or None if the index has no raw sequence for `t`.
    pub fn tx_raw_codes(&self, t: u32, lo: u32, hi: u32) -> Option<Vec<u8>> {
        let (a, b) = (self.tx_raw_pos(t, lo)?, self.tx_raw_pos(t, hi)?);
        let (words, _) = self.raw_words(t)?;
        let mut v = Vec::new();
        decode_into(words, a, b, &mut v);
        Some(v)
    }

    /// Whether the raw bases under parsed `[lo, hi)` of transcript `t` equal those under parsed
    /// `[ulo, uhi)` of `u` (homopolymer run lengths included), compared on the packed sequence.
    /// None if the index lacks either raw sequence.
    pub fn tx_raw_eq(&self, t: u32, lo: u32, hi: u32, u: u32, ulo: u32, uhi: u32) -> Option<bool> {
        let (a, b) = (self.tx_raw_pos(t, lo)?, self.tx_raw_pos(t, hi)?);
        let (ua, ub) = (self.tx_raw_pos(u, ulo)?, self.tx_raw_pos(u, uhi)?);
        if b - a != ub - ua {
            return Some(false);
        }
        let ((tw, _), (uw, _)) = (self.raw_words(t)?, self.raw_words(u)?);
        Some(packed_eq(tw, a, uw, ua, b - a))
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

    /// Ask the caches for the offsets [`Self::occurrences`] reads first.
    #[inline]
    pub fn prefetch_occ_off(&self, u: u32) {
        prefetch(self.occ_off().as_ptr().wrapping_add(u as usize));
    }

    /// Ask the caches for the occurrence lists of `u` (reads the offsets, so call it after
    /// [`Self::prefetch_occ_off`] has had time to land).
    #[inline]
    pub fn prefetch_occurrences(&self, u: u32) {
        let o = self.occ_off()[u as usize] as usize;
        prefetch(self.occ_tx().as_ptr().wrapping_add(o));
        prefetch(self.occ_pos().as_ptr().wrapping_add(o));
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
        #[test]
        fn borrowed_ranges_match_decoded_slices(
            bases in prop::collection::vec(0u8..4, 0..300),
            other in prop::collection::vec(0u8..4, 0..300),
            a in 0usize..300, b in 0usize..300, len in 0usize..300,
        ) {
            let (a, b) = (a.min(bases.len()), b.min(other.len()));
            let (an, bn) = (len.min(bases.len() - a), len.min(other.len() - b));
            let (w, v) = (pack(&bases), pack(&other));
            let x = PackedRange { words: &w, lo: a as u32, hi: (a + an) as u32 };
            let y = PackedRange { words: &v, lo: b as u32, hi: (b + bn) as u32 };
            prop_assert_eq!(x.sequence_eq(&y), bases[a..a+an] == other[b..b+bn]);
            prop_assert_eq!(x.decode(), &bases[a..a+an]);
            prop_assert_eq!(x.count_other(an / 2, 0), bases[a+an/2..a+an].iter().filter(|&&c| c != 0).count());
        }

        #[test]
        fn relocated_ranges_at_all_word_offsets(
            bases in prop::collection::vec(0u8..4, 0..100),
        ) {
            for a in 0..32 {
                for b in 0..32 {
                    let (mut x, mut y) = (vec![3; a], vec![1; b]);
                    x.extend(&bases); y.extend(&bases);
                    // Different padding must not enter the comparison's last partial word.
                    x.push(2); y.push(0);
                    let (w, v) = (pack(&x), pack(&y));
                    let x = PackedRange { words: &w, lo: a as u32, hi: (a+bases.len()) as u32 };
                    let y = PackedRange { words: &v, lo: b as u32, hi: (b+bases.len()) as u32 };
                    prop_assert!(x.sequence_eq(&y));
                }
            }
        }

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
            decode_into(&w, lo, hi, &mut got);
            let want: Vec<u8> = (lo..hi).map(|i| base_at(&w, i)).collect();
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

        /// Packed range comparison equals comparing the bases, for equal and nearly equal
        /// ranges at any alignment.
        #[test]
        fn packed_eq_matches_scalar(
            bases in prop::collection::vec(0u8..4, 1..300),
            a in 0usize..300, b in 0usize..300, len in 0usize..300, flip in prop::option::of(0usize..300),
        ) {
            let n = bases.len();
            let (a, b) = (a % n, b % n);
            let len = len.min(n - a).min(n - b);
            // A second sequence equal to the first except, optionally, one base.
            let mut other = bases.clone();
            if let Some(f) = flip {
                other[f % n] ^= 1;
            }
            let (w, v) = (pack(&bases), pack(&other));
            let want = bases[a..a + len] == other[b..b + len];
            prop_assert_eq!(packed_eq(&w, a as u32, &v, b as u32, len as u32), want);
            // Same range of the same sequence.
            prop_assert!(packed_eq(&w, a as u32, &w, a as u32, len as u32));
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
