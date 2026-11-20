//! The genome index: packed sequence, level-3 seed table, positional level-2 table.
//!
//! All positions are raw (uncompressed) bases on one axis that concatenates the reference
//! sequences in file order, so a position fits 32 bits for genomes under 4.29 Gbp.
//!
//! - `seq`: 2 bits per base (A=0, C=1, G=2, T=3; other symbols are stored as A and listed
//!   in `n_runs`, and never seeded).
//! - Level-3 LCP cores of the homopolymer-compressed sequence: distinct `labels` (sorted),
//!   and per label the raw start and raw length of each occurrence. Labels with more than
//!   `max_occ` occurrences are left out.
//! - Level-2 cores in genome order (`l2_pos`, `l2_id`): the table that exon rescue scans
//!   over one intron. `l2_id` is the rank of the core's label among the genome's distinct
//!   level-2 labels (there are only about 21,000 on the human genome).

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, ensure};
use index::store::{Store, Writer};
use rayon::prelude::*;
use seeds::BaseSeeder;

/// Sequence is seeded in chunks of this many bases, so that long chromosomes parse in
/// parallel and in bounded memory. Cores that would straddle a chunk edge are lost.
const CHUNK: usize = 8 << 20;

pub const L3_SPEC: &str = "lcp:K=3";
pub const L2_SPEC: &str = "lcp:K=2";

#[derive(Clone, Debug)]
pub struct BuildOptions {
    /// Level-3 labels with more occurrences than this are not indexed.
    pub max_occ: u32,
    /// Seed the homopolymer-compressed sequence (longer, sparser cores that survive
    /// homopolymer-length errors) or the sequence as it is (shorter, denser cores; the default:
    /// on simulated and real reads it finds more loci and exons in less time, for a larger
    /// index).
    pub hpc: bool,
}

impl Default for BuildOptions {
    fn default() -> Self {
        BuildOptions {
            max_occ: 1000,
            hpc: false,
        }
    }
}

/// An index held in memory (built, or read from a container).
pub struct GenomeIndex {
    store: Store,
    names: Vec<String>,
    n: Raw,
    /// Known junctions, when the user gave some for this run (not stored in the index).
    pub junctions: Option<crate::junctions::Junctions>,
}

/// Pointers into the container's sections (the `Store` owns the bytes and outlives them).
struct Raw {
    seq: *const [u64],
    offsets: *const [u32],
    labels: *const [u64],
    dir: *const [u32],
    occ_off: *const [u32],
    pos: *const [u32],
    len: *const [u8],
    l2_labels: *const [u64],
    l2_pos: *const [u32],
    l2_id: *const [u16],
    p_labels: *const [u64],
    p_dir: *const [u32],
    p_occ_off: *const [u32],
    p_pos: *const [u32],
    p_len: *const [u8],
    p_bits: u32,
    dir_bits: u32,
    hpc: bool,
    ranked: bool,
}

// SAFETY: the pointers refer to immutable data owned by `store`, which is never moved out
// of the heap/mmap it lives in while the index exists.
unsafe impl Send for GenomeIndex {}
unsafe impl Sync for GenomeIndex {}

#[inline]
pub fn code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => 4,
    }
}

fn read_fasta(path: &Path) -> Result<(Vec<u8>, Vec<String>, Vec<u32>)> {
    let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut seq = Vec::new();
    let mut names = Vec::new();
    let mut offsets = Vec::new();
    for line in BufReader::with_capacity(1 << 24, f).lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            names.push(h.split_whitespace().next().unwrap_or("").to_string());
            offsets.push(seq.len() as u32);
        } else {
            seq.extend(line.bytes().map(|b| b.to_ascii_uppercase()));
        }
        ensure!(
            seq.len() < u32::MAX as usize,
            "reference longer than 4.29 Gbp is not supported yet"
        );
    }
    ensure!(!names.is_empty(), "no sequences in {}", path.display());
    offsets.push(seq.len() as u32);
    Ok((seq, names, offsets))
}

/// Maximal ACGT stretches, cut into chunks of at most `CHUNK` bases.
fn segments(seq: &[u8], offsets: &[u32]) -> Vec<(usize, usize)> {
    let mut v = Vec::new();
    for w in offsets.windows(2) {
        let (mut i, end) = (w[0] as usize, w[1] as usize);
        while i < end {
            if code(seq[i]) > 3 {
                i += 1;
                continue;
            }
            let s = i;
            while i < end && code(seq[i]) <= 3 {
                i += 1;
            }
            let mut a = s;
            while a < i {
                let b = (a + CHUNK).min(i);
                v.push((a, b));
                a = b;
            }
        }
    }
    v
}

/// (label, raw start, raw length) of the top-level cores of every segment.
fn seed_all(seq: &[u8], segs: &[(usize, usize)], spec: &str, hpc: bool) -> Vec<(u64, u32, u32)> {
    let seeder = BaseSeeder::new(spec, hpc).expect("valid seed spec");
    segs.par_iter()
        .flat_map_iter(|&(a, b)| match seeder.parsed(&seq[a..b]) {
            Some(h) => seeder
                .seed_parsed(&h.seq)
                .into_iter()
                .map(|s| {
                    let (rs, re) = h.to_raw(s.start, s.end);
                    (s.label, a as u32 + rs, re - rs)
                })
                .collect::<Vec<_>>(),
            None => seeder
                .seed_parsed(&seq[a..b])
                .into_iter()
                .map(|s| (s.label, a as u32 + s.start, s.end - s.start))
                .collect::<Vec<_>>(),
        })
        .collect()
}

/// The label of two consecutive level-2 cores taken as one seed.
#[inline]
pub fn pair_label(a: u64, b: u64) -> u64 {
    (a.rotate_left(23) ^ b).wrapping_mul(0x9E37_79B9_7F4A_7C15)
}

/// Whether a pair label is in the fallback table (every other one is).
#[inline]
pub fn keep_pair(l: u64) -> bool {
    (l >> 40) & 1 == 0
}

/// A label folded to 16 bits.
#[inline]
fn fold16(l: u64) -> u16 {
    (l ^ (l >> 16) ^ (l >> 32) ^ (l >> 48)) as u16
}

impl GenomeIndex {
    /// Build from a FASTA file (plain text).
    pub fn build(fasta: &Path, opt: &BuildOptions) -> Result<GenomeIndex> {
        let (seq, names, offsets) = read_fasta(fasta)?;
        let segs = segments(&seq, &offsets);

        // Packed sequence and the runs of other symbols.
        let mut words = vec![0u64; seq.len().div_ceil(32) + 1];
        words.par_iter_mut().enumerate().for_each(|(w, out)| {
            let lo = w * 32;
            let hi = (lo + 32).min(seq.len());
            let mut x = 0u64;
            for (k, &b) in seq[lo.min(hi)..hi].iter().enumerate() {
                x |= ((code(b) & 3) as u64) << (2 * k);
            }
            *out = x;
        });
        let mut n_runs: Vec<u32> = Vec::new();
        let mut i = 0;
        while i < seq.len() {
            if code(seq[i]) > 3 {
                let s = i;
                while i < seq.len() && code(seq[i]) > 3 {
                    i += 1;
                }
                n_runs.extend([s as u32, i as u32]);
            } else {
                i += 1;
            }
        }

        // Level 3: sorted by (label, position); frequent labels dropped.
        let mut l3 = seed_all(&seq, &segs, L3_SPEC, opt.hpc);
        l3.par_sort_unstable();
        let mut labels: Vec<u64> = Vec::new();
        let mut occ_off: Vec<u32> = vec![0];
        let mut pos: Vec<u32> = Vec::with_capacity(l3.len());
        let mut len: Vec<u8> = Vec::with_capacity(l3.len());
        let mut dropped = 0u64;
        let mut i = 0;
        while i < l3.len() {
            let mut j = i + 1;
            while j < l3.len() && l3[j].0 == l3[i].0 {
                j += 1;
            }
            if j - i <= opt.max_occ as usize {
                labels.push(l3[i].0);
                pos.extend(l3[i..j].iter().map(|x| x.1));
                len.extend(l3[i..j].iter().map(|x| x.2.min(255) as u8));
                occ_off.push(pos.len() as u32);
            } else {
                dropped += (j - i) as u64;
            }
            i = j;
        }
        drop(l3);
        let dir_bits = (labels.len().max(2).ilog2()).clamp(1, 28);
        let mut dir = vec![0u32; (1usize << dir_bits) + 1];
        for &l in &labels {
            dir[(l >> (64 - dir_bits)) as usize + 1] += 1;
        }
        for b in 1..dir.len() {
            dir[b] += dir[b - 1];
        }

        // Level 2 in genome order.
        let mut l2 = seed_all(&seq, &segs, L2_SPEC, opt.hpc);
        let mut l2_labels: Vec<u64> = l2.par_iter().map(|x| x.0).collect();
        l2_labels.par_sort_unstable();
        l2_labels.dedup();
        // Ids are ranks when the labels fit 16 bits (they do with compression: about 21,000
        // on the human genome); otherwise a 16-bit fold of the label, whose collisions only
        // add chance hits that the gap solver discards.
        let ranked = l2_labels.len() <= u16::MAX as usize;
        l2.par_sort_unstable_by_key(|x| x.1);
        let l2_pos: Vec<u32> = l2.par_iter().map(|x| x.1).collect();
        let l2_id: Vec<u16> = l2
            .par_iter()
            .map(|x| {
                if ranked {
                    l2_labels.binary_search(&x.0).unwrap() as u16
                } else {
                    fold16(x.0)
                }
            })
            .collect();
        if !ranked {
            l2_labels.clear();
        }
        // Fallback seeds for reads that level 3 anchors weakly: pairs of consecutive level-2
        // cores (every other pair, chosen by label), as a label table like level 3's.
        let mut pairs: Vec<(u64, u32, u8)> = l2
            .par_windows(2)
            .filter_map(|w| {
                let l = pair_label(w[0].0, w[1].0);
                (w[1].1 - w[0].1 < 64 && keep_pair(l))
                    .then(|| (l, w[0].1, (w[1].1 + w[1].2 - w[0].1).min(255) as u8))
            })
            .collect();
        pairs.par_sort_unstable();
        let (mut p_labels, mut p_occ_off, mut p_pos, mut p_len) =
            (Vec::new(), vec![0u32], Vec::new(), Vec::new());
        let mut i = 0;
        while i < pairs.len() {
            let mut j = i + 1;
            while j < pairs.len() && pairs[j].0 == pairs[i].0 {
                j += 1;
            }
            if j - i <= opt.max_occ as usize {
                p_labels.push(pairs[i].0);
                p_pos.extend(pairs[i..j].iter().map(|x| x.1));
                p_len.extend(pairs[i..j].iter().map(|x| x.2));
                p_occ_off.push(p_pos.len() as u32);
            }
            i = j;
        }
        drop(pairs);
        let p_bits = (p_labels.len().max(2).ilog2()).clamp(1, 28);
        let mut p_dir = vec![0u32; (1usize << p_bits) + 1];
        for &l in &p_labels {
            p_dir[(l >> (64 - p_bits)) as usize + 1] += 1;
        }
        for b in 1..p_dir.len() {
            p_dir[b] += p_dir[b - 1];
        }
        drop(l2);
        drop(seq);

        let mut name_bytes = Vec::new();
        for n in &names {
            name_bytes.extend_from_slice(n.as_bytes());
            name_bytes.push(b'\n');
        }
        let meta = vec![
            dir_bits as u64,
            opt.max_occ as u64,
            dropped,
            opt.hpc as u64,
            ranked as u64,
            p_bits as u64,
        ];
        let mut w = Writer::default();
        w.add_owned("g_meta", meta);
        w.add_owned("g_names", name_bytes);
        w.add_owned("g_offsets", offsets);
        w.add_owned("g_seq", words);
        w.add_owned("g_n_runs", n_runs);
        w.add_owned("g3_labels", labels);
        w.add_owned("g3_dir", dir);
        w.add_owned("g3_occ_off", occ_off);
        w.add_owned("g3_pos", pos);
        w.add_owned("g3_len", len);
        w.add_owned("g2_labels", l2_labels);
        w.add_owned("g2_pos", l2_pos);
        w.add_owned("g2_id", l2_id);
        w.add_owned("gp_labels", p_labels);
        w.add_owned("gp_dir", p_dir);
        w.add_owned("gp_occ_off", p_occ_off);
        w.add_owned("gp_pos", p_pos);
        w.add_owned("gp_len", p_len);
        Self::from_store(Store::from_bytes(&w.to_bytes())?)
    }

    pub fn open(path: &Path) -> Result<GenomeIndex> {
        Self::from_store(Store::open(path)?)
    }

    pub fn write(&self, path: &Path) -> Result<u64> {
        let mut w = Writer::default();
        let optional = ["gp_labels", "gp_dir", "gp_occ_off", "gp_pos", "gp_len"];
        for name in optional.iter().copied().filter(|n| self.store.has(n)) {
            w.add::<u8>(name, self.store.get(name)?);
        }
        for name in [
            "g_meta",
            "g_names",
            "g_offsets",
            "g_seq",
            "g_n_runs",
            "g3_labels",
            "g3_dir",
            "g3_occ_off",
            "g3_pos",
            "g3_len",
            "g2_labels",
            "g2_pos",
            "g2_id",
        ] {
            w.add::<u8>(name, self.store.get(name)?);
        }
        w.write(path)
    }

    fn from_store(store: Store) -> Result<GenomeIndex> {
        ensure!(store.has("g_meta"), "not a spearfish genome index");
        let meta: &[u64] = store.get("g_meta")?;
        let names = std::str::from_utf8(store.get::<u8>("g_names")?)?
            .lines()
            .map(str::to_string)
            .collect();
        let n = Raw {
            seq: store.get("g_seq")?,
            offsets: store.get("g_offsets")?,
            labels: store.get("g3_labels")?,
            dir: store.get("g3_dir")?,
            occ_off: store.get("g3_occ_off")?,
            pos: store.get("g3_pos")?,
            len: store.get("g3_len")?,
            l2_labels: store.get("g2_labels")?,
            l2_pos: store.get("g2_pos")?,
            l2_id: store.get("g2_id")?,
            p_labels: if store.has("gp_labels") {
                store.get("gp_labels")?
            } else {
                &[]
            },
            p_dir: if store.has("gp_dir") {
                store.get("gp_dir")?
            } else {
                &[]
            },
            p_occ_off: if store.has("gp_occ_off") {
                store.get("gp_occ_off")?
            } else {
                &[]
            },
            p_pos: if store.has("gp_pos") {
                store.get("gp_pos")?
            } else {
                &[]
            },
            p_len: if store.has("gp_len") {
                store.get("gp_len")?
            } else {
                &[]
            },
            p_bits: meta.get(5).copied().unwrap_or(0) as u32,
            dir_bits: meta[0] as u32,
            hpc: meta.get(3).is_none_or(|&x| x != 0),
            ranked: meta.get(4).is_none_or(|&x| x != 0),
        };
        Ok(GenomeIndex {
            store,
            names,
            n,
            junctions: None,
        })
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Start of each sequence on the concatenated axis, plus the total length.
    pub fn offsets(&self) -> &[u32] {
        // SAFETY: see `Raw`.
        unsafe { &*self.n.offsets }
    }

    pub fn total_len(&self) -> u32 {
        *self.offsets().last().unwrap()
    }

    /// The sequence that holds position `p`, and `p` within it.
    pub fn locate(&self, p: u32) -> (usize, u32) {
        let o = self.offsets();
        let i = o.partition_point(|&x| x <= p) - 1;
        (i, p - o[i])
    }

    /// Base code (0..4) at `p`.
    #[inline]
    pub fn base(&self, p: u32) -> u8 {
        // SAFETY: see `Raw`.
        let w = unsafe { &*self.n.seq };
        ((w[(p >> 5) as usize] >> ((p & 31) * 2)) & 3) as u8
    }

    /// Base codes of `[start, end)` into `out` (cleared first).
    pub fn bases(&self, start: u32, end: u32, out: &mut Vec<u8>) {
        out.clear();
        out.extend((start..end).map(|p| self.base(p)));
    }

    /// Calls `f(s)` for every position `s` in `[start, end)` where the genome reads base code
    /// `b0` at `s` and `b1` at `s + 1`, in increasing order. This is how splice-site
    /// candidates are found: 32 positions are tested per machine word of the packed sequence.
    pub fn pair_scan(&self, start: u32, end: u32, b0: u8, b1: u8, mut f: impl FnMut(u32)) {
        // SAFETY: see `Raw`.
        let w = unsafe { &*self.n.seq };
        let end = end.min(self.total_len().saturating_sub(1));
        if start >= end {
            return;
        }
        const LOW: u64 = 0x5555_5555_5555_5555;
        // bit 2i set where base i of the word equals `b`
        let eq = |x: u64, b: u8| {
            let d = x ^ (LOW * b as u64);
            !(d | (d >> 1)) & LOW
        };
        let (first, last) = ((start >> 5) as usize, ((end - 1) >> 5) as usize);
        for wi in first..=last {
            let x = w[wi];
            // base i + 1 of this word, with the next word's first base as base 32
            let next = (x >> 2) | (w.get(wi + 1).copied().unwrap_or(0) << 62);
            let mut hits = eq(x, b0) & eq(next, b1);
            let base = (wi as u32) << 5;
            if wi == first {
                hits &= !0u64 << (2 * (start - base));
            }
            if wi == last && end - base < 32 {
                hits &= (1u64 << (2 * (end - base))) - 1;
            }
            while hits != 0 {
                f(base + hits.trailing_zeros() / 2);
                hits &= hits - 1;
            }
        }
    }

    /// Occurrences (raw start, raw length) of a level-3 label.
    #[inline]
    pub fn occurrences(&self, label: u64) -> (&[u32], &[u8]) {
        // SAFETY: see `Raw`.
        let (labels, dir, off, pos, len) = unsafe {
            (
                &*self.n.labels,
                &*self.n.dir,
                &*self.n.occ_off,
                &*self.n.pos,
                &*self.n.len,
            )
        };
        let b = (label >> (64 - self.n.dir_bits)) as usize;
        let (lo, hi) = (dir[b] as usize, dir[b + 1] as usize);
        match labels[lo..hi].binary_search(&label) {
            Ok(i) => {
                let (s, e) = (off[lo + i] as usize, off[lo + i + 1] as usize);
                (&pos[s..e], &len[s..e])
            }
            Err(_) => (&[], &[]),
        }
    }

    /// Whether the index has the fallback table of level-2 pairs.
    pub fn has_pairs(&self) -> bool {
        // SAFETY: see `Raw`.
        self.n.p_bits > 0 && !unsafe { &*self.n.p_labels }.is_empty()
    }

    /// Occurrences (raw start, raw length) of a level-2 pair label (see [`pair_label`]).
    #[inline]
    pub fn pair_occurrences(&self, label: u64) -> (&[u32], &[u8]) {
        // SAFETY: see `Raw`.
        let (labels, dir, off, pos, len) = unsafe {
            (
                &*self.n.p_labels,
                &*self.n.p_dir,
                &*self.n.p_occ_off,
                &*self.n.p_pos,
                &*self.n.p_len,
            )
        };
        let b = (label >> (64 - self.n.p_bits)) as usize;
        let (lo, hi) = (dir[b] as usize, dir[b + 1] as usize);
        match labels[lo..hi].binary_search(&label) {
            Ok(i) => {
                let (s, e) = (off[lo + i] as usize, off[lo + i + 1] as usize);
                (&pos[s..e], &len[s..e])
            }
            Err(_) => (&[], &[]),
        }
    }

    /// Whether seeds are cores of the homopolymer-compressed sequence.
    pub fn hpc(&self) -> bool {
        self.n.hpc
    }

    /// The id of a level-2 label, if the genome can have it.
    #[inline]
    pub fn l2_id(&self, label: u64) -> Option<u16> {
        if !self.n.ranked {
            return Some(fold16(label));
        }
        // SAFETY: see `Raw`.
        unsafe { &*self.n.l2_labels }
            .binary_search(&label)
            .ok()
            .map(|i| i as u16)
    }

    /// Level-2 cores that start in `[start, end)`: their starts and label ids.
    pub fn l2_window(&self, start: u32, end: u32) -> (&[u32], &[u16]) {
        // SAFETY: see `Raw`.
        let (pos, id) = unsafe { (&*self.n.l2_pos, &*self.n.l2_id) };
        let lo = pos.partition_point(|&p| p < start);
        let hi = lo + pos[lo..].partition_point(|&p| p < end);
        (&pos[lo..hi], &id[lo..hi])
    }

    /// Section sizes in bytes, for reports.
    pub fn sizes(&self) -> Vec<(String, usize)> {
        self.store.sizes()
    }
}
