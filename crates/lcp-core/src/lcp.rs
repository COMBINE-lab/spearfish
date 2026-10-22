//! LCP in the lcptools style: level-1 cores from local rules on the 2-bit alphabet, then
//! repeated deterministic coin tossing (DCT) and the same rules on the reduced alphabet.
//!
//! Core *boundaries* reproduce lcptools (BilkentCompGen/lcptools, `lps.c`) exactly; this is
//! checked by the differential test against `lcp-oracle-sys`. Labels differ by design:
//! lcptools labels are lossy summaries, ours are 64-bit hashes of the full content
//! (level 1: the bases; level k: the child labels and their relative offsets). Short level-1
//! cores (< 32 bases) are labelled by an invertible mix of their packed 2-bit codes, so they
//! never collide with each other; level-k labels fold child labels and offsets with a
//! multiply-rotate mix and a final avalanche.

use xxhash_rust::xxh3::xxh3_64_with_seed;

use crate::bits::{BitArena, dct_value_words, value_width};
use crate::encode::acgt_segments;
use crate::simd::{self, Level as SimdLevel};
use crate::{Core, Hierarchy, Level, Parser, ParserMeta};

/// Bumped whenever a change alters any parse or label.
pub const LCP_PARSER_VERSION: u32 = 3;

/// murmur3 fmix64: an invertible avalanche.
#[inline(always)]
fn fmix64(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    x ^ (x >> 33)
}

/// A level-1 core with 2-bit codes `c`: its label, and its DCT bit string (2 bits per base,
/// last base lowest) appended to `bits` as one item.
#[inline(always)]
fn level1(c: &[u8], seed: u64, bits: &mut BitArena) -> u64 {
    if c.len() < 32 {
        let mut x = 0u64;
        for &b in c {
            x = (x << 2) | b as u64;
        }
        bits.push_item(x, 2 * c.len() as u32);
        // A leading 1 marks the length: injective, and fmix64 is a bijection.
        fmix64((x | (1u64 << (2 * c.len()))) ^ seed)
    } else {
        for &b in c.iter().rev() {
            bits.push_bits(b as u64, 2);
        }
        bits.finish_item();
        xxh3_64_with_seed(c, seed)
    }
}

/// Label of a level-k core from its children's (label, offset from the first child).
#[inline(always)]
fn label_k(children: impl Iterator<Item = (u64, u32)>, seed: u64) -> u64 {
    let mut h = seed;
    for (l, off) in children {
        h = (h ^ l).wrapping_mul(0x9e37_79b9_7f4a_7c15).rotate_left(29)
            ^ (off as u64).wrapping_mul(0xc2b2_ae3d_27d4_eb4f);
    }
    fmix64(h)
}

/// Reusable buffers for `LcpParser::top_cores`.
#[derive(Default)]
pub struct ParseScratch {
    codes: Vec<u8>,
    spans: Vec<(usize, usize)>,
    vals: Vec<u64>,
    bits: [BitArena; 2],
    /// Two levels (ping-pong): labels, starts, ends.
    lv: [(Vec<u64>, Vec<u32>, Vec<u32>); 2],
}

#[derive(Clone, Debug)]
pub struct LcpConfig {
    /// Global label seed, shared by index and reads.
    pub seed: u64,
    /// DCT iterations between levels (lcptools `DCT`, default 1).
    pub dct_iters: u8,
}

impl Default for LcpConfig {
    fn default() -> Self {
        LcpConfig {
            seed: 0x6c63_706d_6170,
            dct_iters: 1,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct LcpParser {
    cfg: LcpConfig,
}

impl LcpParser {
    pub fn new(cfg: LcpConfig) -> Self {
        assert!(cfg.dct_iters >= 1, "at least one DCT iteration is required");
        LcpParser { cfg }
    }

    #[inline]
    fn level_seed(&self, level: u8) -> u64 {
        self.cfg.seed ^ (level as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15)
    }

    /// Parse one run of 2-bit codes. Returns exactly `max_level` levels.
    fn parse_codes(&self, s: &[u8], max_level: u8) -> Vec<Level> {
        let mut levels: Vec<Level> = Vec::with_capacity(max_level as usize);
        if max_level == 0 {
            return levels;
        }

        // Level 1.
        let mut spans = Vec::with_capacity(s.len() / 2);
        parse_rules(s, |lo, hi| spans.push((lo, hi)));
        let seed1 = self.level_seed(1);
        let mut bits = BitArena::with_capacity(spans.len());
        let mut level = Level {
            cores: Vec::with_capacity(spans.len()),
            ctx_start: Vec::with_capacity(spans.len()),
        };
        for &(lo, hi) in &spans {
            level.cores.push(Core {
                label: level1(&s[lo..hi], seed1, &mut bits),
                start: lo as u32,
                end: hi as u32,
                child_lo: lo as u32,
                child_hi: hi as u32,
            });
            level.ctx_start.push(lo as u32);
        }
        levels.push(level);

        let d = self.cfg.dct_iters as usize;
        for k in 2..=max_level {
            let prev = levels.last().unwrap();
            let n = prev.cores.len();
            if n < d + 1 {
                levels.push(Level::default());
                continue;
            }
            let vals = dct(&bits, d);

            spans.clear();
            parse_rules(&vals[d..], |lo, hi| spans.push((lo + d, hi + d)));

            let seed = self.level_seed(k);
            let mut next_bits = BitArena::with_capacity(spans.len());
            let mut level = Level {
                cores: Vec::with_capacity(spans.len()),
                ctx_start: Vec::with_capacity(spans.len()),
            };
            for &(lo, hi) in &spans {
                let children = &prev.cores[lo..hi];
                // Child labels and offsets: the same children can overlap by different amounts
                // in different contexts, so labels alone would not determine the content.
                let label = label_k(
                    children.iter().map(|c| (c.label, c.start - children[0].start)),
                    seed,
                );
                level.cores.push(Core {
                    label,
                    start: children[0].start,
                    end: children.iter().map(|c| c.end).max().unwrap(),
                    child_lo: lo as u32,
                    child_hi: hi as u32,
                });
                // DCT replaced each core's start by that of the core `d` places to its left.
                level.ctx_start.push(prev.ctx_start[lo - d]);
                next_bits.push_values(&vals[lo..hi]);
            }
            levels.push(level);
            bits = next_bits;
        }
        levels
    }
}

impl LcpParser {
    /// The level-`k` cores of `seq` as (label, start, end), appended to `out` in start order;
    /// identical to `parse(seq, k).level(k)` but without building the hierarchy, reusing the
    /// buffers in `sc`.
    pub fn top_cores(&self, seq: &[u8], k: u8, sc: &mut ParseScratch, out: &mut Vec<(u64, u32, u32)>) {
        match simd::active() {
            #[cfg(target_arch = "x86_64")]
            SimdLevel::Avx512 => unsafe { self.top_cores_avx512(seq, k, sc, out) },
            #[cfg(target_arch = "x86_64")]
            SimdLevel::Avx2 => unsafe { self.top_cores_avx2(seq, k, sc, out) },
            // NEON is part of the aarch64 baseline: no entry point of its own.
            #[cfg(target_arch = "aarch64")]
            SimdLevel::Neon => self.top_cores_at(SimdLevel::Neon, seq, k, sc, out),
            _ => self.top_cores_at(SimdLevel::Scalar, seq, k, sc, out),
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx512f,avx512bw,avx512vl,avx2,fma,bmi1,bmi2,popcnt,lzcnt")]
    fn top_cores_avx512(&self, seq: &[u8], k: u8, sc: &mut ParseScratch, out: &mut Vec<(u64, u32, u32)>) {
        self.top_cores_at(SimdLevel::Avx512, seq, k, sc, out)
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2,fma,bmi1,bmi2,popcnt,lzcnt")]
    fn top_cores_avx2(&self, seq: &[u8], k: u8, sc: &mut ParseScratch, out: &mut Vec<(u64, u32, u32)>) {
        self.top_cores_at(SimdLevel::Avx2, seq, k, sc, out)
    }

    #[inline(always)]
    fn top_cores_at(
        &self,
        simd: SimdLevel,
        seq: &[u8],
        k: u8,
        sc: &mut ParseScratch,
        out: &mut Vec<(u64, u32, u32)>,
    ) {
        if k == 0 {
            return;
        }
        let d = self.cfg.dct_iters as usize;
        let mut i = 0;
        while i < seq.len() {
            // Next maximal ACGT segment.
            while i < seq.len() && crate::encode::CODE[seq[i] as usize] == crate::encode::INVALID {
                i += 1;
            }
            let off = i;
            sc.codes.clear();
            while i < seq.len() {
                let c = crate::encode::CODE[seq[i] as usize];
                if c == crate::encode::INVALID {
                    break;
                }
                sc.codes.push(c);
                i += 1;
            }
            if sc.codes.is_empty() {
                continue;
            }
            // Level 1.
            sc.spans.clear();
            parse_spans(simd, &sc.codes, 0, &mut sc.spans);
            let seed1 = self.level_seed(1);
            let cur = &mut sc.lv[0];
            cur.0.clear();
            cur.1.clear();
            cur.2.clear();
            sc.bits[0].clear();
            // Bit strings feed the next level's DCT: none needed for the top level.
            let bits = &mut sc.bits[0];
            let codes = &sc.codes;
            if k > 1 {
                cur.0.extend(
                    sc.spans
                        .iter()
                        .map(|&(lo, hi)| level1(&codes[lo..hi], seed1, bits)),
                );
            } else {
                let mut none = BitArena::default();
                cur.0.extend(sc.spans.iter().map(|&(lo, hi)| {
                    let l = level1(&codes[lo..hi], seed1, &mut none);
                    none.clear();
                    l
                }));
            }
            cur.1.extend(sc.spans.iter().map(|s| s.0 as u32));
            cur.2.extend(sc.spans.iter().map(|s| s.1 as u32));
            let mut ci = 0; // index of the current level's buffers
            for level in 2..=k {
                let n = sc.lv[ci].0.len();
                let ni = 1 - ci;
                if n < d + 1 {
                    sc.lv[ni].0.clear();
                    sc.lv[ni].1.clear();
                    sc.lv[ni].2.clear();
                    ci = ni;
                    break;
                }
                dct_into(&sc.bits[ci], d, &mut sc.vals);
                sc.spans.clear();
                parse_spans(simd, &sc.vals[d..], d, &mut sc.spans);
                let seed = self.level_seed(level);
                let (a, b) = sc.lv.split_at_mut(1);
                let (prev, next) = if ci == 0 {
                    (&a[0], &mut b[0])
                } else {
                    (&b[0], &mut a[0])
                };
                next.0.clear();
                next.1.clear();
                next.2.clear();
                let (ba, bb) = sc.bits.split_at_mut(1);
                let nbits = if ci == 0 { &mut bb[0] } else { &mut ba[0] };
                nbits.clear();
                next.0.extend(sc.spans.iter().map(|&(lo, hi)| {
                    let s0 = prev.1[lo];
                    label_k((lo..hi).map(|j| (prev.0[j], prev.1[j] - s0)), seed)
                }));
                next.1.extend(sc.spans.iter().map(|&(lo, _)| prev.1[lo]));
                // Core ends strictly increase at every level (each rule's core ends past the
                // previous core's end), so a core ends where its last child does.
                next.2.extend(sc.spans.iter().map(|&(_, hi)| prev.2[hi - 1]));
                debug_assert!(next.2.windows(2).all(|w| w[0] < w[1]));
                if level < k {
                    for &(lo, hi) in &sc.spans {
                        nbits.push_values(&sc.vals[lo..hi]);
                    }
                }
                ci = ni;
            }
            let top = &sc.lv[ci];
            let o = off as u32;
            out.extend(
                top.0
                    .iter()
                    .zip(&top.1)
                    .zip(&top.2)
                    .map(|((&l, &s), &e)| (l, s + o, e + o)),
            );
        }
    }
}

impl Parser for LcpParser {
    fn parse(&self, seq: &[u8], max_level: u8) -> Hierarchy {
        let mut h = Hierarchy {
            len: seq.len() as u32,
            levels: vec![Level::default(); max_level as usize],
        };
        for (off, codes) in acgt_segments(seq) {
            let seg = self.parse_codes(&codes, max_level);
            // Child indices of level k point into level k-1, so offset them by the number of
            // level k-1 cores already present.
            let before: Vec<u32> = h.levels.iter().map(|l| l.cores.len() as u32).collect();
            for (k, level) in seg.into_iter().enumerate() {
                let off = off as u32;
                let child_off = if k == 0 { off } else { before[k - 1] };
                let dst = &mut h.levels[k];
                dst.cores.extend(level.cores.into_iter().map(|c| Core {
                    start: c.start + off,
                    end: c.end + off,
                    child_lo: c.child_lo + child_off,
                    child_hi: c.child_hi + child_off,
                    ..c
                }));
                dst.ctx_start.extend(level.ctx_start.into_iter().map(|s| s + off));
            }
        }
        h
    }

    fn meta(&self) -> ParserMeta {
        ParserMeta {
            name: "lcp",
            version: LCP_PARSER_VERSION,
            seed: self.cfg.seed,
            dct_iters: self.cfg.dct_iters,
        }
    }
}

/// DCT values after `d` iterations. Entries `[d..]` are valid; the first `d` cores have no
/// left neighbour at some iteration and are dropped, as in lcptools.
fn dct(bits: &BitArena, d: usize) -> Vec<u64> {
    let mut vals = Vec::new();
    dct_into(bits, d, &mut vals);
    vals
}

fn dct_into(bits: &BitArena, d: usize, vals: &mut Vec<u64>) {
    let n = bits.len();
    vals.clear();
    vals.resize(n, 0);
    if n > 1 {
        bits.dct_all(&mut vals[1..]);
    }
    for t in 1..d {
        // Iteration t compresses cores i > t against their left neighbour's previous value.
        for i in (t + 1..n).rev() {
            let (l, r) = (vals[i - 1], vals[i]);
            vals[i] = dct_value_words(l, value_width(l), r, value_width(r));
        }
    }
}

/// The LCP core rules, shared by all levels (lcptools `parse1`/`parse3`). Calls
/// `emit(lo, hi)` for each core `[lo, hi)` over `s`, in order:
///
/// * RINT: a symbol, a run of length >= 2 of another symbol, and the symbol after the run
///   (skipped if the run reaches the end of `s`);
/// * LMIN: `s[i] > s[i+1] < s[i+2]`;
/// * LMAX: `s[i-1] <= s[i] < s[i+1] > s[i+2] >= s[i+3]`;
/// * SSEQ: fills the gap between consecutive cores, overlapping each by one symbol.
pub(crate) fn parse_rules<T: Ord>(s: &[T], mut emit: impl FnMut(usize, usize)) {
    let n = s.len();
    // End of the previous core; usize::MAX means "none yet" (no SSEQ before the first core).
    let mut prev_end = usize::MAX;
    let mut core = |lo: usize, hi: usize, prev_end: &mut usize| {
        if *prev_end < lo {
            emit(*prev_end - 1, lo + 1);
        }
        emit(lo, hi);
        *prev_end = hi;
    };

    let mut i = 0;
    while i + 2 < n {
        if s[i] == s[i + 1] {
            i += 1;
            continue;
        }
        if s[i + 1] == s[i + 2] {
            let mut t = i + 2;
            while t < n && s[t - 1] == s[t] {
                t += 1;
            }
            if t != n {
                core(i, t + 1, &mut prev_end);
                i += 1;
                continue;
            }
        }
        if s[i] > s[i + 1] && s[i + 1] < s[i + 2] {
            core(i, i + 3, &mut prev_end);
            i += 1;
            continue;
        }
        if i > 0
            && i + 3 < n
            && s[i] < s[i + 1]
            && s[i + 1] > s[i + 2]
            && s[i - 1] <= s[i]
            && s[i + 2] >= s[i + 3]
        {
            core(i, i + 3, &mut prev_end);
        }
        i += 1;
    }
}

/// Symbols [`parse_spans`] scans, with their vector mask kernels: a block of 64 interior
/// positions' LMIN/LMAX and run masks in a few vector compares.
pub(crate) trait Sym: Ord + Copy {
    /// (core, run) masks of positions `blk..blk + 64`, as the scalar pass computes them.
    ///
    /// # Safety
    /// The CPU has AVX-512 F/BW, and `s[blk - 1..blk + 67]` is in bounds.
    #[cfg(target_arch = "x86_64")]
    unsafe fn masks_avx512(s: &[Self], blk: usize) -> (u64, u64);
    /// [`Sym::masks_avx512`] with AVX2.
    ///
    /// # Safety
    /// The CPU has AVX2, and `s[blk - 1..blk + 67]` is in bounds.
    #[cfg(target_arch = "x86_64")]
    unsafe fn masks_avx2(s: &[Self], blk: usize) -> (u64, u64);
    /// Whether [`Sym::masks_neon`] has a kernel (else the NEON level scans these symbols in
    /// the scalar loop).
    #[cfg(target_arch = "aarch64")]
    const NEON: bool;
    /// [`Sym::masks_avx512`] with NEON.
    ///
    /// # Safety
    /// [`Sym::NEON`], and `s[blk - 1..blk + 67]` is in bounds.
    #[cfg(target_arch = "aarch64")]
    unsafe fn masks_neon(s: &[Self], blk: usize) -> (u64, u64);
}

impl Sym for u8 {
    #[cfg(target_arch = "x86_64")]
    #[inline(always)]
    unsafe fn masks_avx512(s: &[u8], blk: usize) -> (u64, u64) {
        // SAFETY: as the caller guarantees.
        unsafe { masks_u8(s, blk) }
    }
    #[cfg(target_arch = "x86_64")]
    #[inline(always)]
    unsafe fn masks_avx2(s: &[u8], blk: usize) -> (u64, u64) {
        // SAFETY: as the caller guarantees.
        unsafe { masks_u8_avx2(s, blk) }
    }
    #[cfg(target_arch = "aarch64")]
    const NEON: bool = true;
    #[cfg(target_arch = "aarch64")]
    #[inline(always)]
    unsafe fn masks_neon(s: &[u8], blk: usize) -> (u64, u64) {
        // SAFETY: as the caller guarantees (NEON is part of the aarch64 baseline).
        unsafe { masks_u8_neon(s, blk) }
    }
}

impl Sym for u64 {
    #[cfg(target_arch = "x86_64")]
    #[inline(always)]
    unsafe fn masks_avx512(s: &[u64], blk: usize) -> (u64, u64) {
        // SAFETY: as the caller guarantees.
        unsafe { masks_u64(s, blk) }
    }
    #[cfg(target_arch = "x86_64")]
    #[inline(always)]
    unsafe fn masks_avx2(s: &[u64], blk: usize) -> (u64, u64) {
        // SAFETY: as the caller guarantees.
        unsafe { masks_u64_avx2(s, blk) }
    }
    // DCT values (level 2 and up) stay in the scalar loop on aarch64.
    #[cfg(target_arch = "aarch64")]
    const NEON: bool = false;
    #[cfg(target_arch = "aarch64")]
    unsafe fn masks_neon(_: &[u64], _: usize) -> (u64, u64) {
        unreachable!()
    }
}

/// How [`spans_core`] computes a whole block's masks: the scalar loop, or a level's kernel.
trait Masks {
    /// Whether the level has a kernel for symbols `T` (a constant per instance).
    fn simd<T: Sym>() -> bool;
    /// # Safety
    /// `Self::simd::<T>()`, the CPU has the level's features, and `s[blk - 1..blk + 67]` is in
    /// bounds.
    unsafe fn masks<T: Sym>(s: &[T], blk: usize) -> (u64, u64);
}

struct ScalarMasks;

#[cfg(target_arch = "x86_64")]
struct Avx2Masks;

#[cfg(target_arch = "x86_64")]
impl Masks for Avx2Masks {
    #[inline(always)]
    fn simd<T: Sym>() -> bool {
        true
    }
    #[inline(always)]
    unsafe fn masks<T: Sym>(s: &[T], blk: usize) -> (u64, u64) {
        // SAFETY: as the caller guarantees.
        unsafe { T::masks_avx2(s, blk) }
    }
}

impl Masks for ScalarMasks {
    #[inline(always)]
    fn simd<T: Sym>() -> bool {
        false
    }
    unsafe fn masks<T: Sym>(_: &[T], _: usize) -> (u64, u64) {
        unreachable!()
    }
}

#[cfg(target_arch = "x86_64")]
struct Avx512Masks;

#[cfg(target_arch = "x86_64")]
impl Masks for Avx512Masks {
    #[inline(always)]
    fn simd<T: Sym>() -> bool {
        true
    }
    #[inline(always)]
    unsafe fn masks<T: Sym>(s: &[T], blk: usize) -> (u64, u64) {
        // SAFETY: as the caller guarantees.
        unsafe { T::masks_avx512(s, blk) }
    }
}

/// NEON kernels where [`Sym::NEON`] (base codes), else the scalar loop.
#[cfg(target_arch = "aarch64")]
struct NeonMasks;

#[cfg(target_arch = "aarch64")]
impl Masks for NeonMasks {
    #[inline(always)]
    fn simd<T: Sym>() -> bool {
        T::NEON
    }
    #[inline(always)]
    unsafe fn masks<T: Sym>(s: &[T], blk: usize) -> (u64, u64) {
        // SAFETY: as the caller guarantees.
        unsafe { T::masks_neon(s, blk) }
    }
}

// Lanes hold positions i; vectors a..e load s[i - 1], s[i], ..., s[i + 3].
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw")]
#[inline]
unsafe fn masks_u8(s: &[u8], blk: usize) -> (u64, u64) {
    use std::arch::x86_64::*;
    // SAFETY: s[blk - 1..blk + 67] is in bounds (caller).
    let ld = |o: usize| unsafe { _mm512_loadu_si512(s.as_ptr().add(blk + o - 1).cast()) };
    let (a, b, c, d, e) = (ld(0), ld(1), ld(2), ld(3), ld(4));
    let lmin = _mm512_cmpgt_epu8_mask(b, c) & _mm512_cmplt_epu8_mask(c, d);
    let lmax = _mm512_cmple_epu8_mask(a, b)
        & _mm512_cmplt_epu8_mask(b, c)
        & _mm512_cmpgt_epu8_mask(c, d)
        & _mm512_cmpge_epu8_mask(d, e);
    let run = _mm512_cmpneq_epu8_mask(b, c) & _mm512_cmpeq_epu8_mask(c, d);
    (lmin | lmax, run)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
#[inline]
unsafe fn masks_u64(s: &[u64], blk: usize) -> (u64, u64) {
    use std::arch::x86_64::*;
    let (mut core, mut run) = (0u64, 0u64);
    for g in 0..8 {
        // SAFETY: s[blk - 1..blk + 67] is in bounds (caller).
        let ld = |o: usize| unsafe { _mm512_loadu_si512(s.as_ptr().add(blk + 8 * g + o - 1).cast()) };
        let (a, b, c, d, e) = (ld(0), ld(1), ld(2), ld(3), ld(4));
        let lmin = _mm512_cmpgt_epu64_mask(b, c) & _mm512_cmplt_epu64_mask(c, d);
        let lmax = _mm512_cmple_epu64_mask(a, b)
            & _mm512_cmplt_epu64_mask(b, c)
            & _mm512_cmpgt_epu64_mask(c, d)
            & _mm512_cmpge_epu64_mask(d, e);
        let r = _mm512_cmpneq_epu64_mask(b, c) & _mm512_cmpeq_epu64_mask(c, d);
        core |= ((lmin | lmax) as u64) << (8 * g);
        run |= (r as u64) << (8 * g);
    }
    (core, run)
}

// [`masks_u8`] in two 32-byte halves: unsigned compares are signed ones after flipping the
// sign bits, and `movemask` gathers the lane bits.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn masks_u8_avx2(s: &[u8], blk: usize) -> (u64, u64) {
    use std::arch::x86_64::*;
    let bias = _mm256_set1_epi8(i8::MIN);
    let half = |h: usize| {
        // SAFETY: s[blk - 1..blk + 67] is in bounds (caller).
        let ld = |o: usize| unsafe {
            _mm256_xor_si256(
                _mm256_loadu_si256(s.as_ptr().add(blk + 32 * h + o - 1).cast()),
                bias,
            )
        };
        let (a, b, c, d, e) = (ld(0), ld(1), ld(2), ld(3), ld(4));
        let gt = |x, y| _mm256_cmpgt_epi8(x, y);
        // LMIN: b > c < d. LMAX: a <= b < c > d >= e. Run: b != c == d.
        let lmin = _mm256_and_si256(gt(b, c), gt(d, c));
        let lmax = _mm256_andnot_si256(
            gt(a, b),
            _mm256_and_si256(gt(c, b), _mm256_andnot_si256(gt(e, d), gt(c, d))),
        );
        let run = _mm256_andnot_si256(_mm256_cmpeq_epi8(b, c), _mm256_cmpeq_epi8(c, d));
        let bits = |v| _mm256_movemask_epi8(v) as u32 as u64;
        (bits(_mm256_or_si256(lmin, lmax)), bits(run))
    };
    let ((c0, r0), (c1, r1)) = (half(0), half(1));
    (c0 | (c1 << 32), r0 | (r1 << 32))
}

// [`masks_u64`] in sixteen groups of four 64-bit lanes, unsigned compares as for
// [`masks_u8_avx2`].
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn masks_u64_avx2(s: &[u64], blk: usize) -> (u64, u64) {
    use std::arch::x86_64::*;
    let bias = _mm256_set1_epi64x(i64::MIN);
    let (mut core, mut run) = (0u64, 0u64);
    for g in 0..16 {
        // SAFETY: s[blk - 1..blk + 67] is in bounds (caller).
        let ld = |o: usize| unsafe {
            _mm256_xor_si256(
                _mm256_loadu_si256(s.as_ptr().add(blk + 4 * g + o - 1).cast()),
                bias,
            )
        };
        let (a, b, c, d, e) = (ld(0), ld(1), ld(2), ld(3), ld(4));
        let gt = |x, y| _mm256_cmpgt_epi64(x, y);
        let lmin = _mm256_and_si256(gt(b, c), gt(d, c));
        let lmax = _mm256_andnot_si256(
            gt(a, b),
            _mm256_and_si256(gt(c, b), _mm256_andnot_si256(gt(e, d), gt(c, d))),
        );
        let r = _mm256_andnot_si256(_mm256_cmpeq_epi64(b, c), _mm256_cmpeq_epi64(c, d));
        let bits = |v| _mm256_movemask_pd(_mm256_castsi256_pd(v)) as u64;
        core |= bits(_mm256_or_si256(lmin, lmax)) << (4 * g);
        run |= bits(r) << (4 * g);
    }
    (core, run)
}

/// [`masks_u8`] in four 16-byte NEON blocks. NEON has no movemask: each compare result
/// (0 or 0xff per byte) keeps its lane's bit within the byte's group of eight, and three
/// rounds of pairwise byte sums pack the 64 lanes into one word, lane i at bit i.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn masks_u8_neon(s: &[u8], blk: usize) -> (u64, u64) {
    use std::arch::aarch64::*;
    let weights: [u8; 16] = [1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128];
    // SAFETY: 16 bytes.
    let w = unsafe { vld1q_u8(weights.as_ptr()) };
    let pack = |m: [uint8x16_t; 4]| {
        let (a, b) = (vandq_u8(m[0], w), vandq_u8(m[1], w));
        let (c, d) = (vandq_u8(m[2], w), vandq_u8(m[3], w));
        let x = vpaddq_u8(vpaddq_u8(a, b), vpaddq_u8(c, d));
        vgetq_lane_u64::<0>(vreinterpretq_u64_u8(vpaddq_u8(x, x)))
    };
    let (mut core, mut run) = ([vdupq_n_u8(0); 4], [vdupq_n_u8(0); 4]);
    for q in 0..4 {
        // SAFETY: s[blk - 1..blk + 67] is in bounds (caller); the last load ends there.
        let ld = |o: usize| unsafe { vld1q_u8(s.as_ptr().add(blk + 16 * q + o - 1)) };
        let (a, b, c, d, e) = (ld(0), ld(1), ld(2), ld(3), ld(4));
        let lmin = vandq_u8(vcgtq_u8(b, c), vcltq_u8(c, d));
        let lmax = vandq_u8(
            vandq_u8(vcleq_u8(a, b), vcltq_u8(b, c)),
            vandq_u8(vcgtq_u8(c, d), vcgeq_u8(d, e)),
        );
        core[q] = vorrq_u8(lmin, lmax);
        // b != c and c == d.
        run[q] = vbicq_u8(vceqq_u8(c, d), vceqq_u8(b, c));
    }
    (pack(core), pack(run))
}

/// [`parse_rules`] in two phases per block of 64 positions: a branch-free pass computes which
/// positions start an LMIN/LMAX core or a run (RINT candidate) as bit masks, then a loop over
/// the set bits emits the cores and the SSEQ gap fillers between them. Positions that start
/// nothing cost no branch. Appends `(lo + add, hi + add)` to `spans` in `parse_rules` order.
///
/// Runs at `level` (which the CPU must support): the whole scan is compiled once per level,
/// so the block kernels inline into it.
pub(crate) fn parse_spans<T: Sym>(level: SimdLevel, s: &[T], add: usize, spans: &mut Vec<(usize, usize)>) {
    assert!(
        simd::detected().includes(level),
        "SIMD level {level} not supported"
    );
    match level {
        // SAFETY: the CPU has the level's features (asserted above).
        #[cfg(target_arch = "x86_64")]
        SimdLevel::Avx512 => unsafe { spans_avx512(s, add, spans) },
        #[cfg(target_arch = "x86_64")]
        SimdLevel::Avx2 => unsafe { spans_avx2(s, add, spans) },
        // NEON is part of the aarch64 baseline: no entry point of its own.
        #[cfg(target_arch = "aarch64")]
        SimdLevel::Neon => spans_core::<T, NeonMasks>(s, add, spans),
        _ => spans_core::<T, ScalarMasks>(s, add, spans),
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vl,avx2,fma,bmi1,bmi2,popcnt,lzcnt")]
fn spans_avx512<T: Sym>(s: &[T], add: usize, spans: &mut Vec<(usize, usize)>) {
    spans_core::<T, Avx512Masks>(s, add, spans)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma,bmi1,bmi2,popcnt,lzcnt")]
fn spans_avx2<T: Sym>(s: &[T], add: usize, spans: &mut Vec<(usize, usize)>) {
    spans_core::<T, Avx2Masks>(s, add, spans)
}

#[inline(always)]
fn spans_core<T: Sym, M: Masks>(s: &[T], add: usize, spans: &mut Vec<(usize, usize)>) {
    let n = s.len();
    if n < 3 {
        return;
    }
    let m = n - 2;
    // At most two spans per position.
    spans.reserve(2 * n);
    let base = spans.len();
    let out = &mut spans.spare_capacity_mut()[..2 * n];
    let mut cnt = 0;
    let mut prev_end = usize::MAX;
    let mut blk = 0;
    while blk < m {
        let end = (blk + 64).min(m);
        let (mut core, mut run) = (0u64, 0u64);
        // Whole interior blocks (every position has a symbol on each side) in SIMD.
        let simd = if M::simd::<T>() && blk > 0 && blk + 67 <= n {
            // SAFETY: the CPU has the level's features (the caller's entry point is compiled
            // for them), and s[blk - 1..blk + 67] is in bounds.
            (core, run) = unsafe { M::masks(s, blk) };
            true
        } else {
            false
        };
        for i in (blk..end).filter(|_| !simd) {
            let (a, b, c) = (s[i], s[i + 1], s[i + 2]);
            let (l, r) = (s[i.saturating_sub(1)], s[(i + 3).min(n - 1)]);
            // LMIN, or LMAX (needs a symbol on each side); exclusive with each other and RINT.
            let lmin = (a > b) & (b < c);
            let lmax = (i > 0) & (i + 3 < n) & (l <= a) & (a < b) & (b > c) & (c >= r);
            core |= ((lmin | lmax) as u64) << (i - blk);
            run |= (((a != b) & (b == c)) as u64) << (i - blk);
        }
        let mut bits = core | run;
        while bits != 0 {
            let j = bits.trailing_zeros();
            bits &= bits - 1;
            let i = blk + j as usize;
            let mut hi = i + 3;
            if (run >> j) & 1 != 0 {
                // RINT: the run of s[i + 1] must end before the end of `s`.
                let mut t = i + 2;
                while t < n && s[t - 1] == s[t] {
                    t += 1;
                }
                if t == n {
                    continue;
                }
                hi = t + 1;
            }
            out[cnt].write((prev_end.wrapping_sub(1).wrapping_add(add), i + 1 + add));
            cnt += (prev_end < i) as usize;
            out[cnt].write((i + add, hi + add));
            cnt += 1;
            prev_end = hi;
        }
        blk = end;
    }
    // SAFETY: the first `cnt` spare slots were written.
    unsafe { spans.set_len(base + cnt) };
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spans(s: &[u8]) -> Vec<(usize, usize)> {
        let codes: Vec<u8> = s.iter().map(|&b| crate::encode::CODE[b as usize]).collect();
        let mut v = Vec::new();
        parse_rules(&codes, |lo, hi| v.push((lo, hi)));
        v
    }

    /// The block masks by `level`'s kernel for `T`, called directly (`None` if the level scans
    /// `T` in the scalar loop).
    fn masks_at<T: Sym>(level: SimdLevel, s: &[T], blk: usize) -> Option<(u64, u64)> {
        assert!(simd::detected().includes(level) && blk > 0 && blk + 67 <= s.len());
        match level {
            // SAFETY: the CPU supports the level (asserted above), and s[blk - 1..blk + 67] is
            // in bounds.
            #[cfg(target_arch = "x86_64")]
            SimdLevel::Avx512 => Some(unsafe { T::masks_avx512(s, blk) }),
            #[cfg(target_arch = "x86_64")]
            SimdLevel::Avx2 => Some(unsafe { T::masks_avx2(s, blk) }),
            #[cfg(target_arch = "aarch64")]
            SimdLevel::Neon if T::NEON => Some(unsafe { T::masks_neon(s, blk) }),
            _ => None,
        }
    }

    /// The scalar scan's masks of the whole block `blk..blk + 64`.
    fn masks_scalar<T: Sym>(s: &[T], blk: usize) -> (u64, u64) {
        let (mut core, mut run) = (0u64, 0u64);
        for i in blk..blk + 64 {
            let (l, a, b, c, r) = (s[i - 1], s[i], s[i + 1], s[i + 2], s[i + 3]);
            core |= ((((a > b) & (b < c)) | ((l <= a) & (a < b) & (b > c) & (c >= r))) as u64) << (i - blk);
            run |= (((a != b) & (b == c)) as u64) << (i - blk);
        }
        (core, run)
    }

    proptest::proptest! {
        /// Every level's vector masks equal the scalar scan's on whole blocks, over all byte
        /// values and across the u64 range.
        #[test]
        fn masks_equal_scalar(
            s in proptest::collection::vec(proptest::prelude::any::<u8>(), 68..300),
            at in 1usize..300,
            hi in proptest::collection::vec(proptest::sample::select(vec![0u64, 1, (1 << 63) - 1, 1 << 63, u64::MAX - 1, u64::MAX]), 300),
        ) {
            let n = s.len();
            let blk = 1 + at % (n - 67);
            // Wide symbols whose order the low byte decides only among equal high parts.
            let w: Vec<u64> = s.iter().zip(&hi).map(|(&b, &h)| (h & !0xff) | b as u64).collect();
            for level in crate::simd::supported() {
                if let Some(got) = masks_at(level, &s, blk) {
                    proptest::prop_assert_eq!(got, masks_scalar(&s, blk), "u8 at {}", level);
                }
                if let Some(got) = masks_at(level, &w, blk) {
                    proptest::prop_assert_eq!(got, masks_scalar(&w, blk), "u64 at {}", level);
                }
            }
        }
    }

    proptest::proptest! {
        /// The branch-free scan emits exactly what the rules do (small alphabets for runs).
        #[test]
        fn parse_spans_equal_rules(
            s in proptest::collection::vec(0u64..4, 0..600),
            add in 0usize..5,
        ) {
            // Every level this CPU supports, each called explicitly.
            for level in crate::simd::supported() {
                let mut want = Vec::new();
                parse_rules(&s, |lo, hi| want.push((lo + add, hi + add)));
                let mut got = vec![(7, 7)];
                parse_spans(level, &s, add, &mut got);
                proptest::prop_assert_eq!(&got[1..], &want[..], "{}", level);
                let b: Vec<u8> = s.iter().map(|&x| (x as u8) >> 1).collect();
                let mut want = Vec::new();
                parse_rules(&b, |lo, hi| want.push((lo, hi)));
                let mut got = Vec::new();
                parse_spans(level, &b, 0, &mut got);
                proptest::prop_assert_eq!(got, want, "{}", level);
            }
        }
    }

    proptest::proptest! {
        /// The same with symbols across the whole unsigned range, where the vector kernels'
        /// signed compares need the sign-bit flip, and with long runs of equal symbols.
        #[test]
        fn parse_spans_equal_rules_wide(
            s in proptest::collection::vec(proptest::sample::select(vec![0u64, 1, 2, (1 << 63) - 1, 1 << 63, (1 << 63) + 1, u64::MAX - 1, u64::MAX]), 0..600),
            run in 1usize..6,
        ) {
            let s: Vec<u64> = s.iter().flat_map(|&x| std::iter::repeat_n(x, if x & 1 == 0 { run } else { 1 })).collect();
            let b: Vec<u8> = s.iter().map(|&x| (x >> 56) as u8 ^ (x as u8)).collect();
            for level in crate::simd::supported() {
                let mut want = Vec::new();
                parse_rules(&s, |lo, hi| want.push((lo, hi)));
                let mut got = Vec::new();
                parse_spans(level, &s, 0, &mut got);
                proptest::prop_assert_eq!(got, want, "u64 at {}", level);
                let mut want = Vec::new();
                parse_rules(&b, |lo, hi| want.push((lo, hi)));
                let mut got = Vec::new();
                parse_spans(level, &b, 0, &mut got);
                proptest::prop_assert_eq!(got, want, "u8 at {}", level);
            }
        }
    }

    /// Isolated block-scan timings per SIMD level (`cargo test --release -p lcp-core --
    /// --ignored --nocapture scan_timings`).
    #[test]
    #[ignore]
    fn scan_timings() {
        use std::hint::black_box;
        let mut x = 5u64;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let n = 1 << 20;
        let codes: Vec<u8> = (0..n).map(|_| (next() & 3) as u8).collect();
        // Level-2 symbols (DCT values) are small integers.
        let vals: Vec<u64> = (0..n).map(|_| next() % 13).collect();
        let mut spans = Vec::with_capacity(2 * n);
        for level in crate::simd::supported() {
            let mut t8 = f64::INFINITY;
            let mut t64 = f64::INFINITY;
            for _ in 0..5 {
                spans.clear();
                let t = std::time::Instant::now();
                parse_spans(level, black_box(&codes), 0, &mut spans);
                t8 = t8.min(t.elapsed().as_nanos() as f64 / n as f64);
                spans.clear();
                let t = std::time::Instant::now();
                parse_spans(level, black_box(&vals), 0, &mut spans);
                t64 = t64.min(t.elapsed().as_nanos() as f64 / n as f64);
            }
            println!(
                "parse_spans {:7}: u8 {t8:.3} ns/symbol, u64 {t64:.3} ns/symbol",
                level.name()
            );
        }
    }

    #[test]
    fn lmin_cores() {
        // A=0 C=1 G=2 T=3. "CAC" is a local minimum at A.
        assert_eq!(spans(b"CACAC"), vec![(0, 3), (2, 5)]);
    }

    #[test]
    fn rint_core_includes_flanks() {
        // C AAAA G: RINT over [0, 6); run reaching the end is not a core.
        assert_eq!(spans(b"CAAAAG"), vec![(0, 6)]);
        assert_eq!(spans(b"CAAAA"), vec![]);
    }

    #[test]
    fn sseq_fills_gaps() {
        // ACGT is strictly increasing: no rule fires until the LMIN "TAC" at 3.
        // Then "ACGTTTC" ... construct: LMIN at 0 ("CAG"), plain ascent, LMIN later.
        let v = spans(b"CAGTTTTGAC");
        // Every gap between consecutive cores is filled.
        for w in v.windows(2) {
            assert!(w[1].0 <= w[0].1, "gap between {:?} and {:?}", w[0], w[1]);
        }
    }

    #[test]
    fn hierarchy_children_are_consistent() {
        let p = LcpParser::default();
        let seq: Vec<u8> = (0..5000u64)
            .map(|i| b"ACGT"[((i.wrapping_mul(0x9e3779b97f4a7c15) >> 61) & 3) as usize])
            .collect();
        let h = p.parse(&seq, 5);
        assert_eq!(h.num_levels(), 5);
        for k in 2..=5u8 {
            let prev = h.level(k - 1);
            for c in h.level(k) {
                let ch = &prev[c.child_lo as usize..c.child_hi as usize];
                assert!(ch.len() >= 2);
                assert_eq!(c.start, ch[0].start);
                assert_eq!(c.end, ch.iter().map(|x| x.end).max().unwrap());
            }
        }
    }

    #[test]
    fn n_splits_parse() {
        let p = LcpParser::default();
        let seq = b"CACACAGTGTCANNNNCACACAGTGTCA";
        let h = p.parse(seq, 1);
        for c in h.level(1) {
            assert!(!seq[c.start as usize..c.end as usize].contains(&b'N'));
        }
        // Same content on both sides of the Ns gives the same labels.
        let left: Vec<u64> = h
            .level(1)
            .iter()
            .filter(|c| c.end <= 12)
            .map(|c| c.label)
            .collect();
        let right: Vec<u64> = h
            .level(1)
            .iter()
            .filter(|c| c.start >= 16)
            .map(|c| c.label)
            .collect();
        assert_eq!(left, right);
    }
}

#[cfg(test)]
mod top_cores_tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(400))]

        /// The allocation-free path equals the hierarchy's top level exactly.
        #[test]
        fn top_cores_equal_parse(
            seq in prop::collection::vec(prop::sample::select(b"AACGTTTGN".to_vec()), 0..3000),
            k in 1u8..=5,
        ) {
            let p = LcpParser::default();
            let want: Vec<(u64, u32, u32)> = p.parse(&seq, k).level(k).iter().map(|c| (c.label, c.start, c.end)).collect();
            let mut sc = ParseScratch::default();
            let mut got = Vec::new();
            p.top_cores(&seq, k, &mut sc, &mut got);
            prop_assert_eq!(&got, &want);
            // Reusing the scratch gives the same result.
            let mut again = Vec::new();
            p.top_cores(&seq, k, &mut sc, &mut again);
            prop_assert_eq!(again, want);
        }
    }
}

#[cfg(test)]
mod label_tests {
    use super::*;

    /// Equal labels at a level imply equal content (spans of equal length and equal bases).
    #[test]
    fn labels_determine_content() {
        let p = LcpParser::default();
        let mut x = 99u64;
        let mut seen: std::collections::HashMap<(usize, u64), Vec<u8>> = std::collections::HashMap::new();
        for _ in 0..200 {
            // Low-entropy sequence: many repeated substrings with varying contexts.
            let seq: Vec<u8> = (0..3000)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    b"AACGT"[(x % 5) as usize]
                })
                .collect();
            let h = p.parse(&seq, 5);
            for (k, l) in h.levels.iter().enumerate() {
                for c in &l.cores {
                    let content = seq[c.start as usize..c.end as usize].to_vec();
                    let e = seen.entry((k, c.label)).or_insert_with(|| content.clone());
                    assert_eq!(*e, content, "level {} label {:x}", k + 1, c.label);
                }
            }
        }
    }
}
