//! Edit distance with Myers' bit-parallel algorithm in multi-word blocks and Ukkonen's
//! cutoff (as in edlib). Inputs are 2-bit codes (0..4). Three modes:
//!
//! - `Infix`: pattern fully aligned, both text ends free (best substring of the text).
//! - `Prefix`: pattern fully aligned starting at text position 0, text end free.
//! - `Global`: pattern and text fully aligned.
//!
//! With a threshold `k`, only blocks of pattern rows that can still score ≤ k are computed
//! (the band), and the result is `None` if the distance exceeds `k`. Cost
//! O(⌈min(m, k)/64⌉ · n) word operations in practice, O(⌈m/64⌉ · n) at worst.

use lcp_core::simd::{self, Level};

/// The SIMD level the kernels here and in `chain` run at: [`simd::active`], resolved once per
/// call (a DP, a batch) by the caller and passed down. Tests force a level on their thread.
#[inline]
pub(crate) fn simd_level() -> Level {
    #[cfg(test)]
    if let Some(l) = FORCE_LEVEL.with(|f| f.get()) {
        return l;
    }
    simd::active()
}

#[cfg(test)]
thread_local! {
    static FORCE_LEVEL: std::cell::Cell<Option<Level>> = const { std::cell::Cell::new(None) };
}

/// Runs `f` with the kernels on this thread at `level` (which the CPU must support).
#[cfg(test)]
pub(crate) fn at_level<R>(level: Level, f: impl FnOnce() -> R) -> R {
    assert!(simd::detected().includes(level));
    let old = FORCE_LEVEL.with(|x| x.replace(Some(level)));
    let r = f();
    FORCE_LEVEL.with(|x| x.set(old));
    r
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Infix,
    Prefix,
    Global,
}

const W: usize = 64;

/// One Myers block step; returns the horizontal delta out of row `bit`.
#[inline(always)]
fn block(pv: &mut u64, mv: &mut u64, eq: u64, hin: i32, bit: usize) -> i32 {
    // The incoming delta's sign as bits, without branches (it is unpredictable).
    let (neg, pos) = ((hin < 0) as u64, (hin > 0) as u64);
    let (p, m) = (*pv, *mv);
    let xv = eq | m;
    let eq = eq | neg;
    let xh = ((eq & p).wrapping_add(p) ^ p) | eq;
    let ph = m | !(xh | p);
    let mh = p & xh;
    let hout = ((ph >> bit) & 1) as i32 - ((mh >> bit) & 1) as i32;
    let ph = (ph << 1) | pos;
    let mh = (mh << 1) | neg;
    *pv = mh | !(xv | ph);
    *mv = ph & xv;
    hout
}

/// A pattern's match masks per 64-row block, built once and reused across texts.
pub struct Pattern {
    m: usize,
    peq: Vec<[u64; 4]>,
}

impl Pattern {
    pub fn new(pattern: &[u8]) -> Self {
        let mut peq = vec![[0u64; 4]; pattern.len().div_ceil(W)];
        for (i, &c) in pattern.iter().enumerate() {
            peq[i / W][(c & 3) as usize] |= 1 << (i % W);
        }
        Pattern {
            m: pattern.len(),
            peq,
        }
    }
}

/// Edit distance in `mode`, or `None` if it exceeds `k`.
pub fn distance(pattern: &[u8], text: &[u8], mode: Mode, k: u32) -> Option<u32> {
    distance_with(&Pattern::new(pattern), text, mode, k)
}

/// [`distance`] with the pattern's masks prebuilt.
pub fn distance_with(pat: &Pattern, text: &[u8], mode: Mode, k: u32) -> Option<u32> {
    let m = pat.m;
    let n = text.len();
    if m == 0 {
        let d = if mode == Mode::Global { n as u32 } else { 0 };
        return (d <= k).then_some(d);
    }
    let k = k.min((m + n) as u32) as i32;
    let blocks = m.div_ceil(W);
    let last_bit = (m - 1) % W;
    let rows = |b: usize| if b + 1 == blocks { last_bit + 1 } else { W };
    let bit = |b: usize| if b + 1 == blocks { last_bit } else { W - 1 };
    let peq = &pat.peq;
    // Horizontal delta entering row 0: D[0][j] is 0 (free text start) or j.
    let top: i32 = if mode == Mode::Infix { 0 } else { 1 };
    let mut pv = vec![u64::MAX; blocks];
    let mut mv = vec![0u64; blocks];
    // score[b] = D[bottom row of block b][current column]; column 0 is D[i][0] = i.
    let mut score: Vec<i32> = (0..blocks).map(|b| (b * W + rows(b)) as i32).collect();
    // Active blocks 0..=last: rows that can still be <= k (D[i][0] = i).
    let mut last: isize = (blocks.min((k as usize + 1).div_ceil(W)) as isize) - 1;
    let mut best: i32 = if last as usize + 1 == blocks && score[blocks - 1] <= k {
        score[blocks - 1]
    } else {
        i32::MAX
    };
    for (j, &c) in text.iter().enumerate() {
        let c = (c & 3) as usize;
        let mut hin = top;
        for b in 0..=last as usize {
            hin = block(&mut pv[b], &mut mv[b], peq[b][c], hin, bit(b));
            score[b] += hin;
        }
        let lb = last as usize;
        // Ukkonen: add the next block if its top can still be <= k, else drop exhausted
        // blocks from the bottom.
        if lb + 1 < blocks && score[lb] - hin <= k && ((peq[lb + 1][c] & 1) != 0 || hin < 0) {
            let nb = lb + 1;
            pv[nb] = u64::MAX;
            mv[nb] = 0;
            let h = block(&mut pv[nb], &mut mv[nb], peq[nb][c], hin, bit(nb));
            score[nb] = score[lb] - hin + rows(nb) as i32 + h;
            last += 1;
        } else {
            while last >= 0 && score[last as usize] >= k + W as i32 {
                last -= 1;
            }
            if last < 0 {
                // No row can be <= k from here on (only possible in the anchored modes):
                // Global needs the last column, so it fails; the others keep their best.
                return if mode == Mode::Global {
                    None
                } else {
                    (best <= k).then_some(best as u32)
                };
            }
        }
        if last as usize + 1 == blocks {
            let s = score[blocks - 1];
            match mode {
                Mode::Infix | Mode::Prefix => best = best.min(s),
                Mode::Global => {
                    if j + 1 == n {
                        best = s;
                    }
                }
            }
        } else if mode == Mode::Global && j + 1 == n {
            best = i32::MAX;
        }
    }
    if mode == Mode::Global && n == 0 {
        best = m as i32;
    }
    (best <= k).then_some(best as u32)
}

/// [`distance_with`] in `mode` for each of `texts` with one threshold `k`: the same results
/// as separate calls. With AVX-512, up to eight texts run at once, one per 64-bit lane,
/// sharing the pattern's blocks; with AVX2, up to four; at the NEON level, four in lockstep
/// (`batch_lanes`).
pub fn batch_distances(pat: &Pattern, texts: &[&[u8]], mode: Mode, k: u32) -> Vec<Option<u32>> {
    match simd_level() {
        #[cfg(target_arch = "x86_64")]
        Level::Avx512 if pat.m > 0 => {
            let mut out = Vec::with_capacity(texts.len());
            // SAFETY: the CPU has AVX-512 (the active level).
            unsafe { batch_avx512(pat, texts, mode, k, &mut out) };
            out
        }
        #[cfg(target_arch = "x86_64")]
        Level::Avx2 if pat.m > 0 => {
            let mut out = Vec::with_capacity(texts.len());
            // SAFETY: the CPU has AVX2 (the active level).
            unsafe { batch_avx2(pat, texts, mode, k, &mut out) };
            out
        }
        #[cfg(target_arch = "aarch64")]
        Level::Neon if pat.m > 0 => {
            let mut out = Vec::with_capacity(texts.len());
            for batch in texts.chunks(BATCH_LANES) {
                batch_lanes::<BATCH_LANES>(pat, batch, mode, k, &mut out);
            }
            out
        }
        _ => texts.iter().map(|t| distance_with(pat, t, mode, k)).collect(),
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vl,avx2,fma,bmi1,bmi2,popcnt,lzcnt")]
fn batch_avx512(pat: &Pattern, texts: &[&[u8]], mode: Mode, k: u32, out: &mut Vec<Option<u32>>) {
    for batch in texts.chunks(8) {
        // SAFETY: AVX-512F is enabled here.
        unsafe { batch_x8(pat, batch, mode, k, out) };
    }
}

/// Distances in `mode` of up to eight texts, one per lane, appended to `out`. Each lane runs
/// [`distance_with`]'s recurrence; the band is the union of the lanes' Ukkonen bands (a block
/// is added if any lane would add it, dropped only if every lane would), and blocks are
/// added with the same upper-bound state. Every lane's band thus contains its own, which
/// already holds every cell scoring <= k, so a lane's distance is exact when <= k, as in
/// `distance_with`, and above k otherwise.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
#[inline]
unsafe fn batch_x8(pat: &Pattern, texts: &[&[u8]], mode: Mode, k: u32, out: &mut Vec<Option<u32>>) {
    use std::arch::x86_64::*;
    debug_assert!((1..=8).contains(&texts.len()) && pat.m > 0);
    let m = pat.m;
    let n_max = texts.iter().map(|t| t.len()).max().unwrap();
    let kk = k.min((m + n_max) as u32) as i64;
    let blocks = m.div_ceil(W);
    let last_bit = (m - 1) % W;
    let rows = |b: usize| if b + 1 == blocks { last_bit + 1 } else { W };
    // Text codes interleaved by lane: column j of lane l at 8j + l.
    let mut tx = vec![0u8; 8 * n_max];
    for (l, t) in texts.iter().enumerate() {
        for (j, &c) in t.iter().enumerate() {
            tx[8 * j + l] = c & 3;
        }
    }
    // Lanes still inside their text, per column: lane l is live while j < len(l).
    let lens = _mm512_set_epi64(
        texts.get(7).map_or(0, |t| t.len() as i64),
        texts.get(6).map_or(0, |t| t.len() as i64),
        texts.get(5).map_or(0, |t| t.len() as i64),
        texts.get(4).map_or(0, |t| t.len() as i64),
        texts.get(3).map_or(0, |t| t.len() as i64),
        texts.get(2).map_or(0, |t| t.len() as i64),
        texts.get(1).map_or(0, |t| t.len() as i64),
        texts[0].len() as i64,
    );
    // Each block's four match masks, twice over, so a lane's code selects its mask.
    let peq: Vec<__m512i> = pat
        .peq
        .iter()
        // SAFETY: four u64s.
        .map(|q| _mm512_broadcast_i64x4(unsafe { _mm256_loadu_si256(q.as_ptr().cast()) }))
        .collect();
    let ones = _mm512_set1_epi64(-1);
    let one = _mm512_set1_epi64(1);
    let zero = _mm512_setzero_si512();
    let mut pv = vec![ones; blocks];
    let mut mv = vec![zero; blocks];
    let mut score: Vec<__m512i> = (0..blocks)
        .map(|b| _mm512_set1_epi64((b * W + rows(b)) as i64))
        .collect();
    let kv = _mm512_set1_epi64(kk);
    let kw = _mm512_set1_epi64(kk + W as i64);
    let mut last: isize = (blocks.min((kk as usize + 1).div_ceil(W)) as isize) - 1;
    // Column 0: D[m][0] = m (Global: only for an empty text).
    let mut best = if mode == Mode::Global {
        _mm512_mask_mov_epi64(
            _mm512_set1_epi64(i64::MAX),
            _mm512_cmpeq_epi64_mask(lens, zero),
            _mm512_set1_epi64(m as i64),
        )
    } else if last as usize + 1 == blocks && m as i64 <= kk {
        _mm512_set1_epi64(m as i64)
    } else {
        _mm512_set1_epi64(i64::MAX)
    };
    // Horizontal delta entering row 0 (as -1/+1 bits): D[0][j] is 0 (free text start) or j.
    let h0 = if mode == Mode::Infix {
        (zero, zero)
    } else {
        (zero, one)
    };
    // One block step on all lanes (`distance_with`'s `block`): (h out as +1/-1 bits).
    #[inline]
    #[target_feature(enable = "avx512f")]
    fn step(
        pv: &mut __m512i,
        mv: &mut __m512i,
        eq: __m512i,
        (neg, pos): (__m512i, __m512i),
        bitv: __m128i,
        live: __mmask8,
    ) -> (__m512i, __m512i) {
        let (p, m) = (*pv, *mv);
        let one = _mm512_set1_epi64(1);
        let xv = _mm512_or_si512(eq, m);
        let eq = _mm512_or_si512(eq, neg);
        let xh = _mm512_or_si512(
            _mm512_xor_si512(_mm512_add_epi64(_mm512_and_si512(eq, p), p), p),
            eq,
        );
        let ph = _mm512_or_si512(m, _mm512_ternarylogic_epi64::<0x01>(xh, p, p));
        let mh = _mm512_and_si512(p, xh);
        let hp = _mm512_and_si512(_mm512_srl_epi64(ph, bitv), one);
        let hm = _mm512_and_si512(_mm512_srl_epi64(mh, bitv), one);
        let ph = _mm512_or_si512(_mm512_slli_epi64::<1>(ph), pos);
        let mh = _mm512_or_si512(_mm512_slli_epi64::<1>(mh), neg);
        let npv = _mm512_or_si512(mh, _mm512_ternarylogic_epi64::<0x01>(xv, ph, ph));
        let nmv = _mm512_and_si512(ph, xv);
        *pv = _mm512_mask_mov_epi64(p, live, npv);
        *mv = _mm512_mask_mov_epi64(m, live, nmv);
        (hm, hp)
    }
    let bitw = _mm_cvtsi64_si128((W - 1) as i64);
    let bitl = _mm_cvtsi64_si128(last_bit as i64);
    let bitv = |b: usize| if b + 1 == blocks { bitl } else { bitw };
    for j in 0..n_max {
        let live: __mmask8 = _mm512_cmpgt_epi64_mask(lens, _mm512_set1_epi64(j as i64));
        // SAFETY: tx holds 8 codes per column.
        let c = _mm512_cvtepu8_epi64(unsafe { _mm_loadl_epi64(tx.as_ptr().add(8 * j).cast()) });
        let mut h = h0;
        for b in 0..=last as usize {
            let eq = _mm512_permutexvar_epi64(c, peq[b]);
            h = step(&mut pv[b], &mut mv[b], eq, h, bitv(b), live);
            score[b] = _mm512_mask_add_epi64(score[b], live, score[b], _mm512_sub_epi64(h.1, h.0));
        }
        let lb = last as usize;
        let hin = _mm512_sub_epi64(h.1, h.0);
        let mut added = false;
        if lb + 1 < blocks {
            let eqn = _mm512_permutexvar_epi64(c, peq[lb + 1]);
            let top = _mm512_sub_epi64(score[lb], hin);
            let want = _mm512_mask_cmple_epi64_mask(live, top, kv)
                & (_mm512_test_epi64_mask(eqn, one) | _mm512_cmplt_epi64_mask(hin, zero));
            if want != 0 {
                let nb = lb + 1;
                pv[nb] = ones;
                mv[nb] = zero;
                let hn = step(&mut pv[nb], &mut mv[nb], eqn, h, bitv(nb), live);
                score[nb] = _mm512_add_epi64(
                    _mm512_add_epi64(top, _mm512_set1_epi64(rows(nb) as i64)),
                    _mm512_sub_epi64(hn.1, hn.0),
                );
                last += 1;
                added = true;
            }
        }
        if !added {
            while last >= 0 && _mm512_mask_cmplt_epi64_mask(live, score[last as usize], kw) == 0 {
                last -= 1;
            }
            if last < 0 {
                break;
            }
        }
        if last as usize + 1 == blocks {
            if mode == Mode::Global {
                let end = _mm512_cmpeq_epi64_mask(lens, _mm512_set1_epi64(j as i64 + 1));
                best = _mm512_mask_mov_epi64(best, end, score[blocks - 1]);
            } else {
                best = _mm512_mask_min_epi64(best, live, best, score[blocks - 1]);
            }
        }
    }
    let mut bs = [0i64; 8];
    // SAFETY: eight i64s.
    unsafe { _mm512_storeu_si512(bs.as_mut_ptr().cast(), best) };
    for (l, t) in texts.iter().enumerate() {
        let kl = k.min((m + t.len()) as u32) as i64;
        out.push((bs[l] <= kl).then_some(bs[l] as u32));
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma,bmi1,bmi2,popcnt,lzcnt")]
fn batch_avx2(pat: &Pattern, texts: &[&[u8]], mode: Mode, k: u32, out: &mut Vec<Option<u32>>) {
    for batch in texts.chunks(4) {
        // SAFETY: AVX2 is enabled here.
        unsafe { batch_x4(pat, batch, mode, k, out) };
    }
}

/// [`batch_x8`] for up to four texts in AVX2, one per 64-bit lane: the same recurrence and
/// band union, with blends for the per-lane masked updates and vector compares for the mask
/// tests. A lane's code selects its block mask through a dword permute.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn batch_x4(pat: &Pattern, texts: &[&[u8]], mode: Mode, k: u32, out: &mut Vec<Option<u32>>) {
    use std::arch::x86_64::*;
    debug_assert!((1..=4).contains(&texts.len()) && pat.m > 0);
    let m = pat.m;
    let n_max = texts.iter().map(|t| t.len()).max().unwrap();
    let kk = k.min((m + n_max) as u32) as i64;
    let blocks = m.div_ceil(W);
    let last_bit = (m - 1) % W;
    let rows = |b: usize| if b + 1 == blocks { last_bit + 1 } else { W };
    // Column j of lane l at 4j + l: the dword indices (2c, 2c + 1) of the lane's code c in a
    // block's four masks (code 0 past the lane's text, whose lane is not updated there).
    let mut tx = vec![1u64 << 32; 4 * n_max];
    for (l, t) in texts.iter().enumerate() {
        for (j, &c) in t.iter().enumerate() {
            let c = (c & 3) as u64;
            tx[4 * j + l] = (2 * c) | ((2 * c + 1) << 32);
        }
    }
    let len = |l: usize| texts.get(l).map_or(0, |t| t.len() as i64);
    // Lanes still inside their text, per column: lane l is live while j < len(l).
    let lens = _mm256_set_epi64x(len(3), len(2), len(1), len(0));
    // SAFETY: four u64s.
    let peq: Vec<__m256i> = pat
        .peq
        .iter()
        .map(|q| unsafe { _mm256_loadu_si256(q.as_ptr().cast()) })
        .collect();
    let ones = _mm256_set1_epi64x(-1);
    let one = _mm256_set1_epi64x(1);
    let zero = _mm256_setzero_si256();
    let mut pv = vec![ones; blocks];
    let mut mv = vec![zero; blocks];
    let mut score: Vec<__m256i> = (0..blocks)
        .map(|b| _mm256_set1_epi64x((b * W + rows(b)) as i64))
        .collect();
    let kv = _mm256_set1_epi64x(kk);
    let kw = _mm256_set1_epi64x(kk + W as i64);
    let mut last: isize = (blocks.min((kk as usize + 1).div_ceil(W)) as isize) - 1;
    // Column 0: D[m][0] = m (Global: only for an empty text).
    let mut best = if mode == Mode::Global {
        _mm256_blendv_epi8(
            _mm256_set1_epi64x(i64::MAX),
            _mm256_set1_epi64x(m as i64),
            _mm256_cmpeq_epi64(lens, zero),
        )
    } else if last as usize + 1 == blocks && m as i64 <= kk {
        _mm256_set1_epi64x(m as i64)
    } else {
        _mm256_set1_epi64x(i64::MAX)
    };
    // Horizontal delta entering row 0 (as -1/+1 bits): D[0][j] is 0 (free text start) or j.
    let h0 = if mode == Mode::Infix {
        (zero, zero)
    } else {
        (zero, one)
    };
    // One block step on all lanes (`distance_with`'s `block`): (h out as +1/-1 bits).
    #[inline]
    #[target_feature(enable = "avx2")]
    fn step(
        pv: &mut __m256i,
        mv: &mut __m256i,
        eq: __m256i,
        (neg, pos): (__m256i, __m256i),
        bitv: __m128i,
        live: __m256i,
    ) -> (__m256i, __m256i) {
        let (p, m) = (*pv, *mv);
        let one = _mm256_set1_epi64x(1);
        let ones = _mm256_set1_epi64x(-1);
        let xv = _mm256_or_si256(eq, m);
        let eq = _mm256_or_si256(eq, neg);
        let xh = _mm256_or_si256(
            _mm256_xor_si256(_mm256_add_epi64(_mm256_and_si256(eq, p), p), p),
            eq,
        );
        let ph = _mm256_or_si256(m, _mm256_xor_si256(_mm256_or_si256(xh, p), ones));
        let mh = _mm256_and_si256(p, xh);
        let hp = _mm256_and_si256(_mm256_srl_epi64(ph, bitv), one);
        let hm = _mm256_and_si256(_mm256_srl_epi64(mh, bitv), one);
        let ph = _mm256_or_si256(_mm256_slli_epi64::<1>(ph), pos);
        let mh = _mm256_or_si256(_mm256_slli_epi64::<1>(mh), neg);
        let npv = _mm256_or_si256(mh, _mm256_xor_si256(_mm256_or_si256(xv, ph), ones));
        let nmv = _mm256_and_si256(ph, xv);
        *pv = _mm256_blendv_epi8(p, npv, live);
        *mv = _mm256_blendv_epi8(m, nmv, live);
        (hm, hp)
    }
    let bitw = _mm_cvtsi64_si128((W - 1) as i64);
    let bitl = _mm_cvtsi64_si128(last_bit as i64);
    let bitv = |b: usize| if b + 1 == blocks { bitl } else { bitw };
    for j in 0..n_max {
        let live = _mm256_cmpgt_epi64(lens, _mm256_set1_epi64x(j as i64));
        // SAFETY: tx holds 4 indices per column.
        let idx = unsafe { _mm256_loadu_si256(tx.as_ptr().add(4 * j).cast()) };
        let mut h = h0;
        for b in 0..=last as usize {
            let eq = _mm256_permutevar8x32_epi32(peq[b], idx);
            h = step(&mut pv[b], &mut mv[b], eq, h, bitv(b), live);
            score[b] = _mm256_add_epi64(score[b], _mm256_and_si256(_mm256_sub_epi64(h.1, h.0), live));
        }
        let lb = last as usize;
        let hin = _mm256_sub_epi64(h.1, h.0);
        let mut added = false;
        if lb + 1 < blocks {
            let eqn = _mm256_permutevar8x32_epi32(peq[lb + 1], idx);
            let top = _mm256_sub_epi64(score[lb], hin);
            // Live lanes with top <= k whose next block's top row can still be <= k.
            let want = _mm256_and_si256(
                _mm256_andnot_si256(_mm256_cmpgt_epi64(top, kv), live),
                _mm256_or_si256(
                    _mm256_cmpeq_epi64(_mm256_and_si256(eqn, one), one),
                    _mm256_cmpgt_epi64(zero, hin),
                ),
            );
            if _mm256_testz_si256(want, want) == 0 {
                let nb = lb + 1;
                pv[nb] = ones;
                mv[nb] = zero;
                let hn = step(&mut pv[nb], &mut mv[nb], eqn, h, bitv(nb), live);
                score[nb] = _mm256_add_epi64(
                    _mm256_add_epi64(top, _mm256_set1_epi64x(rows(nb) as i64)),
                    _mm256_sub_epi64(hn.1, hn.0),
                );
                last += 1;
                added = true;
            }
        }
        if !added {
            // Drop the bottom block while no live lane scores below k + W in it.
            while last >= 0 && {
                let below = _mm256_and_si256(live, _mm256_cmpgt_epi64(kw, score[last as usize]));
                _mm256_testz_si256(below, below) != 0
            } {
                last -= 1;
            }
            if last < 0 {
                break;
            }
        }
        if last as usize + 1 == blocks {
            let s = score[blocks - 1];
            let lt = if mode == Mode::Global {
                _mm256_cmpeq_epi64(lens, _mm256_set1_epi64x(j as i64 + 1))
            } else {
                _mm256_and_si256(live, _mm256_cmpgt_epi64(best, s))
            };
            best = _mm256_blendv_epi8(best, s, lt);
        }
    }
    let mut bs = [0i64; 4];
    // SAFETY: four i64s.
    unsafe { _mm256_storeu_si256(bs.as_mut_ptr().cast(), best) };
    for (l, t) in texts.iter().enumerate() {
        let kl = k.min((m + t.len()) as u32) as i64;
        out.push((bs[l] <= kl).then_some(bs[l] as u32));
    }
}

/// Lanes of [`batch_lanes`] on aarch64.
#[cfg(target_arch = "aarch64")]
const BATCH_LANES: usize = 4;

/// [`batch_x8`] in scalar code, `L` texts in lockstep: the lanes' block steps are
/// independent, so a wide core overlaps their carried horizontal deltas, the critical path
/// of one [`distance_with`] call. Same band rules as `batch_x8`, so the same results. A lane
/// past its text's end keeps stepping (on code 0) but no longer counts in the band or its
/// best, so it changes no output.
#[cfg_attr(not(target_arch = "aarch64"), allow(dead_code))]
fn batch_lanes<const L: usize>(
    pat: &Pattern,
    texts: &[&[u8]],
    mode: Mode,
    k: u32,
    out: &mut Vec<Option<u32>>,
) {
    debug_assert!((1..=L).contains(&texts.len()) && pat.m > 0);
    let m = pat.m;
    let n_max = texts.iter().map(|t| t.len()).max().unwrap();
    let kk = k.min((m + n_max) as u32) as i32;
    let blocks = m.div_ceil(W);
    let last_bit = (m - 1) % W;
    let rows = |b: usize| if b + 1 == blocks { last_bit + 1 } else { W };
    let bit = |b: usize| if b + 1 == blocks { last_bit } else { W - 1 };
    // Text codes interleaved by lane: column j of lane l at L * j + l.
    let mut tx = vec![0u8; L * n_max];
    let mut lens = [0usize; L];
    for (l, t) in texts.iter().enumerate() {
        lens[l] = t.len();
        for (j, &c) in t.iter().enumerate() {
            tx[L * j + l] = c & 3;
        }
    }
    let peq = &pat.peq;
    let mut pv = vec![[u64::MAX; L]; blocks];
    let mut mv = vec![[0u64; L]; blocks];
    let mut score: Vec<[i32; L]> = (0..blocks).map(|b| [(b * W + rows(b)) as i32; L]).collect();
    let kw = kk + W as i32;
    let mut last: isize = (blocks.min((kk as usize + 1).div_ceil(W)) as isize) - 1;
    // Column 0: D[m][0] = m (Global: only for an empty text).
    let best0 = if last as usize + 1 == blocks && m as i32 <= kk {
        m as i32
    } else {
        i32::MAX
    };
    let mut best: [i32; L] = std::array::from_fn(|l| match mode {
        Mode::Global if lens[l] == 0 => m as i32,
        Mode::Global => i32::MAX,
        _ => best0,
    });
    // Horizontal delta entering row 0: D[0][j] is 0 (free text start) or j.
    let h0: i32 = if mode == Mode::Infix { 0 } else { 1 };
    for (j, c) in tx.as_chunks::<L>().0.iter().enumerate() {
        let c: [usize; L] = std::array::from_fn(|l| c[l] as usize);
        let live: [bool; L] = std::array::from_fn(|l| j < lens[l]);
        let mut h = [h0; L];
        for b in 0..=last as usize {
            let (p, mm, s, q, bt) = (&mut pv[b], &mut mv[b], &mut score[b], &peq[b], bit(b));
            for l in 0..L {
                h[l] = block(&mut p[l], &mut mm[l], q[c[l]], h[l], bt);
                s[l] += h[l];
            }
        }
        let lb = last as usize;
        let mut added = false;
        if lb + 1 < blocks {
            let nb = lb + 1;
            let top: [i32; L] = std::array::from_fn(|l| score[lb][l] - h[l]);
            let want = (0..L).any(|l| live[l] && top[l] <= kk && ((peq[nb][c[l]] & 1) != 0 || h[l] < 0));
            if want {
                pv[nb] = [u64::MAX; L];
                mv[nb] = [0; L];
                for l in 0..L {
                    let hn = block(&mut pv[nb][l], &mut mv[nb][l], peq[nb][c[l]], h[l], bit(nb));
                    score[nb][l] = top[l] + rows(nb) as i32 + hn;
                }
                last += 1;
                added = true;
            }
        }
        if !added {
            while last >= 0 && !(0..L).any(|l| live[l] && score[last as usize][l] < kw) {
                last -= 1;
            }
            if last < 0 {
                break;
            }
        }
        if last as usize + 1 == blocks {
            for l in 0..L {
                if mode == Mode::Global {
                    if j + 1 == lens[l] {
                        best[l] = score[blocks - 1][l];
                    }
                } else if live[l] {
                    best[l] = best[l].min(score[blocks - 1][l]);
                }
            }
        }
    }
    for (l, t) in texts.iter().enumerate() {
        let kl = k.min((m + t.len()) as u32) as i32;
        out.push((best[l] <= kl).then_some(best[l] as u32));
    }
}

/// Minimum edit distance between `pattern` and any substring of `text` (no cutoff).
pub fn infix_distance(pattern: &[u8], text: &[u8]) -> u32 {
    distance(pattern, text, Mode::Infix, u32::MAX).unwrap()
}

/// Affine-gap scoring (maximized): `+mat` per match, `-mis` per mismatch, `-(open + ext·l)`
/// per gap of length `l`. Defaults are minimap2's map-ont single-affine part.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Affine {
    pub mat: i32,
    pub mis: i32,
    pub open: i32,
    pub ext: i32,
}

impl Default for Affine {
    fn default() -> Self {
        Self::ONT
    }
}

impl Affine {
    /// minimap2 map-ont (single-affine part of -A2 -B4 -O4,24 -E2,1).
    pub const ONT: Affine = Affine {
        mat: 2,
        mis: 4,
        open: 4,
        ext: 2,
    };
    /// minimap2 map-hifi (single-affine part of -A1 -B4 -O6,26 -E2,1).
    pub const HIFI: Affine = Affine {
        mat: 1,
        mis: 4,
        open: 6,
        ext: 2,
    };
}

const NEG: i32 = i32::MIN / 4;

/// Best affine score aligning all of `pattern` against `text` (Gotoh). `Global`: all of the
/// text too; `Prefix`: a prefix of the text (text end free); `Infix`: any substring. Cells
/// more than `band` off the diagonals spanning the two lengths are not computed, so the
/// result is a lower bound (pass `usize::MAX` for the exact score).
pub fn affine(pattern: &[u8], text: &[u8], mode: Mode, sc: &Affine, band: usize) -> i32 {
    affine_dp(pattern, text, mode, sc, band, None, text.len()).score
}

/// Best affine score of an extension from the start of both sequences, with both far ends
/// free: the maximum over all DP cells (as minimap2's extension, without z-drop). Read tails
/// that don't align (adapters, junk) are clipped instead of forced through the alignment.
///
/// With `xdrop`, the extension stops once a whole row of the band scores more than `xdrop`
/// below the best cell so far (minimap2-style X-drop): the rest of the read is junk for
/// this text. Pass `i32::MAX` for none.
pub fn affine_extend(
    pattern: &[u8],
    text: &[u8],
    sc: &Affine,
    band: usize,
    end_bonus: i32,
    xdrop: i32,
) -> i32 {
    // Reaching the pattern's end earns `end_bonus` (minimap2's --end-bonus).
    let d = affine_dp(
        pattern,
        text,
        Mode::Prefix,
        sc,
        band,
        Some((xdrop, false)),
        text.len(),
    );
    d.clip_best.max(d.score.saturating_add(end_bonus))
}

/// Best affine local score of `read` in `text` (2-bit codes; a read base of 4 or more matches
/// nothing) over the diagonals `dlo..=dhi` (text position minus read position), both read
/// ends free to clip: the best global alignment of a piece of the read with a piece of the
/// text, plus `end_bonus` for each read end the piece reaches (0 if nothing scores above it).
pub fn local_band(read: &[u8], text: &[u8], sc: &Affine, dlo: i64, dhi: i64, end_bonus: i32) -> i32 {
    debug_assert!(text.iter().all(|&b| b < 4));
    match simd_level() {
        // SAFETY: the CPU has the active level's features.
        #[cfg(target_arch = "x86_64")]
        Level::Avx512 => unsafe { local_band_avx512(read, text, sc, dlo, dhi, end_bonus) },
        #[cfg(target_arch = "x86_64")]
        Level::Avx2 => unsafe { local_band_avx2(read, text, sc, dlo, dhi, end_bonus) },
        #[cfg(target_arch = "aarch64")]
        Level::Neon => local_band_core::<NeonRows>(read, text, sc, dlo, dhi, end_bonus),
        _ => local_band_core::<ScalarRows>(read, text, sc, dlo, dhi, end_bonus),
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vl,avx2,fma,bmi1,bmi2,popcnt,lzcnt")]
fn local_band_avx512(read: &[u8], text: &[u8], sc: &Affine, dlo: i64, dhi: i64, end_bonus: i32) -> i32 {
    local_band_core::<Avx512Rows>(read, text, sc, dlo, dhi, end_bonus)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma,bmi1,bmi2,popcnt,lzcnt")]
fn local_band_avx2(read: &[u8], text: &[u8], sc: &Affine, dlo: i64, dhi: i64, end_bonus: i32) -> i32 {
    local_band_core::<Avx2Rows>(read, text, sc, dlo, dhi, end_bonus)
}

/// [`local_band`] by the row kernel `R`, whose features the caller is compiled for.
#[inline(always)]
fn local_band_core<R: Rows>(
    read: &[u8],
    text: &[u8],
    sc: &Affine,
    dlo: i64,
    dhi: i64,
    end_bonus: i32,
) -> i32 {
    let (m, n) = (read.len(), text.len());
    if m == 0 || dlo > dhi {
        return 0;
    }
    let (go, ext) = (sc.open + sc.ext, sc.ext);
    // Columns 0..=n and the band's right sentinel.
    let w = n + 2;
    let (mut hp, mut fp, mut h, mut f) = (vec![NEG; w], vec![NEG; w], vec![NEG; w], vec![NEG; w]);
    let (mut hpv, mut gsv) = (Vec::new(), Vec::new());
    // Row `i`'s columns (first, one past the last): from one row to the next, each end stays
    // or moves one column right, so a row reads only what the row above wrote (its cells and
    // its sentinel).
    let rng = |i: usize| {
        let lo = (i as i64 + dlo).clamp(0, n as i64 + 1) as usize;
        let hi = (i as i64 + dhi + 1).clamp(0, n as i64 + 1) as usize;
        (lo, hi.max(lo))
    };
    let (lo, hi) = rng(0);
    hp[lo..hi].fill(end_bonus);
    let mut best = 0;
    for i in 1..=m {
        let (lo, hi) = rng(i);
        // An alignment may start at any cell, but not in the last row: it would hold nothing
        // of the read.
        let floor = if i < m { 0 } else { NEG };
        let (mut start, mut e, mut rowmax) = (lo, NEG, NEG);
        if lo == 0 && hi > 0 {
            // Column 0: only a vertical gap reaches it.
            let fj = (hp[0] - go).max(fp[0] - ext).max(NEG);
            let v = fj.max(floor);
            (h[0], f[0]) = (v, fj);
            (start, e, rowmax) = (1, v - go, v);
        }
        if start < hi {
            // SAFETY: the entry point is compiled for the kernel's features, which the CPU has.
            let row = unsafe {
                R::row(
                    &hp[start - 1..hi],
                    &fp[start..hi],
                    &text[start - 1..hi - 1],
                    &mut h[start..=hi],
                    &mut f[start..=hi],
                    read[i - 1],
                    (sc.mat, sc.mis, go, ext),
                    e,
                    floor,
                    (&mut hpv, &mut gsv),
                )
            };
            rowmax = rowmax.max(row);
        } else {
            (h[hi], f[hi]) = (NEG, NEG);
        }
        best = best.max(rowmax);
        if i == m {
            best = best.max(rowmax.saturating_add(end_bonus));
        }
        std::mem::swap(&mut hp, &mut h);
        std::mem::swap(&mut fp, &mut f);
    }
    best
}

/// [`affine`] for each of `texts`, which must be sorted: equal to separate calls, but rows
/// that read only a prefix shared with the previous text are not recomputed.
pub fn affine_sorted(pattern: &[u8], texts: &[&[u8]], mode: Mode, sc: &Affine, band: usize) -> Vec<i32> {
    affine_dp_sorted(pattern, texts, mode, sc, band, None, |n| n, usize::MAX)
        .into_iter()
        .map(|d| d.score)
        .collect()
}

/// [`affine_extend`] for each of `texts`, which must be sorted, sharing rows as
/// [`affine_sorted`] does, with each text's score also over its first `cut` bases only:
/// (whole text, text cut at `cut`, text bases the whole text's best extension covers), from
/// one pass (`cut = usize::MAX` for none).
///
/// With `diagonal`, the band stays on the diagonal from the start (`2 * band + 1` cells per
/// row) however much longer the text is: for accurate reads, whose extensions stay near it. A
/// text shorter than the pattern (a transcript ending inside the read's flank) keeps the band
/// of [`affine_extend`], which can place the read's extra bases. Without `diagonal`, the band is
/// laid out for a text of at least `window` bases, as [`affine_extend_window`] does: a text cut
/// short where its transcript ends gets the band of the full window, so a read whose far end
/// doesn't align (a primer, an adapter) is placed as it would be against a longer transcript.
///
/// Without X-drop, the cut score equals the pass over the cut text (same band layout: `window`
/// at least as long as the texts); with it, the cut score may extend further, since the pass
/// stops only when the whole text's does.
#[allow(clippy::too_many_arguments)]
pub fn affine_extend_sorted_cut(
    pattern: &[u8],
    texts: &[&[u8]],
    sc: &Affine,
    band: usize,
    end_bonus: i32,
    xdrop: i32,
    diagonal: bool,
    window: usize,
    cut: usize,
) -> Vec<(i32, i32, usize)> {
    let m = pattern.len();
    let geom = |n: usize| if diagonal { m.min(n) } else { n.max(window) };
    affine_dp_sorted(
        pattern,
        texts,
        Mode::Prefix,
        sc,
        band,
        Some((xdrop, cut != usize::MAX)),
        geom,
        cut,
    )
    .into_iter()
    .map(|d| {
        let end = d.score.saturating_add(end_bonus);
        let whole = d.clip_best.max(end);
        let cols = if d.clip_best >= end {
            d.clip_cell.1
        } else {
            d.score_col
        };
        let cut = if cut == usize::MAX {
            whole
        } else {
            d.cut_best.max(d.cut_score.saturating_add(end_bonus))
        };
        (whole, cut, cols)
    })
    .collect()
}

/// [`affine_extend`] with where the best extension ends, (score, pattern bases, text bases,
/// text bases read), on a text cut from a window of `window` bases (`text.len() <= window`;
/// `text.len()` for a plain extension): the band is laid out for the whole window, so the
/// result depends on the text only
/// through the bases read and, if it read them all, the text's length. Texts cut from windows
/// of one length share results on that basis.
pub fn affine_extend_window(
    pattern: &[u8],
    text: &[u8],
    window: usize,
    sc: &Affine,
    band: usize,
    end_bonus: i32,
    xdrop: i32,
) -> (i32, usize, usize, usize) {
    // Reaching the pattern's end earns `end_bonus` (minimap2's --end-bonus).
    let d = affine_dp(
        pattern,
        text,
        Mode::Prefix,
        sc,
        band,
        Some((xdrop, true)),
        window.max(text.len()),
    );
    let whole = d.score.saturating_add(end_bonus);
    if d.clip_best >= whole {
        (d.clip_best, d.clip_cell.0, d.clip_cell.1, d.text_read)
    } else {
        (whole, pattern.len(), d.score_col, d.text_read)
    }
}

/// Exclusive prefix maxima of `xs` in place, starting from `init`: `xs[k]` becomes the
/// maximum of `init` and `xs[..k]` (the scalar rows' horizontal-gap pass; the vector row
/// kernels scan in registers instead).
#[inline(always)]
fn excl_prefix_max(xs: &mut [i32], init: i32) {
    let mut g = init;
    for x in xs.iter_mut() {
        let next = g.max(*x);
        *x = g;
        g = next;
    }
}

/// A DP row kernel for `run_core`: [`row_scalar`], or a fused vector row.
trait Rows {
    /// Row `i` from row `i - 1`, as [`row_scalar`] computes it (the arguments are its).
    ///
    /// # Safety
    /// The CPU has the kernel's features.
    #[allow(clippy::too_many_arguments)]
    unsafe fn row(
        hprev: &[i32],
        fu: &[i32],
        tx: &[u8],
        nh: &mut [i32],
        nf: &mut [i32],
        c: u8,
        sc: (i32, i32, i32, i32),
        e: i32,
        floor: i32,
        scratch: (&mut Vec<i32>, &mut Vec<i32>),
    ) -> i32;
}

/// One DP row over `w` cells by three passes over the row (the first and last autovectorize
/// for the entry point's level): `hprev` is H[i-1][start-1..=hi] (`w + 1` cells), `fu`
/// F[i-1][start..=hi], `tx` the text bases of columns start..=hi, `c` the pattern base, `e`
/// the horizontal gap entering column `start`, `floor` (NEG or above) the least a cell's H is
/// (NEG: none; 0: an alignment may start at the cell); writes H and F of row i to `nh` and `nf`
/// (`w + 1` cells: the last gets NEG, the band's right sentinel) and returns the row's maximum.
/// `scratch` holds two rows of temporaries.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn row_scalar(
    hprev: &[i32],
    fu: &[i32],
    tx: &[u8],
    nh: &mut [i32],
    nf: &mut [i32],
    c: u8,
    (mat, mis, go, ext): (i32, i32, i32, i32),
    e: i32,
    floor: i32,
    (hpv, gsv): (&mut Vec<i32>, &mut Vec<i32>),
) -> i32 {
    let w = fu.len();
    if hpv.len() < w {
        hpv.resize(w, NEG);
        gsv.resize(w, NEG);
    }
    let (hp, gs) = (&mut hpv[..w], &mut gsv[..w]);
    let hd = &hprev[..w]; // H[i-1][j-1]
    let hu = &hprev[1..=w]; // H[i-1][j]
    // The band's right sentinel.
    (nh[w], nf[w]) = (NEG, NEG);
    let (tx, nf, nh) = (&tx[..w], &mut nf[..w], &mut nh[..w]);
    // Pass 1: vertical gaps and diagonal moves (H' = H without horizontal gaps), from the
    // previous row only.
    for k in 0..w {
        let fj = (hu[k] - go).max(fu[k] - ext).max(NEG);
        let s = if c == tx[k] & 3 { mat } else { -mis };
        let x = (hd[k] + s).max(fj).max(floor);
        hp[k] = x;
        gs[k] = x - go + ext * (k as i32 + 1);
        nf[k] = fj;
    }
    // Horizontal gaps. With open >= 0, E[k] = max(H'[k-1] - go, E[k-1] - ext): opening from
    // H = max(H', E) adds nothing. With g[k] = E[k] + ext*k this is a running max (pass 2:
    // exclusive prefix maxima, one op per cell on the critical path); pass 3 takes
    // H = max(H', E).
    excl_prefix_max(gs, e);
    let mut m = NEG;
    for k in 0..w {
        let v = hp[k].max(gs[k] - ext * k as i32);
        nh[k] = v;
        m = m.max(v);
    }
    m
}

/// The scalar passes.
struct ScalarRows;

impl Rows for ScalarRows {
    #[inline(always)]
    unsafe fn row(
        hprev: &[i32],
        fu: &[i32],
        tx: &[u8],
        nh: &mut [i32],
        nf: &mut [i32],
        c: u8,
        sc: (i32, i32, i32, i32),
        e: i32,
        floor: i32,
        scratch: (&mut Vec<i32>, &mut Vec<i32>),
    ) -> i32 {
        row_scalar(hprev, fu, tx, nh, nf, c, sc, e, floor, scratch)
    }
}

/// Cells `k0..w` of a row one at a time: [`row_scalar`]'s steps fused per cell (no row
/// buffers), with `carry` the horizontal-gap maximum entering cell `k0` and `m` the row's
/// maximum so far, which it returns updated. The vector kernels' short rows and tails.
#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn cells_scalar(
    hprev: &[i32],
    fu: &[i32],
    tx: &[u8],
    nh: &mut [i32],
    nf: &mut [i32],
    c: u8,
    (mat, mis, go, ext): (i32, i32, i32, i32),
    floor: i32,
    k0: usize,
    mut carry: i32,
    mut m: i32,
) -> i32 {
    for k in k0..fu.len() {
        let fj = (hprev[k + 1] - go).max(fu[k] - ext).max(NEG);
        let s = if c == tx[k] & 3 { mat } else { -mis };
        let x = (hprev[k] + s).max(fj).max(floor);
        let v = x.max(carry - ext * k as i32);
        carry = carry.max(x - go + ext * (k as i32 + 1));
        nh[k] = v;
        nf[k] = fj;
        m = m.max(v);
    }
    m
}

#[cfg(target_arch = "x86_64")]
struct Avx512Rows;

#[cfg(target_arch = "x86_64")]
impl Rows for Avx512Rows {
    #[inline(always)]
    unsafe fn row(
        hprev: &[i32],
        fu: &[i32],
        tx: &[u8],
        nh: &mut [i32],
        nf: &mut [i32],
        c: u8,
        sc: (i32, i32, i32, i32),
        e: i32,
        floor: i32,
        _: (&mut Vec<i32>, &mut Vec<i32>),
    ) -> i32 {
        // SAFETY: as the caller guarantees.
        unsafe { row_avx512(hprev, fu, tx, nh, nf, c, sc, e, floor) }
    }
}

/// [`row_scalar`] with its three passes fused per 16-lane block: diagonal and vertical moves,
/// then horizontal gaps as an in-register prefix max carried across blocks, with no
/// intermediate row written. Same integer operations as the scalar passes, so the same results.
/// As in [`row_avx2`], whole vectors only: a row's last `w % 16` cells are covered by one more
/// vector ending at the row's end (rows under 16 cells take one masked block).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vl")]
#[inline]
#[allow(clippy::too_many_arguments)]
unsafe fn row_avx512(
    hprev: &[i32],
    fu: &[i32],
    tx: &[u8],
    nh: &mut [i32],
    nf: &mut [i32],
    c: u8,
    (mat, mis, go, ext): (i32, i32, i32, i32),
    e: i32,
    floor: i32,
) -> i32 {
    use std::arch::x86_64::*;
    let w = fu.len();
    assert!(hprev.len() == w + 1 && tx.len() == w && nh.len() == w + 1 && nf.len() == w + 1);
    (nh[w], nf[w]) = (NEG, NEG);
    let (neg, vfloor) = (_mm512_set1_epi32(NEG), _mm512_set1_epi32(floor));
    let (vgo, vext) = (_mm512_set1_epi32(go), _mm512_set1_epi32(ext));
    let (vmat, vmis) = (_mm512_set1_epi32(mat), _mm512_set1_epi32(-mis));
    let vc = _mm512_set1_epi32(c as i32);
    let three = _mm512_set1_epi32(3);
    let last = _mm512_set1_epi32(15);
    let lane = _mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
    let lext = _mm512_mullo_epi32(lane, vext);
    let mut mx = neg;
    // Cells k..k + 16 (the lanes in `mask`) with `cv` the maximum entering cell k (broadcast)
    // and `ek` ext times each lane's cell index: stores H and F and returns the exclusive
    // maxima entering each cell and the maximum after the last.
    let mut cells = |k: usize, mask: __mmask16, cv: __m512i, ek: __m512i| {
        // SAFETY: lanes in `mask` are in bounds of every slice (all 16 unless the row is
        // shorter); the rest are masked off.
        let (hd, hu, fv, t) = unsafe {
            if mask == !0 {
                (
                    _mm512_loadu_si512(hprev.as_ptr().add(k).cast()),
                    _mm512_loadu_si512(hprev.as_ptr().add(k + 1).cast()),
                    _mm512_loadu_si512(fu.as_ptr().add(k).cast()),
                    _mm_loadu_si128(tx.as_ptr().add(k).cast()),
                )
            } else {
                (
                    _mm512_maskz_loadu_epi32(mask, hprev.as_ptr().add(k)),
                    _mm512_maskz_loadu_epi32(mask, hprev.as_ptr().add(k + 1)),
                    _mm512_maskz_loadu_epi32(mask, fu.as_ptr().add(k)),
                    _mm_maskz_loadu_epi8(mask, tx.as_ptr().add(k).cast()),
                )
            }
        };
        let fj = _mm512_max_epi32(
            _mm512_max_epi32(_mm512_sub_epi32(hu, vgo), _mm512_sub_epi32(fv, vext)),
            neg,
        );
        let t = _mm512_and_si512(_mm512_cvtepu8_epi32(t), three);
        let s = _mm512_mask_blend_epi32(_mm512_cmpeq_epi32_mask(t, vc), vmis, vmat);
        let x = _mm512_max_epi32(_mm512_max_epi32(_mm512_add_epi32(hd, s), fj), vfloor);
        // gs[k] = x - go + ext * (k + 1); lanes past the row are NEG (they come after every
        // real lane, so they change no output).
        let g = _mm512_mask_mov_epi32(
            neg,
            mask,
            _mm512_add_epi32(_mm512_sub_epi32(x, vgo), _mm512_add_epi32(ek, vext)),
        );
        // Inclusive prefix max within the block, then the carry; exclusive by one lane.
        let mut v = _mm512_max_epi32(g, _mm512_alignr_epi32::<15>(g, neg));
        v = _mm512_max_epi32(v, _mm512_alignr_epi32::<14>(v, neg));
        v = _mm512_max_epi32(v, _mm512_alignr_epi32::<12>(v, neg));
        v = _mm512_max_epi32(v, _mm512_alignr_epi32::<8>(v, neg));
        v = _mm512_max_epi32(v, cv);
        let excl = _mm512_alignr_epi32::<15>(v, cv);
        let hv = _mm512_max_epi32(x, _mm512_sub_epi32(excl, ek));
        // SAFETY: as for the loads.
        unsafe {
            if mask == !0 {
                _mm512_storeu_si512(nh.as_mut_ptr().add(k).cast(), hv);
                _mm512_storeu_si512(nf.as_mut_ptr().add(k).cast(), fj);
            } else {
                _mm512_mask_storeu_epi32(nh.as_mut_ptr().add(k), mask, hv);
                _mm512_mask_storeu_epi32(nf.as_mut_ptr().add(k), mask, fj);
            }
        }
        mx = _mm512_mask_max_epi32(mx, mask, mx, hv);
        (excl, _mm512_permutexvar_epi32(last, v))
    };
    let mut cv = _mm512_set1_epi32(e);
    if w < 16 {
        cells(0, (1u16 << w) - 1, cv, lext);
        return _mm512_reduce_max_epi32(mx);
    }
    let mut excl = cv;
    let mut ek = lext;
    let step = _mm512_set1_epi32(16 * ext);
    let mut k = 0;
    while k + 16 <= w {
        (excl, cv) = cells(k, !0, cv, ek);
        ek = _mm512_add_epi32(ek, step);
        k += 16;
    }
    if k < w {
        // The last vector, at w - 16: lane w - k of the previous one is the maximum entering
        // its first cell.
        let cv = _mm512_permutexvar_epi32(_mm512_set1_epi32((w - k) as i32), excl);
        cells(
            w - 16,
            !0,
            cv,
            _mm512_add_epi32(lext, _mm512_set1_epi32(ext * (w - 16) as i32)),
        );
    }
    _mm512_reduce_max_epi32(mx)
}

/// A cache line of 16-bit cells: the 16-bit rows are stored and loaded as whole lines (see
/// [`run16_avx512`]).
#[derive(Clone, Copy)]
#[repr(C, align(64))]
struct Block16([i16; 32]);

/// Lane indices for [`run16_avx512`]: `LANE16` counts its 32 lanes, and `DOWN16` shifts them by
/// one (lane 31 of a second vector, then lanes ..31 of the first).
#[cfg(target_arch = "x86_64")]
static LANE16: [i16; 32] = {
    let mut a = [0i16; 32];
    let mut i = 0;
    while i < 32 {
        a[i] = i as i16;
        i += 1;
    }
    a
};
#[cfg(target_arch = "x86_64")]
static DOWN16: [i16; 32] = {
    let mut a = LANE16;
    let mut i = 1;
    while i < 32 {
        a[i] -= 1;
        i += 1;
    }
    a[0] = 63;
    a
};

/// For each 16-bit lane of a 128-bit lane, the bytes of its last one.
#[cfg(target_arch = "x86_64")]
static LAST16: [u8; 64] = {
    let mut a = [14u8; 64];
    let mut i = 1;
    while i < 64 {
        a[i] = 15;
        i += 2;
    }
    a
};

/// The maximum of 32 16-bit lanes.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vl")]
#[inline]
fn max16_avx512(v: std::arch::x86_64::__m512i) -> i16 {
    use std::arch::x86_64::*;
    let a = _mm256_max_epi16(_mm512_castsi512_si256(v), _mm512_extracti64x4_epi64::<1>(v));
    let b = _mm_max_epi16(_mm256_castsi256_si128(a), _mm256_extracti128_si256::<1>(a));
    // The unsigned minimum of eight lanes is one instruction: x ^ 0x7fff turns the signed
    // order around into the unsigned one.
    let flip = _mm_set1_epi16(0x7fff);
    (_mm_cvtsi128_si32(_mm_minpos_epu16(_mm_xor_si128(b, flip))) as i16) ^ 0x7fff
}

/// The least of the 16-bit lanes of `col` whose lane of `max` is `m` (0xffff if none is).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vl")]
#[inline]
fn col16_avx512(max: std::arch::x86_64::__m512i, col: std::arch::x86_64::__m512i, m: i16) -> usize {
    use std::arch::x86_64::*;
    let at = _mm512_mask_mov_epi16(
        _mm512_set1_epi16(-1),
        _mm512_cmpeq_epi16_mask(max, _mm512_set1_epi16(m)),
        col,
    );
    let a = _mm256_min_epu16(_mm512_castsi512_si256(at), _mm512_extracti64x4_epi64::<1>(at));
    let a = _mm_min_epu16(_mm256_castsi256_si128(a), _mm256_extracti128_si256::<1>(a));
    _mm_cvtsi128_si32(_mm_minpos_epu16(a)) as u16 as usize
}

/// What [`run16_avx512`] keeps of a row: each lane's maximum over the row's lines, and the
/// column of the first line holding it.
#[cfg(target_arch = "x86_64")]
#[derive(Clone, Copy)]
#[repr(C, align(64))]
struct Seen16([std::arch::x86_64::__m512i; 2]);

/// [`run_core`] on 16-bit cells in 32-lane vectors, for the passes whose scores fit them
/// ([`narrow`]): the same pass (the same cells, and so the same result), written as one loop
/// with its row kernel, since at these widths what a row costs besides its cells matters as
/// much as they do.
///
/// The cells' NEG is `i16::MIN` and the arithmetic saturates. Every cell of a band is
/// reachable in the modes this serves, with a score the cells hold; an unreachable neighbour
/// outside the band stays near `i16::MIN` through what is added to it, below every real score.
///
/// Rows are laid out for the store buffer: a row reads what the row before it just wrote, and
/// a load is served from a store still in flight only if it reads exactly the bytes the store
/// wrote, within one cache line. A load that isn't waits until the store reaches the cache.
/// So the rows are whole cache lines (columns 32q to 32q + 31), written and read as such
/// whatever the band's columns are: a line's cells outside the band are written as NEG, the
/// line after the band's last is too (the sentinel on the right), and the diagonal neighbours
/// (the row above, one column left) are made in registers from each line and the one before
/// it. A line of NEG before column 0 stands for the one before the first.
///
/// Column 0 is a cell like the others: only a vertical gap reaches it (its diagonal neighbour
/// is that line of NEG), and the horizontal gap leaving it is the one the scan below gives.
///
/// A horizontal gap into cell k from cell j < k of its line scores `x[j] - open - ext * (k -
/// j)`, that is `g[j] - ext * k` for `g[j] = x[j] - open + ext * j`: the best is a prefix
/// maximum of `g`, taken within each 128-bit lane by three shifts and then across them. The
/// maximum may include `g[k]` itself, which gives `x[k] - open`, no more than `x[k]`.
///
/// What the driver needs of a row stays in vectors: each lane's maximum over the row's lines
/// and the first line holding it, reduced only for a row that raises the best score (its
/// column only once, at the end), and the maximum up to the cut over all rows.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vl,avx2,fma,bmi1,bmi2,popcnt,lzcnt")]
#[allow(clippy::too_many_arguments)]
fn run16_avx512<const TRACK: bool>(
    pattern: &[u8],
    text: &[u8],
    mode: Mode,
    sc: &Affine,
    (band, geom): (usize, usize),
    clip: Option<(i32, bool)>,
    cut: usize,
    b: &mut DpBuf,
    resume: Option<&Snap>,
    snaps: &mut [Snap],
) -> Dp {
    use std::arch::x86_64::*;
    const L: usize = 32;
    // The most lines of a row that are held in registers.
    const HELD: usize = 4;
    const MIN: i16 = i16::MIN;
    let (m, n) = (pattern.len(), text.len());
    debug_assert!(narrow(m, mode, sc, (band, geom), clip));
    debug_assert!(TRACK || (cut == usize::MAX && !clip.is_some_and(|(_, cell)| cell)));
    let (below, above) = (
        band.saturating_add(m.saturating_sub(geom)),
        band.saturating_add(geom.saturating_sub(m)),
    );
    let range = |i: usize| (i.saturating_sub(below), i.saturating_add(above).min(n));
    // A row is a line of NEG, then H and F over `w` columns each: the lines of columns 0..=n
    // and one more for the sentinel.
    let w = (n / L + 2) * L;
    let row = L + 2 * w;
    if b.narrow.len() * L < 2 * row {
        b.narrow.resize((2 * row).div_ceil(L), Block16([MIN; L]));
    }
    // The text's bases as the kernel compares them, at their columns (column j reads base
    // j - 1; nothing equals column 0's).
    if b.bases.len() < w {
        b.bases.resize(w, 4);
    }
    b.bases[0] = 4;
    for (x, &t) in b.bases[1..].iter_mut().zip(text) {
        *x = t & 3;
    }
    let bases = b.bases.as_ptr();
    let rows: *mut i16 = b.narrow.as_mut_ptr().cast();
    // SAFETY: the two rows are disjoint parts of `b.narrow`, which nothing else uses below.
    let (mut cur, mut nxt) = unsafe { (rows, rows.add(row)) };
    // Columns `lo..=hi` of H and F of the row at `p`.
    // SAFETY (of its uses): `lo <= hi < w`, and the row is not written while they live.
    let cells = |p: *mut i16, lo: usize, hi: usize| unsafe {
        (
            std::slice::from_raw_parts_mut(p.add(L + lo), hi + 1 - lo),
            std::slice::from_raw_parts_mut(p.add(L + w + lo), hi + 1 - lo),
        )
    };
    let (mat, mis, open, ext) = (sc.mat as i16, sc.mis as i16, sc.open as i16, sc.ext as i16);
    let mut clip_best = 0i32;
    let mut clip_cell = (0usize, 0usize);
    let mut cut_best = 0i32;
    // The first row's band, with the lines it lies in (and the sentinel's) NEG around it.
    let (mut prev, first) = {
        let (lo, hi) = resume.map_or(range(0), |s| (s.lo, s.hi));
        let (h, f) = cells(cur, lo / L * L, ((hi + 1) / L + 1) * L - 1);
        h.fill(MIN);
        f.fill(MIN);
        let (h, f) = cells(cur, lo, hi);
        match resume {
            Some(s) => {
                h.copy_from_slice(&s.h16);
                f.copy_from_slice(&s.f16);
                clip_best = s.clip_best;
                clip_cell = s.clip_cell;
                cut_best = s.cut_best;
            }
            None => {
                for (x, j) in h.iter_mut().zip(lo..) {
                    *x = if j == 0 {
                        0
                    } else {
                        -(sc.open + sc.ext * j as i32) as i16
                    };
                }
            }
        }
        ((lo, hi), resume.map_or(1, |s| s.row + 1))
    };
    let neg = _mm512_set1_epi16(MIN);
    // SAFETY: both rows have the line.
    unsafe {
        _mm512_store_si512(cur.cast(), neg);
        _mm512_store_si512(nxt.cast(), neg);
    }
    let (vgo, vopen, vnext) = (
        _mm512_set1_epi16(open + ext),
        _mm512_set1_epi16(open),
        _mm512_set1_epi16(-ext),
    );
    let (vmat, vmis) = (_mm512_set1_epi16(mat), _mm512_set1_epi16(-mis));
    // SAFETY: the statics hold a vector each.
    let (lane, down, last) = unsafe {
        (
            _mm512_loadu_si512(LANE16.as_ptr().cast()),
            _mm512_loadu_si512(DOWN16.as_ptr().cast()),
            _mm512_loadu_si512(LAST16.as_ptr().cast()),
        )
    };
    let top = _mm512_set1_epi16(L as i16 - 1);
    let (cols, ek0, eks) = (
        _mm512_set1_epi16(L as i16),
        _mm512_mullo_epi16(lane, _mm512_set1_epi16(ext)),
        _mm512_set1_epi16(ext * L as i16),
    );
    let (clipping, xdrop, cell) = clip.map_or((false, 0, false), |(x, c)| (true, x, c));
    // A row stops the pass if no cell of it reaches `clip_best - xdrop`.
    let stop_at = |best: i32| best.saturating_sub(xdrop).max(i32::from(MIN)) as i16;
    let (mut best16, mut stop16) = (clip_best as i16, stop_at(clip_best));
    // The row that last raised the best score: its lanes' maxima and columns (from line
    // `bbase`'s first), and the row (0: none has). These are kept in memory, written whole:
    // a vector the compiler spills to a stack it hasn't aligned for it straddles cache lines.
    let mut bseen = Seen16([neg; 2]);
    let (mut bbase, mut brow) = (0usize, 0usize);
    // Each lane's maximum over the columns up to the cut, over the rows so far.
    let mut cutv = neg;
    let mut next_snap = snaps.partition_point(|s| s.row < first);
    let mut next_row = snaps.get(next_snap).map_or(usize::MAX, |s| s.row);
    // The row the pass stopped at, with the text it had read.
    let mut stopped = None;
    // The last row computed.
    let mut done = first - 1;
    // For rows held in registers: ext, a match's and a mismatch's scores with it, and the
    // last lane of each 128-bit lane (kept from the compiler, which splits the shuffles).
    let (vext, tmat, tmis) = (
        _mm512_set1_epi16(ext),
        _mm512_set1_epi16(mat + ext),
        _mm512_set1_epi16(ext - mis),
    );
    // The lanes after the first, second and third 128-bit lanes.
    let lanes: [__mmask32; 3] = std::hint::black_box([0xffff_ff00, 0xffff_0000, 0xff00_0000]);
    let ends: [__m512i; 4] = std::hint::black_box([
        _mm512_set1_epi16(7),
        _mm512_set1_epi16(15),
        _mm512_set1_epi16(23),
        _mm512_set1_epi16(31),
    ]);
    // The lanes' columns and ext times those, by the line's place among a row's first few.
    let (mut colb, mut ekb) = ([lane; HELD], [ek0; HELD]);
    for b in 1..HELD {
        colb[b] = _mm512_add_epi16(colb[b - 1], cols);
        ekb[b] = _mm512_add_epi16(ekb[b - 1], eks);
    }
    // A line of a row: H and F of its cells from those of the line in the row before
    // (`$hup`, `$fup`) and H of the line before that one (`$hprev`), for the lanes `$real`
    // (NEG in the others) whose bases equal the row's in `$eq`; `$ek` is ext times the lanes'
    // columns counted from any one column, and `$cv` the horizontal-gap maximum entering the
    // line in g's terms (left as the one leaving it).
    macro_rules! line {
        ($hup:expr, $fup:expr, $hprev:expr, $real:expr, $eq:expr, $ek:expr, $cv:ident) => {{
            let (hup, real, ek): (__m512i, __mmask32, __m512i) = ($hup, $real, $ek);
            let fj = _mm512_max_epi16(_mm512_subs_epi16(hup, vgo), _mm512_adds_epi16($fup, vnext));
            let hdiag = _mm512_permutex2var_epi16(hup, down, $hprev);
            let x = _mm512_max_epi16(
                _mm512_adds_epi16(hdiag, _mm512_mask_blend_epi16($eq, vmis, vmat)),
                fj,
            );
            let g = _mm512_mask_adds_epi16(neg, real, _mm512_subs_epi16(x, vopen), ek);
            let p = _mm512_max_epi16(g, _mm512_alignr_epi8::<14>(g, neg));
            let p = _mm512_max_epi16(p, _mm512_alignr_epi8::<12>(p, neg));
            let p = _mm512_max_epi16(p, _mm512_alignr_epi8::<8>(p, neg));
            // Each 128-bit lane's maximum, to the lanes after it.
            let l = _mm512_shuffle_epi8(p, last);
            let p = _mm512_max_epi16(
                _mm512_max_epi16(p, _mm512_alignr_epi64::<6>(l, neg)),
                _mm512_max_epi16(_mm512_alignr_epi64::<4>(l, neg), _mm512_alignr_epi64::<2>(l, neg)),
            );
            let hv = _mm512_mask_max_epi16(neg, real, x, _mm512_subs_epi16(_mm512_max_epi16(p, $cv), ek));
            // Only this max is carried from line to line.
            $cv = _mm512_max_epi16($cv, _mm512_permutexvar_epi16(top, p));
            (hv, _mm512_mask_blend_epi16(real, neg, fj))
        }};
    }
    // Row `$i`'s part in the clip, from its lanes' maxima and their columns (from line
    // `$base`'s first): whether it stops the pass.
    macro_rules! clipped {
        ($mrow:expr, $mcol:expr, $base:expr, $i:expr, $hi:expr) => {{
            let mut stop = false;
            if clipping {
                if _mm512_cmpgt_epi16_mask($mrow, _mm512_set1_epi16(best16)) != 0 {
                    best16 = max16_avx512($mrow);
                    clip_best = i32::from(best16);
                    stop16 = stop_at(clip_best);
                    if TRACK {
                        // SAFETY: a local.
                        unsafe { std::ptr::write_volatile(&mut bseen, Seen16([$mrow, $mcol])) };
                        (bbase, brow) = ($base, $i);
                    }
                } else if _mm512_cmpge_epi16_mask($mrow, _mm512_set1_epi16(stop16)) == 0 {
                    stopped = Some(($i, $hi));
                    stop = true;
                }
            }
            stop
        }};
    }
    // Rows `$i..=$last`, which lie in the same `$nb` lines from `$base`'s: the lines stay in
    // registers from row to row (a row read back from where the one before was just stored
    // waits for those stores), and go back to `cur` after the last. There they are tilted, a
    // cell holding ext times its column (from `$base`) more than its score: a horizontal gap
    // then costs the same wherever it ends, so a row's are a prefix maximum with nothing
    // added before or taken off after, and that is the longest step from one row to the
    // next.
    macro_rules! held {
        ($nb:literal, $i:ident, $last:ident, $base:ident, $end:ident, $within:ident) => {{
            const NB: usize = $nb;
            let (mut h, mut f, mut t) = ([neg; NB], [neg; NB], [_mm256_setzero_si256(); NB]);
            // H of each lane's left neighbour.
            let mut hd = [neg; NB];
            // The lanes of the row before.
            let was = if prev.1 < $base {
                0
            } else {
                (u128::MAX << (prev.0.max($base) - $base)) & (u128::MAX >> (4 * L - 1 - (prev.1 - $base)))
            };
            // SAFETY: the row and the bases have these lines, and the row the line before
            // them (the NEG one, before the first).
            unsafe {
                let mut hprev = _mm512_subs_epi16(_mm512_load_si512(cur.add($base).cast()), vext);
                for b in 0..NB {
                    let real = (was >> (L * b)) as u32;
                    let at = cur.add(L + $base + b * L);
                    h[b] = _mm512_mask_adds_epi16(neg, real, _mm512_load_si512(at.cast()), ekb[b]);
                    f[b] = _mm512_mask_adds_epi16(neg, real, _mm512_load_si512(at.add(w).cast()), ekb[b]);
                    t[b] = _mm256_loadu_si256(bases.add($base + b * L).cast());
                    hd[b] = _mm512_permutex2var_epi16(h[b], down, hprev);
                    hprev = h[b];
                }
            }
            for r in $i..=$last {
                let (lo, hi) = range(r);
                let vc = _mm256_set1_epi8((pattern[r - 1] & 3) as i8);
                // The row's lanes among the lines', and the lanes after one of those.
                let among = (u128::MAX << (lo - $base)) & (u128::MAX >> (4 * L - 1 - (hi - $base)));
                let after = among << 1;
                // The horizontal-gap maximum entering a line.
                let mut cv = neg;
                let (mut mrow, mut mcol) = (neg, neg);
                // The line before's cells, before their horizontal gaps.
                let mut xprev = neg;
                for b in 0..NB {
                    let (real, sreal) = ((among >> (L * b)) as u32, (after >> (L * b)) as u32);
                    let eq = _mm256_cmpeq_epi8_mask(vc, t[b]);
                    let fj = _mm512_max_epi16(_mm512_subs_epi16(h[b], vgo), _mm512_adds_epi16(f[b], vnext));
                    let x = _mm512_max_epi16(
                        _mm512_adds_epi16(hd[b], _mm512_mask_blend_epi16(eq, tmis, tmat)),
                        fj,
                    );
                    let g = _mm512_mask_subs_epi16(neg, real, x, vopen);
                    // The maximum up to each lane within its 128-bit lane: over four lanes,
                    // then eight.
                    let p = _mm512_max_epi16(
                        _mm512_max_epi16(g, _mm512_alignr_epi8::<14>(g, neg)),
                        _mm512_max_epi16(_mm512_alignr_epi8::<12>(g, neg), _mm512_alignr_epi8::<10>(g, neg)),
                    );
                    let p = _mm512_max_epi16(p, _mm512_alignr_epi8::<8>(p, neg));
                    // The 128-bit lanes' maxima, which the lanes after each take.
                    let l = [
                        _mm512_permutexvar_epi16(ends[0], p),
                        _mm512_permutexvar_epi16(ends[1], p),
                        _mm512_permutexvar_epi16(ends[2], p),
                        _mm512_permutexvar_epi16(ends[3], p),
                    ];
                    let l12 = _mm512_mask_max_epi16(l[1], lanes[2], l[1], l[2]);
                    let y = _mm512_max_epi16(x, p);
                    let y = _mm512_mask_max_epi16(y, lanes[0], y, l[0]);
                    let y = _mm512_mask_max_epi16(y, lanes[1], y, l12);
                    let hv = _mm512_mask_max_epi16(neg, real, y, cv);
                    // The same of each lane's left neighbour, for the next row's diagonal:
                    // from the parts as they come, where shifting H would follow it.
                    let y = _mm512_max_epi16(
                        _mm512_permutex2var_epi16(x, down, xprev),
                        _mm512_alignr_epi8::<14>(p, neg),
                    );
                    xprev = x;
                    let y = _mm512_mask_max_epi16(y, lanes[0], y, l[0]);
                    let y = _mm512_mask_max_epi16(y, lanes[1], y, l12);
                    hd[b] = _mm512_mask_max_epi16(neg, sreal, y, cv);
                    cv = _mm512_max_epi16(
                        cv,
                        _mm512_max_epi16(_mm512_max_epi16(l[0], l[1]), _mm512_max_epi16(l[2], l[3])),
                    );
                    (h[b], f[b]) = (hv, _mm512_mask_blend_epi16(real, neg, fj));
                    // The scores.
                    let hv = _mm512_subs_epi16(hv, ekb[b]);
                    if TRACK {
                        mcol = _mm512_mask_mov_epi16(mcol, _mm512_cmpgt_epi16_mask(hv, mrow), colb[b]);
                        cutv =
                            _mm512_mask_max_epi16(cutv, _mm512_cmple_epi16_mask(colb[b], $within), cutv, hv);
                    }
                    mrow = _mm512_max_epi16(mrow, hv);
                }
                prev = (lo, hi);
                done = r;
                if clipped!(mrow, mcol, $base, r, hi) {
                    break;
                }
            }
            if stopped.is_none() {
                // SAFETY: as for the loads, and the rows have a line after column n's.
                unsafe {
                    for b in 0..NB {
                        let at = cur.add(L + $base + b * L);
                        _mm512_store_si512(at.cast(), _mm512_subs_epi16(h[b], ekb[b]));
                        _mm512_store_si512(at.add(w).cast(), _mm512_subs_epi16(f[b], ekb[b]));
                    }
                    _mm512_store_si512(cur.add(L + $end + L).cast(), neg);
                    _mm512_store_si512(cur.add(L + w + $end + L).cast(), neg);
                }
            }
        }};
    }
    let mut i = first;
    while i <= m {
        let (lo, hi) = range(i);
        debug_assert!((lo == 0 || lo == prev.0 + 1) && hi <= prev.1 + 1);
        if lo > hi {
            // Past the text's end by more than the band (only when geom > n): no
            // cell is reachable, so the rest of the pattern can't align.
            stopped = Some((i, n));
            break;
        }
        // The row's lines: columns `base..` to the one holding `hi`.
        let (base, end) = (lo / L * L, hi / L * L);
        // The cut's last column, from `base`.
        let within = if TRACK {
            _mm512_set1_epi16(if cut < base {
                -1
            } else {
                (cut - base).min(i16::MAX as usize) as i16
            })
        } else {
            neg
        };
        // The last of the rows after this one that lie in the same lines (so not past the
        // text's end by more than the band), up to one with a snapshot to take.
        let last = m
            .min(next_row.max(i))
            .min(n + below)
            .min(base + L - 1 + below)
            .min(if n < end + L {
                usize::MAX
            } else {
                end + L - 1 - above
            });
        match (end - base) / L {
            0 => held!(1, i, last, base, end, within),
            1 => held!(2, i, last, base, end, within),
            2 => held!(3, i, last, base, end, within),
            3 => held!(4, i, last, base, end, within),
            _ => {
                let vc = _mm256_set1_epi8((pattern[i - 1] & 3) as i8);
                // The lanes' columns from `base`, the row's among them, and ext times those.
                let (from, to) = (
                    _mm512_set1_epi16((lo - base) as i16),
                    _mm512_set1_epi16((hi - base) as i16),
                );
                let (mut col, mut ek) = (lane, ek0);
                let mut cv = neg;
                let (mut mrow, mut mcol) = (neg, neg);
                // The previous row's line before this one (its last cell is this line's first
                // diagonal neighbour).
                // SAFETY: the line before `base`'s (the NEG one, before the first) is in the row.
                let mut hprev = unsafe { _mm512_load_si512(cur.add(base).cast()) };
                let mut k = base;
                loop {
                    let real = _mm512_mask_cmple_epi16_mask(_mm512_cmpge_epi16_mask(col, from), col, to);
                    // SAFETY: `k <= end <= n`, and the rows and the bases have the lines of
                    // columns 0..=n and one more.
                    let (hup, fup, eq) = unsafe {
                        (
                            _mm512_load_si512(cur.add(L + k).cast()),
                            _mm512_load_si512(cur.add(L + w + k).cast()),
                            _mm256_cmpeq_epi8_mask(vc, _mm256_loadu_si256(bases.add(k).cast())),
                        )
                    };
                    let (hv, fj) = line!(hup, fup, hprev, real, eq, ek, cv);
                    hprev = hup;
                    if TRACK {
                        mcol = _mm512_mask_mov_epi16(mcol, _mm512_cmpgt_epi16_mask(hv, mrow), col);
                        cutv = _mm512_mask_max_epi16(cutv, _mm512_cmple_epi16_mask(col, within), cutv, hv);
                    }
                    mrow = _mm512_max_epi16(mrow, hv);
                    // SAFETY: as for the loads.
                    unsafe {
                        _mm512_store_si512(nxt.add(L + k).cast(), hv);
                        _mm512_store_si512(nxt.add(L + w + k).cast(), fj);
                    }
                    if k == end {
                        break;
                    }
                    (k, col, ek) = (k + L, _mm512_add_epi16(col, cols), _mm512_add_epi16(ek, eks));
                }
                // The sentinel, where it is the next line's first cell (and harmless where
                // not: the row after this one reads that line only then).
                // SAFETY: the rows have a line after column n's.
                unsafe {
                    _mm512_store_si512(nxt.add(L + end + L).cast(), neg);
                    _mm512_store_si512(nxt.add(L + w + end + L).cast(), neg);
                }
                (cur, nxt) = (nxt, cur);
                prev = (lo, hi);
                done = i;
                clipped!(mrow, mcol, base, i, hi);
            }
        }
        if stopped.is_some() {
            break;
        }
        i = done + 1;
        if done == next_row {
            let (lo, hi) = prev;
            let s = &mut snaps[next_snap];
            (s.lo, s.hi) = (lo, hi);
            let (h, f) = cells(cur, lo, hi);
            s.h16.clear();
            s.h16.extend_from_slice(h);
            s.f16.clear();
            s.f16.extend_from_slice(f);
            if cell && brow != 0 {
                // SAFETY: a local.
                let Seen16([bmax, bcol]) = unsafe { std::ptr::read_volatile(&bseen) };
                clip_cell = (brow, bbase + col16_avx512(bmax, bcol, best16));
                brow = 0;
            }
            if cut != usize::MAX {
                cut_best = cut_best.max(i32::from(max16_avx512(cutv)));
            }
            (s.clip_best, s.clip_cell, s.cut_best) = (clip_best, clip_cell, cut_best);
            next_snap += 1;
            next_row = snaps.get(next_snap).map_or(usize::MAX, |s| s.row);
        }
    }
    if cell && brow != 0 {
        // SAFETY: a local.
        let Seen16([bmax, bcol]) = unsafe { std::ptr::read_volatile(&bseen) };
        clip_cell = (brow, bbase + col16_avx512(bmax, bcol, best16));
    }
    if cut != usize::MAX {
        cut_best = cut_best.max(i32::from(max16_avx512(cutv)));
    }
    if crate::stats::on() {
        let cells: usize = (first..=done).map(|i| range(i).1 + 1 - range(i).0).sum();
        DP_WORK.with(|x| {
            let (c, r) = x.get();
            x.set((c + cells as u64, r + (done + 1 - first) as u64));
        });
    }
    if let Some((last_row, text_read)) = stopped {
        return Dp {
            score: NEG,
            score_col: 0,
            clip_best,
            clip_cell,
            text_read,
            last_row,
            cut_best,
            cut_score: NEG,
        };
    }
    let h = &*cells(cur, prev.0, prev.1).0;
    let most = |v: &[i16]| v.iter().copied().max().map_or(NEG, i32::from);
    let (score, score_col) = match mode {
        // A band laid out for another text's length may end short of column n.
        Mode::Global => (
            if prev.1 == n {
                i32::from(h[n - prev.0])
            } else {
                NEG
            },
            n,
        ),
        Mode::Prefix | Mode::Infix => {
            let x = most(h);
            (x, prev.0 + h.iter().position(|&y| i32::from(y) == x).unwrap_or(0))
        }
    };
    let cut_score = if cut != usize::MAX && prev.0 <= cut {
        most(&h[..=prev.1.min(cut) - prev.0])
    } else {
        NEG
    };
    Dp {
        score,
        score_col,
        clip_best,
        clip_cell,
        text_read: prev.1,
        last_row: m,
        cut_best,
        cut_score,
    }
}

/// [`row_avx512`] in 8-lane AVX2 vectors. The horizontal-gap prefix max is a three-step
/// shift-and-max scan within each vector (shifts cross the two 128-bit halves through a
/// permute that also fills the low lanes with NEG), with the running maximum carried across
/// vectors in a register. A row's last `w % 8` cells are covered by one more whole vector
/// ending at the row's end, which recomputes the cells it shares with the previous vector
/// from the maximum entering its first cell (no masked loads or stores: those of cells the
/// previous row just stored stall on store forwarding, and masked stores are slow on some
/// CPUs). Rows under 8 cells take the scalar steps.
///
/// Same results as the scalar passes: the integer operations per cell are theirs, and the scan
/// differs from the exact prefix maximum only where both lie below NEG, where neither can win
/// against the cell's diagonal/vertical score (which is at least NEG).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
#[allow(clippy::too_many_arguments)]
unsafe fn row_avx2(
    hprev: &[i32],
    fu: &[i32],
    tx: &[u8],
    nh: &mut [i32],
    nf: &mut [i32],
    c: u8,
    (mat, mis, go, ext): (i32, i32, i32, i32),
    e: i32,
    floor: i32,
) -> i32 {
    use std::arch::x86_64::*;
    let w = fu.len();
    assert!(hprev.len() == w + 1 && tx.len() == w && nh.len() == w + 1 && nf.len() == w + 1);
    (nh[w], nf[w]) = (NEG, NEG);
    if w < 8 {
        return cells_scalar(hprev, fu, tx, nh, nf, c, (mat, mis, go, ext), floor, 0, e, NEG);
    }
    let (neg, vfloor) = (_mm256_set1_epi32(NEG), _mm256_set1_epi32(floor));
    let (vgo, vext) = (_mm256_set1_epi32(go), _mm256_set1_epi32(ext));
    let (vmat, vmis) = (_mm256_set1_epi32(mat), _mm256_set1_epi32(-mis));
    let vc = _mm256_set1_epi32(c as i32);
    let three = _mm256_set1_epi32(3);
    let seven = _mm256_set1_epi32(7);
    let lane = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    let lext = _mm256_mullo_epi32(lane, vext);
    // [NEG x 4, v's low half]: `v` shifted up by four lanes, and the source of the lanes that
    // cross into the high half when shifting by fewer.
    let up4 = |v: __m256i| _mm256_permute2x128_si256::<0x02>(v, neg);
    let mut mx = neg;
    // Cells k..k + 8 with `cv` the maximum entering cell k (broadcast) and `ek` ext times
    // each lane's cell index: stores H and F and returns the exclusive maxima entering each
    // cell and the maximum after the last.
    let mut cells = |k: usize, cv: __m256i, ek: __m256i| {
        // SAFETY: k + 8 <= w, so cells k..k + 8 are in bounds of every slice (`hprev` has
        // w + 1).
        let (hd, hu, fv, t) = unsafe {
            (
                _mm256_loadu_si256(hprev.as_ptr().add(k).cast()),
                _mm256_loadu_si256(hprev.as_ptr().add(k + 1).cast()),
                _mm256_loadu_si256(fu.as_ptr().add(k).cast()),
                _mm_loadl_epi64(tx.as_ptr().add(k).cast()),
            )
        };
        let fj = _mm256_max_epi32(
            _mm256_max_epi32(_mm256_sub_epi32(hu, vgo), _mm256_sub_epi32(fv, vext)),
            neg,
        );
        let t = _mm256_and_si256(_mm256_cvtepu8_epi32(t), three);
        let s = _mm256_blendv_epi8(vmis, vmat, _mm256_cmpeq_epi32(t, vc));
        let x = _mm256_max_epi32(_mm256_max_epi32(_mm256_add_epi32(hd, s), fj), vfloor);
        // gs[k] = x - go + ext * (k + 1).
        let g = _mm256_add_epi32(_mm256_sub_epi32(x, vgo), _mm256_add_epi32(ek, vext));
        // Inclusive prefix max within the vector: shifts by one, two, four lanes.
        let v = _mm256_max_epi32(g, _mm256_alignr_epi8::<12>(g, up4(g)));
        let v = _mm256_max_epi32(v, _mm256_alignr_epi8::<8>(v, up4(v)));
        let v4 = up4(v);
        let incl = _mm256_max_epi32(v, v4);
        // Exclusive: shifted up one lane, with the carry. (`v4` stands in for `up4(incl)`:
        // they differ only in lanes below NEG.)
        let excl = _mm256_max_epi32(_mm256_alignr_epi8::<12>(incl, v4), cv);
        let hv = _mm256_max_epi32(x, _mm256_sub_epi32(excl, ek));
        // SAFETY: as for the loads.
        unsafe {
            _mm256_storeu_si256(nh.as_mut_ptr().add(k).cast(), hv);
            _mm256_storeu_si256(nf.as_mut_ptr().add(k).cast(), fj);
        }
        mx = _mm256_max_epi32(mx, hv);
        (
            excl,
            _mm256_max_epi32(cv, _mm256_permutevar8x32_epi32(incl, seven)),
        )
    };
    let mut cv = _mm256_set1_epi32(e);
    let mut excl = cv;
    let mut ek = lext;
    let step = _mm256_set1_epi32(8 * ext);
    let mut k = 0;
    while k + 8 <= w {
        (excl, cv) = cells(k, cv, ek);
        ek = _mm256_add_epi32(ek, step);
        k += 8;
    }
    if k < w {
        // The last vector, at w - 8: lane w - k of the previous one is the maximum entering
        // its first cell.
        let cv = _mm256_permutevar8x32_epi32(excl, _mm256_set1_epi32((w - k) as i32));
        cells(
            w - 8,
            cv,
            _mm256_add_epi32(lext, _mm256_set1_epi32(ext * (w - 8) as i32)),
        );
    }
    let m = _mm_max_epi32(_mm256_castsi256_si128(mx), _mm256_extracti128_si256::<1>(mx));
    let m = _mm_max_epi32(m, _mm_shuffle_epi32::<0b01_00_11_10>(m));
    _mm_cvtsi128_si32(_mm_max_epi32(m, _mm_shuffle_epi32::<0b10_11_00_01>(m)))
}

#[cfg(target_arch = "x86_64")]
struct Avx2Rows;

#[cfg(target_arch = "x86_64")]
impl Rows for Avx2Rows {
    #[inline(always)]
    unsafe fn row(
        hprev: &[i32],
        fu: &[i32],
        tx: &[u8],
        nh: &mut [i32],
        nf: &mut [i32],
        c: u8,
        sc: (i32, i32, i32, i32),
        e: i32,
        floor: i32,
        _: (&mut Vec<i32>, &mut Vec<i32>),
    ) -> i32 {
        // SAFETY: as the caller guarantees.
        unsafe { row_avx2(hprev, fu, tx, nh, nf, c, sc, e, floor) }
    }
}

/// [`row_avx512`] in 4-lane NEON blocks, the rest scalar: the horizontal-gap scan is two
/// shift-and-max steps per block, and only the block maximum is carried, so the carried
/// dependency is one max per four cells (the scalar scan's is one per cell, and its rows
/// round-trip through memory between passes).
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
#[allow(clippy::too_many_arguments)]
unsafe fn row_neon(
    hprev: &[i32],
    fu: &[i32],
    tx: &[u8],
    nh: &mut [i32],
    nf: &mut [i32],
    c: u8,
    (mat, mis, go, ext): (i32, i32, i32, i32),
    e: i32,
    floor: i32,
) -> i32 {
    use std::arch::aarch64::*;
    let w = fu.len();
    assert!(hprev.len() == w + 1 && tx.len() == w && nh.len() == w + 1 && nf.len() == w + 1);
    (nh[w], nf[w]) = (NEG, NEG);
    let (neg, vfloor) = (vdupq_n_s32(NEG), vdupq_n_s32(floor));
    let (vgo, vext) = (vdupq_n_s32(go), vdupq_n_s32(ext));
    let (vmat, vmis) = (vdupq_n_s32(mat), vdupq_n_s32(-mis));
    let (vc, three) = (vdupq_n_u32(c as u32), vdupq_n_u32(3));
    // ext * k for the block's lanes, advanced by 4 * ext per block.
    let lane: [i32; 4] = [0, 1, 2, 3];
    // SAFETY: `lane` holds four i32s.
    let mut ek = vmulq_s32(unsafe { vld1q_s32(lane.as_ptr()) }, vext);
    let step = vdupq_n_s32(4 * ext);
    // max(e, gs[..k]) in every lane.
    let mut cv = vdupq_n_s32(e);
    let mut mx = neg;
    let mut k = 0;
    while k + 4 <= w {
        // SAFETY: k + 4 <= w, and hprev holds w + 1 cells.
        let (hd, hu, fv, t4) = unsafe {
            (
                vld1q_s32(hprev.as_ptr().add(k)),
                vld1q_s32(hprev.as_ptr().add(k + 1)),
                vld1q_s32(fu.as_ptr().add(k)),
                tx.as_ptr().add(k).cast::<u32>().read_unaligned(),
            )
        };
        let t = vmovl_u16(vget_low_u16(vmovl_u8(vreinterpret_u8_u32(vdup_n_u32(t4)))));
        let t = vandq_u32(t, three);
        let fj = vmaxq_s32(vmaxq_s32(vsubq_s32(hu, vgo), vsubq_s32(fv, vext)), neg);
        let s = vbslq_s32(vceqq_u32(t, vc), vmat, vmis);
        let x = vmaxq_s32(vmaxq_s32(vaddq_s32(hd, s), fj), vfloor);
        // gs[k] = x - go + ext * (k + 1).
        let g = vaddq_s32(vsubq_s32(x, vgo), vaddq_s32(ek, vext));
        // Inclusive prefix max within the block; then the carry, exclusive by one lane.
        let mut v = vmaxq_s32(g, vextq_s32::<3>(neg, g));
        v = vmaxq_s32(v, vextq_s32::<2>(neg, v));
        let excl = vextq_s32::<3>(cv, vmaxq_s32(v, cv));
        cv = vmaxq_s32(cv, vdupq_laneq_s32::<3>(v));
        let hv = vmaxq_s32(x, vsubq_s32(excl, ek));
        // SAFETY: as for the loads.
        unsafe {
            vst1q_s32(nh.as_mut_ptr().add(k), hv);
            vst1q_s32(nf.as_mut_ptr().add(k), fj);
        }
        mx = vmaxq_s32(mx, hv);
        ek = vaddq_s32(ek, step);
        k += 4;
    }
    let (m, g) = (vmaxvq_s32(mx), vgetq_lane_s32::<0>(cv));
    cells_scalar(hprev, fu, tx, nh, nf, c, (mat, mis, go, ext), floor, k, g, m)
}

/// [`row_neon`] (NEON is part of the aarch64 baseline).
#[cfg(target_arch = "aarch64")]
struct NeonRows;

#[cfg(target_arch = "aarch64")]
impl Rows for NeonRows {
    #[inline(always)]
    unsafe fn row(
        hprev: &[i32],
        fu: &[i32],
        tx: &[u8],
        nh: &mut [i32],
        nf: &mut [i32],
        c: u8,
        sc: (i32, i32, i32, i32),
        e: i32,
        floor: i32,
        _: (&mut Vec<i32>, &mut Vec<i32>),
    ) -> i32 {
        // SAFETY: as the caller guarantees.
        unsafe { row_neon(hprev, fu, tx, nh, nf, c, sc, e, floor) }
    }
}

/// What one DP pass found.
#[derive(Clone, Copy)]
struct Dp {
    /// Score in the pass's mode, and the text column where it ends (Prefix mode).
    score: i32,
    score_col: usize,
    /// With `clip`, the best cell anywhere and its (pattern, text) position.
    clip_best: i32,
    clip_cell: (usize, usize),
    /// Text bases the pass read (a prefix of the text).
    text_read: usize,
    /// The last row computed: the pattern length, or the row where the pass stopped.
    last_row: usize,
    /// With a cut: `clip_best` and the Prefix-mode score over text columns up to the cut only
    /// (the pass over the text cut there, X-drop aside).
    cut_best: i32,
    cut_score: i32,
}

/// A computed DP row (`h` and `f` over its band `[lo, hi]`) and the clip state after it, from
/// which a pass over a text sharing the bases the rows up to it read can resume.
#[derive(Clone, Default)]
struct Snap {
    row: usize,
    lo: usize,
    hi: usize,
    h: Vec<i32>,
    f: Vec<i32>,
    /// The row, if the pass is on 16-bit cells.
    #[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
    h16: Vec<i16>,
    #[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
    f16: Vec<i16>,
    clip_best: i32,
    clip_cell: (usize, usize),
    cut_best: i32,
}

/// Reused per-thread DP rows.
#[derive(Default)]
struct DpBuf {
    h: Vec<i32>,
    f: Vec<i32>,
    nh: Vec<i32>,
    nf: Vec<i32>,
    hp: Vec<i32>,
    gs: Vec<i32>,
    /// The 16-bit passes' two rows, and their text's bases.
    #[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
    narrow: Vec<Block16>,
    #[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
    bases: Vec<u8>,
    snaps: Vec<Snap>,
}

thread_local! {
    static DP_BUF: std::cell::RefCell<DpBuf> = std::cell::RefCell::new(DpBuf::default());
    /// DP cells and rows computed on this thread (counted only with `SPEARFISH_STATS`).
    static DP_WORK: std::cell::Cell<(u64, u64)> = const { std::cell::Cell::new((0, 0)) };
}

/// Cells and rows `affine_dp` has computed on this thread so far (profiling).
pub fn dp_work() -> (u64, u64) {
    DP_WORK.with(|w| w.get())
}

/// The score in `mode`; with `clip = Some((xdrop, cell))`, also the best score anywhere (0 if
/// none) and, if `cell`, where it is, and stops early (returning NEG for the score in `mode`)
/// once a row falls more than `xdrop` below it. Rows live in per-thread buffers, and each pass runs over
/// explicit slices so it vectorizes without bounds checks.
/// The band is laid out for a text of `geom` bases (usually `text.len()`).
fn affine_dp(
    pattern: &[u8],
    text: &[u8],
    mode: Mode,
    sc: &Affine,
    band: usize,
    clip: Option<(i32, bool)>,
    geom: usize,
) -> Dp {
    let level = simd_level();
    DP_BUF.with(|b| {
        run(
            level,
            pattern,
            text,
            mode,
            sc,
            (band, geom),
            clip,
            usize::MAX,
            &mut b.borrow_mut(),
            None,
            &mut [],
        )
    })
}

/// [`affine_dp`] over several texts in sorted order with one pattern, band and clip: a pass
/// resumes from the previous text's rows wherever they read only bases the two texts share
/// and the bands have the same shape, so isoforms that agree near the anchor are aligned
/// there once. Results equal separate [`affine_dp`] calls.
#[allow(clippy::too_many_arguments)]
fn affine_dp_sorted(
    pattern: &[u8],
    texts: &[&[u8]],
    mode: Mode,
    sc: &Affine,
    band: usize,
    clip: Option<(i32, bool)>,
    geom: impl Fn(usize) -> usize,
    cut: usize,
) -> Vec<Dp> {
    let m = pattern.len();
    let shape = |n: usize| {
        let g = geom(n);
        (
            band.saturating_add(m.saturating_sub(g)),
            band.saturating_add(g.saturating_sub(m)),
        )
    };
    let full = mode == Mode::Infix || band == usize::MAX;
    // Row up to which text k's pass equals text k - 1's: rows i read text bases below i + above
    // (capped at the text's end), so rows with i + above <= the common prefix agree. The last
    // row is always computed (its score depends on the mode and the text's length).
    let resume: Vec<usize> = (0..texts.len())
        .map(|k| {
            if k == 0 || full || shape(texts[k].len()) != shape(texts[k - 1].len()) {
                return 0;
            }
            let l = texts[k]
                .iter()
                .zip(texts[k - 1])
                .take_while(|(a, b)| a == b)
                .count();
            l.saturating_sub(shape(texts[k].len()).1).min(m.saturating_sub(1))
        })
        .collect();
    let mut want: Vec<usize> = resume.iter().copied().filter(|&r| r > 0).collect();
    let level = simd_level();
    want.sort_unstable();
    want.dedup();
    DP_BUF.with(|b| {
        let b = &mut *b.borrow_mut();
        let mut snaps = std::mem::take(&mut b.snaps);
        snaps.resize_with(want.len(), Snap::default);
        for (s, &r) in snaps.iter_mut().zip(&want) {
            s.row = r;
        }
        let mut out: Vec<Dp> = Vec::with_capacity(texts.len());
        for (k, t) in texts.iter().enumerate() {
            let r = resume[k];
            // A pass depends on its text only through the bases it read (below the band's
            // clamp at the text's end): a text agreeing with the previous one there, with a
            // band of the same shape, has the same result. With clipping the pass usually
            // reads far less than the window (X-drop).
            let same = k > 0 && clip.is_some() && !full && {
                let (p, tr) = (texts[k - 1], out[k - 1].text_read);
                shape(t.len()) == shape(p.len()) && tr < p.len() && tr < t.len() && t[..tr] == p[..tr]
            };
            let d = if same {
                out[k - 1]
            } else if r == 0 {
                run(
                    level,
                    pattern,
                    t,
                    mode,
                    sc,
                    (band, geom(t.len())),
                    clip,
                    cut,
                    b,
                    None,
                    &mut snaps,
                )
            } else if out[k - 1].last_row <= r {
                // The previous pass stopped within the shared rows: so does this one.
                out[k - 1]
            } else {
                let i = want.binary_search(&r).unwrap();
                let s = std::mem::take(&mut snaps[i]);
                let d = run(
                    level,
                    pattern,
                    t,
                    mode,
                    sc,
                    (band, geom(t.len())),
                    clip,
                    cut,
                    b,
                    Some(&s),
                    &mut snaps,
                );
                snaps[i] = s;
                d
            };
            out.push(d);
        }
        b.snaps = snaps;
        out
    })
}

/// One DP pass at `level` (which the CPU must support): the whole pass is compiled once per
/// level, so the row kernel inlines into the row loop and nothing is dispatched per row.
#[allow(clippy::too_many_arguments)]
fn run(
    level: Level,
    pattern: &[u8],
    text: &[u8],
    mode: Mode,
    sc: &Affine,
    bg: (usize, usize),
    clip: Option<(i32, bool)>,
    cut: usize,
    b: &mut DpBuf,
    resume: Option<&Snap>,
    snaps: &mut [Snap],
) -> Dp {
    debug_assert!(simd::detected().includes(level));
    match level {
        // SAFETY: the CPU has the level's features (callers pass the active level).
        #[cfg(target_arch = "x86_64")]
        Level::Avx512 => unsafe {
            if !narrow(pattern.len(), mode, sc, bg, clip) {
                run_avx512(pattern, text, mode, sc, bg, clip, cut, b, resume, snaps)
            } else if cut != usize::MAX || clip.is_some_and(|(_, cell)| cell) {
                run16_avx512::<true>(pattern, text, mode, sc, bg, clip, cut, b, resume, snaps)
            } else {
                run16_avx512::<false>(pattern, text, mode, sc, bg, clip, cut, b, resume, snaps)
            }
        },
        #[cfg(target_arch = "x86_64")]
        Level::Avx2 => unsafe { run_avx2(pattern, text, mode, sc, bg, clip, cut, b, resume, snaps) },
        // NEON is part of the aarch64 baseline: no entry point of its own.
        #[cfg(target_arch = "aarch64")]
        Level::Neon => run_core::<NeonRows>(pattern, text, mode, sc, bg, clip, cut, b, resume, snaps),
        _ => run_core::<ScalarRows>(pattern, text, mode, sc, bg, clip, cut, b, resume, snaps),
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vl,avx2,fma,bmi1,bmi2,popcnt,lzcnt")]
#[allow(clippy::too_many_arguments)]
fn run_avx512(
    pattern: &[u8],
    text: &[u8],
    mode: Mode,
    sc: &Affine,
    bg: (usize, usize),
    clip: Option<(i32, bool)>,
    cut: usize,
    b: &mut DpBuf,
    resume: Option<&Snap>,
    snaps: &mut [Snap],
) -> Dp {
    run_core::<Avx512Rows>(pattern, text, mode, sc, bg, clip, cut, b, resume, snaps)
}

#[cfg(all(test, target_arch = "x86_64"))]
thread_local! {
    /// Set by tests to keep every pass on 32-bit cells.
    static WIDE_ONLY: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Whether a banded pass's scores fit 16-bit cells with room for the row kernel's
/// temporaries, so that [`run16_avx512`] computes what 32-bit rows would. Every cell of the
/// band is reached from the origin by diagonal steps and then a gap no longer than the band's
/// reach to either side: its score is at least `-(mis * m + open + ext * reach)`, and it is at
/// most `mat * m`. The kernel adds gap costs over a row's width to those, and a match's score
/// to an unreachable cell's.
#[cfg(target_arch = "x86_64")]
fn narrow(
    m: usize,
    mode: Mode,
    sc: &Affine,
    (band, geom): (usize, usize),
    clip: Option<(i32, bool)>,
) -> bool {
    #[cfg(test)]
    if WIDE_ONLY.with(|w| w.get()) {
        return false;
    }
    // The pass relies on a row that raises the best score not stopping it.
    if mode == Mode::Infix || band == usize::MAX || clip.is_some_and(|(xdrop, _)| xdrop < 0) {
        return false;
    }
    let (Ok(mat), Ok(mis), Ok(open), Ok(ext)) = (
        u64::try_from(sc.mat),
        u64::try_from(sc.mis),
        u64::try_from(sc.open),
        u64::try_from(sc.ext),
    ) else {
        return false;
    };
    let below = band.saturating_add(m.saturating_sub(geom)) as u64;
    let above = band.saturating_add(geom.saturating_sub(m)) as u64;
    // A row's cells, with the lanes of its lines before and after them.
    let width = below.saturating_add(above).saturating_add(66);
    let gaps = below.max(above).saturating_add(width).saturating_add(1);
    let low = mis
        .saturating_mul(m as u64)
        .saturating_add(open.saturating_mul(2))
        .saturating_add(ext.saturating_mul(gaps))
        .saturating_add(mat);
    let high = mat
        .saturating_mul(m as u64)
        .saturating_add(ext.saturating_mul(width));
    // The kernel also keeps columns within a row as 16-bit lanes.
    low.max(high) <= 32_000 && width <= 30_000
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma,bmi1,bmi2,popcnt,lzcnt")]
#[allow(clippy::too_many_arguments)]
fn run_avx2(
    pattern: &[u8],
    text: &[u8],
    mode: Mode,
    sc: &Affine,
    bg: (usize, usize),
    clip: Option<(i32, bool)>,
    cut: usize,
    b: &mut DpBuf,
    resume: Option<&Snap>,
    snaps: &mut [Snap],
) -> Dp {
    run_core::<Avx2Rows>(pattern, text, mode, sc, bg, clip, cut, b, resume, snaps)
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn run_core<R: Rows>(
    pattern: &[u8],
    text: &[u8],
    mode: Mode,
    sc: &Affine,
    (band, geom): (usize, usize),
    clip: Option<(i32, bool)>,
    cut: usize,
    b: &mut DpBuf,
    resume: Option<&Snap>,
    snaps: &mut [Snap],
) -> Dp {
    let (m, n) = (pattern.len(), text.len());
    // With `clip`, the best cell anywhere (0 aligns nothing), and among text columns up to
    // `cut` (usize::MAX: no cut).
    let mut clip_best = 0i32;
    let mut clip_cell = (0usize, 0usize);
    let mut cut_best = 0i32;
    let (mat, mis, ext) = (sc.mat, sc.mis, sc.ext);
    let go = sc.open + sc.ext;
    // Row i (pattern) covers text columns [i - band - (m - g)⁺, i + band +
    // (g - m)⁺] up to n, for the band's text length g; Infix is unbanded.
    let full = mode == Mode::Infix || band == usize::MAX;
    let (below, above) = (
        band.saturating_add(m.saturating_sub(geom)),
        band.saturating_add(geom.saturating_sub(m)),
    );
    let range = |i: usize| {
        if full {
            (0, n)
        } else {
            (i.saturating_sub(below), i.saturating_add(above).min(n))
        }
    };
    // h: best score ending at (i, j); f: ending in a gap that consumes pattern
    // (vertical). Rows alternate between two buffers.
    // Rows are not cleared: a row's band is written in full with a NEG sentinel just right of
    // it (stored by the row kernel), and a row reads the previous row only within that row's
    // band and sentinel. None is needed on the left: a banded row starting at lo > 0 reads
    // from column lo - 1, where the previous row starts.
    for v in [&mut b.h, &mut b.f, &mut b.nh, &mut b.nf] {
        if v.len() < n + 2 {
            v.resize(n + 2, NEG);
        }
    }
    // Rows are swapped as local slices (in registers), not as the buffers' Vec headers in
    // `b`: a store and reload of those per row stalls depending on where `b` lies.
    let DpBuf {
        h: hv,
        f: fv,
        nh: nhv,
        nf: nfv,
        hp: hpv,
        gs: gsv,
        ..
    } = &mut *b;
    let (mut h, mut f, mut nh, mut nf) = (&mut hv[..], &mut fv[..], &mut nhv[..], &mut nfv[..]);
    let edges = |v: &mut [i32], lo: usize, hi: usize| {
        if lo > 0 {
            v[lo - 1] = NEG;
        }
        v[hi + 1] = NEG;
    };
    debug_assert!(sc.open >= 0 && sc.ext >= 0);
    let (mut prev, first) = match resume {
        Some(s) => {
            h[s.lo..=s.hi].copy_from_slice(&s.h);
            f[s.lo..=s.hi].copy_from_slice(&s.f);
            edges(h, s.lo, s.hi);
            edges(f, s.lo, s.hi);
            clip_best = s.clip_best;
            clip_cell = s.clip_cell;
            cut_best = s.cut_best;
            ((s.lo, s.hi), s.row + 1)
        }
        None => {
            let (lo0, hi0) = range(0);
            for (j, x) in h.iter_mut().enumerate().take(hi0 + 1).skip(lo0) {
                *x = if j == 0 || mode == Mode::Infix {
                    0
                } else {
                    -(sc.open + sc.ext * j as i32)
                };
            }
            f[lo0..=hi0].fill(NEG);
            edges(h, lo0, hi0);
            edges(f, lo0, hi0);
            ((lo0, hi0), 1)
        }
    }; // band of the row in h
    let mut next_snap = snaps.partition_point(|s| s.row < first);
    let stats = crate::stats::on();
    for i in first..=m {
        let (lo, hi) = range(i);
        debug_assert!((lo == 0 || lo == prev.0 + 1) && hi <= prev.1 + 1);
        if lo > hi {
            // Past the text's end by more than the band (only when geom > n): no
            // cell is reachable, so the rest of the pattern can't align.
            return Dp {
                score: NEG,
                score_col: 0,
                clip_best,
                clip_cell,
                text_read: n,
                last_row: i,
                cut_best,
                cut_score: NEG,
            };
        }
        let c = pattern[i - 1] & 3;
        let mut start = lo;
        // Horizontal gap (consumes text) ending at (i, start - 1): from column
        // 0 if it is in the band, else unreachable.
        let mut e = NEG;
        // The row's maximum, kept as the row is written (for clipping).
        let mut row_max = NEG;
        if lo == 0 {
            // Column 0: only a vertical gap reaches it.
            let fj = (h[0] - go).max(f[0] - ext);
            nh[0] = fj;
            nf[0] = fj;
            e = fj - go;
            start = 1;
            row_max = fj;
        }
        if stats {
            DP_WORK.with(|x| {
                let (c, r) = x.get();
                x.set((c + (hi + 1 - lo) as u64, r + 1));
            });
        }
        if start <= hi {
            // SAFETY: the entry point is compiled for the kernel's features, which the CPU has.
            let m = unsafe {
                R::row(
                    &h[start - 1..=hi],
                    &f[start..=hi],
                    &text[start - 1..hi],
                    &mut nh[start..=hi + 1],
                    &mut nf[start..=hi + 1],
                    c,
                    (mat, mis, go, ext),
                    e,
                    NEG,
                    (hpv, gsv),
                )
            };
            row_max = row_max.max(m);
        } else {
            (nh[hi + 1], nf[hi + 1]) = (NEG, NEG);
        }
        if cut != usize::MAX && lo <= cut {
            let row = if hi <= cut {
                row_max
            } else {
                nh[lo..=cut].iter().copied().max().unwrap_or(NEG)
            };
            cut_best = cut_best.max(row);
        }
        if let Some((xdrop, cell)) = clip {
            let row = row_max;
            if row > clip_best {
                clip_best = row;
                if cell {
                    let j = nh[lo..=hi].iter().position(|&x| x == row).unwrap();
                    clip_cell = (i, lo + j);
                }
            }
            if row < clip_best.saturating_sub(xdrop) {
                return Dp {
                    score: NEG,
                    score_col: 0,
                    clip_best,
                    clip_cell,
                    text_read: hi,
                    last_row: i,
                    cut_best,
                    cut_score: NEG,
                };
            }
        }
        std::mem::swap(&mut h, &mut nh);
        std::mem::swap(&mut f, &mut nf);
        prev = (lo, hi);
        if next_snap < snaps.len() && snaps[next_snap].row == i {
            let s = &mut snaps[next_snap];
            (s.lo, s.hi) = (lo, hi);
            s.h.clear();
            s.h.extend_from_slice(&h[lo..=hi]);
            s.f.clear();
            s.f.extend_from_slice(&f[lo..=hi]);
            (s.clip_best, s.clip_cell, s.cut_best) = (clip_best, clip_cell, cut_best);
            next_snap += 1;
        }
    }
    let (score, score_col) = match mode {
        Mode::Global => (h[n], n),
        Mode::Prefix | Mode::Infix => {
            let row = &h[prev.0..=prev.1];
            let x = row.iter().copied().max().unwrap_or(NEG);
            (x, prev.0 + row.iter().position(|&y| y == x).unwrap_or(0))
        }
    };
    let cut_score = if cut != usize::MAX && prev.0 <= cut {
        h[prev.0..=prev.1.min(cut)].iter().copied().max().unwrap_or(NEG)
    } else {
        NEG
    };
    Dp {
        score,
        score_col,
        clip_best,
        clip_cell,
        text_read: prev.1,
        last_row: m,
        cut_best,
        cut_score,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    /// [`affine_extend_sorted_cut`]'s whole-text scores.
    #[allow(clippy::too_many_arguments)]
    fn affine_extend_sorted(
        pattern: &[u8],
        texts: &[&[u8]],
        sc: &Affine,
        band: usize,
        end_bonus: i32,
        xdrop: i32,
        diagonal: bool,
        window: usize,
    ) -> Vec<i32> {
        affine_extend_sorted_cut(
            pattern,
            texts,
            sc,
            band,
            end_bonus,
            xdrop,
            diagonal,
            window,
            usize::MAX,
        )
        .into_iter()
        .map(|(x, _, _)| x)
        .collect()
    }

    /// One text's score as [`affine_extend_sorted_cut`] computes it with `diagonal`.
    fn affine_extend_diagonal(
        pattern: &[u8],
        text: &[u8],
        sc: &Affine,
        band: usize,
        end_bonus: i32,
        xdrop: i32,
    ) -> i32 {
        let geom = pattern.len().min(text.len());
        let d = affine_dp(pattern, text, Mode::Prefix, sc, band, Some((xdrop, false)), geom);
        d.clip_best.max(d.score.saturating_add(end_bonus))
    }

    /// [`affine_extend_window`] on a whole text.
    fn affine_extend_end(
        pattern: &[u8],
        text: &[u8],
        sc: &Affine,
        band: usize,
        end_bonus: i32,
        xdrop: i32,
    ) -> (i32, usize, usize, usize) {
        affine_extend_window(pattern, text, text.len(), sc, band, end_bonus, xdrop)
    }

    /// Row kernel timing (not a test): `cargo test --release -p spearfish-map bench_rows --
    /// --ignored --nocapture`.
    #[cfg(target_arch = "x86_64")]
    #[test]
    #[ignore]
    fn bench_rows() {
        let mut x = 0x9e3779b97f4a7c15u64;
        let mut rnd = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let only: Option<Vec<usize>> = std::env::var("KB")
            .ok()
            .map(|v| v.split(',').map(|x| x.parse().unwrap()).collect());
        for (m, band) in [
            (80usize, 16usize),
            (80, 32),
            (80, 48),
            (400, 32),
            (400, 64),
            (400, 80),
            (2000, 100),
        ] {
            let pattern: Vec<u8> = (0..m).map(|_| (rnd() & 3) as u8).collect();
            let mut text = Vec::new();
            for &c in &pattern {
                match rnd() % 100 {
                    0..=2 => {}
                    3..=5 => text.extend([c, (rnd() & 3) as u8]),
                    6..=8 => text.push((rnd() & 3) as u8),
                    _ => text.push(c),
                }
            }
            let reps = 40_000_000 / m;
            for wide in [false, true] {
                if only
                    .as_ref()
                    .is_some_and(|o| (o[0], o[1], o[2] != 0) != (m, band, wide))
                {
                    continue;
                }
                WIDE_ONLY.with(|w| w.set(wide));
                let mut best = f64::MAX;
                let mut sum = 0i64;
                for _ in 0..5 {
                    let t0 = std::time::Instant::now();
                    for _ in 0..reps {
                        let r = affine_extend_sorted_cut(
                            std::hint::black_box(&pattern),
                            &[&text[..]],
                            &Affine::ONT,
                            band,
                            10,
                            200,
                            true,
                            0,
                            usize::MAX,
                        );
                        sum += i64::from(r[0].0);
                    }
                    best = best.min(t0.elapsed().as_secs_f64());
                }
                WIDE_ONLY.with(|w| w.set(false));
                println!(
                    "m {m} band {band} {} {:.2} ns/row ({sum})",
                    if wide { "wide" } else { "narrow" },
                    best * 1e9 / (reps * m) as f64
                );
            }
        }
    }

    fn naive(p: &[u8], t: &[u8], mode: Mode) -> u32 {
        let (m, n) = (p.len(), t.len());
        let mut prev: Vec<u32> = (0..=m as u32).collect();
        let mut best = if mode == Mode::Global && n > 0 {
            u32::MAX
        } else {
            prev[m]
        };
        for j in 1..=n {
            let mut cur = vec![0u32; m + 1];
            cur[0] = if mode == Mode::Infix { 0 } else { j as u32 };
            for i in 1..=m {
                let sub = prev[i - 1] + (p[i - 1] != t[j - 1]) as u32;
                cur[i] = sub.min(prev[i] + 1).min(cur[i - 1] + 1);
            }
            match mode {
                Mode::Global => {
                    if j == n {
                        best = cur[m]
                    }
                }
                _ => best = best.min(cur[m]),
            }
            prev = cur;
        }
        best
    }

    /// Runs `f` at every SIMD level this CPU supports, each forced explicitly.
    fn each_level(mut f: impl FnMut(Level) -> Result<(), TestCaseError>) -> Result<(), TestCaseError> {
        for l in simd::supported() {
            at_level(l, || f(l))?;
        }
        Ok(())
    }

    /// One row by the row kernel of `level` (which the CPU must support), called directly.
    #[allow(clippy::too_many_arguments)]
    fn row_at(
        level: Level,
        hprev: &[i32],
        fu: &[i32],
        tx: &[u8],
        nh: &mut [i32],
        nf: &mut [i32],
        c: u8,
        sc: (i32, i32, i32, i32),
        e: i32,
        floor: i32,
        scratch: (&mut Vec<i32>, &mut Vec<i32>),
    ) -> i32 {
        assert!(simd::detected().includes(level));
        match level {
            // SAFETY: the CPU supports the level (asserted above).
            #[cfg(target_arch = "x86_64")]
            Level::Avx512 => unsafe { row_avx512(hprev, fu, tx, nh, nf, c, sc, e, floor) },
            #[cfg(target_arch = "x86_64")]
            Level::Avx2 => unsafe { row_avx2(hprev, fu, tx, nh, nf, c, sc, e, floor) },
            #[cfg(target_arch = "aarch64")]
            Level::Neon => unsafe { row_neon(hprev, fu, tx, nh, nf, c, sc, e, floor) },
            _ => row_scalar(hprev, fu, tx, nh, nf, c, sc, e, floor, scratch),
        }
    }

    /// One row through every row kernel this CPU supports: (level, H, F, max).
    fn rows_at_each_level(
        hprev: &[i32],
        fu: &[i32],
        tx: &[u8],
        c: u8,
        sc: (i32, i32, i32, i32),
        e: i32,
        floor: i32,
    ) -> Vec<(Level, Vec<i32>, Vec<i32>, i32)> {
        let w = fu.len();
        let mut out = Vec::new();
        for l in simd::supported() {
            let (mut nh, mut nf) = (vec![7; w + 1], vec![7; w + 1]);
            let (mut a, mut b) = (Vec::new(), Vec::new());
            let m = row_at(
                l,
                hprev,
                fu,
                tx,
                &mut nh,
                &mut nf,
                c,
                sc,
                e,
                floor,
                (&mut a, &mut b),
            );
            assert_eq!((nh[w], nf[w]), (NEG, NEG), "sentinel at {l}");
            out.push((l, nh, nf, m));
        }
        out
    }

    fn check(p: &[u8], t: &[u8], k: u32) -> Result<(), TestCaseError> {
        for mode in [Mode::Infix, Mode::Prefix, Mode::Global] {
            let want = naive(p, t, mode);
            let got = distance(p, t, mode, k);
            prop_assert_eq!(got, (want <= k).then_some(want), "mode {:?} k {}", mode, k);
        }
        Ok(())
    }

    #[test]
    fn examples() {
        assert_eq!(infix_distance(&[0, 1, 2, 3], &[3, 3, 0, 1, 2, 3, 3]), 0);
        assert_eq!(infix_distance(&[0, 1, 2, 3], &[0, 1, 3]), 1);
        assert_eq!(infix_distance(&[0, 1, 2, 3], &[]), 4);
        assert_eq!(
            distance(&[0, 1, 2, 3], &[3, 3, 0, 1, 2, 3, 3], Mode::Prefix, 10),
            Some(2)
        );
        assert_eq!(
            distance(&[0, 1, 2, 3], &[3, 3, 0, 1, 2, 3, 3], Mode::Global, 10),
            Some(3)
        );
        assert_eq!(
            distance(&[0, 1, 2, 3], &[3, 3, 0, 1, 2, 3, 3], Mode::Global, 2),
            None
        );
        assert_eq!(distance(&[], &[0, 1], Mode::Global, 5), Some(2));
    }

    fn naive_affine(p: &[u8], t: &[u8], mode: Mode, sc: &Affine) -> i32 {
        let (m, n) = (p.len(), t.len());
        let go = sc.open + sc.ext;
        let mut h = vec![vec![NEG; n + 1]; m + 1];
        let mut e = vec![vec![NEG; n + 1]; m + 1];
        let mut f = vec![vec![NEG; n + 1]; m + 1];
        h[0][0] = 0;
        for j in 1..=n {
            h[0][j] = if mode == Mode::Infix {
                0
            } else {
                -(sc.open + sc.ext * j as i32)
            };
            e[0][j] = h[0][j];
        }
        for i in 1..=m {
            h[i][0] = -(sc.open + sc.ext * i as i32);
            f[i][0] = h[i][0];
            for j in 1..=n {
                e[i][j] = (h[i][j - 1] - go).max(e[i][j - 1] - sc.ext);
                f[i][j] = (h[i - 1][j] - go).max(f[i - 1][j] - sc.ext);
                let s = if p[i - 1] == t[j - 1] { sc.mat } else { -sc.mis };
                h[i][j] = (h[i - 1][j - 1] + s).max(e[i][j]).max(f[i][j]);
            }
        }
        match mode {
            Mode::Global => h[m][n],
            _ => *h[m].iter().max().unwrap(),
        }
    }

    fn naive_extend(p: &[u8], t: &[u8], sc: &Affine) -> i32 {
        (0..=p.len())
            .flat_map(|i| (0..=t.len()).map(move |j| (i, j)))
            .map(|(i, j)| naive_affine(&p[..i], &t[..j], Mode::Global, sc))
            .max()
            .unwrap()
    }

    /// [`local_band`] by its definition: the best global alignment of a piece of the read
    /// with a piece of the text, plus the bonus per read end reached (no band).
    fn naive_local(p: &[u8], t: &[u8], sc: &Affine, bonus: i32) -> i32 {
        let mut best = 0;
        for i0 in 0..p.len() {
            for i1 in i0 + 1..=p.len() {
                let ends = (i0 == 0) as i32 + (i1 == p.len()) as i32;
                for j0 in 0..=t.len() {
                    for j1 in j0..=t.len() {
                        let x = naive_affine(&p[i0..i1], &t[j0..j1], Mode::Global, sc);
                        best = best.max(x + bonus * ends);
                    }
                }
            }
        }
        best
    }

    /// [`local_band`] over whole matrices: cells off the diagonals `dlo..=dhi` do not exist.
    fn naive_local_band(p: &[u8], t: &[u8], sc: &Affine, dlo: i64, dhi: i64, bonus: i32) -> i32 {
        let (m, n) = (p.len(), t.len());
        let go = sc.open + sc.ext;
        let inside = |i: usize, j: usize| (dlo..=dhi).contains(&(j as i64 - i as i64));
        let mut h = vec![vec![NEG; n + 1]; m + 1];
        let mut e = vec![vec![NEG; n + 1]; m + 1];
        let mut f = vec![vec![NEG; n + 1]; m + 1];
        for (j, x) in h[0].iter_mut().enumerate() {
            if inside(0, j) {
                *x = bonus;
            }
        }
        let mut best = 0;
        for i in 1..=m {
            for j in (0..=n).filter(|&j| inside(i, j)) {
                f[i][j] = (h[i - 1][j] - go).max(f[i - 1][j] - sc.ext);
                let mut x = f[i][j];
                if j > 0 {
                    e[i][j] = (h[i][j - 1] - go).max(e[i][j - 1] - sc.ext);
                    let s = if p[i - 1] == t[j - 1] { sc.mat } else { -sc.mis };
                    x = x.max(e[i][j]).max(h[i - 1][j - 1] + s);
                }
                h[i][j] = if i < m { x.max(0) } else { x };
                best = best.max(h[i][j] + if i == m { bonus } else { 0 });
            }
        }
        best
    }

    #[test]
    fn local_band_examples() {
        let sc = Affine::default();
        let t = [3, 3, 0, 1, 2, 3, 0, 1, 2, 2, 2];
        let p = [0, 1, 2, 3, 0, 1];
        // The read lies at text position 2: found on that diagonal, with both end bonuses.
        assert_eq!(local_band(&p, &t, &sc, 2, 2, 5), 12 + 10);
        assert_eq!(local_band(&p, &t, &sc, -6, 11, 5), 12 + 10);
        // Off it, only what matches there: nothing on the main diagonal, three bases at 6.
        assert_eq!(local_band(&p, &t, &sc, 0, 0, 0), 0);
        assert_eq!(local_band(&p, &t, &sc, 6, 6, 0), 6);
        // Junk at both read ends is clipped, and earns no bonus.
        let junk = [2, 2, 2, 0, 1, 2, 3, 0, 1, 3, 3, 3];
        assert_eq!(local_band(&junk, &t, &sc, -12, 11, 5), 12);
        // An ambiguous read base matches nothing.
        assert_eq!(local_band(&[4, 4, 4], &[0, 0, 0], &sc, -3, 3, 0), 0);
        assert_eq!(local_band(&[], &t, &sc, -3, 3, 5), 0);
        assert_eq!(local_band(&p, &t, &sc, 3, 2, 5), 0);
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(300))]

        /// With every diagonal in the band, `local_band` is the best piece-to-piece alignment.
        #[test]
        fn local_band_matches_definition(
            p in prop::collection::vec(0u8..4, 0..8), t in prop::collection::vec(0u8..4, 0..8),
            bonus in 0i32..8,
        ) {
            each_level(|l| {
                for sc in [Affine::default(), Affine::HIFI] {
                    let (m, n) = (p.len() as i64, t.len() as i64);
                    prop_assert_eq!(local_band(&p, &t, &sc, -m, n, bonus), naive_local(&p, &t, &sc, bonus), "at {}", l);
                }
                Ok(())
            })?;
        }

        /// Within a band, it is the same recurrence over the cells the band keeps, wherever
        /// the band lies (also off the matrix); a wider band never scores less.
        #[test]
        fn local_band_matches_matrix(
            base in prop::collection::vec(0u8..4, 1..60),
            edits in prop::collection::vec((0usize..60, 0u8..4, 0u8..3), 0..10),
            cut in (0usize..60, 0usize..60),
            dlo in -70i64..70, width in 0i64..40, bonus in 0i32..8,
        ) {
            let mut p = base[cut.0 % base.len()..].to_vec();
            p.truncate(p.len().saturating_sub(cut.1 % 20));
            for (pos, b, kind) in edits {
                if p.is_empty() { break; }
                let i = pos % p.len();
                match kind { 0 => p[i] = b, 1 => p.insert(i, b), _ => { p.remove(i); } }
            }
            each_level(|l| {
                for sc in [Affine::default(), Affine::HIFI] {
                    let got = local_band(&p, &base, &sc, dlo, dlo + width, bonus);
                    prop_assert_eq!(got, naive_local_band(&p, &base, &sc, dlo, dlo + width, bonus), "at {}", l);
                    prop_assert!(got <= local_band(&p, &base, &sc, dlo - 3, dlo + width + 5, bonus));
                }
                Ok(())
            })?;
        }

        /// The same over bands of several vectors, on reads with ambiguous bases and junk ends.
        #[test]
        fn local_band_matches_matrix_wide(
            base in prop::collection::vec(0u8..4, 1..300),
            edits in prop::collection::vec((0usize..300, 0u8..5, 0u8..3), 0..30),
            cut in (0usize..300, 0usize..100),
            junk in (prop::collection::vec(0u8..4, 0..30), prop::collection::vec(0u8..4, 0..30)),
            dlo in -320i64..320, width in 0i64..200, bonus in 0i32..8,
        ) {
            let mut p = base[cut.0 % base.len()..].to_vec();
            p.truncate(p.len().saturating_sub(cut.1));
            for (pos, b, kind) in edits {
                if p.is_empty() { break; }
                let i = pos % p.len();
                match kind { 0 => p[i] = b, 1 => p.insert(i, b), _ => { p.remove(i); } }
            }
            let p = [junk.0, p, junk.1].concat();
            each_level(|l| {
                for sc in [Affine::default(), Affine::HIFI] {
                    let got = local_band(&p, &base, &sc, dlo, dlo + width, bonus);
                    prop_assert_eq!(got, naive_local_band(&p, &base, &sc, dlo, dlo + width, bonus), "at {}", l);
                }
                Ok(())
            })?;
        }
    }

    #[test]
    fn affine_extend_examples() {
        let sc = Affine::default();
        // Read with a junk tail: the matching prefix scores, the tail is clipped.
        let t = [0, 1, 2, 3, 0, 1, 2, 3];
        let p = [0, 1, 2, 3, 0, 1, 3, 3, 3, 3, 3];
        assert_eq!(
            affine_extend(&p, &t, &sc, usize::MAX, 0, i32::MAX),
            naive_extend(&p, &t, &sc)
        );
        assert_eq!(affine_extend(&[], &t, &sc, usize::MAX, 0, i32::MAX), 0);
    }

    #[test]
    fn affine_extend_end_examples() {
        let sc = Affine::default();
        let end = |p: &[u8], t: &[u8], bonus| {
            let (x, i, j, _) = affine_extend_end(p, t, &sc, 16, bonus, 50);
            (x, i, j)
        };
        // The read's last 5 bases are junk: the extension ends after 6 read and text bases.
        let t = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        assert_eq!(end(&[0, 1, 2, 3, 0, 1, 3, 3, 3, 3, 3], &t, 0), (12, 6, 6));
        // A whole read with one deleted text base reaches its end at text base 8.
        assert_eq!(end(&[0, 1, 2, 0, 1, 2, 3], &t, 5), (14 - 6 + 5, 7, 8));
        // Nothing aligns: no bases.
        assert_eq!(end(&[3, 3, 3], &[0, 0, 0], 0), (0, 0, 0));
    }

    /// In a window, a text cut anywhere after the bases the DP read gives the same result, and
    /// a text it read to the end is matched only at that length.
    #[test]
    fn affine_extend_window_cuts() {
        let sc = Affine::default();
        let t: Vec<u8> = (0..140).map(|i| ((i * 5 + i / 4) % 4) as u8).collect();
        let mut p = t[..70].to_vec();
        p.extend((0..40).map(|i| [3, 3, 0, 0][i % 4]));
        let r = affine_extend_window(&p, &t, t.len(), &sc, 16, 5, 30);
        assert!(r.3 < t.len(), "{r:?}");
        for cut in [r.3, r.3 + 1, t.len() - 1] {
            assert_eq!(
                affine_extend_window(&p, &t[..cut], t.len(), &sc, 16, 5, 30).0,
                r.0
            );
        }
        // Cut inside the alignment: the DP reads the whole (shorter) text and stops there.
        let short = affine_extend_window(&p, &t[..40], t.len(), &sc, 16, 5, 30);
        assert_eq!((short.2, short.3), (40, 40));
    }

    /// The result depends on the text only through its length and the prefix the DP read.
    #[test]
    fn affine_extend_end_reads_a_prefix() {
        let sc = Affine::default();
        let t: Vec<u8> = (0..120).map(|i| ((i * 7 + i / 3) % 4) as u8).collect();
        let mut p = t[..60].to_vec();
        // A 40-base junk tail: the X-drop stops the DP before it reads the whole text.
        p.extend((0..40).map(|i| [3, 3, 0, 0][i % 4]));
        let r = affine_extend_end(&p, &t, &sc, 16, 5, 30);
        assert!(r.3 < t.len(), "{r:?}");
        let mut u = t.clone();
        for x in &mut u[r.3..] {
            *x = 3 - *x;
        }
        assert_eq!(affine_extend_end(&p, &u, &sc, 16, 5, 30), r);
    }

    #[test]
    fn affine_examples() {
        let sc = Affine::default();
        let p = [0, 1, 2, 3, 0, 1, 2, 3];
        assert_eq!(affine(&p, &p, Mode::Global, &sc, usize::MAX), 16);
        // One 2-base deletion from the text: 6 matches, one gap of 2.
        let t = [0, 1, 2, 3, 0, 1, 2, 3, 2, 2];
        assert_eq!(affine(&p, &t, Mode::Global, &sc, usize::MAX), 16 - 8);
        assert_eq!(affine(&p, &t, Mode::Prefix, &sc, usize::MAX), 16);
        assert_eq!(affine(&[], &[0, 1], Mode::Global, &sc, usize::MAX), -8);
        assert_eq!(affine(&[0, 1], &[], Mode::Global, &sc, usize::MAX), -8);
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(600))]

        /// The scalar passes' prefix maxima equal the definition.
        #[test]
        fn prefix_max_matches_scalar(
            xs in prop::collection::vec(-100000i32..100000, 0..200), init in -100000i32..100000,
        ) {
            let want: Vec<i32> = (0..xs.len()).map(|k| xs[..k].iter().fold(init, |a, &x| a.max(x))).collect();
            let mut got = xs.clone();
            excl_prefix_max(&mut got, init);
            prop_assert_eq!(got, want);
        }

        /// Every row kernel computes the scalar passes' row exactly, on rows of every length
        /// (whole vectors and partial tails), with cells and the entering gap down at and below
        /// NEG (unreachable cells and column 0 next to them), and with a floor under the cells
        /// (local alignment).
        #[test]
        fn row_kernels_match_scalar(
            cells in prop::collection::vec((prop_oneof![-3000i32..3000, (NEG - 300)..(NEG + 300)], prop_oneof![-3000i32..3000, (NEG - 300)..(NEG + 300)]), 1..90),
            tx in prop::collection::vec(any::<u8>(), 90),
            c in 0u8..4,
            sc in prop_oneof![Just((2, 4, 4, 2)), Just((1, 4, 6, 2)), (0i32..5, 0i32..9, 0i32..30, 0i32..6)],
            e in prop_oneof![Just(NEG), (NEG - 300)..(NEG + 300), -3000i32..3000],
            floor in prop_oneof![Just(NEG), Just(0), -50i32..50],
        ) {
            let hprev: Vec<i32> = cells.iter().map(|x| x.0).collect();
            let w = hprev.len() - 1;
            let fu: Vec<i32> = cells[1..].iter().map(|x| x.1).collect();
            let (mat, mis, open, ext) = sc;
            let sc = (mat, mis, open + ext, ext);
            let all = rows_at_each_level(&hprev, &fu, &tx[..w], c, sc, e, floor);
            let want = &all[0];
            prop_assert_eq!(want.0, Level::Scalar);
            for got in &all[1..] {
                prop_assert_eq!(&got.1, &want.1, "H at {}", got.0);
                prop_assert_eq!(&got.2, &want.2, "F at {}", got.0);
                prop_assert_eq!(got.3, want.3, "max at {}", got.0);
            }
        }

        /// Sorted texts sharing prefixes (isoform flanks) give the same scores through the
        /// row-sharing batch as one at a time, with and without X-drop stops and bands.
        #[test]
        fn sorted_batch_matches_single(
            base in prop::collection::vec(0u8..4, 1..160),
            edits in prop::collection::vec((0usize..160, 0u8..4, 0u8..3), 0..12),
            variants in prop::collection::vec((0usize..160, prop::collection::vec(0u8..4, 0..40), 0usize..60), 1..8),
            band in 2usize..24, xdrop in 5i32..60, bonus in 0i32..8, window_extra in 0usize..40,
        ) {
            let mut p = base.clone();
            for (pos, b, kind) in edits {
                if p.is_empty() { break; }
                let i = pos % p.len();
                match kind { 0 => p[i] = b, 1 => p.insert(i, b), _ => { p.remove(i); } }
            }
            // Each variant keeps a prefix of the base text, then diverges and is cut to a length.
            let mut texts: Vec<Vec<u8>> = vec![base.clone()];
            for (keep, tail, cut) in variants {
                let mut t = base[..keep.min(base.len())].to_vec();
                t.extend(tail);
                t.extend_from_slice(&base[keep.min(base.len())..]);
                t.truncate(t.len().saturating_sub(cut % 20));
                texts.push(t);
            }
            texts.sort();
            texts.dedup();
            let refs: Vec<&[u8]> = texts.iter().map(|t| t.as_slice()).collect();
            each_level(|_| {
            for sc in [Affine::default(), Affine::HIFI] {
                for mode in [Mode::Global, Mode::Prefix] {
                    let got = affine_sorted(&p, &refs, mode, &sc, band);
                    let want: Vec<i32> = refs.iter().map(|t| affine(&p, t, mode, &sc, band)).collect();
                    prop_assert_eq!(got, want, "mode {:?}", mode);
                }
                for diag in [false, true] {
                    let got = affine_extend_sorted(&p, &refs, &sc, band, bonus, xdrop, diag, 0);
                    let want: Vec<i32> = refs
                        .iter()
                        .map(|t| if diag {
                            affine_extend_diagonal(&p, t, &sc, band, bonus, xdrop)
                        } else {
                            affine_extend(&p, t, &sc, band, bonus, xdrop)
                        })
                        .collect();
                    prop_assert_eq!(got, want, "diag {}", diag);
                }
                // A band laid out for a window at least as long as each text, as in
                // affine_extend_window.
                let w = refs.iter().map(|t| t.len()).max().unwrap_or(0) + window_extra;
                let got = affine_extend_sorted(&p, &refs, &sc, band, bonus, xdrop, false, w);
                let want: Vec<i32> = refs
                    .iter()
                    .map(|t| affine_extend_window(&p, t, w, &sc, band, bonus, xdrop).0)
                    .collect();
                prop_assert_eq!(got, want, "window {}", w);
                // Scores over each text cut at `c`, from the same pass, equal passes over the cut
                // texts (no X-drop, so neither pass stops early).
                let c = refs.iter().map(|t| t.len()).min().unwrap_or(0) + window_extra / 2;
                let both = affine_extend_sorted_cut(&p, &refs, &sc, band, bonus, i32::MAX, false, w, c);
                let want: Vec<i32> = refs
                    .iter()
                    .map(|t| affine_extend_window(&p, &t[..t.len().min(c)], w, &sc, band, bonus, i32::MAX).0)
                    .collect();
                prop_assert_eq!(both.iter().map(|x| x.1).collect::<Vec<_>>(), want, "cut {}", c);
                // Tracking the cut changes nothing about the whole texts' scores.
                let whole = affine_extend_sorted(&p, &refs, &sc, band, bonus, i32::MAX, false, w);
                prop_assert_eq!(both.iter().map(|x| x.0).collect::<Vec<_>>(), whole, "whole with cut");
            }
            Ok(())
            })?;
        }

        #[test]
        fn affine_vs_naive(
            t in prop::collection::vec(0u8..4, 0..200),
            edits in prop::collection::vec((0usize..200, 0u8..4, 0u8..3), 0..20),
            a in 0usize..30, len in 0usize..180,
        ) {
            let end = (a + len).min(t.len());
            let mut p: Vec<u8> = t[a.min(end)..end].to_vec();
            for (pos, base, kind) in edits {
                if p.is_empty() { break; }
                let i = pos % p.len();
                match kind { 0 => p[i] = base, 1 => p.insert(i, base), _ => { p.remove(i); } }
            }
            let sc = Affine::default();
            let wants: Vec<i32> = [Mode::Global, Mode::Prefix, Mode::Infix].iter().map(|&m| naive_affine(&p, &t, m, &sc)).collect();
            let small = p.len() <= 40 && t.len() <= 60;
            let want_ext = if small { naive_extend(&p, &t, &sc) } else { 0 };
            each_level(|level| {
                for (mode, &want) in [Mode::Global, Mode::Prefix, Mode::Infix].into_iter().zip(&wants) {
                    prop_assert_eq!(affine(&p, &t, mode, &sc, usize::MAX), want, "mode {:?} at {}", mode, level);
                    // A band never improves the score, and a band covering everything is exact.
                    prop_assert!(affine(&p, &t, mode, &sc, 8) <= want);
                    prop_assert_eq!(affine(&p, &t, mode, &sc, 400), want);
                }
                if small {
                    prop_assert_eq!(affine_extend(&p, &t, &sc, usize::MAX, 0, i32::MAX), want_ext, "at {}", level);
                }
                Ok(())
            })?;
        }

        #[test]
        fn random_vs_naive(
            p in prop::collection::vec(0u8..4, 0..300),
            t in prop::collection::vec(0u8..4, 0..300),
            k in prop_oneof![Just(u32::MAX), 0u32..400],
        ) {
            check(&p, &t, k)?;
        }

        /// Related sequences (the realistic case): a mutated copy of a substring, with
        /// thresholds around the true distance.
        #[test]
        fn related_vs_naive(
            t in prop::collection::vec(0u8..4, 100..500),
            edits in prop::collection::vec((0usize..500, 0u8..4, 0u8..3), 0..30),
            a in 0usize..50, len in 60usize..450, slack in 0u32..20,
        ) {
            let end = (a + len).min(t.len());
            let mut p: Vec<u8> = t[a.min(end)..end].to_vec();
            for (pos, base, kind) in edits {
                if p.is_empty() { break; }
                let i = pos % p.len();
                match kind { 0 => p[i] = base, 1 => p.insert(i, base), _ => { p.remove(i); } }
            }
            let d = naive(&p, &t, Mode::Infix);
            for k in [d.saturating_sub(1), d, d + slack, u32::MAX] {
                check(&p, &t, k)?;
            }
        }

        /// Batched distances equal one `distance_with` call per text, in every mode: texts are
        /// mutated copies of one sequence (as isoform regions are) of differing lengths,
        /// some empty or unrelated, in batches of any size.
        #[test]
        fn batch_distances_match_single(
            base in prop::collection::vec(0u8..4, 0..700),
            texts in prop::collection::vec(
                (0usize..700, 0usize..700, prop::collection::vec((0usize..700, 0u8..4, 0u8..3), 0..25), any::<bool>()),
                1..20,
            ),
            pa in 0usize..300, plen in 0usize..500,
            pedits in prop::collection::vec((0usize..700, 0u8..4, 0u8..3), 0..40),
            kdelta in prop_oneof![Just(None), (0u32..60).prop_map(Some)],
            rnd in prop::collection::vec(0u8..4, 0..300),
        ) {
            let mutate = |mut v: Vec<u8>, edits: &[(usize, u8, u8)]| {
                for &(pos, b, kind) in edits {
                    if v.is_empty() { break; }
                    let i = pos % v.len();
                    match kind { 0 => v[i] = b, 1 => v.insert(i, b), _ => { v.remove(i); } }
                }
                v
            };
            let end = (pa + plen).min(base.len());
            let p = mutate(base[pa.min(end)..end].to_vec(), &pedits);
            let ts: Vec<Vec<u8>> = texts
                .iter()
                .map(|(a, l, e, unrelated)| {
                    if *unrelated {
                        rnd[..(*l).min(rnd.len())].to_vec()
                    } else {
                        let end = (a + l).min(base.len());
                        mutate(base[(*a).min(end)..end].to_vec(), e)
                    }
                })
                .collect();
            let pat = Pattern::new(&p);
            let refs: Vec<&[u8]> = ts.iter().map(|t| t.as_slice()).collect();
            // The lockstep lanes at any width (on every target, not only where they are used).
            fn lanes<const L: usize>(pat: &Pattern, refs: &[&[u8]], mode: Mode, k: u32) -> Vec<Option<u32>> {
                let mut out = Vec::new();
                refs.chunks(L).for_each(|c| batch_lanes::<L>(pat, c, mode, k, &mut out));
                out
            }
            for mode in [Mode::Infix, Mode::Prefix, Mode::Global] {
                let exact: Vec<u32> = ts.iter().map(|t| distance_with(&pat, t, mode, u32::MAX).unwrap()).collect();
                let k = match kdelta { None => u32::MAX, Some(x) => exact.iter().min().unwrap() + x };
                let want: Vec<Option<u32>> = ts.iter().map(|t| distance_with(&pat, t, mode, k)).collect();
                each_level(|level| {
                    prop_assert_eq!(batch_distances(&pat, &refs, mode, k), want.clone(), "at {} {:?}", level, mode);
                    Ok(())
                })?;
                if !p.is_empty() {
                    prop_assert_eq!(&lanes::<1>(&pat, &refs, mode, k), &want, "{:?}", mode);
                    prop_assert_eq!(&lanes::<2>(&pat, &refs, mode, k), &want, "{:?}", mode);
                    prop_assert_eq!(&lanes::<4>(&pat, &refs, mode, k), &want, "{:?}", mode);
                    prop_assert_eq!(&lanes::<8>(&pat, &refs, mode, k), &want, "{:?}", mode);
                }
            }
        }
    }

    #[cfg(target_arch = "x86_64")]
    proptest! {
        #![proptest_config(ProptestConfig::with_cases(400))]

        /// A pass on 16-bit cells finds what the pass on 32-bit cells does: over bands of one
        /// line of cells to several, texts sharing rows, cuts, X-drop stops with and without
        /// the best cell's position, and scorings up to what 16-bit cells hold.
        #[test]
        fn narrow_cells_match_wide(
            seed in 1u64.., n in 1usize..700, div in 2u64..30,
            variants in prop::collection::vec((0usize..700, prop::collection::vec(0u8..4, 0..60), 0usize..80), 0..4),
            band in 1usize..140, xdrop in 0i32..300, window_extra in 0usize..200, cut in 0usize..700,
            scale in 1i32..10,
        ) {
            if !simd::supported().contains(&Level::Avx512) {
                return Ok(());
            }
            let (p, t) = rand_pair(n, div, seed);
            let mut texts = vec![t.clone()];
            for (keep, tail, drop) in variants {
                let mut v = t[..keep.min(t.len())].to_vec();
                v.extend(tail);
                v.extend_from_slice(&t[keep.min(t.len())..]);
                v.truncate(v.len().saturating_sub(drop));
                texts.push(v);
            }
            texts.sort();
            texts.dedup();
            let refs: Vec<&[u8]> = texts.iter().map(|t| t.as_slice()).collect();
            let sc = Affine { mat: 2 * scale, mis: 4 * scale, open: 4 * scale, ext: 2 * scale };
            let (m, w) = (p.len(), refs.iter().map(|t| t.len()).max().unwrap_or(0) + window_extra);
            let all = |wide: bool| {
                WIDE_ONLY.with(|x| x.set(wide));
                let mut out = Vec::new();
                let mut keep = |ds: Vec<Dp>| {
                    out.extend(ds.into_iter().map(|d| {
                        (d.score, d.score_col, d.clip_best, d.clip_cell, d.text_read, d.last_row, d.cut_best, d.cut_score)
                    }))
                };
                for mode in [Mode::Global, Mode::Prefix] {
                    keep(affine_dp_sorted(&p, &refs, mode, &sc, band, None, |n| n, usize::MAX));
                }
                for (xdrop, cell) in [(xdrop * scale, false), (xdrop * scale, true), (i32::MAX, false)] {
                    for cut in [usize::MAX, cut] {
                        let clip = Some((xdrop, cell));
                        keep(affine_dp_sorted(&p, &refs, Mode::Prefix, &sc, band, clip, |n| m.min(n), cut));
                        keep(affine_dp_sorted(&p, &refs, Mode::Prefix, &sc, band, clip, |n| n.max(w), cut));
                    }
                }
                WIDE_ONLY.with(|x| x.set(false));
                out
            };
            let (narrow, wide) = at_level(Level::Avx512, || (all(false), all(true)));
            prop_assert_eq!(narrow, wide);
        }
    }

    /// Pseudo-random 2-bit codes, and a copy with about one edit in `div`.
    fn rand_pair(n: usize, div: u64, seed: u64) -> (Vec<u8>, Vec<u8>) {
        let mut x = seed;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let t: Vec<u8> = (0..n).map(|_| (next() & 3) as u8).collect();
        let mut p = Vec::with_capacity(n);
        for &b in &t {
            match next() % (3 * div) {
                0 => p.push((b + 1) & 3),
                1 => {}
                2 => p.extend([b, (next() & 3) as u8]),
                _ => p.push(b),
            }
        }
        (p, t)
    }

    /// Times `f` (best of 5) in ns per `per` units.
    fn time(per: f64, mut f: impl FnMut()) -> f64 {
        (0..5)
            .map(|_| {
                let t = std::time::Instant::now();
                f();
                t.elapsed().as_nanos() as f64 / per
            })
            .fold(f64::INFINITY, f64::min)
    }

    /// Isolated kernel timings per SIMD level (`cargo test --release -p spearfish-map -- --ignored
    /// --nocapture kernel_timings`).
    #[test]
    #[ignore]
    fn kernel_timings() {
        use std::hint::black_box;
        // One DP row kernel, called directly on rows of a few widths.
        for w in [24usize, 63, 129, 400] {
            let (_, t) = rand_pair(w + 1, 10, w as u64);
            let hprev: Vec<i32> = (0..=w as i32).map(|j| 3 * j - 200).collect();
            let fu: Vec<i32> = (0..w as i32).map(|j| 2 * j - 300).collect();
            let sc = (2, 4, 6, 2);
            let reps = 2_000_000 / w;
            for l in simd::supported() {
                let (mut nh, mut nf) = (vec![0; w + 1], vec![0; w + 1]);
                let (mut a, mut b) = (Vec::new(), Vec::new());
                let ns = time((reps * w) as f64, || {
                    for _ in 0..reps {
                        let m = row_at(
                            l,
                            black_box(&hprev),
                            &fu,
                            &t[..w],
                            &mut nh,
                            &mut nf,
                            1,
                            sc,
                            NEG,
                            NEG,
                            (&mut a, &mut b),
                        );
                        black_box(m);
                    }
                });
                println!("row w={w:4} {:7}: {ns:.3} ns/cell", l.name());
            }
        }
        // Whole banded extension DPs (the row loop around the kernel).
        let (p, t) = rand_pair(3000, 12, 7);
        for l in simd::supported() {
            let ns = time(10.0, || {
                at_level(l, || {
                    for _ in 0..10 {
                        black_box(affine_extend(black_box(&p), &t, &Affine::ONT, 64, 0, 400));
                    }
                })
            });
            println!(
                "affine_extend 3 kb band 64 {:7}: {:.1} us/DP",
                l.name(),
                ns / 1000.0
            );
        }
        // Batched infix distances: a 600-base pattern against 16 related 900-base texts.
        let (p, base) = rand_pair(600, 15, 11);
        let texts: Vec<Vec<u8>> = (0..16u64)
            .map(|s| {
                let mut t = rand_pair(150, 4, s + 100).1;
                t.extend(
                    rand_pair(600, 12, s + 200)
                        .0
                        .iter()
                        .zip(&base)
                        .map(|(&a, &b)| if s % 3 == 0 { a } else { b }),
                );
                t.extend(rand_pair(150, 4, s + 300).1);
                t
            })
            .collect();
        let refs: Vec<&[u8]> = texts.iter().map(|t| t.as_slice()).collect();
        let pat = Pattern::new(&p);
        for k in [60u32, u32::MAX] {
            for l in simd::supported() {
                let ns = time(20.0 * 16.0, || {
                    at_level(l, || {
                        for _ in 0..20 {
                            black_box(batch_distances(&pat, black_box(&refs), Mode::Infix, k));
                        }
                    })
                });
                println!("infix k={k:10} {:7}: {:.1} us/text", l.name(), ns / 1000.0);
            }
        }
    }
}
