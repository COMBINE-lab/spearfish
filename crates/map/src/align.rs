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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Infix,
    Prefix,
    Global,
}

const W: usize = 64;

/// One Myers block step; returns the horizontal delta out of row `bit`.
#[inline(always)]
fn block(pv: &mut u64, mv: &mut u64, mut eq: u64, hin: i32, bit: usize) -> i32 {
    let (p, m) = (*pv, *mv);
    let xv = eq | m;
    if hin < 0 {
        eq |= 1;
    }
    let xh = ((eq & p).wrapping_add(p) ^ p) | eq;
    let mut ph = m | !(xh | p);
    let mut mh = p & xh;
    let hout = ((ph >> bit) & 1) as i32 - ((mh >> bit) & 1) as i32;
    ph <<= 1;
    mh <<= 1;
    if hin < 0 {
        mh |= 1;
    } else if hin > 0 {
        ph |= 1;
    }
    *pv = mh | !(xv | ph);
    *mv = ph & xv;
    hout
}

/// Edit distance in `mode`, or `None` if it exceeds `k`.
pub fn distance(pattern: &[u8], text: &[u8], mode: Mode, k: u32) -> Option<u32> {
    let m = pattern.len();
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
    let mut peq = vec![[0u64; 4]; blocks];
    for (i, &c) in pattern.iter().enumerate() {
        peq[i / W][(c & 3) as usize] |= 1 << (i % W);
    }
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

/// Minimum edit distance between `pattern` and any substring of `text` (no cutoff).
pub fn infix_distance(pattern: &[u8], text: &[u8]) -> u32 {
    distance(pattern, text, Mode::Infix, u32::MAX).unwrap()
}

/// Affine-gap scoring (maximized): `+mat` per match, `-mis` per mismatch, `-(open + ext·l)`
/// per gap of length `l`. Defaults are minimap2's map-ont single-affine part.
#[derive(Clone, Copy, Debug)]
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

/// [`affine`] for each of `texts`, which must be sorted: equal to separate calls, but rows
/// that read only a prefix shared with the previous text are not recomputed.
pub fn affine_sorted(pattern: &[u8], texts: &[&[u8]], mode: Mode, sc: &Affine, band: usize) -> Vec<i32> {
    affine_dp_sorted(pattern, texts, mode, sc, band, None, |n| n)
        .into_iter()
        .map(|d| d.score)
        .collect()
}

/// [`affine_extend`] (or, with `diagonal`, [`affine_extend_diagonal`]) for each of `texts`,
/// which must be sorted, sharing rows as [`affine_sorted`] does.
pub fn affine_extend_sorted(
    pattern: &[u8],
    texts: &[&[u8]],
    sc: &Affine,
    band: usize,
    end_bonus: i32,
    xdrop: i32,
    diagonal: bool,
) -> Vec<i32> {
    let m = pattern.len();
    let geom = |n: usize| if diagonal { m.min(n) } else { n };
    affine_dp_sorted(pattern, texts, Mode::Prefix, sc, band, Some((xdrop, false)), geom)
        .into_iter()
        .map(|d| d.clip_best.max(d.score.saturating_add(end_bonus)))
        .collect()
}

/// [`affine_extend`] with the band on the diagonal from the start (`2 * band + 1` cells per
/// row) however much longer the text is: for accurate reads, whose extensions stay near it.
/// A text shorter than the pattern (a transcript ending inside the read's flank) keeps the
/// band of [`affine_extend`], which can place the read's extra bases.
pub fn affine_extend_diagonal(
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

/// [`affine_extend`] with where the best extension ends: (score, pattern bases, text bases,
/// text bases read).
pub fn affine_extend_end(
    pattern: &[u8],
    text: &[u8],
    sc: &Affine,
    band: usize,
    end_bonus: i32,
    xdrop: i32,
) -> (i32, usize, usize, usize) {
    affine_extend_window(pattern, text, text.len(), sc, band, end_bonus, xdrop)
}

/// [`affine_extend_end`] on a text cut from a window of `window` bases (`text.len() <=
/// window`): the band is laid out for the whole window, so the result depends on the text only
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
/// maximum of `init` and `xs[..k]`. The scalar loop is one dependent max per element; with
/// AVX-512 a 16-lane block is scanned in four shift-and-max steps and the block maximum carried.
#[inline]
fn excl_prefix_max(xs: &mut [i32], init: i32) {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx512f") {
        // SAFETY: AVX-512F is available (checked above).
        unsafe { excl_prefix_max_avx512(xs, init) };
        return;
    }
    let mut g = init;
    for x in xs.iter_mut() {
        let next = g.max(*x);
        *x = g;
        g = next;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn excl_prefix_max_avx512(xs: &mut [i32], init: i32) {
    use std::arch::x86_64::*;
    let neg = _mm512_set1_epi32(NEG);
    let last = _mm512_set1_epi32(15);
    // The running maximum, broadcast to every lane (kept in a register between blocks).
    let mut cv = _mm512_set1_epi32(init);
    let (blocks, rest) = xs.as_chunks_mut::<16>();
    for c in blocks {
        // SAFETY: `c` holds 16 i32s; unaligned loads and stores.
        let mut v = unsafe { _mm512_loadu_si512(c.as_ptr().cast()) };
        // Inclusive prefix max within the block: lane k takes the max of lanes <= k.
        v = _mm512_max_epi32(v, _mm512_alignr_epi32::<15>(v, neg));
        v = _mm512_max_epi32(v, _mm512_alignr_epi32::<14>(v, neg));
        v = _mm512_max_epi32(v, _mm512_alignr_epi32::<12>(v, neg));
        v = _mm512_max_epi32(v, _mm512_alignr_epi32::<8>(v, neg));
        v = _mm512_max_epi32(v, cv);
        // Exclusive: shift up one lane, the carry entering lane 0.
        let out = _mm512_alignr_epi32::<15>(v, cv);
        unsafe { _mm512_storeu_si512(c.as_mut_ptr().cast(), out) };
        // Lane 15 holds the block's (inclusive) maximum.
        cv = _mm512_permutexvar_epi32(last, v);
    }
    let mut carry = _mm512_cvtsi512_si32(cv);
    for x in rest {
        let next = carry.max(*x);
        *x = carry;
        carry = next;
    }
}

/// First and last columns of `h` within `band` scoring at least `thr` ((1, 0) if none).
#[inline]
fn live_cols(h: &[i32], band: (usize, usize), thr: i32) -> (usize, usize) {
    let row = &h[band.0..=band.1];
    match (
        row.iter().position(|&x| x >= thr),
        row.iter().rposition(|&x| x >= thr),
    ) {
        (Some(a), Some(z)) => (band.0 + a, band.0 + z),
        _ => (1, 0),
    }
}

/// Experimental X-drop band pruning, on with `SPEARFISH_XDROP_PRUNE` in the environment.
fn xdrop_prune() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("SPEARFISH_XDROP_PRUNE").is_some())
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
    clip_best: i32,
    clip_cell: (usize, usize),
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
    DP_BUF.with(|b| {
        run(
            pattern,
            text,
            mode,
            sc,
            (band, geom),
            clip,
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
fn affine_dp_sorted(
    pattern: &[u8],
    texts: &[&[u8]],
    mode: Mode,
    sc: &Affine,
    band: usize,
    clip: Option<(i32, bool)>,
    geom: impl Fn(usize) -> usize,
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
                    pattern,
                    t,
                    mode,
                    sc,
                    (band, geom(t.len())),
                    clip,
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
                    pattern,
                    t,
                    mode,
                    sc,
                    (band, geom(t.len())),
                    clip,
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

#[allow(clippy::too_many_arguments)]
#[multiversion::multiversion(targets(
    "x86_64+avx512f+avx512bw+avx512vl+avx2+fma",
    "x86_64+avx2+fma",
    "x86_64+avx",
    "aarch64+neon"
))]
fn run(
    pattern: &[u8],
    text: &[u8],
    mode: Mode,
    sc: &Affine,
    (band, geom): (usize, usize),
    clip: Option<(i32, bool)>,
    b: &mut DpBuf,
    resume: Option<&Snap>,
    snaps: &mut [Snap],
) -> Dp {
    let (m, n) = (pattern.len(), text.len());
    // With `clip`, the best cell anywhere (0 aligns nothing).
    let mut clip_best = 0i32;
    let mut clip_cell = (0usize, 0usize);
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
    // (vertical). Rows alternate between two buffers; cells outside a row's
    // band stay NEG.
    for v in [&mut b.h, &mut b.f, &mut b.nh, &mut b.nf] {
        v.clear();
        v.resize(n + 1, NEG);
    }
    debug_assert!(sc.open >= 0 && sc.ext >= 0);
    let (mut prev, first) = match resume {
        Some(s) => {
            b.h[s.lo..=s.hi].copy_from_slice(&s.h);
            b.f[s.lo..=s.hi].copy_from_slice(&s.f);
            clip_best = s.clip_best;
            clip_cell = s.clip_cell;
            ((s.lo, s.hi), s.row + 1)
        }
        None => {
            let (lo0, hi0) = range(0);
            for j in lo0..=hi0 {
                b.h[j] = if j == 0 || mode == Mode::Infix {
                    0
                } else {
                    -(sc.open + sc.ext * j as i32)
                };
            }
            ((lo0, hi0), 1)
        }
    }; // band of the row in h
    let mut stale = (1, 0); // band of the row left in nh (empty)
    let mut next_snap = snaps.partition_point(|s| s.row < first);
    // X-drop pruning (experimental): cells more than `xdrop` below the best are dead, and a
    // row's band shrinks to the columns its predecessor's live cells can reach (one diagonal
    // step plus as far as a horizontal gap stays above the threshold).
    let prune = !full && xdrop_prune() && clip.is_some();
    let (xd, reach) = match clip {
        Some((x, _)) => (x, ((x - go).max(0) / ext.max(1) + 1) as usize),
        None => (0, 0),
    };
    let mut live = if prune {
        live_cols(&b.h, prev, clip_best.saturating_sub(xd))
    } else {
        prev
    };
    for i in first..=m {
        let (mut lo, mut hi) = range(i);
        if prune && lo <= hi {
            if live.0 > live.1 {
                return Dp {
                    score: NEG,
                    score_col: 0,
                    clip_best,
                    clip_cell,
                    text_read: prev.1,
                    last_row: i,
                };
            }
            // Dead edges of the previous row propagate nothing.
            b.h[prev.0..live.0].fill(NEG);
            b.f[prev.0..live.0].fill(NEG);
            if live.1 < prev.1 {
                b.h[live.1 + 1..=prev.1].fill(NEG);
                b.f[live.1 + 1..=prev.1].fill(NEG);
            }
            lo = lo.max(live.0);
            hi = hi.min(live.1 + 1 + reach);
            if lo > hi {
                return Dp {
                    score: NEG,
                    score_col: 0,
                    clip_best,
                    clip_cell,
                    text_read: prev.1,
                    last_row: i,
                };
            }
        }
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
            };
        }
        b.nh[stale.0..=stale.1.min(n)].fill(NEG);
        b.nf[stale.0..=stale.1.min(n)].fill(NEG);
        let c = pattern[i - 1] & 3;
        let mut start = lo;
        // Horizontal gap (consumes text) ending at (i, start - 1): from column
        // 0 if it is in the band, else unreachable.
        let mut e = NEG;
        // The row's maximum, kept as the row is written (for clipping).
        let mut row_max = NEG;
        if lo == 0 {
            // Column 0: only a vertical gap reaches it.
            let fj = (b.h[0] - go).max(b.f[0] - ext);
            b.nh[0] = fj;
            b.nf[0] = fj;
            e = fj - go;
            start = 1;
            row_max = fj;
        }
        if crate::stats::on() {
            DP_WORK.with(|x| {
                let (c, r) = x.get();
                x.set((c + (hi + 1 - lo) as u64, r + 1));
            });
        }
        if start <= hi {
            // Pass 1 (vectorizes): vertical gaps and diagonal moves (H' = H
            // without horizontal gaps), from the previous row only.
            let w = hi + 1 - start;
            if b.hp.len() < w {
                b.hp.resize(w, NEG);
                b.gs.resize(w, NEG);
            }
            let DpBuf {
                h, f, nh, nf, hp, gs, ..
            } = &mut *b;
            let (hp, gs) = (&mut hp[..w], &mut gs[..w]);
            let hd = &h[start - 1..hi]; // H[i-1][j-1]
            let hu = &h[start..=hi]; // H[i-1][j]
            let fu = &f[start..=hi];
            let tx = &text[start - 1..hi];
            let nfr = &mut nf[start..=hi];
            for k in 0..w {
                let fj = (hu[k] - go).max(fu[k] - ext).max(NEG);
                let s = if c == tx[k] & 3 { mat } else { -mis };
                let x = (hd[k] + s).max(fj).max(NEG);
                hp[k] = x;
                gs[k] = x - go + ext * (k as i32 + 1);
                nfr[k] = fj;
            }
            // Horizontal gaps. With open >= 0, E[k] = max(H'[k-1] - go, E[k-1]
            // - ext): opening from H = max(H', E) adds nothing. With g[k] =
            // E[k] + ext*k this is a running max (pass 2: exclusive prefix
            // maxima, one op per cell on the critical path); pass 3
            // (vectorizes) takes H = max(H', E).
            excl_prefix_max(gs, e);
            let nhr = &mut nh[start..=hi];
            let mut m = NEG;
            for k in 0..w {
                let v = hp[k].max(gs[k] - ext * k as i32);
                nhr[k] = v;
                m = m.max(v);
            }
            row_max = row_max.max(m);
        }
        if let Some((xdrop, cell)) = clip {
            let row = row_max;
            if row > clip_best {
                clip_best = row;
                if cell {
                    let j = b.nh[lo..=hi].iter().position(|&x| x == row).unwrap();
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
                };
            }
            if prune {
                live = live_cols(&b.nh, (lo, hi), clip_best.saturating_sub(xdrop));
            }
        }
        std::mem::swap(&mut b.h, &mut b.nh);
        std::mem::swap(&mut b.f, &mut b.nf);
        stale = prev;
        prev = (lo, hi);
        if next_snap < snaps.len() && snaps[next_snap].row == i {
            let s = &mut snaps[next_snap];
            (s.lo, s.hi) = (lo, hi);
            s.h.clear();
            s.h.extend_from_slice(&b.h[lo..=hi]);
            s.f.clear();
            s.f.extend_from_slice(&b.f[lo..=hi]);
            (s.clip_best, s.clip_cell) = (clip_best, clip_cell);
            next_snap += 1;
        }
    }
    let (score, score_col) = match mode {
        Mode::Global => (b.h[n], n),
        Mode::Prefix | Mode::Infix => {
            let row = &b.h[prev.0..=prev.1];
            let x = row.iter().copied().max().unwrap_or(NEG);
            (x, prev.0 + row.iter().position(|&y| y == x).unwrap_or(0))
        }
    };
    Dp {
        score,
        score_col,
        clip_best,
        clip_cell,
        text_read: prev.1,
        last_row: m,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

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

        /// The vectorized prefix maxima equal the scalar loop.
        #[test]
        fn prefix_max_matches_scalar(
            xs in prop::collection::vec(-100000i32..100000, 0..200), init in -100000i32..100000,
        ) {
            let mut want = xs.clone();
            let mut g = init;
            for x in want.iter_mut() {
                let next = g.max(*x);
                *x = g;
                g = next;
            }
            let mut got = xs.clone();
            excl_prefix_max(&mut got, init);
            prop_assert_eq!(got, want);
        }

        /// Sorted texts sharing prefixes (isoform flanks) give the same scores through the
        /// row-sharing batch as one at a time, with and without X-drop stops and bands.
        #[test]
        fn sorted_batch_matches_single(
            base in prop::collection::vec(0u8..4, 1..160),
            edits in prop::collection::vec((0usize..160, 0u8..4, 0u8..3), 0..12),
            variants in prop::collection::vec((0usize..160, prop::collection::vec(0u8..4, 0..40), 0usize..60), 1..8),
            band in 2usize..24, xdrop in 5i32..60, bonus in 0i32..8,
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
            for sc in [Affine::default(), Affine::HIFI] {
                for mode in [Mode::Global, Mode::Prefix] {
                    let got = affine_sorted(&p, &refs, mode, &sc, band);
                    let want: Vec<i32> = refs.iter().map(|t| affine(&p, t, mode, &sc, band)).collect();
                    prop_assert_eq!(got, want, "mode {:?}", mode);
                }
                for diag in [false, true] {
                    let got = affine_extend_sorted(&p, &refs, &sc, band, bonus, xdrop, diag);
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
            }
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
            for mode in [Mode::Global, Mode::Prefix, Mode::Infix] {
                let want = naive_affine(&p, &t, mode, &sc);
                prop_assert_eq!(affine(&p, &t, mode, &sc, usize::MAX), want, "mode {:?}", mode);
                // A band never improves the score, and a band covering everything is exact.
                prop_assert!(affine(&p, &t, mode, &sc, 8) <= want);
                prop_assert_eq!(affine(&p, &t, mode, &sc, 400), want);
            }
            if p.len() <= 40 && t.len() <= 60 {
                prop_assert_eq!(affine_extend(&p, &t[..t.len().min(60)], &sc, usize::MAX, 0, i32::MAX), naive_extend(&p, &t[..t.len().min(60)], &sc));
            }
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
    }
}
