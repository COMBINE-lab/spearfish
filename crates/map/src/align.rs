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

/// What one DP pass found.
struct Dp {
    /// Score in the pass's mode, and the text column where it ends (Prefix mode).
    score: i32,
    score_col: usize,
    /// With `clip`, the best cell anywhere and its (pattern, text) position.
    clip_best: i32,
    clip_cell: (usize, usize),
    /// Text bases the pass read (a prefix of the text).
    text_read: usize,
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
}

thread_local! {
    static DP_BUF: std::cell::RefCell<DpBuf> = std::cell::RefCell::new(DpBuf::default());
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
        let b = &mut *b.borrow_mut();
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
            let (lo0, hi0) = range(0);
            for j in lo0..=hi0 {
                b.h[j] = if j == 0 || mode == Mode::Infix {
                    0
                } else {
                    -(sc.open + sc.ext * j as i32)
                };
            }
            let mut prev = (lo0, hi0); // band of the row in h
            let mut stale = (1, 0); // band of the row left in nh (empty)
            for i in 1..=m {
                let (lo, hi) = range(i);
                if lo > hi {
                    // Past the text's end by more than the band (only when geom > n): no
                    // cell is reachable, so the rest of the pattern can't align.
                    return Dp {
                        score: NEG,
                        score_col: 0,
                        clip_best,
                        clip_cell,
                        text_read: n,
                    };
                }
                b.nh[stale.0..=stale.1.min(n)].fill(NEG);
                b.nf[stale.0..=stale.1.min(n)].fill(NEG);
                let c = pattern[i - 1] & 3;
                let mut start = lo;
                // Horizontal gap (consumes text) ending at (i, start - 1): from column
                // 0 if it is in the band, else unreachable.
                let mut e = NEG;
                if lo == 0 {
                    // Column 0: only a vertical gap reaches it.
                    let fj = (b.h[0] - go).max(b.f[0] - ext);
                    b.nh[0] = fj;
                    b.nf[0] = fj;
                    e = fj - go;
                    start = 1;
                }
                if start <= hi {
                    // Pass 1 (vectorizes): vertical gaps and diagonal moves (H' = H
                    // without horizontal gaps), from the previous row only.
                    let w = hi + 1 - start;
                    if b.hp.len() < w {
                        b.hp.resize(w, NEG);
                        b.gs.resize(w, NEG);
                    }
                    let DpBuf { h, f, nh, nf, hp, gs } = &mut *b;
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
                    let mut g = e;
                    for x in gs.iter_mut() {
                        let next = g.max(*x);
                        *x = g;
                        g = next;
                    }
                    let nhr = &mut nh[start..=hi];
                    for k in 0..w {
                        nhr[k] = hp[k].max(gs[k] - ext * k as i32);
                    }
                }
                if let Some((xdrop, cell)) = clip {
                    let row = b.nh[lo..=hi].iter().copied().max().unwrap_or(NEG);
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
                        };
                    }
                }
                std::mem::swap(&mut b.h, &mut b.nh);
                std::mem::swap(&mut b.f, &mut b.nf);
                stale = prev;
                prev = (lo, hi);
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
            }
        }
        run(pattern, text, mode, sc, (band, geom), clip, b)
    })
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
