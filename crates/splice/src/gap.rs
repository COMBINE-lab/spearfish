//! Explaining read bases between two anchors, or past the last one, with zero or more exons.
//!
//! The read bases start at a fixed point that corresponds to a genome position. They are
//! explained by continuing on the genome from there (the upstream exon), then optionally
//! passing through candidate exon windows, each entered after an intron and left before the
//! next one, and finally either arriving at a far fixed point (a gap between two anchors) or
//! stopping with the rest of the read clipped (a read end). Every read base is paid for:
//! aligned (banded affine costs), or clipped. Each
//! window used costs `EXTRA_INTRON` on top of its splice sites, so an exon is added only when
//! it explains more than it costs.
//!
//! A `View` reads the genome forwards or backwards from an origin, so that the read's left end
//! is solved by the same code as its right end.

use crate::index::GenomeIndex;
use crate::map::{
    GAP_EXT, GAP_OPEN, INF, INS_EXT, MISMATCH, acceptor_motif, acceptor_pen, band_prefix, donor_motif,
    donor_pen, site_pens,
};

/// Price of each exon added between two fixed points (besides its splice sites).
pub(crate) const EXTRA_INTRON: u32 = 20;
/// Price of a read base left unaligned at a read end (the match it forgoes).
pub(crate) const CLIP: u32 = 2;
const BIG: u32 = 1 << 28;

/// The genome read from `origin`, forwards or backwards. Position `p` of the view is the
/// boundary after its first `p` bases.
#[derive(Clone, Copy, Debug)]
pub(crate) struct View {
    pub origin: u32,
    pub rev: bool,
}

impl View {
    #[inline]
    pub fn base(&self, ix: &GenomeIndex, k: u32) -> u8 {
        if self.rev {
            ix.base(self.origin - 1 - k)
        } else {
            ix.base(self.origin + k)
        }
    }
    #[inline]
    pub fn real(&self, p: u32) -> u32 {
        if self.rev {
            self.origin - p
        } else {
            self.origin + p
        }
    }
    /// The 32 bases from view position `p` on, in the direction of travel, packed two bits
    /// each with the first base lowest; `None` too close to the genome's start to read them.
    #[inline]
    pub fn word(&self, ix: &GenomeIndex, p: u32) -> Option<u64> {
        if !self.rev {
            return Some(ix.word(self.origin + p));
        }
        let x = self.origin - p; // the bases x-1, x-2, ... in that order
        if x < 32 {
            return None;
        }
        // reverse the order of the 32 two-bit groups
        let mut w = ix.word(x - 32).swap_bytes();
        w = ((w >> 4) & 0x0f0f_0f0f_0f0f_0f0f) | ((w & 0x0f0f_0f0f_0f0f_0f0f) << 4);
        w = ((w >> 2) & 0x3333_3333_3333_3333) | ((w & 0x3333_3333_3333_3333) << 2);
        Some(w)
    }
    pub fn room(&self, ix: &GenomeIndex) -> u32 {
        if self.rev {
            self.origin
        } else {
            ix.total_len() - self.origin
        }
    }
    /// Cost of an intron starting at boundary `p` in the direction of travel.
    #[inline]
    fn leave(&self, ix: &GenomeIndex, p: u32, minus: bool) -> u32 {
        let x = self.real(p);
        (if self.rev {
            acceptor_pen(ix, x, minus)
        } else {
            donor_pen(ix, x, minus)
        }) as u32
    }
    /// Cost of an intron ending at boundary `p` in the direction of travel.
    /// The motif part of [`Self::enter`] alone.
    #[inline]
    pub fn enter_motif(&self, ix: &GenomeIndex, p: u32, minus: bool) -> u32 {
        let x = self.real(p);
        (if self.rev {
            donor_motif(ix, x, minus)
        } else {
            acceptor_motif(ix, x, minus)
        }) as u32
    }
    /// For the `len + 1` boundaries from `vs` on: the cost of an intron ending there plus
    /// `add`, the cost of one starting there, and (when wanted) the base before each (4
    /// before the first).
    #[allow(clippy::too_many_arguments)]
    fn sites(
        &self,
        ix: &GenomeIndex,
        vs: u32,
        len: usize,
        minus: bool,
        add: u32,
        s: &mut Sites,
        want_bases: bool,
    ) {
        let l = len as u32;
        let (lo, hi) = if self.rev {
            ((self.origin - vs).wrapping_sub(l), self.origin - vs)
        } else {
            (self.origin + vs, self.origin + vs + l)
        };
        s.enter.clear();
        s.leave.clear();
        s.gb.clear();
        let fits = !self.rev || self.origin - vs >= l;
        if fits && site_pens(ix, lo, hi, minus, &mut s.bases, &mut s.donors, &mut s.acceptors) {
            if self.rev {
                s.enter.extend(s.donors.iter().rev().map(|&c| c + add));
                s.leave.extend(s.acceptors.iter().rev());
            } else {
                s.enter.extend(s.acceptors.iter().map(|&c| c + add));
                s.leave.extend_from_slice(&s.donors);
            }
            if want_bases {
                s.gb.push(4);
                if self.rev {
                    s.gb.extend(s.bases[3..len + 3].iter().rev());
                } else {
                    s.gb.extend_from_slice(&s.bases[3..len + 3]);
                }
            }
            return;
        }
        s.enter
            .extend((0..=l).map(|j| add + self.enter(ix, vs + j, minus)));
        s.leave.extend((0..=l).map(|j| self.leave(ix, vs + j, minus)));
        if want_bases {
            s.gb.push(4);
            s.gb.extend((0..l).map(|j| self.base(ix, vs + j)));
        }
    }
    #[inline]
    pub fn enter(&self, ix: &GenomeIndex, p: u32, minus: bool) -> u32 {
        let x = self.real(p);
        (if self.rev {
            donor_pen(ix, x, minus)
        } else {
            acceptor_pen(ix, x, minus)
        }) as u32
    }
}

/// Splice-site costs and bases of a stretch of a view (see [`View::sites`]), and the buffers
/// they are worked out in.
#[derive(Default)]
pub(crate) struct Sites {
    enter: Vec<u32>,
    leave: Vec<u32>,
    gb: Vec<u8>,
    bases: Vec<u8>,
    donors: Vec<u32>,
    acceptors: Vec<u32>,
}

/// The read from position `q`, forwards or backwards.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ReadView {
    pub q: u32,
    pub rev: bool,
}

impl ReadView {
    #[inline]
    pub fn base(&self, read: &[u8], i: usize) -> u8 {
        if self.rev {
            read[self.q as usize - 1 - i]
        } else {
            read[self.q as usize + i]
        }
    }
    #[inline]
    pub fn real(&self, i: u32) -> u32 {
        if self.rev { self.q - i } else { self.q + i }
    }
}

/// One intron of a solution, in the direction of travel and in real coordinates: the read
/// boundary where it falls, where the genome is left and where it is entered again.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Hop {
    pub q: u32,
    pub leave: u32,
    pub enter: u32,
}

#[derive(Clone, Debug, Default)]
pub(crate) struct Solved {
    pub cost: u32,
    pub hops: Vec<Hop>,
    /// Read ends only: where the alignment stops (read boundary, genome boundary), real.
    pub end: (u32, u32),
}

/// The solver's working storage, kept between solves. `cur`, `curp`: per read position, the
/// cost of having left the genome there and where, for the level being built. `lev_p`: `curp`
/// as it was after each level, one after another. `via`: per level and read position, where
/// the level's window was entered (read position, view position) by the path that left the
/// genome there, or `NO_VIA`. `rows`: the rows of a window's table.
#[derive(Default)]
struct Work {
    cur: Vec<u32>,
    curp: Vec<u32>,
    lev_p: Vec<u32>,
    via: Vec<(u32, u32)>,
    rows: [Vec<u64>; 4],
}
const NO_VIA: (u32, u32) = (u32::MAX, 0);

/// What a pair of prefix tables was computed for: where the read and the genome start, the
/// direction, the length and the far point (`None` before any).
type Key = Option<(u32, bool, usize, u32, Option<u32>)>;

#[derive(Default)]
pub(crate) struct Bufs {
    e1: Vec<u16>,
    e2: Vec<u16>,
    tmp: Vec<u16>,
    /// What `e1` (and `e2`) currently hold: the tables depend only on where the read and the
    /// genome start, the direction, the length and the far point, not on the strand or the
    /// windows, so repeated solves of one gap reuse them.
    key1: Key,
    key2: Key,
    /// Level 0 and the far side, cached per tables and strand.
    l0_key: Option<(Key, bool)>,
    l0_cur: Vec<u32>,
    l0_p: Vec<u32>,
    l0_end: (u32, usize, u32),
    r_key: Option<(Key, bool)>,
    r_cost: Vec<u32>,
    r_j: Vec<u32>,
    r_leave: Vec<u32>,
    sites: Sites,
    work: Work,
    /// Work counters: solves, prefix-table cells filled, window cells filled.
    /// Set while the windows given are known exons (and the answer is kept only if every
    /// junction in it is a known one): an exon then costs nothing on top of its sites, so
    /// that a known exon of three bases beats writing them as an insertion.
    pub known_only: bool,
    pub calls: u64,
    pub prefix_cells: u64,
    pub window_cells: u64,
}

/// A copy of the prefix tables of one gap, to solve it again later without recomputing them.
#[derive(Clone, Default)]
pub(crate) struct Tables {
    key: Key,
    e1: Vec<u16>,
    e2: Vec<u16>,
}

impl Bufs {
    /// The tables of the gap solved last (both must be present: a gap between two anchors).
    pub fn stash(&self) -> Tables {
        debug_assert!(self.key1 == self.key2);
        Tables {
            key: self.key1,
            e1: self.e1.clone(),
            e2: self.e2.clone(),
        }
    }

    /// Put stashed tables back, unless they are the current ones.
    pub fn restore(&mut self, t: &Tables) {
        if t.key.is_some() && self.key1 != t.key {
            self.e1.clear();
            self.e1.extend_from_slice(&t.e1);
            self.e2.clear();
            self.e2.extend_from_slice(&t.e2);
            self.key1 = t.key;
            self.key2 = t.key;
            self.l0_key = None;
            self.r_key = None;
        }
    }

    /// Forget the cached tables (the read changed).
    pub fn reset(&mut self) {
        self.key1 = None;
        self.key2 = None;
        self.l0_key = None;
        self.r_key = None;
    }
}

/// Solve `m` read bases from `rv` against the genome from `gv`.
///
/// `far` is the real genome position of the far fixed point (the read bases then end exactly
/// there; `gv` and `rv` must be forward); without it the read end is free and unexplained
/// bases are clipped. `windows` are candidate exons as (view position, length, read offset
/// that goes with the window's start, half-width of the band around that diagonal or 0 for the
/// usual one), in travel order; they may overlap. `minus` is the transcript strand.
const SH: u32 = 32;
const BIGP: u64 = (BIG as u64) << SH;

/// One row of a candidate exon's table. A cell is one word: its cost above, and below it where
/// its path entered the window. `pd`, `pu` are the row above at the same and the next k, `ur`
/// its cells that end in read bases left out (next k); `hr`, `nu` receive this row. Returns the
/// cheapest way to leave the window in this row and the cheapest cell to end the read in, each
/// as cost above and offset in the row below (`u64::MAX`: none), and the row's lowest cost.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn window_row_body<const WIDE: bool>(
    pd: &[u64],
    pu: &[u64],
    ur: &[u64],
    hr: &mut [u64],
    nu: &mut [u64],
    en: &[u32],
    lv: &[u32],
    g: &[u8],
    ri: u8,
    entry: u64,
    j_lo: usize,
    j_min: usize,
) -> (u64, u64, u64) {
    let open_ins = ((GAP_OPEN + INS_EXT) as u64) << SH;
    let ins_ext = (INS_EXT as u64) << SH;
    let mismatch = (MISMATCH as u64) << SH;
    let len = hr.len();
    let (pd, pu, ur, nu) = (&pd[..len], &pu[..len], &ur[..len], &mut nu[..len]);
    let (en, lv, g) = (&en[..len], &lv[..len], &g[..len]);
    // every way into a cell but along the row
    for x in 0..len {
        let j = j_lo + x;
        // (selections are written as masks: these branches are not predictable)
        let may = ((j >= j_min) as u64).wrapping_neg();
        let e = ((entry + (((en[x] as u64) << SH) | j as u64)) & may) | (BIGP & !may);
        // (i - 1, j) is k + 1 in the row above
        let f = (ur[x] + ins_ext).min(pu[x] + open_ins).min(BIGP);
        // (i - 1, j - 1) is the same k in the row above
        let d = pd[x] + ((ri != g[x]) as u64) * mismatch;
        nu[x] = f;
        hr[x] = e.min(d).min(f).min(BIGP);
    }
    // along the row (`left` never exceeds BIGP + open_del, since a cell is at most BIGP)
    #[cfg(target_arch = "x86_64")]
    if WIDE {
        // SAFETY: WIDE is set only by the copy built (and called) with AVX-512.
        unsafe { row_gap_avx512(hr) };
    } else {
        row_gap(hr);
    }
    #[cfg(not(target_arch = "x86_64"))]
    row_gap(hr);
    // something of the window is used: it can be left here, or the read can end
    let (mut bl, mut be, mut low) = (u64::MAX, u64::MAX, u64::MAX);
    for x in 0..len {
        let c = hr[x];
        low = low.min(c);
        let unused = (((c < BIGP) & ((j_lo + x) as u64 > (c & 0xffff))) as u64).wrapping_sub(1);
        let cost = c >> SH;
        let l = (((cost + lv[x] as u64) << 16) | x as u64) | unused;
        let t = ((cost << 16) | x as u64) | unused;
        bl = bl.min(l);
        be = be.min(t);
    }
    (bl, be, low >> SH)
}

/// The gap along a row: every cell becomes the cheaper of itself and a gap from a cell before
/// it (opening costs GAP_OPEN + GAP_EXT, each further base GAP_EXT).
#[inline(always)]
fn row_gap(hr: &mut [u64]) {
    let open_del = ((GAP_OPEN + GAP_EXT) as u64) << SH;
    let del_ext = (GAP_EXT as u64) << SH;
    // (`left` never exceeds BIGP + open_del, since a cell is at most BIGP; a gap opened from
    // a cell that a gap just lowered never beats that gap extended, so the cell as it was is
    // enough here and the chain from cell to cell is one addition and one comparison)
    let mut left = BIGP;
    for c in hr.iter_mut() {
        let open = *c + open_del;
        *c = (*c).min(left);
        left = (left + del_ext).min(open);
    }
}

/// `row_gap`, eight cells at a time. The gap that reaches cell x from cell y < x costs
/// open + (x - 1 - y) * ext, so with a[y] = cell y + open + (len - 1 - y) * ext the cheapest
/// gap into x is the least a[y] before x, less (len - x) * ext: a running minimum, taken by
/// three doubling steps inside the eight and one carried word between them. The results are
/// those of `row_gap` exactly.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn row_gap_avx512(hr: &mut [u64]) {
    use std::arch::x86_64::*;
    let open_del = ((GAP_OPEN + GAP_EXT) as u64) << SH;
    let del_ext = (GAP_EXT as u64) << SH;
    let len = hr.len();
    // SAFETY (all of it): loads and stores are of whole eights inside `hr`, or masked to
    // the cells that remain.
    unsafe {
        let ramp = _mm512_mullo_epi64(
            _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0),
            _mm512_set1_epi64(del_ext as i64),
        );
        let none = _mm512_set1_epi64(-1);
        let last = _mm512_set1_epi64(7);
        let mut carry = _mm512_set1_epi64((BIGP + len as u64 * del_ext) as i64);
        let mut x = 0;
        while x < len {
            let rem = len - x;
            let k: __mmask8 = if rem >= 8 { 0xff } else { (1u8 << rem) - 1 };
            let at = hr.as_mut_ptr().add(x) as *mut i64;
            // (cells past the end read as "none"; what they turn into is never stored and
            // never reaches a cell before them)
            let c = _mm512_mask_loadu_epi64(none, k, at);
            let a0 = (open_del + (len - 1 - x) as u64 * del_ext) as i64;
            let mut p = _mm512_sub_epi64(_mm512_add_epi64(c, _mm512_set1_epi64(a0)), ramp);
            p = _mm512_min_epu64(p, _mm512_alignr_epi64::<7>(p, none));
            p = _mm512_min_epu64(p, _mm512_alignr_epi64::<6>(p, none));
            p = _mm512_min_epu64(p, _mm512_alignr_epi64::<4>(p, none));
            // (the running minimum includes the cell itself, as a gap of no length that
            // costs GAP_OPEN: never cheaper than the cell)
            let before = _mm512_min_epu64(carry, p);
            let back = _mm512_set1_epi64(((len - x) as u64 * del_ext) as i64);
            let left = _mm512_add_epi64(_mm512_sub_epi64(before, back), ramp);
            _mm512_mask_storeu_epi64(at, k, _mm512_min_epu64(c, left));
            carry = _mm512_min_epu64(carry, _mm512_permutexvar_epi64(last, p));
            x += 8;
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[allow(clippy::too_many_arguments)]
mod window_row_x86 {
    use super::*;
    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn avx2(
        pd: &[u64],
        pu: &[u64],
        ur: &[u64],
        hr: &mut [u64],
        nu: &mut [u64],
        en: &[u32],
        lv: &[u32],
        g: &[u8],
        ri: u8,
        entry: u64,
        j_lo: usize,
        j_min: usize,
    ) -> (u64, u64, u64) {
        window_row_body::<false>(pd, pu, ur, hr, nu, en, lv, g, ri, entry, j_lo, j_min)
    }
    #[target_feature(enable = "avx512f,avx512bw,avx512vl,avx512dq")]
    pub(super) unsafe fn avx512(
        pd: &[u64],
        pu: &[u64],
        ur: &[u64],
        hr: &mut [u64],
        nu: &mut [u64],
        en: &[u32],
        lv: &[u32],
        g: &[u8],
        ri: u8,
        entry: u64,
        j_lo: usize,
        j_min: usize,
    ) -> (u64, u64, u64) {
        window_row_body::<true>(pd, pu, ur, hr, nu, en, lv, g, ri, entry, j_lo, j_min)
    }
}

/// `window_row_body` built for the instruction set found when the program runs.
#[allow(clippy::too_many_arguments)]
#[inline]
fn window_row(
    pd: &[u64],
    pu: &[u64],
    ur: &[u64],
    hr: &mut [u64],
    nu: &mut [u64],
    en: &[u32],
    lv: &[u32],
    g: &[u8],
    ri: u8,
    entry: u64,
    j_lo: usize,
    j_min: usize,
) -> (u64, u64, u64) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx512bw")
            && std::is_x86_feature_detected!("avx512vl")
            && std::is_x86_feature_detected!("avx512dq")
        {
            // SAFETY: the features were just detected.
            return unsafe { window_row_x86::avx512(pd, pu, ur, hr, nu, en, lv, g, ri, entry, j_lo, j_min) };
        }
        if std::is_x86_feature_detected!("avx2") {
            // SAFETY: AVX2 was just detected.
            return unsafe { window_row_x86::avx2(pd, pu, ur, hr, nu, en, lv, g, ri, entry, j_lo, j_min) };
        }
    }
    window_row_body::<false>(pd, pu, ur, hr, nu, en, lv, g, ri, entry, j_lo, j_min)
}

/// The cheapest cell of a prefix-table row once the cost of leaving the genome there is
/// added: cost above, offset in the row in the low 8 bits (the first of equals);
/// `u32::MAX` when no cell is reachable.
#[inline(always)]
fn row_best_body(row: &[u16], lv: &[u32]) -> u32 {
    let len = row.len();
    let lv = &lv[..len];
    let mut best = u32::MAX;
    for x in 0..len {
        let e = row[x] as u32;
        let c = (((e + lv[x]) << 8) | x as u32) | ((e < INF as u32) as u32).wrapping_sub(1);
        best = best.min(c);
    }
    best
}
avx2_variant! {
    fn row_best = row_best_body, row_best_avx2(row: &[u16], lv: &[u32]) -> u32
}

/// The cheapest cell of a prefix-table row, the nearest to the diagonal of equals (`off0` is
/// the first cell's signed distance from it): cost above, then distance, then offset.
#[inline(always)]
fn row_end_body(row: &[u16], off0: i32) -> u32 {
    let mut end = u32::MAX;
    for (x, &e) in row.iter().enumerate() {
        let e = e as u32;
        let off = (off0 + x as i32).unsigned_abs();
        let c = ((e << 16) | (off << 8) | x as u32) | ((e < INF as u32) as u32).wrapping_sub(1);
        end = end.min(c);
    }
    end
}
avx2_variant! {
    pub(crate) fn row_end = row_end_body, row_end_avx2(row: &[u16], off0: i32) -> u32
}

/// The cheapest pairing of a level (`cur`, `curp`: per read position, the cost of having left
/// the genome and where) with the far side: cost, read position, where the genome is entered
/// (real). The far side's cheapest arrival per number of read bases is cached, so a pairing is
/// one addition; only when that arrival would overlap the place the genome was left is the
/// row searched.
#[allow(clippy::too_many_arguments)]
fn pair_far(
    b: &Bufs,
    cur: &[u32],
    curp: &[u32],
    gv: View,
    t1: u32,
    m: usize,
    n: usize,
    band: usize,
    min_gap: u32,
) -> (u32, usize, u32) {
    let w = 2 * band + 1;
    let mut best = (BIG, 0usize, 0u32);
    for (x, &c0) in cur.iter().enumerate() {
        let rest = m - x;
        if c0 >= BIG || b.r_cost[rest] >= BIG || c0 + b.r_cost[rest] >= best.0 {
            continue;
        }
        let enter = t1 - b.r_j[rest];
        if gv.real(curp[x]) + min_gap <= enter {
            best = (c0 + b.r_cost[rest], x, enter);
            continue;
        }
        let (j_lo, j_hi) = (rest.saturating_sub(band), (rest + band).min(n));
        for j in j_lo..=j_hi {
            let e = b.e2[rest * w + j + band - rest] as u32;
            if e >= INF as u32 {
                continue;
            }
            let c = c0 + e + b.r_leave[j];
            let enter = t1 - j as u32;
            if c < best.0 && gv.real(curp[x]) + min_gap <= enter {
                best = (c, x, enter);
            }
        }
    }
    best
}

#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
pub(crate) fn solve(
    ix: &GenomeIndex,
    read: &[u8],
    rv: ReadView,
    m: usize,
    gv: View,
    far: Option<u32>,
    windows: &[(u32, u32, i64, i64)],
    minus: bool,
    min_gap: u32,
    b: &mut Bufs,
) -> Solved {
    let band = (16 + m / 8).min(40);
    let w = 2 * band + 1;
    let limit = match far {
        Some(t1) => (t1 - gv.origin) as usize,
        None => gv.room(ix) as usize,
    };
    let n = (m + band).min(limit);
    let key = Some((rv.q, rv.rev, m, gv.origin, far));
    b.calls += 1;
    if b.key1 != key {
        b.prefix_cells += ((m + 1) * w) as u64;
        band_prefix(
            |i| rv.base(read, i),
            m,
            |j| gv.base(ix, j as u32),
            n,
            band,
            &mut b.e1,
            &mut b.tmp,
            0,
        );
        b.key1 = key;
    }
    // Level 0: the upstream exon alone. cur[i] = cost of read[..i] and of leaving at p[i].
    // It depends on the tables and the strand only, so repeated solves of one gap reuse it.
    if b.l0_key != Some((key, minus)) {
        b.l0_cur.clear();
        b.l0_cur.resize(m + 1, BIG);
        b.l0_p.clear();
        b.l0_p.resize(m + 1, 0);
        b.l0_end = (BIG, 0, 0);
        gv.sites(ix, 0, n, minus, 0, &mut b.sites, false);
        let leave0 = &b.sites.leave;
        // (per row: the cheapest cell with its leaving cost, the first of equals; and, for a
        // read end, the cheapest cell, the nearest to the diagonal of equals)
        let want_end = far.is_none();
        for i in 0..=m {
            let (j_lo, j_hi) = (i.saturating_sub(band), (i + band).min(n));
            if j_hi < j_lo {
                continue;
            }
            let len = j_hi - j_lo + 1;
            let row = &b.e1[i * w + j_lo + band - i..][..len];
            let lv = &leave0[j_lo..][..len];
            let best = row_best(row, lv);
            if best != u32::MAX {
                b.l0_cur[i] = best >> 8;
                b.l0_p[i] = (j_lo + (best & 0xff) as usize) as u32;
            }
            if want_end {
                let end = row_end(row, j_lo as i32 - i as i32);
                if end != u32::MAX {
                    let c = (end >> 16) + CLIP * (m - i) as u32;
                    if c < b.l0_end.0 {
                        b.l0_end = (c, i, (j_lo + (end & 0xff) as usize) as u32);
                    }
                }
            }
        }
        b.l0_key = Some((key, minus));
    }
    let mut wk = std::mem::take(&mut b.work);
    let Work {
        cur,
        curp,
        lev_p,
        via,
        rows: [prev, here, up, up2],
    } = &mut wk;
    cur.clear();
    cur.extend_from_slice(&b.l0_cur);
    curp.clear();
    curp.extend_from_slice(&b.l0_p);
    let m1 = m + 1;
    let mut best_end = (BIG, 0usize, 0usize, 0u32, None::<(u32, u32)>); // cost, level, i, end p, entry
    if far.is_none() && b.l0_end.0 < BIG {
        best_end = (b.l0_end.0, 0, b.l0_end.1, b.l0_end.2, None);
    }
    // The far side (its tables, cached like level 0), before the windows: the answer without
    // any window bounds what a window's cells may cost and still matter.
    let mut bound = BIG;
    if let Some(t1) = far {
        let fv = View {
            origin: t1,
            rev: true,
        };
        let rr = ReadView {
            q: rv.q + m as u32,
            rev: true,
        };
        if b.key2 != key {
            band_prefix(
                |x| rr.base(read, x),
                m,
                |j| fv.base(ix, j as u32),
                n,
                band,
                &mut b.e2,
                &mut b.tmp,
                0,
            );
            b.key2 = key;
        }
        // Per number of read bases it takes: its cheapest arrival, cached like
        // level 0. The pairing with level 0 then costs one addition per read position; only
        // when that arrival would overlap the place the genome was left is the row searched.
        if b.r_key != Some((key, minus)) {
            fv.sites(ix, 0, n, minus, 0, &mut b.sites, false);
            b.r_cost.clear();
            b.r_cost.resize(m + 1, BIG);
            b.r_j.clear();
            b.r_j.resize(m + 1, 0);
            b.r_leave.clear();
            b.r_leave.extend_from_slice(&b.sites.leave);
            for rest in 0..=m {
                let (j_lo, j_hi) = (rest.saturating_sub(band), (rest + band).min(n));
                if j_hi < j_lo {
                    continue;
                }
                let len = j_hi - j_lo + 1;
                let row = &b.e2[rest * w + j_lo + band - rest..][..len];
                let lv = &b.r_leave[j_lo..][..len];
                let best = row_best(row, lv);
                if best != u32::MAX {
                    b.r_cost[rest] = best >> 8;
                    b.r_j[rest] = (j_lo + (best & 0xff) as usize) as u32;
                }
            }
            b.r_key = Some((key, minus));
        }
        if !windows.is_empty() {
            bound = pair_far(b, cur, curp, gv, t1, m, n, band, min_gap).0;
        }
    }
    // Candidate exons, one level each. A level keeps, per read position, the cheapest way to
    // have left the genome there, through this window or not.
    lev_p.clear();
    lev_p.extend_from_slice(curp);
    via.clear();
    let mut levels = 0usize;
    for &(vs, vl, i0, band) in windows {
        let vl = vl as usize;
        // Only cells near the diagonal the window's seed hits lie on: window position j goes
        // with read position i0 + j, give or take `wb`.
        let wb = if band > 0 { band } else { 20 + vl as i64 / 16 };
        let ww = (2 * wb + 1) as usize;
        let i_lo = (i0 - wb).clamp(0, m as i64) as usize;
        let i_hi = (i0 + vl as i64 + wb).clamp(0, m as i64) as usize;
        // A cell that costs more than an answer already in hand cannot be part of a better
        // one, and costs only grow along a path. So a row matters only if the row above has
        // such a cell or the window can be entered in it; the rows before the first place
        // it can be entered and after the last, once nothing is left, are not filled.
        let ub = if far.is_some() { bound } else { best_end.0 };
        let extra = if b.known_only { 0 } else { EXTRA_INTRON };
        let enters = |c: u32, ub: u32| c < BIG && (ub >= BIG || c + extra <= ub);
        let first = (i_lo..=i_hi).find(|&i| enters(cur[i], ub));
        let Some(first) = first else {
            lev_p.extend_from_slice(curp);
            via.resize(via.len() + m1, NO_VIA);
            levels += 1;
            continue;
        };
        let last = (first..=i_hi).rev().find(|&i| enters(cur[i], ub)).unwrap();
        let rows = i_hi - first + 1;
        b.window_cells += (rows * ww) as u64;
        // (gb[j] is the genome base before window position j)
        gv.sites(ix, vs, vl, minus, extra, &mut b.sites, true);
        let (enter, leave, gb) = (&b.sites.enter, &b.sites.leave, &b.sites.gb);
        // (a row reads its own entry of `cur` and `curp` before it replaces them, and no
        // other row's, so the level is built in place)
        let v0 = via.len();
        via.resize(v0 + m1, NO_VIA);
        // Affine gaps, as outside the windows. Only two rows are kept. A cell is one word: its
        // cost above, and below it where its path entered the window (read position and window
        // position), so the cheapest of several ways is one comparison and no table is walked
        // back. `up` holds, per k, the best cell of the row above that ends in read bases left
        // out. All three rows have one spare cell at the end.
        for r in [&mut *prev, &mut *here, &mut *up, &mut *up2] {
            r.clear();
            r.resize(ww + 1, BIGP);
        }
        let mut alive = false;
        for i in first..=i_hi {
            let ub = if far.is_some() { bound } else { best_end.0 };
            let may_enter = enters(cur[i], ub);
            if !alive && !may_enter {
                // (the rows above are empty already)
                if i > last {
                    break;
                }
                continue;
            }
            let ri = if i > 0 { rv.base(read, i - 1) } else { 0 };
            // window position j = base_j + k
            let base_j = i as i64 - i0 - wb;
            let k_lo = (-base_j).max(0) as usize;
            let k_hi = (vl as i64 - base_j).min(ww as i64 - 1);
            here.fill(BIGP);
            if k_hi < k_lo as i64 {
                up.fill(BIGP);
                std::mem::swap(prev, here);
                alive = false;
                continue;
            }
            let k_hi = k_hi as usize;
            let j_lo = (base_j + k_lo as i64) as usize;
            let len = k_hi - k_lo + 1;
            // entering here is allowed where the genome was left far enough before
            let j_min = if may_enter {
                (curp[i] + min_gap).saturating_sub(vs) as usize
            } else {
                usize::MAX
            };
            let entry = ((cur[i] as u64) << SH) | ((i as u64) << 16);
            let (bl, be, low) = window_row(
                &prev[k_lo..k_lo + len],
                &prev[k_lo + 1..k_lo + 1 + len],
                &up[k_lo + 1..k_lo + 1 + len],
                &mut here[k_lo..k_lo + len],
                &mut up2[k_lo..k_lo + len],
                &enter[j_lo..j_lo + len],
                &leave[j_lo..j_lo + len],
                &gb[j_lo..j_lo + len],
                ri,
                entry,
                j_lo,
                j_min,
            );
            up2[..k_lo].fill(BIGP);
            up2[k_hi + 1..].fill(BIGP);
            std::mem::swap(up, up2);
            if bl != u64::MAX && ((bl >> 16) as u32) < cur[i] {
                let x = (bl & 0xffff) as usize;
                let o = here[k_lo + x] as u32;
                cur[i] = (bl >> 16) as u32;
                curp[i] = vs + (j_lo + x) as u32;
                via[v0 + i] = ((o >> 16) & 0xffff, vs + (o & 0xffff));
            }
            if far.is_none() && be != u64::MAX {
                let x = (be & 0xffff) as usize;
                let o = here[k_lo + x] as u32;
                let e = (be >> 16) as u32 + CLIP * (m - i) as u32;
                if e < best_end.0 {
                    best_end = (
                        e,
                        levels + 1,
                        i,
                        vs + (j_lo + x) as u32,
                        Some(((o >> 16) & 0xffff, vs + (o & 0xffff))),
                    );
                }
            }
            std::mem::swap(prev, here);
            alive = ub >= BIG || low <= ub as u64;
            if !alive {
                prev.fill(BIGP);
                up.fill(BIGP);
            }
        }
        lev_p.extend_from_slice(curp);
        levels += 1;
    }
    // Finish, then walk back through the levels.
    let mut out = Solved {
        cost: BIG,
        ..Default::default()
    };
    'done: {
        let (mut level, mut i);
        let mut hops: Vec<(u32, u32, u32)> = Vec::new(); // read i, leave p, enter (real)
        if let Some(t1) = far {
            let best = pair_far(b, cur, curp, gv, t1, m, n, band, min_gap);
            if best.0 >= BIG {
                break 'done;
            }
            out.cost = best.0;
            level = levels;
            i = best.1;
            hops.push((i as u32, gv.real(lev_p[level * m1 + i]), best.2));
        } else {
            if best_end.0 >= BIG {
                break 'done;
            }
            out.cost = best_end.0;
            out.end = (rv.real(best_end.2 as u32), gv.real(best_end.3));
            level = best_end.1;
            i = best_end.2;
            if let Some((i1, ep)) = best_end.4 {
                level -= 1;
                hops.push((i1, gv.real(lev_p[level * m1 + i1 as usize]), gv.real(ep)));
                i = i1 as usize;
            }
        }
        while level > 0 {
            let (i1, ep) = via[(level - 1) * m1 + i];
            if (i1, ep) != NO_VIA {
                hops.push((i1, gv.real(lev_p[(level - 1) * m1 + i1 as usize]), gv.real(ep)));
                i = i1 as usize;
            }
            level -= 1;
        }
        hops.reverse();
        out.hops = hops
            .into_iter()
            .map(|(i, leave, enter)| Hop {
                q: rv.real(i),
                leave,
                enter,
            })
            .collect();
    }
    b.work = wk;
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The gap pass of a window's row, eight cells at a time, against the cell-by-cell one.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn row_gaps_match() {
        if !std::is_x86_feature_detected!("avx512f") {
            return;
        }
        let mut st = 0x9E37_79B9_7F4A_7C15u64;
        let mut rng = move || {
            st ^= st << 13;
            st ^= st >> 7;
            st ^= st << 17;
            st
        };
        for round in 0..2000 {
            let len = 1 + (rng() % 70) as usize;
            let row: Vec<u64> = (0..len)
                .map(|_| match rng() % 8 {
                    0 => BIGP,
                    1 => ((rng() % (BIG as u64)) << SH) | (rng() & 0xffff_ffff),
                    _ => ((rng() % 40) << SH) | (rng() & 0xffff_ffff),
                })
                .collect();
            let (mut x, mut y) = (row.clone(), row.clone());
            row_gap(&mut x);
            // SAFETY: the feature was just detected.
            unsafe { row_gap_avx512(&mut y) };
            assert_eq!(x, y, "round {round} {row:?}");
        }
    }
}
