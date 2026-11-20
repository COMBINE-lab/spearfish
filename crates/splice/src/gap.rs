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
    donor_pen,
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
    /// Work counters: solves, prefix-table cells filled, window cells filled.
    pub calls: u64,
    pub prefix_cells: u64,
    pub window_cells: u64,
}

impl Bufs {
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
    let band = (4 + m / 10).min(24);
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
        let leave0: Vec<u32> = (0..=n).map(|j| gv.leave(ix, j as u32, minus)).collect();
        for i in 0..=m {
            let mut end = (INF as u32, 0u32);
            let (j_lo, j_hi) = (i.saturating_sub(band), (i + band).min(n));
            for j in j_lo..=j_hi {
                let e = b.e1[i * w + j + band - i] as u32;
                if e >= INF as u32 {
                    continue;
                }
                let c = e + leave0[j];
                if c < b.l0_cur[i] {
                    b.l0_cur[i] = c;
                    b.l0_p[i] = j as u32;
                }
                if e < end.0 || (e == end.0 && j.abs_diff(i) < (end.1 as usize).abs_diff(i)) {
                    end = (e, j as u32);
                }
            }
            if end.0 < INF as u32 {
                let c = end.0 + CLIP * (m - i) as u32;
                if c < b.l0_end.0 {
                    b.l0_end = (c, i, end.1);
                }
            }
        }
        b.l0_key = Some((key, minus));
    }
    let mut cur = b.l0_cur.clone();
    let mut curp = b.l0_p.clone();
    let mut best_end = (BIG, 0usize, 0usize, 0u32, None::<(u32, u32)>); // cost, level, i, end p, entry
    if far.is_none() && b.l0_end.0 < BIG {
        best_end = (b.l0_end.0, 0, b.l0_end.1, b.l0_end.2, None);
    }
    // Candidate exons, one level each. A level keeps, per read position, the cheapest way to
    // have left the genome there, through this window or not.
    let mut lev_p: Vec<Vec<u32>> = vec![curp.clone()];
    let mut via: Vec<Vec<Option<(u32, u32)>>> = Vec::new();
    for &(vs, vl, i0, band) in windows {
        let vl = vl as usize;
        // Only cells near the diagonal the window's seed hits lie on: window position j goes
        // with read position i0 + j, give or take `wb`.
        let wb = if band > 0 { band } else { 14 + vl as i64 / 24 };
        let ww = (2 * wb + 1) as usize;
        let i_lo = (i0 - wb).clamp(0, m as i64) as usize;
        let i_hi = (i0 + vl as i64 + wb).clamp(0, m as i64) as usize;
        let rows = i_hi - i_lo + 1;
        b.window_cells += (rows * ww) as u64;
        let enter: Vec<u32> = (0..=vl)
            .map(|j| EXTRA_INTRON + gv.enter(ix, vs + j as u32, minus))
            .collect();
        let leave: Vec<u32> = (0..=vl).map(|j| gv.leave(ix, vs + j as u32, minus)).collect();
        let gb: Vec<u8> = (0..vl).map(|j| gv.base(ix, vs + j as u32)).collect();
        let mut next = cur.clone();
        let mut nextp = curp.clone();
        let mut v: Vec<Option<(u32, u32)>> = vec![None; m + 1];
        // Affine gaps, as outside the windows. Only two rows are kept: every cell carries where
        // its path entered the window (read position and window position, packed), so no table
        // is walked back. `uc`/`uo` hold, per k, the best cost of a cell of the row above that
        // ends in read bases left out.
        let open_del = (GAP_OPEN + GAP_EXT) as u32;
        let open_ins = (GAP_OPEN + INS_EXT) as u32;
        let (ins_ext, del_ext, mismatch) = (INS_EXT as u32, GAP_EXT as u32, MISMATCH as u32);
        let (mut pc, mut po) = (vec![BIG; ww], vec![0u32; ww]);
        let (mut cc, mut co) = (vec![BIG; ww], vec![0u32; ww]);
        let (mut uc, mut uo) = (vec![BIG; ww], vec![0u32; ww]);
        for i in i_lo..=i_hi {
            let ri = if i > 0 { rv.base(read, i - 1) } else { 0 };
            // window position j = base_j + k
            let base_j = i as i64 - i0 - wb;
            let k_lo = (-base_j).max(0) as usize;
            let k_hi = (vl as i64 - base_j).min(ww as i64 - 1);
            cc.fill(BIG);
            if k_hi < k_lo as i64 {
                uc.fill(BIG);
                std::mem::swap(&mut pc, &mut cc);
                std::mem::swap(&mut po, &mut co);
                continue;
            }
            let k_hi = k_hi as usize;
            // entering here is allowed where the genome was left far enough before
            let entry_ok = cur[i] < BIG;
            let j_min = (curp[i] + min_gap).saturating_sub(vs) as usize;
            let first_row = i == i_lo;
            let (mut left_c, mut left_o) = (BIG, 0u32);
            let (mut bl, mut blj, mut blo) = (next[i], usize::MAX, 0u32);
            let (mut be, mut bej, mut beo) = (BIG, 0usize, 0u32);
            for k in k_lo..=k_hi {
                let j = (base_j + k as i64) as usize;
                let (mut c, mut o) = if entry_ok && j >= j_min {
                    (cur[i] + enter[j], ((i as u32) << 16) | j as u32)
                } else {
                    (BIG, 0)
                };
                // (i - 1, j) is k + 1 in the row above
                let (fc, fo) = if !first_row && k + 1 < ww {
                    let (ext, opn) = (uc[k + 1] + ins_ext, pc[k + 1] + open_ins);
                    if ext <= opn {
                        (ext, uo[k + 1])
                    } else {
                        (opn, po[k + 1])
                    }
                } else {
                    (BIG, 0)
                };
                if k > k_lo {
                    let (ext, opn) = (left_c + del_ext, cc[k - 1] + open_del);
                    if ext <= opn {
                        left_c = ext;
                    } else {
                        left_c = opn;
                        left_o = co[k - 1];
                    }
                } else {
                    left_c = BIG;
                }
                // (i - 1, j - 1) is the same k in the row above
                if !first_row && j > 0 {
                    let d = pc[k] + if ri != gb[j - 1] { mismatch } else { 0 };
                    if d < c {
                        c = d;
                        o = po[k];
                    }
                }
                if fc < c {
                    c = fc;
                    o = fo;
                }
                if left_c < c {
                    c = left_c;
                    o = left_o;
                }
                let c = c.min(BIG);
                left_c = left_c.min(BIG);
                cc[k] = c;
                co[k] = o;
                uc[k] = fc.min(BIG);
                uo[k] = fo;
                if c < BIG && j as u32 > (o & 0xffff) {
                    // something of the window is used: it can be left here, or the read can end
                    let l = c + leave[j];
                    if l < bl {
                        (bl, blj, blo) = (l, j, o);
                    }
                    if c < be {
                        (be, bej, beo) = (c, j, o);
                    }
                }
            }
            uc[..k_lo].fill(BIG);
            uc[k_hi + 1..].fill(BIG);
            if blj != usize::MAX {
                next[i] = bl;
                nextp[i] = vs + blj as u32;
                v[i] = Some((blo >> 16, vs + (blo & 0xffff)));
            }
            if far.is_none() && be < BIG {
                let e = be + CLIP * (m - i) as u32;
                if e < best_end.0 {
                    best_end = (
                        e,
                        via.len() + 1,
                        i,
                        vs + bej as u32,
                        Some((beo >> 16, vs + (beo & 0xffff))),
                    );
                }
            }
            std::mem::swap(&mut pc, &mut cc);
            std::mem::swap(&mut po, &mut co);
        }
        cur = next;
        curp = nextp;
        lev_p.push(curp.clone());
        via.push(v);
    }
    // Finish, then walk back through the levels.
    let mut out = Solved {
        cost: BIG,
        ..Default::default()
    };
    let (mut level, mut i);
    let mut hops: Vec<(u32, u32, u32)> = Vec::new(); // read i, leave p, enter (real)
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
        // The far side per number of read bases it takes: its cheapest arrival, cached like
        // level 0. The pairing with level 0 then costs one addition per read position; only
        // when that arrival would overlap the place the genome was left is the row searched.
        if b.r_key != Some((key, minus)) {
            let far_leave: Vec<u32> = (0..=n).map(|j| fv.leave(ix, j as u32, minus)).collect();
            b.r_cost.clear();
            b.r_cost.resize(m + 1, BIG);
            b.r_j.clear();
            b.r_j.resize(m + 1, 0);
            b.r_leave = far_leave;
            for rest in 0..=m {
                let (j_lo, j_hi) = (rest.saturating_sub(band), (rest + band).min(n));
                for j in j_lo..=j_hi {
                    let e = b.e2[rest * w + j + band - rest] as u32;
                    if e >= INF as u32 {
                        continue;
                    }
                    let c = e + b.r_leave[j];
                    if c < b.r_cost[rest] {
                        b.r_cost[rest] = c;
                        b.r_j[rest] = j as u32;
                    }
                }
            }
            b.r_key = Some((key, minus));
        }
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
        if best.0 >= BIG {
            return out;
        }
        out.cost = best.0;
        level = via.len();
        i = best.1;
        hops.push((i as u32, gv.real(lev_p[level][i]), best.2));
    } else {
        if best_end.0 >= BIG {
            return out;
        }
        out.cost = best_end.0;
        out.end = (rv.real(best_end.2 as u32), gv.real(best_end.3));
        level = best_end.1;
        i = best_end.2;
        if let Some((i1, ep)) = best_end.4 {
            level -= 1;
            hops.push((i1, gv.real(lev_p[level][i1 as usize]), gv.real(ep)));
            i = i1 as usize;
        }
    }
    while level > 0 {
        if let Some((i1, ep)) = via[level - 1][i] {
            hops.push((i1, gv.real(lev_p[level - 1][i1 as usize]), gv.real(ep)));
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
    out
}
