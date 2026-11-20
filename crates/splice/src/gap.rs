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
use crate::map::{GAP_EXT, GAP_OPEN, INF, INS_EXT, MISMATCH, acceptor_pen, band_prefix, donor_pen};

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
    #[inline]
    fn enter(&self, ix: &GenomeIndex, p: u32, minus: bool) -> u32 {
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
    fn base(&self, read: &[u8], i: usize) -> u8 {
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

#[derive(Default)]
pub(crate) struct Bufs {
    e1: Vec<u16>,
    e2: Vec<u16>,
    tmp: Vec<u16>,
    c: Vec<u32>,
    o: Vec<(u16, u16)>,
}

/// Solve `m` read bases from `rv` against the genome from `gv`.
///
/// `far` is the real genome position of the far fixed point (the read bases then end exactly
/// there; `gv` and `rv` must be forward); without it the read end is free and unexplained
/// bases are clipped. `windows` are candidate exons as (view position, length, read offset
/// that goes with the window's start), in travel order and not overlapping. `minus` is the transcript strand.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve(
    ix: &GenomeIndex,
    read: &[u8],
    rv: ReadView,
    m: usize,
    gv: View,
    far: Option<u32>,
    windows: &[(u32, u32, i64)],
    minus: bool,
    min_gap: u32,
    b: &mut Bufs,
) -> Solved {
    let band = 4 + m / 10;
    let w = 2 * band + 1;
    let limit = match far {
        Some(t1) => (t1 - gv.origin) as usize,
        None => gv.room(ix) as usize,
    };
    let n = (m + band).min(limit);
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
    // Level 0: the upstream exon alone. cur[i] = cost of read[..i] and of leaving at p[i].
    let mut cur = vec![BIG; m + 1];
    let mut curp = vec![0u32; m + 1];
    let mut best_end = (BIG, 0usize, 0usize, 0u32, None::<(u32, u32)>); // cost, level, i, end p, entry
    for i in 0..=m {
        let mut end = (INF as u32, 0u32);
        for k in 0..w {
            let j = (i + k) as isize - band as isize;
            if j < 0 || j as usize > n {
                continue;
            }
            let e = b.e1[i * w + k] as u32;
            if e >= INF as u32 {
                continue;
            }
            let c = e + gv.leave(ix, j as u32, minus);
            if c < cur[i] {
                cur[i] = c;
                curp[i] = j as u32;
            }
            if e < end.0 || (e == end.0 && (j as usize).abs_diff(i) < (end.1 as usize).abs_diff(i)) {
                end = (e, j as u32);
            }
        }
        if far.is_none() && end.0 < INF as u32 {
            let c = end.0 + CLIP * (m - i) as u32;
            if c < best_end.0 {
                best_end = (c, 0, i, end.1, None);
            }
        }
    }
    // Candidate exons, one level each. A level keeps, per read position, the cheapest way to
    // have left the genome there, through this window or not.
    let mut lev_p: Vec<Vec<u32>> = vec![curp.clone()];
    let mut via: Vec<Vec<Option<(u32, u32)>>> = Vec::new();
    for &(vs, vl, i0) in windows {
        let vl = vl as usize;
        // Only cells near the diagonal the window's seed hits lie on: window position j goes
        // with read position i0 + j, give or take `wb`.
        let wb = 16 + vl as i64 / 16;
        let ww = (2 * wb + 1) as usize;
        let i_lo = (i0 - wb).clamp(0, m as i64) as usize;
        let i_hi = (i0 + vl as i64 + wb).clamp(0, m as i64) as usize;
        let rows = i_hi - i_lo + 1;
        b.c.clear();
        b.c.resize(rows * ww, BIG);
        b.o.clear();
        b.o.resize(rows * ww, (0, 0));
        let enter: Vec<u32> = (0..=vl)
            .map(|j| EXTRA_INTRON + gv.enter(ix, vs + j as u32, minus))
            .collect();
        let leave: Vec<u32> = (0..=vl).map(|j| gv.leave(ix, vs + j as u32, minus)).collect();
        let gb: Vec<u8> = (0..vl).map(|j| gv.base(ix, vs + j as u32)).collect();
        let mut next = cur.clone();
        let mut nextp = curp.clone();
        let mut v: Vec<Option<(u32, u32)>> = vec![None; m + 1];
        // Affine gaps, as outside the windows: `up` holds, per k, the best cost of a cell that
        // ends in read bases left out (and where that path entered the window).
        let mut up: Vec<(u32, (u16, u16))> = vec![(BIG, (0, 0)); ww];
        let open_del = (GAP_OPEN + GAP_EXT) as u32;
        let open_ins = (GAP_OPEN + INS_EXT) as u32;
        for i in i_lo..=i_hi {
            let ri = if i > 0 { rv.base(read, i - 1) } else { 0 };
            let row = (i - i_lo) * ww;
            let mut left = (BIG, (0u16, 0u16)); // ends in genome bases skipped
            for k in 0..ww {
                let j = i as i64 - i0 + k as i64 - wb;
                if j < 0 || j > vl as i64 {
                    left = (BIG, (0, 0));
                    up[k] = (BIG, (0, 0));
                    continue;
                }
                let j = j as usize;
                let at = row + k;
                // enter here, if the genome was left far enough before this point
                let (mut c, mut o) = if cur[i] < BIG && vs + j as u32 >= curp[i] + min_gap {
                    (cur[i] + enter[j], (i as u16, j as u16))
                } else {
                    (BIG, (0, 0))
                };
                // (i - 1, j - 1) is the same k one row up; (i - 1, j) is k + 1 there
                let f = if i > i_lo && k + 1 < ww {
                    let (ext, opn) = (up[k + 1].0 + INS_EXT as u32, b.c[at - ww + 1] + open_ins);
                    if ext <= opn {
                        (ext, up[k + 1].1)
                    } else {
                        (opn, b.o[at - ww + 1])
                    }
                } else {
                    (BIG, (0, 0))
                };
                left = if k > 0 && j > 0 {
                    let (ext, opn) = (left.0 + GAP_EXT as u32, b.c[at - 1] + open_del);
                    if ext <= opn {
                        (ext, left.1)
                    } else {
                        (opn, b.o[at - 1])
                    }
                } else {
                    (BIG, (0, 0))
                };
                if i > i_lo && j > 0 {
                    let d = b.c[at - ww] + if ri != gb[j - 1] { MISMATCH as u32 } else { 0 };
                    if d < c {
                        c = d;
                        o = b.o[at - ww];
                    }
                }
                if f.0 < c {
                    c = f.0;
                    o = f.1;
                }
                if left.0 < c {
                    c = left.0;
                    o = left.1;
                }
                b.c[at] = c.min(BIG);
                b.o[at] = o;
                up[k] = (f.0.min(BIG), f.1);
                if c >= BIG || j as u16 <= o.1 {
                    continue; // nothing of the window used yet
                }
                let entry = Some((o.0 as u32, vs + o.1 as u32));
                if c + leave[j] < next[i] {
                    next[i] = c + leave[j];
                    nextp[i] = vs + j as u32;
                    v[i] = entry;
                }
                if far.is_none() {
                    let e = c + CLIP * (m - i) as u32;
                    if e < best_end.0 {
                        best_end = (e, via.len() + 1, i, vs + j as u32, entry);
                    }
                }
            }
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
        let mut best = (BIG, 0usize, 0u32);
        for (x, &c0) in cur.iter().enumerate() {
            if c0 >= BIG {
                continue;
            }
            let rest = m - x;
            for k in 0..w {
                let j = (rest + k) as isize - band as isize;
                if j < 0 || j as usize > n {
                    continue;
                }
                let e = b.e2[rest * w + k] as u32;
                if e >= INF as u32 {
                    continue;
                }
                let c = c0 + e + fv.leave(ix, j as u32, minus);
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
