//! Multi-transcript chaining with lazy color refinement.
//!
//! One DP over a candidate's anchors computes the best chain score of every transcript at once.
//! A state is (anchor, D, f): D is the set of transcript placements of the anchor for which a
//! chain ending there has score f. Extending a state to a later anchor splits D by the path
//! length between the two unitigs in each transcript (the only quantity the transition cost
//! depends on), so transcripts share states until their isoform paths diverge. Each placement
//! is then kept only in its best-scoring state (which merges identical sets and drops dominated
//! states), and an optional beam keeps the B best. With no beam the scores equal [`OracleChainer`](crate::OracleChainer)'s
//! bit for bit: the window, validity rule and transition arithmetic are the same.
//!
//! Sets are bitsets over the candidate's transcripts, one per occurrence rank: rank r holds the
//! transcripts whose r-th occurrence of the unitig is in the set (repeats are rare, so almost
//! every set has one rank).

use rustc_hash::FxHashMap;

use index::Index;

use crate::anchor::Candidate;
use crate::chain::{ChainParams, Chainer, TxHit, hit_from_path};

pub struct RefineChainer {
    pub p: ChainParams,
    /// Keep at most this many states per anchor (None: unbounded, exact).
    pub beam: Option<usize>,
}

const NONE: u32 = u32::MAX;

/// One distinct unitig of the candidate: occurrence positions per (rank, local transcript)
/// and the full set of its placements.
struct Unitig {
    ranks: usize,
    /// `pos[pos_at + r * ntx + j]`: position of the r-th occurrence in local transcript j, or
    /// NONE.
    pos_at: u32,
    /// Offset of the all-placements set (ranks × nw words) in `words`.
    all: u32,
}

/// Transcripts whose unitig path from `u` (rank r) to `u′` (rank r′) has length `dp`.
struct Part {
    r: u32,
    r2: u32,
    bits: u32,
}

struct Group {
    dp: i64,
    parts: std::ops::Range<u32>,
}

#[derive(Clone, Copy)]
struct State {
    anchor: u32,
    f: f32,
    /// Offset of the set (ranks(u) × nw words) in `state_words`.
    set: u32,
    back: u32,
    /// Path length from the predecessor's unitig, for traceback.
    dp: i64,
    /// Gap discrepancy of the incoming transition (0 for a chain start).
    d: f32,
}

/// A transition candidate before its set is computed: predecessor state, path-length group,
/// discrepancy and score.
#[derive(Clone, Copy)]
struct Cand {
    s: u32,
    g: u32,
    d: f32,
    f: f32,
}

#[derive(Default)]
struct Scratch {
    local: FxHashMap<u32, u32>,
    tx: Vec<u32>,
    uid: FxHashMap<u32, u32>,
    units: Vec<Unitig>,
    pos: Vec<u32>,
    anchor_unit: Vec<u32>,
    pairs: FxHashMap<(u32, u32), std::ops::Range<u32>>,
    items: Vec<(i64, u32, u32, u32)>,
    groups: Vec<Group>,
    parts: Vec<Part>,
    /// Bitset words: unitig sets and pair parts.
    words: Vec<u64>,
    states: Vec<State>,
    state_words: Vec<u64>,
    /// States of each anchor: `first[i]..first[i + 1]`.
    first: Vec<u32>,
    cands: Vec<Cand>,
    keys: Vec<u64>,
    covered: Vec<u64>,
    set: Vec<u64>,
}

/// Whether to collect state statistics (`REFINE_STATS` set in the environment).
fn stats_on() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("REFINE_STATS").is_some())
}

thread_local! {
    /// Per-thread state statistics when `REFINE_STATS` is set: (anchors, candidates, states,
    /// max states at an anchor, max transcripts in a candidate).
    pub static STATS: std::cell::RefCell<(usize, usize, usize, usize, usize)> = const { std::cell::RefCell::new((0, 0, 0, 0, 0)) };
    static SC: std::cell::RefCell<Scratch> = std::cell::RefCell::new(Scratch::default());
}

/// Sort key: decreasing f, then increasing candidate index.
fn key(f: f32, i: usize) -> u64 {
    let b = f.to_bits();
    let ord = if b & 0x8000_0000 != 0 { !b } else { b | 0x8000_0000 };
    ((!ord as u64) << 32) | i as u64
}

impl Chainer for RefineChainer {
    fn chain(&self, idx: &Index, cand: &Candidate) -> Vec<TxHit> {
        SC.with(|sc| self.chain_with(idx, cand, &mut sc.borrow_mut()))
    }
}

impl RefineChainer {
    fn chain_with(&self, idx: &Index, cand: &Candidate, sc: &mut Scratch) -> Vec<TxHit> {
        let p = &self.p;
        let unitig = idx.node_unitig();
        let offset = idx.node_offset();
        let anchors = &cand.anchors;
        let n = anchors.len();

        // Local transcript ids, in first-seen order.
        sc.local.clear();
        sc.tx.clear();
        for a in anchors {
            for &t in idx.occurrences(unitig[a.node as usize]).0 {
                let next = sc.tx.len() as u32;
                if *sc.local.entry(t).or_insert(next) == next {
                    sc.tx.push(t);
                }
            }
        }
        let ntx = sc.tx.len();
        let nw = ntx.div_ceil(64).max(1);

        // Distinct unitigs: positions per (rank, transcript) and their full placement sets.
        sc.uid.clear();
        sc.units.clear();
        sc.pos.clear();
        sc.words.clear();
        sc.anchor_unit.clear();
        for a in anchors {
            let u = unitig[a.node as usize];
            let next = sc.units.len() as u32;
            let id = *sc.uid.entry(u).or_insert(next);
            if id == next {
                let (ts, ps) = idx.occurrences(u);
                let mut ranks = 0;
                let mut k = 0;
                while k < ts.len() {
                    let e = k + ts[k..].iter().take_while(|&&t| t == ts[k]).count();
                    ranks = ranks.max(e - k);
                    k = e;
                }
                let pos_at = sc.pos.len();
                sc.pos.resize(pos_at + ranks * ntx, NONE);
                let all = sc.words.len() as u32;
                sc.words.resize(sc.words.len() + ranks * nw, 0);
                let mut k = 0;
                while k < ts.len() {
                    let j = sc.local[&ts[k]] as usize;
                    let mut r = 0;
                    while k < ts.len() && sc.local[&ts[k]] as usize == j {
                        sc.pos[pos_at + r * ntx + j] = ps[k];
                        sc.words[all as usize + r * nw + j / 64] |= 1 << (j % 64);
                        r += 1;
                        k += 1;
                    }
                }
                sc.units.push(Unitig {
                    ranks,
                    pos_at: pos_at as u32,
                    all,
                });
            }
            sc.anchor_unit.push(id);
        }
        sc.pairs.clear();
        sc.groups.clear();
        sc.parts.clear();

        sc.states.clear();
        sc.state_words.clear();
        sc.first.clear();
        sc.first.push(0);
        let beam = self.beam.unwrap_or(usize::MAX).max(1);
        for k in 0..n {
            let y = &anchors[k];
            let uk = sc.anchor_unit[k] as usize;
            let sw = sc.units[uk].ranks * nw;
            let (yqs, yqe) = (y.qs as i32, y.qe as i32);
            let y_off = offset[y.node as usize] as i64;
            let qlen = (yqe - yqs) as f32;
            // Score every transition first; a candidate can win placements only if it beats
            // the fresh start (f = w), which holds every placement.
            sc.cands.clear();
            let mut memo: Option<(usize, std::ops::Range<usize>)> = None;
            let lo = k.saturating_sub(p.h);
            for (i, x) in anchors.iter().enumerate().take(k).skip(lo) {
                let (xqs, xqe) = (x.qs as i32, x.qe as i32);
                let dq = yqs - xqe;
                if yqs <= xqs || dq < -(p.o_max as i32) || dq > p.g_max as i32 {
                    continue;
                }
                let (s0, s1) = (sc.first[i] as usize, sc.first[i + 1] as usize);
                if s0 == s1 {
                    continue;
                }
                let ui = sc.anchor_unit[i] as usize;
                let groups = match &memo {
                    Some((u, r)) if *u == ui => r.clone(),
                    _ => {
                        let r = self.pair(sc, ui, uk, ntx, nw);
                        memo = Some((ui, r.clone()));
                        r
                    }
                };
                let x_off = offset[x.node as usize] as i64;
                let x_len = (x.qe - x.qs) as i64;
                let overlap = (xqe - yqs).max(0).min(yqe - yqs) as f32;
                let gain = y.weight * (1.0 - overlap / qlen);
                for g in groups {
                    // ts' - ts = dp + off' - off; dt = ts' - (ts + len).
                    let dts = sc.groups[g].dp + y_off - x_off;
                    let dt = dts - x_len;
                    if dts <= 0 || dt < -p.o_max || dt > p.g_max {
                        continue;
                    }
                    let dt = dt as i32;
                    let d = (dq as f32 - p.rho * dt as f32).abs();
                    let gamma = if p.beta == 0.0 {
                        p.alpha * d
                    } else {
                        p.alpha * d + p.beta * (1.0 + dq.max(dt).max(0) as f32).log2()
                    };
                    let delta = gain - gamma;
                    for s in s0..s1 {
                        let f = sc.states[s].f + delta;
                        if f > y.weight {
                            sc.cands.push(Cand {
                                s: s as u32,
                                g: g as u32,
                                d,
                                f,
                            });
                        }
                    }
                }
            }
            sc.keys.clear();
            sc.keys
                .extend(sc.cands.iter().enumerate().map(|(i, c)| key(c.f, i)));
            sc.keys.sort_unstable();
            // In decreasing f, each candidate takes the placements no better candidate has.
            // This merges identical sets and drops dominated states (D ⊆ D″, f ≤ f″) exactly,
            // since a placement's future depends only on its own coordinates: the states
            // partition the anchor's placements by optimal score.
            let full = sc.units[uk].all as usize;
            sc.covered.clear();
            sc.covered.resize(sw, 0);
            let mut made = 0;
            for ki in 0..sc.keys.len() {
                if made == beam {
                    break;
                }
                let c = sc.cands[(sc.keys[ki] & 0xffff_ffff) as usize];
                let st = sc.states[c.s as usize];
                let grp = &sc.groups[c.g as usize];
                sc.set.clear();
                sc.set.resize(sw, 0);
                for pi in grp.parts.clone() {
                    let pt = &sc.parts[pi as usize];
                    let src = &sc.state_words[st.set as usize + pt.r as usize * nw..][..nw];
                    let bits = &sc.words[pt.bits as usize..][..nw];
                    let dst = &mut sc.set[pt.r2 as usize * nw..][..nw];
                    for w in 0..nw {
                        dst[w] |= src[w] & bits[w];
                    }
                }
                let mut any = 0u64;
                for w in 0..sw {
                    let v = sc.set[w] & !sc.covered[w];
                    sc.set[w] = v;
                    sc.covered[w] |= v;
                    any |= v;
                }
                if any == 0 {
                    continue;
                }
                let set = sc.state_words.len() as u32;
                sc.state_words.extend_from_slice(&sc.set);
                sc.states.push(State {
                    anchor: k as u32,
                    f: c.f,
                    set,
                    back: c.s,
                    dp: grp.dp,
                    d: c.d,
                });
                made += 1;
                if (0..sw).all(|w| sc.covered[w] == sc.words[full + w]) {
                    break;
                }
            }
            // The fresh start keeps the rest.
            let set = sc.state_words.len() as u32;
            let mut any = 0u64;
            for w in 0..sw {
                let v = sc.words[full + w] & !sc.covered[w];
                sc.state_words.push(v);
                any |= v;
            }
            if any != 0 {
                sc.states.push(State {
                    anchor: k as u32,
                    f: y.weight,
                    set,
                    back: NONE,
                    dp: 0,
                    d: 0.0,
                });
            } else {
                sc.state_words.truncate(set as usize);
            }
            if stats_on() {
                let kept = sc.states.len() - sc.first[k] as usize;
                STATS.with(|x| {
                    let mut x = x.borrow_mut();
                    x.0 += 1;
                    x.1 += sc.cands.len();
                    x.2 += kept;
                    x.3 = x.3.max(kept);
                    x.4 = x.4.max(ntx);
                });
            }
            sc.first.push(sc.states.len() as u32);
        }
        self.extract(idx, cand, sc, ntx, nw)
    }

    /// Path-length groups between unitigs `ui` and `uk` (cached per candidate), sorted by dp.
    fn pair(&self, sc: &mut Scratch, ui: usize, uk: usize, ntx: usize, nw: usize) -> std::ops::Range<usize> {
        if let Some(r) = sc.pairs.get(&(ui as u32, uk as u32)) {
            return r.start as usize..r.end as usize;
        }
        // (dp, r, r2, transcript) for every pair of occurrences in the same transcript.
        let mut items = std::mem::take(&mut sc.items);
        items.clear();
        let (a, b) = (&sc.units[ui], &sc.units[uk]);
        let (pa_at, pb_at) = (a.pos_at as usize, b.pos_at as usize);
        for j in 0..ntx {
            for r in 0..a.ranks {
                let pa = sc.pos[pa_at + r * ntx + j];
                if pa == NONE {
                    continue;
                }
                for r2 in 0..b.ranks {
                    let pb = sc.pos[pb_at + r2 * ntx + j];
                    if pb != NONE {
                        items.push((pb as i64 - pa as i64, r as u32, r2 as u32, j as u32));
                    }
                }
            }
        }
        items.sort_unstable();
        let g0 = sc.groups.len();
        let mut k = 0;
        while k < items.len() {
            let dp = items[k].0;
            let p0 = sc.parts.len() as u32;
            while k < items.len() && items[k].0 == dp {
                let (r, r2) = (items[k].1, items[k].2);
                let bits = sc.words.len() as u32;
                sc.words.resize(sc.words.len() + nw, 0);
                while k < items.len() && items[k].0 == dp && items[k].1 == r && items[k].2 == r2 {
                    let j = items[k].3 as usize;
                    sc.words[bits as usize + j / 64] |= 1 << (j % 64);
                    k += 1;
                }
                sc.parts.push(Part { r, r2, bits });
            }
            sc.groups.push(Group {
                dp,
                parts: p0..sc.parts.len() as u32,
            });
        }
        sc.items = items;
        let r = g0 as u32..sc.groups.len() as u32;
        sc.pairs.insert((ui as u32, uk as u32), r.clone());
        r.start as usize..r.end as usize
    }

    /// Best state per transcript (earliest anchor and rank on ties, as the oracle), then the
    /// chain by following back pointers.
    fn extract(&self, idx: &Index, cand: &Candidate, sc: &Scratch, ntx: usize, nw: usize) -> Vec<TxHit> {
        let state_words = &sc.state_words;
        let offset = idx.node_offset();
        // Per transcript: (f, state, rank).
        let mut best: Vec<(f32, u32, u32)> = vec![(f32::NEG_INFINITY, NONE, 0); ntx];
        for (s, st) in sc.states.iter().enumerate() {
            let u = &sc.units[sc.anchor_unit[st.anchor as usize] as usize];
            for r in 0..u.ranks {
                let set = &state_words[st.set as usize + r * nw..][..nw];
                for (w, &word) in set.iter().enumerate() {
                    let mut bits = word;
                    while bits != 0 {
                        let j = w * 64 + bits.trailing_zeros() as usize;
                        bits &= bits - 1;
                        let b = &mut best[j];
                        // Earlier states come from earlier anchors; same anchor: lower rank.
                        if st.f > b.0
                            || (st.f == b.0
                                && sc.states[b.1 as usize].anchor == st.anchor
                                && (r as u32) < b.2)
                        {
                            *b = (st.f, s as u32, r as u32);
                        }
                    }
                }
            }
        }
        let mut hits = Vec::with_capacity(ntx);
        let mut path: Vec<(u32, u32, u32, u32)> = Vec::new();
        for (j, &(f, s0, r0)) in best.iter().enumerate() {
            if s0 == NONE {
                continue;
            }
            path.clear();
            let (mut s, mut r) = (s0 as usize, r0 as usize);
            let mut maxd = 0f32;
            loop {
                let st = &sc.states[s];
                let a = &cand.anchors[st.anchor as usize];
                let u = &sc.units[sc.anchor_unit[st.anchor as usize] as usize];
                let pos = sc.pos[u.pos_at as usize + r * ntx + j];
                path.push((a.qs, a.qe, pos + offset[a.node as usize], a.qe - a.qs));
                maxd = maxd.max(st.d);
                if st.back == NONE {
                    break;
                }
                // The predecessor's occurrence of j is `dp` before this one.
                let prev = &sc.states[st.back as usize];
                let pu = &sc.units[sc.anchor_unit[prev.anchor as usize] as usize];
                let want = pos as i64 - st.dp;
                r = (0..pu.ranks)
                    .find(|&q| {
                        sc.pos[pu.pos_at as usize + q * ntx + j] as i64 == want
                            && state_words[prev.set as usize + q * nw + j / 64] >> (j % 64) & 1 == 1
                    })
                    .expect("traceback: predecessor holds the transcript");
                s = st.back as usize;
            }
            path.reverse();
            hits.push(hit_from_path(sc.tx[j], f, maxd, path.iter().copied()));
        }
        hits
    }
}
