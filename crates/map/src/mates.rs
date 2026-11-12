//! Antisense mates: per transcript, the transcripts that hold stretches of its reverse
//! complement, kept as a few bytes per stretch. A mapper reads them from the index
//! ([`build_into`], `Mates::load`) or makes a transcript's on first use.

use std::sync::{Arc, OnceLock};

use index::Index;
use seeds::BaseSeeder;

use crate::tables::{Arr, put, varint};

/// A stretch of a transcript that occurs reverse-complemented in another one, its antisense
/// mate: the base-level seeds of the transcript's reverse complement found in the mate on one
/// diagonal.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct MateSeg {
    /// The mate.
    pub t: u32,
    /// The stretch in the transcript (forward coordinates).
    pub x: (u32, u32),
    /// A seed's first base in the transcript + its first base in the mate + its length: the
    /// same for every seed of one antisense match.
    pub c: u32,
    /// Bases of the stretch's seeds.
    pub w: u32,
}

/// Seeds with more occurrences than this name no antisense mates: they are repeats, which
/// every transcript shares with hundreds of unrelated ones.
pub const MATE_MAX_OCC: usize = 1024;

/// Transcript occurrences (transcript, first base) of the seed at index nodes `nodes`, until
/// `f` returns true; whether it did.
fn node_occurrences(idx: &Index, nodes: std::ops::Range<u32>, mut f: impl FnMut(u32, u32) -> bool) -> bool {
    let flags = idx.node_flags();
    for node in nodes {
        if flags[node as usize] & index::FLAG_MASKED != 0 {
            continue;
        }
        let off = idx.node_offset()[node as usize];
        let (ts, ps) = idx.occurrences(idx.node_unitig()[node as usize]);
        if ts.iter().zip(ps).any(|(&t, &p)| f(t, p + off)) {
            return true;
        }
    }
    false
}

/// The [`MateSeg`]s of transcript `x`, by mate: of the mates sharing at least two seeds with
/// it. Seeds of more than `h_max` index nodes are skipped, as in anchoring.
pub(crate) fn segs(idx: &Index, seeder: &BaseSeeder, h_max: usize, x: u32) -> Vec<MateSeg> {
    let tlen = idx.tx_plen()[x as usize];
    let mut codes = idx.tx_codes(x, 0, tlen);
    codes.reverse();
    codes.iter_mut().for_each(|c| *c = 3 - *c);
    let mut seeds = Vec::new();
    seeder.seed_codes_into(&codes, &mut seeds, &mut Vec::new());
    // Per occurrence in another transcript: (mate, c, first base here, length).
    let mut occ: Vec<(u32, u32, u32, u32)> = Vec::new();
    for s in &seeds {
        let r = idx.lookup(s.label);
        if r.is_empty() || r.len() > h_max {
            continue;
        }
        let (from, at, len) = (occ.len(), tlen - s.end, s.end - s.start);
        let over = node_occurrences(idx, r, |t, p| {
            if t != x {
                occ.push((t, at + p + len, at, len));
            }
            occ.len() - from > MATE_MAX_OCC
        });
        if over {
            occ.truncate(from);
        }
    }
    occ.sort_unstable();
    let mut out: Vec<MateSeg> = Vec::new();
    let mut i = 0;
    while i < occ.len() {
        let t = occ[i].0;
        let (from, mut n) = (out.len(), 0);
        while i < occ.len() && occ[i].0 == t {
            let (_, c, at, len) = occ[i];
            match out[from..].last_mut() {
                Some(m) if m.c == c => {
                    m.x.1 = m.x.1.max(at + len);
                    m.w += len;
                }
                _ => out.push(MateSeg {
                    t,
                    x: (at, at + len),
                    c,
                    w: len,
                }),
            }
            n += 1;
            i += 1;
        }
        if n < 2 {
            out.truncate(from);
        }
    }
    out
}

#[inline]
fn zig(x: u32) -> u32 {
    (x << 1) ^ ((x as i32) >> 31) as u32
}

#[inline]
fn zag(x: u32) -> u32 {
    (x >> 1) ^ (x & 1).wrapping_neg()
}

/// `segs` (grouped by mate, mates ascending) as bytes ([`Coded`] reads them back; nothing for
/// no segments). Isoforms of one antisense gene share most stretches, so the distinct
/// stretches `(x, w)` come first, once, and each mate then lists stretch numbers and `c`; a
/// mate whose list is the previous mate's with every `c` shifted alike holds the shift only.
/// Varints throughout: the number of stretches; per stretch its start (less the previous
/// one's), length and `w`; per mate `t` (less the previous mate's), its number of segments
/// (0: the previous list again), its first `c` (zigzag, less the previous mate's) and, unless
/// repeated, per segment the stretch number (less the previous one's) and `c` (zigzag, less
/// the previous one's). Segments of a mate are in order of stretch number, then `c`.
pub(crate) fn encode(segs: &[MateSeg]) -> Vec<u8> {
    let mut out = Vec::new();
    if segs.is_empty() {
        return out;
    }
    let mut st: Vec<(u32, u32, u32)> = segs.iter().map(|m| (m.x.0, m.x.1, m.w)).collect();
    st.sort_unstable();
    st.dedup();
    put(&mut out, st.len() as u32);
    let mut x0 = 0;
    for s in &st {
        put(&mut out, s.0 - x0);
        put(&mut out, s.1 - s.0);
        put(&mut out, s.2);
        x0 = s.0;
    }
    let (mut i, mut pt, mut pc) = (0, 0u32, 0u32);
    // Per segment of the previous mate and of this one: (stretch number, c less the first).
    let (mut last, mut cur) = (Vec::<(u32, u32)>::new(), Vec::<(u32, u32)>::new());
    while i < segs.len() {
        let t = segs[i].t;
        cur.clear();
        while i < segs.len() && segs[i].t == t {
            let m = &segs[i];
            let id = st.binary_search(&(m.x.0, m.x.1, m.w)).unwrap_or(0);
            cur.push((id as u32, m.c));
            i += 1;
        }
        cur.sort_unstable();
        let c0 = cur[0].1;
        cur.iter_mut().for_each(|e| e.1 = e.1.wrapping_sub(c0));
        put(&mut out, t - pt);
        put(&mut out, if cur == last { 0 } else { cur.len() as u32 });
        put(&mut out, zig(c0.wrapping_sub(pc)));
        (pt, pc) = (t, c0);
        if cur != last {
            let (mut id, mut c) = (0u32, 0u32);
            for e in &cur {
                put(&mut out, e.0 - id);
                put(&mut out, zig(e.1.wrapping_sub(c)));
                (id, c) = *e;
            }
            std::mem::swap(&mut last, &mut cur);
        }
    }
    out
}

/// The mates of one transcript as [`encode`] wrote them.
#[derive(Clone, Copy)]
pub(crate) struct Coded<'a>(pub &'a [u8]);

impl Coded<'_> {
    pub(crate) fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// The distinct stretches as (start, end, `w`), into `out`; where the mates start.
    pub(crate) fn stretches(&self, out: &mut Vec<(u32, u32, u32)>) -> usize {
        out.clear();
        let b = self.0;
        let (n, mut p) = varint(b, 0);
        let mut x0 = 0u32;
        for _ in 0..n {
            if p > b.len() {
                break;
            }
            let (d, q) = varint(b, p);
            let (l, q) = varint(b, q);
            let (w, q) = varint(b, q);
            x0 = x0.wrapping_add(d);
            out.push((x0, x0.wrapping_add(l), w));
            p = q;
        }
        p
    }

    /// Each segment as (mate, stretch number, `c`), from `p` ([`Self::stretches`]): mates
    /// ascending, a mate's segments together.
    #[inline]
    pub(crate) fn for_each(&self, mut p: usize, mut f: impl FnMut(u32, usize, u32)) {
        let b = self.0;
        let (mut t, mut c0, mut list, mut len) = (0u32, 0u32, 0usize, 0u32);
        while p < b.len() {
            let (dt, q) = varint(b, p);
            let (n, q) = varint(b, q);
            let (dc, q) = varint(b, q);
            t = t.wrapping_add(dt);
            c0 = c0.wrapping_add(zag(dc));
            if n != 0 {
                (list, len) = (q, n);
            }
            let (mut r, mut id, mut c) = (list, 0u32, c0);
            for _ in 0..len {
                if r >= b.len() {
                    break;
                }
                let (di, s) = varint(b, r);
                let (dc, s) = varint(b, s);
                id = id.wrapping_add(di);
                c = c.wrapping_add(zag(dc));
                f(t, id as usize, c);
                r = s;
            }
            p = if n != 0 { r } else { q };
        }
    }

    /// The segments, by mate (the inverse of [`encode`] up to the order within a mate).
    #[cfg(test)]
    fn decode(&self) -> Vec<MateSeg> {
        let mut st = Vec::new();
        let p = self.stretches(&mut st);
        let mut out = Vec::new();
        self.for_each(p, |t, id, c| {
            let (a, b, w) = st[id];
            out.push(MateSeg { t, x: (a, b), c, w });
        });
        out
    }
}

/// Every transcript's mates: those kept in the index, or each made on first use.
pub(crate) enum Mates {
    Kept {
        /// Where each transcript's bytes start in `bytes`, and where the last one's end.
        off: Arr<u64>,
        bytes: Arr<u8>,
        _idx: Arc<Index>,
    },
    Lazy(Vec<OnceLock<Box<[u8]>>>),
}

/// The layout of the mates in an index: a reader takes only mates of its own.
const FORMAT: &str = "1";

/// The names of the mates' sections in an index.
mod section {
    pub const META: &str = "mates_meta";
    pub const OFF: &str = "mates_off";
    pub const BYTES: &str = "mates_bytes";
}

/// Every transcript's mates, coded and laid end to end, and where each one's start (with the
/// end of the last), on `nt` threads.
fn build(idx: &Index, seeder: &BaseSeeder, h_max: usize, nt: usize) -> (Vec<u64>, Vec<u8>) {
    let n = idx.n_transcripts();
    let nt = nt.clamp(1, n.max(1));
    // Transcripts are dealt in blocks, since neighbours (isoforms) cost alike.
    const BLOCK: usize = 64;
    let next = std::sync::atomic::AtomicUsize::new(0);
    let mut parts: Vec<(usize, Vec<u32>, Vec<u8>)> = std::thread::scope(|sc| {
        let workers: Vec<_> = (0..nt)
            .map(|_| {
                sc.spawn(|| {
                    let mut mine = Vec::new();
                    loop {
                        let from = next.fetch_add(BLOCK, std::sync::atomic::Ordering::Relaxed);
                        if from >= n {
                            break;
                        }
                        let (mut lens, mut bytes) = (Vec::new(), Vec::new());
                        for x in from..(from + BLOCK).min(n) {
                            let b = encode(&segs(idx, seeder, h_max, x as u32));
                            lens.push(b.len() as u32);
                            bytes.extend_from_slice(&b);
                        }
                        mine.push((from, lens, bytes));
                    }
                    mine
                })
            })
            .collect();
        workers
            .into_iter()
            .flat_map(|w| w.join().expect("a mates worker panicked"))
            .collect()
    });
    parts.sort_unstable_by_key(|p| p.0);
    let mut off = Vec::with_capacity(n + 1);
    let mut bytes = Vec::with_capacity(parts.iter().map(|p| p.2.len()).sum());
    for (_, lens, b) in &parts {
        let mut o = bytes.len() as u64;
        for &l in lens {
            off.push(o);
            o += u64::from(l);
        }
        bytes.extend_from_slice(b);
    }
    off.push(bytes.len() as u64);
    (off, bytes)
}

/// Make every transcript's mates of `idx` (`data` laid out), for anchoring with `h_max`, and
/// keep them in `data` as sections of the index: a mapper anchoring so then reads them in
/// place (`Mates::load`) instead of making each on first use.
pub fn build_into(data: &mut index::IndexData, h_max: usize, nt: usize) -> anyhow::Result<()> {
    let (off, bytes) = {
        let idx = Index::from_data(data)?;
        let seeder = seeder_of(&idx)?;
        build(&idx, &seeder, h_max, nt)
    };
    let meta = format!("format={FORMAT}\nh_max={h_max}\nmax_occ={MATE_MAX_OCC}\n");
    data.add_section(section::META, meta.as_bytes());
    data.add_section(section::OFF, &off);
    data.add_section(section::BYTES, &bytes);
    Ok(())
}

/// The seeder an index was built with, for sequence already parsed.
fn seeder_of(idx: &Index) -> anyhow::Result<BaseSeeder> {
    use anyhow::Context;
    let spec = idx.meta("spec").context("index lacks its seed spec")?;
    let hpc = idx.meta("hpc").context("index lacks its HPC flag")? == "true";
    Ok(BaseSeeder::new(spec, hpc)?)
}

impl Mates {
    /// The mates kept in `idx` ([`build_into`]), read in place, if they were made for
    /// anchoring with `h_max` and their sizes agree; otherwise none yet, each transcript's
    /// to be made on first use. Kept values are not checked (see [`crate::tables::load`]).
    pub(crate) fn load(idx: &Arc<Index>, h_max: usize) -> Self {
        Self::kept(idx, h_max)
            .unwrap_or_else(|| Mates::Lazy((0..idx.n_transcripts()).map(|_| Default::default()).collect()))
    }

    fn kept(idx: &Arc<Index>, h_max: usize) -> Option<Self> {
        let meta = std::str::from_utf8(idx.section::<u8>(section::META)?).ok()?;
        let field = |k: &str| {
            meta.lines()
                .filter_map(|l| l.split_once('='))
                .find(|x| x.0 == k)
                .map(|x| x.1)
        };
        if field("format")? != FORMAT
            || field("h_max")?.parse() != Ok(h_max)
            || field("max_occ")?.parse() != Ok(MATE_MAX_OCC)
        {
            return None;
        }
        let off: &[u64] = idx.section(section::OFF)?;
        let bytes: &[u8] = idx.section(section::BYTES)?;
        let whole = off.len() == idx.n_transcripts() + 1
            && off[0] == 0
            && off.windows(2).all(|w| w[0] <= w[1])
            && off.last().is_some_and(|&x| x as usize == bytes.len());
        whole.then(|| Mates::Kept {
            off: Arr::of(off),
            bytes: Arr::of(bytes),
            _idx: idx.clone(),
        })
    }

    /// The mates of transcript `x`.
    pub(crate) fn of(&self, idx: &Index, seeder: &BaseSeeder, h_max: usize, x: u32) -> Coded<'_> {
        match self {
            Mates::Kept { off, bytes, .. } => {
                Coded(&bytes[off[x as usize] as usize..off[x as usize + 1] as usize])
            }
            Mates::Lazy(v) => {
                Coded(v[x as usize].get_or_init(|| encode(&segs(idx, seeder, h_max, x)).into_boxed_slice()))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seg(t: u32, x: (u32, u32), c: u32, w: u32) -> MateSeg {
        MateSeg { t, x, c, w }
    }

    /// Coded segments read back as they were, a mate's in order of stretch then `c`.
    fn round_trip(segs: &[MateSeg]) {
        let bytes = encode(segs);
        let mut want = segs.to_vec();
        want.sort_unstable_by_key(|m| (m.t, m.x, m.w, m.c));
        assert_eq!(Coded(&bytes).decode(), want);
    }

    #[test]
    fn no_segments_take_no_bytes() {
        assert!(encode(&[]).is_empty());
        assert!(Coded(&[]).decode().is_empty());
    }

    #[test]
    fn segments_round_trip() {
        round_trip(&[seg(7, (10, 40), 500, 25)]);
        // Isoforms sharing stretches at shifted, equal and smaller `c`; a stretch twice in
        // one mate; a mate with a list of its own between two that repeat.
        round_trip(&[
            seg(3, (0, 30), 1000, 30),
            seg(3, (50, 90), 1200, 31),
            seg(3, (50, 90), 1300, 31),
            seg(4, (0, 30), 1010, 30),
            seg(4, (50, 90), 1210, 31),
            seg(4, (50, 90), 1310, 31),
            seg(5, (0, 30), 700, 30),
            seg(5, (50, 90), 900, 31),
            seg(5, (50, 90), 1000, 31),
            seg(9, (50, 90), 100, 31),
            seg(9, (0, 30), 400, 30),
            seg(200_000, (0, 30), 5, 30),
            seg(200_000, (50, 90), 205, 31),
            seg(200_000, (50, 90), 305, 31),
            seg(200_001, (70_000, 70_300), u32::MAX - 3, 290),
        ]);
    }

    #[test]
    fn a_repeated_list_is_not_written_again() {
        let one = [seg(3, (0, 30), 1000, 30), seg(3, (50, 90), 1200, 31)];
        let mut two = one.to_vec();
        two.extend(one.iter().map(|m| MateSeg {
            t: 4,
            c: m.c + 17,
            ..*m
        }));
        // Mate, a zero and the shift.
        assert_eq!(encode(&two).len(), encode(&one).len() + 3);
    }

    #[test]
    fn cut_bytes_end_the_reading() {
        let bytes = encode(&[
            seg(3, (0, 30), 1000, 30),
            seg(3, (50, 90), 1200, 31),
            seg(4, (0, 30), 1010, 30),
            seg(4, (50, 90), 1210, 31),
        ]);
        for n in 0..bytes.len() {
            let c = Coded(&bytes[..n]);
            let mut st = Vec::new();
            let p = c.stretches(&mut st);
            c.for_each(p, |_, _, _| {});
        }
    }
}
