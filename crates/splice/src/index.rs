//! The genome index: packed sequence, level-3 seed table, positional level-2 table.
//!
//! All positions are raw (uncompressed) bases on one axis that concatenates the reference
//! sequences in file order, so a position fits 32 bits for genomes under 4.29 Gbp.
//!
//! - `seq`: 2 bits per base (A=0, C=1, G=2, T=3; other symbols are stored as a base that
//!   depends on the position, so that they match nothing, are listed in `n_runs`, and are
//!   never seeded).
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
use index::store::{FileWriter, SectionSink, Store, Writer};
use rayon::prelude::*;
use seeds::BaseSeeder;

/// Sequence is seeded in chunks of this many bases, so that long chromosomes parse in
/// parallel and in bounded memory. Cores that would straddle a chunk edge are lost.
const CHUNK: usize = 8 << 20;

/// Segments parsed at one time, at most. Parsing a segment of `CHUNK` bases holds about a third
/// of a gigabyte, and parsing is a small part of the build, so more threads than this buy
/// little time for much memory.
const PARSE_THREADS: usize = 8;

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
    n_runs: *const [u32],
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

/// The genome as it is read: packed sequence, sequence names and starts, runs of other
/// symbols. The sequence is never held unpacked; a FASTA file is packed block by block.
struct Packed {
    words: Vec<u64>,
    names: Vec<String>,
    /// Start of each sequence, plus the total length.
    offsets: Vec<u32>,
    /// `[start, end)` of each maximal run of symbols other than A, C, G, T, flattened.
    n_runs: Vec<u32>,
}

/// Bases packed in one step while reading (a multiple of 32).
const BLOCK: usize = 1 << 26;

impl Packed {
    /// Pack `block`, whose first base is at `at` (a multiple of 32), and note its runs of
    /// other symbols.
    fn push_block(&mut self, at: usize, block: &[u8]) {
        const PIECE: usize = 1 << 16;
        let parts: Vec<(Vec<u64>, Vec<u32>)> = block
            .par_chunks(PIECE)
            .enumerate()
            .map(|(k, piece)| {
                let lo = at + k * PIECE;
                let mut runs = Vec::new();
                let words = piece
                    .chunks(32)
                    .enumerate()
                    .map(|(w, bases)| {
                        let mut x = 0u64;
                        for (i, &b) in bases.iter().enumerate() {
                            let p = lo + w * 32 + i;
                            // A symbol other than A, C, G, T must not look like any base: it
                            // is stored as a base that depends on its position, so that a
                            // run of N matches nothing (a read's poly(A) tail would otherwise
                            // match an N spacer stored as A).
                            let c = match code(b) {
                                c @ 0..=3 => c as u64,
                                _ => {
                                    match runs.last_mut() {
                                        Some(end) if *end == p as u32 => *end += 1,
                                        _ => runs.extend([p as u32, p as u32 + 1]),
                                    }
                                    (p as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 62
                                }
                            };
                            x |= c << (2 * i);
                        }
                        x
                    })
                    .collect();
                (words, runs)
            })
            .collect();
        for (words, runs) in parts {
            self.words.extend(words);
            let mut runs = &runs[..];
            // A run that continues the last one of the piece before.
            if let (Some(end), Some(&start)) = (self.n_runs.last_mut(), runs.first())
                && *end == start
            {
                *end = runs[1];
                runs = &runs[2..];
            }
            self.n_runs.extend(runs);
        }
    }

    fn read(mut input: impl BufRead) -> Result<Packed> {
        let mut g = Packed {
            words: Vec::new(),
            names: Vec::new(),
            offsets: Vec::new(),
            n_runs: Vec::new(),
        };
        let mut seen = std::collections::HashSet::new();
        // Bases read and not yet packed, and the number packed.
        let (mut block, mut done) = (Vec::with_capacity(BLOCK + (1 << 16)), 0usize);
        let mut line = Vec::new();
        loop {
            line.clear();
            if input.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            let text = line.trim_ascii_end();
            if let Some(h) = text.strip_prefix(b">") {
                let h = std::str::from_utf8(h).context("a sequence name is not valid UTF-8")?;
                let name = h.split_whitespace().next().unwrap_or("").to_string();
                ensure!(!name.is_empty(), "sequence {} has no name", g.names.len() + 1);
                ensure!(seen.insert(name.clone()), "two sequences are named {name}");
                g.names.push(name);
                g.offsets.push((done + block.len()) as u32);
            } else {
                ensure!(
                    !g.names.is_empty() || text.is_empty(),
                    "not a FASTA file: sequence before the first `>` line"
                );
                block.extend_from_slice(text);
            }
            ensure!(
                done + block.len() < u32::MAX as usize,
                "reference longer than 4.29 Gbp is not supported yet"
            );
            if block.len() >= BLOCK {
                let full = block.len() & !31;
                g.push_block(done, &block[..full]);
                block.drain(..full);
                done += full;
            }
        }
        ensure!(!g.names.is_empty(), "no sequences in the FASTA input");
        g.push_block(done, &block);
        let total = done + block.len();
        g.offsets.push(total as u32);
        g.words.resize(total.div_ceil(32) + 1, 0);
        Ok(g)
    }

    /// Maximal ACGT stretches of each sequence, cut into chunks of at most `CHUNK` bases.
    fn segments(&self) -> Vec<(usize, usize)> {
        let mut v = Vec::new();
        let mut runs = self.n_runs.chunks(2).peekable();
        for w in self.offsets.windows(2) {
            let (mut i, end) = (w[0] as usize, w[1] as usize);
            while i < end {
                // The run of other symbols at or after `i`, if it starts in this sequence.
                while runs.peek().is_some_and(|r| r[1] as usize <= i) {
                    runs.next();
                }
                let (rs, re) = runs.peek().map_or((end, end), |r| (r[0] as usize, r[1] as usize));
                let stop = rs.clamp(i, end);
                let mut a = i;
                while a < stop {
                    let b = (a + CHUNK).min(stop);
                    v.push((a, b));
                    a = b;
                }
                i = if rs < end { re.min(end).max(stop) } else { end };
            }
        }
        v
    }

    /// Bases `[a, b)` as letters.
    fn bases(&self, a: usize, b: usize) -> Vec<u8> {
        (a..b)
            .map(|p| b"ACGT"[((self.words[p / 32] >> (2 * (p % 32))) & 3) as usize])
            .collect()
    }
}

/// A seed occurrence in a table: (label, raw start, raw length up to 255).
type Entry = (u64, u32, u8);

/// The cores of `spec` in bases `[a, b)` of the genome as (label, raw start, raw length), in
/// genome order.
fn cores(seeder: &BaseSeeder, g: &Packed, a: usize, b: usize, mut each: impl FnMut(u64, u32, u32)) {
    let bases = g.bases(a, b);
    let h = seeder.parsed(&bases);
    let parsed: &[u8] = h.as_ref().map_or(&bases, |h| &h.seq);
    for s in seeder.seed_parsed(parsed) {
        let (rs, re) = h.as_ref().map_or((s.start, s.end), |h| h.to_raw(s.start, s.end));
        each(s.label, a as u32 + rs, re - rs);
    }
}

/// Run `work` on every segment and hand the results to `sink` in genome order.
///
/// Segments are parsed a batch at a time into buffers that are reused, on at most
/// [`PARSE_THREADS`] threads of their own, so the memory held besides what `sink` keeps is that
/// of one batch and of those threads' parse buffers, which end with the threads.
fn each_segment<B: Default + Send>(
    segs: &[(usize, usize)],
    work: impl Fn(usize, usize, &mut B) + Sync,
    mut sink: impl FnMut(&B),
) {
    let threads = PARSE_THREADS.min(rayon::current_num_threads().max(1));
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("thread pool");
    let mut bufs: Vec<B> = Vec::new();
    for batch in segs.chunks(threads) {
        bufs.resize_with(bufs.len().max(batch.len()), B::default);
        pool.install(|| {
            bufs.par_iter_mut()
                .zip(batch)
                .for_each(|(buf, &(a, b))| work(a, b, buf))
        });
        bufs[..batch.len()].iter().for_each(&mut sink);
    }
}

/// What the level-2 pass makes of one segment.
#[derive(Default)]
struct Level2 {
    pos: Vec<u32>,
    /// Labels folded to 16 bits.
    folded: Vec<u16>,
    /// The labels, and the distinct ones in order if they are few, while ranks may be wanted.
    labels: Vec<u64>,
    distinct: Option<Vec<u64>>,
    /// Pairs of consecutive cores within the segment.
    pairs: Vec<Entry>,
    /// The first core (label, start, length) and the last (label, start).
    first: Option<(u64, u32, u32)>,
    last: Option<(u64, u32)>,
}

/// The fallback seed of two consecutive level-2 cores, if it is one that is kept.
#[inline]
fn pair_entry(l0: u64, p0: u32, l1: u64, p1: u32, len1: u32) -> Option<Entry> {
    let l = pair_label(l0, l1);
    (p1 - p0 < 64 && keep_pair(l)).then(|| (l, p0, (p1 + len1 - p0).min(255) as u8))
}

/// `x` and `y` (sorted, distinct) merged, unless that is more than `limit` values.
fn merge_up_to(x: &[u64], y: &[u64], limit: usize) -> Option<Vec<u64>> {
    let mut v = Vec::with_capacity(x.len() + y.len());
    let (mut i, mut j) = (0, 0);
    while i < x.len() || j < y.len() {
        let l = if j == y.len() || (i < x.len() && x[i] <= y[j]) {
            x[i]
        } else {
            y[j]
        };
        i += (i < x.len() && x[i] == l) as usize;
        j += (j < y.len() && y[j] == l) as usize;
        v.push(l);
    }
    (v.len() <= limit).then_some(v)
}

/// Write the label table of `entries` (sorted) as sections `<prefix>_labels`, `_dir`,
/// `_occ_off`, `_pos` and `_len`, leaving out labels with more than `max_occ` occurrences.
/// Returns the directory's bits and the number of occurrences left out.
fn put_table(
    sink: &mut impl SectionSink,
    prefix: &str,
    mut entries: Vec<Entry>,
    max_occ: u32,
) -> Result<(u32, u64)> {
    let mut labels: Vec<u64> = Vec::new();
    let mut occ_off: Vec<u32> = vec![0];
    // Kept occurrences are moved to the front, so the table needs no second copy of them.
    let (mut kept, mut i) = (0, 0);
    while i < entries.len() {
        let mut j = i + 1;
        while j < entries.len() && entries[j].0 == entries[i].0 {
            j += 1;
        }
        if j - i <= max_occ as usize {
            labels.push(entries[i].0);
            entries.copy_within(i..j, kept);
            kept += j - i;
            occ_off.push(kept as u32);
        }
        i = j;
    }
    let dropped = (entries.len() - kept) as u64;
    entries.truncate(kept);
    let bits = (labels.len().max(2).ilog2()).clamp(1, 28);
    let mut dir = vec![0u32; (1usize << bits) + 1];
    for &l in &labels {
        dir[(l >> (64 - bits)) as usize + 1] += 1;
    }
    for b in 1..dir.len() {
        dir[b] += dir[b - 1];
    }
    sink.put(&format!("{prefix}_labels"), labels)?;
    sink.put(&format!("{prefix}_dir"), dir)?;
    sink.put(&format!("{prefix}_occ_off"), occ_off)?;
    sink.put(
        &format!("{prefix}_pos"),
        entries.par_iter().map(|x| x.1).collect::<Vec<u32>>(),
    )?;
    sink.put(
        &format!("{prefix}_len"),
        entries.par_iter().map(|x| x.2).collect::<Vec<u8>>(),
    )?;
    Ok((bits, dropped))
}

/// Sections of a genome index, besides a manifest.
const SECTIONS: usize = 18;

/// Build the index of `fasta` into `sink`; returns the sequence names and their starts (with
/// the total length last).
///
/// The build is in stages, each of which hands its finished tables to `sink` before the next
/// begins: the packed genome is read, then parsed once for the level-3 table and once more for
/// the level-2 table and the table of pairs. With a sink that writes sections out, the most
/// held at one time is the packed genome, one stage's table as it is sorted, and the parse
/// buffers.
fn build_into(
    fasta: impl BufRead,
    opt: &BuildOptions,
    sink: &mut impl SectionSink,
) -> Result<(Vec<String>, Vec<u32>)> {
    let g = Packed::read(fasta)?;
    let segs = g.segments();
    let total = *g.offsets.last().unwrap() as usize;
    let l3_seeder = BaseSeeder::new(L3_SPEC, opt.hpc).expect("valid seed spec");
    let l2_seeder = BaseSeeder::new(L2_SPEC, opt.hpc).expect("valid seed spec");

    // Level 3: sorted by (label, position); frequent labels dropped. The space reserved here
    // and below is an estimate from the human genome (0.08 level-3 cores, 0.18 level-2 cores
    // and 0.09 pairs per base without compression); it is not touched until it is used.
    let mut l3: Vec<Entry> = Vec::with_capacity(total / 8);
    each_segment(
        &segs,
        |a, b, buf: &mut Vec<Entry>| {
            buf.clear();
            cores(&l3_seeder, &g, a, b, |l, p, n| buf.push((l, p, n.min(255) as u8)));
        },
        |buf| l3.extend_from_slice(buf),
    );
    l3.par_sort_unstable();
    let (dir_bits, dropped) = put_table(sink, "g3", l3, opt.max_occ)?;

    // Level 2 in genome order: positions and 16-bit ids. Ids are ranks when the labels fit
    // 16 bits (they do for small genomes); otherwise a 16-bit fold of the label, whose
    // collisions only add chance hits that the gap solver discards. The labels themselves are
    // kept only while ranks are still possible.
    //
    // The same pass makes the fallback seeds for reads that level 3 anchors weakly: pairs of
    // consecutive level-2 cores (every other pair, chosen by label), as a label table like
    // level 3's.
    let limit = u16::MAX as usize;
    let few = std::sync::atomic::AtomicBool::new(true);
    let (mut l2_pos, mut l2_id) = (
        Vec::<u32>::with_capacity(total / 4),
        Vec::<u16>::with_capacity(total / 4),
    );
    let mut pairs: Vec<Entry> = Vec::with_capacity(total / 8);
    let mut ranks: Option<(Vec<u64>, Vec<u64>)> = Some((Vec::new(), Vec::new()));
    let mut prev: Option<(u64, u32)> = None;
    each_segment(
        &segs,
        |a, b, s: &mut Level2| {
            let few = few.load(std::sync::atomic::Ordering::Relaxed);
            let (mut labels, mut lens) = (std::mem::take(&mut s.labels), Vec::new());
            labels.clear();
            s.pos.clear();
            cores(&l2_seeder, &g, a, b, |l, p, n| {
                labels.push(l);
                s.pos.push(p);
                lens.push(n);
            });
            assert!(s.pos.is_sorted(), "level-2 cores are not in genome order");
            s.folded.clear();
            s.folded.extend(labels.iter().map(|&l| fold16(l)));
            s.pairs.clear();
            s.pairs.extend(
                (1..labels.len())
                    .filter_map(|i| pair_entry(labels[i - 1], s.pos[i - 1], labels[i], s.pos[i], lens[i])),
            );
            s.first = labels.first().map(|&l| (l, s.pos[0], lens[0]));
            s.last = labels.last().map(|&l| (l, *s.pos.last().unwrap()));
            s.distinct = few
                .then(|| {
                    let mut d = labels.clone();
                    d.sort_unstable();
                    d.dedup();
                    (d.len() <= limit).then_some(d)
                })
                .flatten();
            s.labels = labels;
        },
        |s| {
            if let Some((l, p, n)) = s.first {
                if let Some((l0, p0)) = prev {
                    assert!(p0 <= p, "level-2 cores are not in genome order");
                    pairs.extend(pair_entry(l0, p0, l, p, n));
                }
                prev = s.last;
            }
            l2_pos.extend_from_slice(&s.pos);
            l2_id.extend_from_slice(&s.folded);
            pairs.extend_from_slice(&s.pairs);
            ranks = ranks.take().and_then(|(mut all, distinct)| {
                let distinct = merge_up_to(&distinct, s.distinct.as_ref()?, limit)?;
                all.extend_from_slice(&s.labels);
                Some((all, distinct))
            });
            if ranks.is_none() {
                few.store(false, std::sync::atomic::Ordering::Relaxed);
            }
        },
    );
    let ranked = ranks.is_some();
    let l2_labels = match ranks {
        Some((all, distinct)) => {
            l2_id = all
                .par_iter()
                .map(|l| distinct.binary_search(l).unwrap() as u16)
                .collect();
            distinct
        }
        None => Vec::new(),
    };
    let Packed {
        words,
        names,
        offsets,
        n_runs,
    } = g;
    sink.put("g_seq", words)?;
    sink.put("g2_labels", l2_labels)?;
    sink.put("g2_pos", l2_pos)?;
    sink.put("g2_id", l2_id)?;
    pairs.par_sort_unstable();
    // No label of a pair is left out for being frequent without a count: `dropped` is level 3's.
    let (p_bits, _) = put_table(sink, "gp", pairs, opt.max_occ)?;

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
    sink.put("g_meta", meta)?;
    sink.put("g_names", name_bytes)?;
    sink.put("g_offsets", offsets.clone())?;
    sink.put("g_n_runs", n_runs)?;
    Ok((names, offsets))
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
    /// Build in memory from a FASTA file (plain text; see [`GenomeIndex::build_from_reader`]
    /// for any other source). [`GenomeIndex::build_to_file`] needs much less memory.
    pub fn build(fasta: &Path, opt: &BuildOptions) -> Result<GenomeIndex> {
        let f = File::open(fasta).with_context(|| format!("opening {}", fasta.display()))?;
        Self::build_from_reader(BufReader::with_capacity(1 << 24, f), opt)
            .with_context(|| format!("reading {}", fasta.display()))
    }

    /// Build in memory from FASTA text. Sequence names are the header up to the first white
    /// space.
    pub fn build_from_reader(fasta: impl BufRead, opt: &BuildOptions) -> Result<GenomeIndex> {
        let mut w = Writer::default();
        build_into(fasta, opt, &mut w)?;
        Self::from_store(w.into_store()?)
    }

    /// Build from FASTA text straight into the file `path`, and open it. Each table is written
    /// as soon as it is made, so the build holds little more than the largest of them.
    ///
    /// `manifest` is given the sequence names and their starts (with the total length last)
    /// once the index is made, and returns the record to store in it, if any.
    pub fn build_to_file(
        fasta: impl BufRead,
        opt: &BuildOptions,
        path: &Path,
        manifest: impl FnOnce(&[String], &[u32]) -> Result<Option<Vec<u8>>>,
    ) -> Result<(GenomeIndex, u64)> {
        let mut w = FileWriter::create(path, SECTIONS + 1)?;
        let (names, offsets) = build_into(fasta, opt, &mut w)?;
        if let Some(m) = manifest(&names, &offsets)? {
            w.put("g_manifest", m)?;
        }
        let bytes = w.finish()?;
        Ok((Self::open(path)?, bytes))
    }

    pub fn open(path: &Path) -> Result<GenomeIndex> {
        Self::from_store(Store::open(path)?)
    }

    /// Whether the file at `path` is a genome index (and not, say, a transcriptome index,
    /// which shares the container format). An error if it is no spearfish index at all.
    pub fn is_genome_index(path: &Path) -> Result<bool> {
        Ok(Store::open(path)?.has("g_meta"))
    }

    /// Write the index to `path`; returns the bytes written.
    pub fn write(&self, path: &Path) -> Result<u64> {
        self.write_with_manifest(path, None)
    }

    /// Write the index with a record of how it was built (any text; the `spearfish` library
    /// stores JSON), replacing the one it holds, if any.
    pub fn write_with_manifest(&self, path: &Path, manifest: Option<&[u8]>) -> Result<u64> {
        let mut w = Writer::default();
        match manifest {
            Some(m) => w.add::<u8>("g_manifest", m),
            None if self.store.has("g_manifest") => w.add::<u8>("g_manifest", self.store.get("g_manifest")?),
            None => {}
        }
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
            n_runs: store.get("g_n_runs")?,
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

    /// The record of how the index was built, when it was written with one.
    pub fn manifest(&self) -> Option<&[u8]> {
        self.store.get::<u8>("g_manifest").ok()
    }

    /// Level-3 labels with more occurrences than this were left out of the index, and how
    /// many labels that was.
    pub fn max_occ(&self) -> (u32, u64) {
        let meta: &[u64] = self.store.get("g_meta").expect("checked when opened");
        (meta[1] as u32, meta[2])
    }

    /// Sequence names, in file order.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Length of sequence `i`.
    pub fn seq_len(&self, i: usize) -> u32 {
        self.offsets()[i + 1] - self.offsets()[i]
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

    /// The runs of symbols other than A, C, G, T, as sorted `start, end` pairs one after the
    /// other (the sequence holds some base there).
    pub fn n_runs(&self) -> &[u32] {
        // SAFETY: see `Raw`.
        unsafe { &*self.n.n_runs }
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

    /// The 32 bases from `p` on, packed two bits each with the first base lowest (bases past
    /// the end of the genome read as A).
    #[inline]
    pub fn word(&self, p: u32) -> u64 {
        // SAFETY: see `Raw`.
        let w = unsafe { &*self.n.seq };
        let (i, sh) = ((p >> 5) as usize, (p & 31) * 2);
        let lo = w.get(i).copied().unwrap_or(0) >> sh;
        if sh == 0 {
            lo
        } else {
            lo | (w.get(i + 1).copied().unwrap_or(0) << (64 - sh))
        }
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

    /// The tables of level-3 labels, or those of level-2 pairs.
    #[inline]
    #[allow(clippy::type_complexity)]
    fn tables(&self, pairs: bool) -> (&[u64], &[u32], &[u32], &[u32], &[u8], u32) {
        // SAFETY: see `Raw`.
        unsafe {
            if pairs {
                (
                    &*self.n.p_labels,
                    &*self.n.p_dir,
                    &*self.n.p_occ_off,
                    &*self.n.p_pos,
                    &*self.n.p_len,
                    self.n.p_bits,
                )
            } else {
                (
                    &*self.n.labels,
                    &*self.n.dir,
                    &*self.n.occ_off,
                    &*self.n.pos,
                    &*self.n.len,
                    self.n.dir_bits,
                )
            }
        }
    }

    /// Looking a label up touches four places in memory, each known only from the one before.
    /// A read's labels are therefore looked up together, a step at a time: this asks for a
    /// label's directory entry, [`Self::prefetch_bucket`] for its bucket of labels, and
    /// [`Self::occ_range`] finds its occurrence list and asks for that.
    #[inline]
    pub fn prefetch_dir(&self, label: u64, pairs: bool) {
        let (_, dir, _, _, _, bits) = self.tables(pairs);
        prefetch(dir.as_ptr().wrapping_add((label >> (64 - bits)) as usize));
    }

    /// See [`Self::prefetch_dir`].
    #[inline]
    pub fn prefetch_bucket(&self, label: u64, pairs: bool) {
        let (labels, dir, off, _, _, bits) = self.tables(pairs);
        let lo = dir[(label >> (64 - bits)) as usize] as usize;
        prefetch(labels.as_ptr().wrapping_add(lo));
        prefetch(off.as_ptr().wrapping_add(lo));
    }

    /// Where a label's occurrences are, for [`Self::occ`] (an empty range: none).
    #[inline]
    pub fn occ_range(&self, label: u64, pairs: bool) -> (u32, u32) {
        let (labels, dir, off, pos, len, bits) = self.tables(pairs);
        let b = (label >> (64 - bits)) as usize;
        let (lo, hi) = (dir[b] as usize, dir[b + 1] as usize);
        match labels[lo..hi].binary_search(&label) {
            Ok(i) => {
                let (s, e) = (off[lo + i], off[lo + i + 1]);
                prefetch(pos.as_ptr().wrapping_add(s as usize));
                prefetch(len.as_ptr().wrapping_add(s as usize));
                (s, e)
            }
            Err(_) => (0, 0),
        }
    }

    /// Occurrences (raw start, raw length) in a range from [`Self::occ_range`].
    #[inline]
    pub fn occ(&self, range: (u32, u32), pairs: bool) -> (&[u32], &[u8]) {
        let (_, _, _, pos, len, _) = self.tables(pairs);
        let (s, e) = (range.0 as usize, range.1 as usize);
        (&pos[s..e], &len[s..e])
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
        // (cores start at different bases, so at most `end - start` of them start in the range)
        let most = (lo + (end.saturating_sub(start)) as usize).min(pos.len());
        let hi = lo + pos[lo..most].partition_point(|&p| p < end);
        (&pos[lo..hi], &id[lo..hi])
    }

    /// Section sizes in bytes, for reports.
    pub fn sizes(&self) -> Vec<(String, usize)> {
        self.store.sizes()
    }
}

/// A cache hint for the line holding `p` (x86-64 only; never faults).
#[inline(always)]
fn prefetch<T>(p: *const T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: prefetching is a hint: it reads nothing and cannot fault, whatever the address.
    unsafe {
        core::arch::x86_64::_mm_prefetch(p as *const i8, core::arch::x86_64::_MM_HINT_T0);
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = p;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A FASTA text of three sequences with runs of N at the places that matter: across the
    /// end of a sequence, across the pieces packed in parallel, and at both ends.
    fn fasta() -> (String, Vec<u8>, Vec<u32>) {
        let mut x = 0x9e37_79b9_7f4a_7c15u64;
        let mut base = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGTacgt"[(x % 8) as usize]
        };
        let mut seqs: Vec<Vec<u8>> = [150_000usize, 70_000, 33]
            .iter()
            .map(|&n| (0..n).map(|_| base()).collect())
            .collect();
        for (s, range) in [
            (0, 0..40),
            (0, 65_530..65_545),
            (0, 149_990..150_000),
            (1, 0..7),
            (1, 500..501),
            (2, 30..33),
        ] {
            seqs[s][range].fill(b'N');
        }
        let (mut text, mut all, mut offsets) = (String::new(), Vec::new(), vec![0u32]);
        for (i, s) in seqs.iter().enumerate() {
            text += &format!(">s{i} description\n");
            for line in s.chunks(61) {
                text += std::str::from_utf8(line).unwrap();
                text += "\r\n";
            }
            all.extend_from_slice(s);
            offsets.push(all.len() as u32);
        }
        (text, all, offsets)
    }

    #[test]
    fn packed_genome() {
        let (text, all, offsets) = fasta();
        let g = Packed::read(text.as_bytes()).unwrap();
        assert_eq!(g.names, ["s0", "s1", "s2"]);
        assert_eq!(g.offsets, offsets);
        assert_eq!(g.words.len(), all.len().div_ceil(32) + 1);
        // Runs of N, found one symbol at a time; a run goes on across the end of a sequence.
        let (mut runs, mut segs) = (Vec::new(), Vec::new());
        for (p, &b) in all.iter().enumerate() {
            let p = p as u32;
            if code(b) > 3 {
                match runs.last_mut() {
                    Some(end) if *end == p => *end += 1,
                    _ => runs.extend([p, p + 1]),
                }
            } else {
                assert_eq!(g.bases(p as usize, p as usize + 1), [b.to_ascii_uppercase()]);
                match segs.last_mut() {
                    Some((_, end)) if *end == p as usize && !offsets.contains(&p) => *end += 1,
                    _ => segs.push((p as usize, p as usize + 1)),
                }
            }
        }
        assert_eq!(g.n_runs, runs);
        assert_eq!(g.segments(), segs);
    }

    #[test]
    fn file_build_is_the_memory_build() {
        let (text, ..) = fasta();
        let dir = std::env::temp_dir().join(format!("spearfish-gidx-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        for hpc in [false, true] {
            let opt = BuildOptions { max_occ: 1000, hpc };
            let mem = GenomeIndex::build_from_reader(text.as_bytes(), &opt).unwrap();
            let path = dir.join("g.gidx");
            let (file, bytes) = GenomeIndex::build_to_file(text.as_bytes(), &opt, &path, |names, _| {
                Ok(Some(names.join(",").into_bytes()))
            })
            .unwrap();
            assert_eq!(bytes, std::fs::metadata(&path).unwrap().len());
            assert_eq!(file.manifest(), Some(&b"s0,s1,s2"[..]));
            assert_eq!(mem.manifest(), None);
            let mut sizes = file.store.sizes();
            sizes.retain(|(name, _)| name != "g_manifest");
            sizes.sort();
            let mut other = mem.store.sizes();
            other.sort();
            assert_eq!(sizes, other);
            for (name, _) in &sizes {
                assert_eq!(
                    file.store.get::<u8>(name).unwrap(),
                    mem.store.get::<u8>(name).unwrap(),
                    "{name}"
                );
            }
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
