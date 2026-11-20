//! Stage 0 probe: can LCP cores seed spliced genome mapping?
//!
//! `genome_probe <genome.fa> <cap> (<reads.fa> <gtruth.tsv>)...`
//!
//! Seeds the genome (HPC, LCP levels 3 and 2), then for simulated reads with known genome
//! blocks reports: seed counts and occurrence spread on the genome; the share of true exon
//! blocks, by length, that keep at least one level-3 or level-2 seed; whether the true locus
//! collects more seed hits than any other 64 kb bin; and how often a short exon is the only
//! match of its length between splice-site dinucleotides inside its intron window.
//!
//! All positions are raw (uncompressed) genome coordinates on one concatenated axis.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

use rayon::prelude::*;
use rustc_hash::FxHashMap;
use seeds::BaseSeeder;

const CHUNK: usize = 8 << 20;
const BINS: [(u32, u32, &str); 5] = [
    (0, 10, "<10"),
    (10, 20, "10-19"),
    (20, 30, "20-29"),
    (30, 60, "30-59"),
    (60, u32::MAX, ">=60"),
];

fn bin_of(n: u32) -> usize {
    BINS.iter().position(|&(lo, hi, _)| lo <= n && n < hi).unwrap()
}

struct Genome {
    seq: Vec<u8>,
    offset: FxHashMap<String, u32>,
}

fn read_genome(path: &str) -> anyhow::Result<Genome> {
    let mut seq = Vec::with_capacity(3_200_000_000);
    let mut offset = FxHashMap::default();
    for line in BufReader::with_capacity(1 << 24, File::open(path)?).lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            offset.insert(h.split_whitespace().next().unwrap().to_string(), seq.len() as u32);
        } else {
            seq.extend(line.bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    anyhow::ensure!(
        seq.len() < u32::MAX as usize,
        "genome too long for 32-bit positions"
    );
    Ok(Genome { seq, offset })
}

/// Maximal ACGT stretches, cut into chunks of at most `CHUNK` bases.
fn segments(seq: &[u8]) -> Vec<(usize, usize)> {
    let mut v = Vec::new();
    let mut i = 0;
    while i < seq.len() {
        if !matches!(seq[i], b'A' | b'C' | b'G' | b'T') {
            i += 1;
            continue;
        }
        let s = i;
        while i < seq.len() && matches!(seq[i], b'A' | b'C' | b'G' | b'T') {
            i += 1;
        }
        let mut a = s;
        while a < i {
            let b = (a + CHUNK).min(i);
            v.push((a, b));
            a = b;
        }
    }
    v
}

/// How seeds are formed from the LCP cores of one sequence.
#[derive(Clone, Copy, PartialEq)]
enum Kind {
    L3,
    L2,
    /// Two or three consecutive level-2 cores as one seed, placed at the first core.
    L2Run(usize),
}

impl Kind {
    fn name(self) -> &'static str {
        match self {
            Kind::L3 => "L3",
            Kind::L2 => "L2",
            Kind::L2Run(2) => "L2x2",
            Kind::L2Run(_) => "L2x3",
        }
    }
}

fn mix(a: u64, b: u64) -> u64 {
    (a.rotate_left(23) ^ b).wrapping_mul(0x9E37_79B9_7F4A_7C15)
}

/// (label, raw start) of the seeds of `raw`, in start order. `s3`/`s2` are the level-3 and
/// level-2 seeders.
fn seeds_raw(kind: Kind, s3: &BaseSeeder, s2: &BaseSeeder, raw: &[u8], base: u32, out: &mut Vec<(u64, u32)>) {
    let sd = if kind == Kind::L3 { s3 } else { s2 };
    let h = sd.parsed(raw).unwrap();
    let v = sd.seed_parsed(&h.seq);
    match kind {
        Kind::L3 | Kind::L2 => out.extend(v.iter().map(|s| (s.label, base + h.to_raw(s.start, s.end).0))),
        Kind::L2Run(n) => out.extend(v.windows(n).map(|w| {
            let l = w.iter().fold(0u64, |a, s| mix(a, s.label));
            (l, base + h.to_raw(w[0].start, w[0].end).0)
        })),
    }
}

/// Every seed of the genome, sorted by (label, start).
fn seed_genome(g: &Genome, kind: Kind) -> Vec<(u64, u32)> {
    let s3 = BaseSeeder::new("lcp:K=3", true).unwrap();
    let s2 = BaseSeeder::new("lcp:K=2", true).unwrap();
    let segs = segments(&g.seq);
    let mut all: Vec<(u64, u32)> = segs
        .par_iter()
        .flat_map_iter(|&(a, b)| {
            let mut v = Vec::new();
            seeds_raw(kind, &s3, &s2, &g.seq[a..b], a as u32, &mut v);
            v
        })
        .collect();
    all.par_sort_unstable();
    all
}

fn occ(t: &[(u64, u32)], label: u64) -> &[(u64, u32)] {
    let lo = t.partition_point(|x| x.0 < label);
    let hi = lo + t[lo..].partition_point(|x| x.0 == label);
    &t[lo..hi]
}

fn occupancy_report(name: &str, t: &[(u64, u32)]) {
    let mut distinct = 0u64;
    let caps = [1usize, 4, 16, 64, 256, 1024, 4096];
    let mut over = [0u64; 7];
    let mut labels_over = [0u64; 7];
    let mut i = 0;
    while i < t.len() {
        let mut j = i + 1;
        while j < t.len() && t[j].0 == t[i].0 {
            j += 1;
        }
        distinct += 1;
        for (k, &c) in caps.iter().enumerate() {
            if j - i > c {
                over[k] += (j - i) as u64;
                labels_over[k] += 1;
            }
        }
        i = j;
    }
    println!("genome\t{name}\tseeds\t{}", t.len());
    println!("genome\t{name}\tdistinct_labels\t{distinct}");
    for (k, &c) in caps.iter().enumerate() {
        println!(
            "genome\t{name}\tocc>{c}\tlabels={}\tshare_of_seeds={:.4}",
            labels_over[k],
            over[k] as f64 / t.len() as f64
        );
    }
}

struct Truth {
    gene_minus: bool,
    blocks: Vec<(u32, u32)>,
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        })
        .collect()
}

const KINDS: [Kind; 4] = [Kind::L3, Kind::L2, Kind::L2Run(2), Kind::L2Run(3)];

#[derive(Default, Clone)]
struct Acc {
    reads: u64,
    read_bases: u64,
    seeds: [u64; 4],
    cand: [u64; 4],       // occurrences under the cap, both orientations
    any: [[u64; 4]; 4],   // per kind: reads with >= 1, 2, 3, 5 true anchors
    win: [u64; 4],        // per kind: true anchors (weighted) beat the best other bin
    win_combo: u64,       // L3 + L2x3 together
    any_combo: [u64; 4],  // L3 + L2x3: reads with >= 1, 2, 3, 5 true anchors
    blk_n: [[u64; 5]; 2], // [internal, read end][length bin]
    blk_hit: [[[u64; 5]; 2]; 4],
    blk_hit_l3_or: [[[u64; 5]; 2]; 4],
    no_l3_short: u64, // reads without a true level-3 anchor that are under 300 bases
    no_l3: u64,
}

impl Acc {
    fn add(mut self, o: &Acc) -> Acc {
        self.reads += o.reads;
        self.read_bases += o.read_bases;
        self.win_combo += o.win_combo;
        self.no_l3 += o.no_l3;
        self.no_l3_short += o.no_l3_short;
        for k in 0..4 {
            self.seeds[k] += o.seeds[k];
            self.cand[k] += o.cand[k];
            self.win[k] += o.win[k];
            self.any_combo[k] += o.any_combo[k];
            for i in 0..4 {
                self.any[k][i] += o.any[k][i];
            }
            for e in 0..2 {
                for b in 0..5 {
                    self.blk_hit[k][e][b] += o.blk_hit[k][e][b];
                    self.blk_hit_l3_or[k][e][b] += o.blk_hit_l3_or[k][e][b];
                }
            }
        }
        for e in 0..2 {
            for b in 0..5 {
                self.blk_n[e][b] += o.blk_n[e][b];
            }
        }
        self
    }
}

/// One orientation of one read against one seed table.
struct Side {
    per: Vec<u32>, // true anchors per block
    tot: u32,
    wtrue: f64,                 // sum of 1/occ over true anchors
    votes: FxHashMap<u32, f64>, // 64 kb bin -> sum of 1/occ, true bins excluded
    cand: u64,
    n: u64,
}

fn side(t: &[(u64, u32)], seeds: &[(u64, u32)], tr: &Truth, cap: usize) -> Side {
    let nb = tr.blocks.len();
    let (ts, te) = (tr.blocks[0].0, tr.blocks[nb - 1].1);
    let mut pre = Vec::with_capacity(nb);
    let mut c = 0u32;
    for &(x, y) in &tr.blocks {
        pre.push(c);
        c += y - x;
    }
    let mut s = Side {
        per: vec![0; nb],
        tot: 0,
        wtrue: 0.0,
        votes: FxHashMap::default(),
        cand: 0,
        n: seeds.len() as u64,
    };
    for &(label, q) in seeds {
        let o = occ(t, label);
        let capped = o.len() <= cap;
        let w = 1.0 / o.len().max(1) as f64;
        if capped {
            s.cand += o.len() as u64;
            for &(_, p) in o {
                if !(p + (1 << 16) >= ts && p < te + (1 << 16)) {
                    *s.votes.entry(p >> 16).or_default() += w;
                }
            }
        }
        let lo = o.partition_point(|x| x.1 < ts);
        for &(_, p) in &o[lo..] {
            if p >= te {
                break;
            }
            let b = tr.blocks.partition_point(|x| x.1 <= p);
            if b < nb && tr.blocks[b].0 <= p {
                // the read position this genome position should have, give or take indels
                let expect = (pre[b] + p - tr.blocks[b].0) as i64;
                if (q as i64 - expect).abs() <= 50 + expect / 20 {
                    s.per[b] += 1;
                    s.tot += 1;
                    if capped {
                        s.wtrue += w;
                    }
                }
            }
        }
    }
    s
}

fn main() -> anyhow::Result<()> {
    let a: Vec<String> = std::env::args().collect();
    anyhow::ensure!(
        a.len() >= 4,
        "usage: genome_probe <genome.fa> <cap> (<reads.fa> <gtruth.tsv>)..."
    );
    let cap: usize = a[2].parse()?;
    let t0 = Instant::now();
    let g = read_genome(&a[1])?;
    eprintln!(
        "genome: {} bases, {:.1}s",
        g.seq.len(),
        t0.elapsed().as_secs_f64()
    );
    println!("genome\tthreads\t{}", rayon::current_num_threads());
    let tables: Vec<Vec<(u64, u32)>> = KINDS
        .iter()
        .map(|&k| {
            let t = Instant::now();
            let v = seed_genome(&g, k);
            println!(
                "genome\t{}\tseed_and_sort_seconds\t{:.1}",
                k.name(),
                t.elapsed().as_secs_f64()
            );
            occupancy_report(k.name(), &v);
            v
        })
        .collect();
    for pair in a[3..].chunks(2) {
        run_set(&g, &tables, &pair[0], &pair[1], cap)?;
    }
    Ok(())
}

fn run_set(
    g: &Genome,
    tables: &[Vec<(u64, u32)>],
    reads_fa: &str,
    gtruth: &str,
    cap: usize,
) -> anyhow::Result<()> {
    let set = std::path::Path::new(reads_fa)
        .file_stem()
        .unwrap()
        .to_string_lossy()
        .to_string();
    // Truth: read -> blocks on the concatenated axis.
    let mut truth: FxHashMap<String, Truth> = FxHashMap::default();
    for line in BufReader::new(File::open(gtruth)?).lines() {
        let line = line?;
        let p: Vec<&str> = line.split('\t').collect();
        let Some(&off) = g.offset.get(p[1]) else { continue };
        let blocks = p[4]
            .split(',')
            .map(|b| {
                let (x, y) = b.split_once('-').unwrap();
                (off + x.parse::<u32>().unwrap(), off + y.parse::<u32>().unwrap())
            })
            .collect();
        truth.insert(
            p[0].to_string(),
            Truth {
                gene_minus: p[2] == "-",
                blocks,
            },
        );
    }
    let mut reads: Vec<(String, Vec<u8>)> = Vec::new();
    let mut name = String::new();
    for line in BufReader::new(File::open(reads_fa)?).lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            name = h.split_whitespace().next().unwrap().to_string();
        } else if truth.contains_key(&name) {
            reads.push((name.clone(), line.into_bytes()));
        }
    }
    eprintln!("{set}: reads with truth: {}", reads.len());

    let s3 = BaseSeeder::new("lcp:K=3", true).unwrap();
    let s2 = BaseSeeder::new("lcp:K=2", true).unwrap();
    let acc = reads
        .par_iter()
        .fold(Acc::default, |mut acc, (name, seq)| {
            let tr = &truth[name];
            let nb = tr.blocks.len();
            let rc = revcomp(seq);
            let mut sides: Vec<[Side; 4]> = Vec::new();
            for s in [seq.as_slice(), rc.as_slice()] {
                let mut v = Vec::new();
                sides.push(std::array::from_fn(|k| {
                    v.clear();
                    seeds_raw(KINDS[k], &s3, &s2, s, 0, &mut v);
                    side(&tables[k], &v, tr, cap)
                }));
            }
            // the orientation that reads along the genome's forward strand has the true anchors
            let f = if sides[0][0].tot + sides[0][3].tot >= sides[1][0].tot + sides[1][3].tot {
                0
            } else {
                1
            };
            let (fw, rv) = (&sides[f], &sides[1 - f]);
            acc.reads += 1;
            acc.read_bases += seq.len() as u64;
            let best_other = |k: usize| {
                let m = |s: &Side| s.votes.values().cloned().fold(0.0f64, f64::max);
                // the other orientation's hits near the true locus also compete
                m(&fw[k]).max(m(&rv[k])).max(rv[k].wtrue)
            };
            for k in 0..4 {
                acc.seeds[k] += fw[k].n;
                acc.cand[k] += fw[k].cand + rv[k].cand;
                for (i, &m) in [1u32, 2, 3, 5].iter().enumerate() {
                    acc.any[k][i] += (fw[k].tot >= m) as u64;
                }
                acc.win[k] += (fw[k].wtrue > best_other(k)) as u64;
            }
            // level 3 and level-2 triples in one vote
            let mut other: FxHashMap<u32, f64> = FxHashMap::default();
            let mut om = rv[0].wtrue + rv[3].wtrue;
            for s in [fw, rv] {
                other.clear();
                for k in [0, 3] {
                    for (&b, &w) in &s[k].votes {
                        *other.entry(b).or_default() += w;
                    }
                }
                om = other.values().cloned().fold(om, f64::max);
            }
            acc.win_combo += (fw[0].wtrue + fw[3].wtrue > om) as u64;
            for (i, &m) in [1u32, 2, 3, 5].iter().enumerate() {
                acc.any_combo[i] += (fw[0].tot + fw[3].tot >= m) as u64;
            }
            if fw[0].tot == 0 {
                acc.no_l3 += 1;
                acc.no_l3_short += (seq.len() < 300) as u64;
            }
            for (i, &(x, y)) in tr.blocks.iter().enumerate() {
                let b = bin_of(y - x);
                let e = (i == 0 || i + 1 == nb) as usize;
                acc.blk_n[e][b] += 1;
                for k in 0..4 {
                    acc.blk_hit[k][e][b] += (fw[k].per[i] > 0) as u64;
                    acc.blk_hit_l3_or[k][e][b] += (fw[k].per[i] > 0 || fw[0].per[i] > 0) as u64;
                }
            }
            acc
        })
        .reduce(Acc::default, |x, y| x.add(&y));

    let n = acc.reads as f64;
    println!("{set}\treads\tn\t{}", acc.reads);
    println!(
        "{set}\treads\twithout_true_L3_anchor\t{}\tof_which_under_300_bases\t{}",
        acc.no_l3, acc.no_l3_short
    );
    for (k, kind) in KINDS.iter().enumerate() {
        let nm = kind.name();
        println!(
            "{set}\treads\t{nm}\tseeds_per_kb\t{:.2}",
            1000.0 * acc.seeds[k] as f64 / acc.read_bases as f64
        );
        println!(
            "{set}\treads\t{nm}\tcandidate_hits_per_read_cap{cap}\t{:.1}",
            acc.cand[k] as f64 / n
        );
        println!(
            "{set}\treads\t{nm}\ttrue_anchors>=1,2,3,5\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
            acc.any[k][0] as f64 / n,
            acc.any[k][1] as f64 / n,
            acc.any[k][2] as f64 / n,
            acc.any[k][3] as f64 / n
        );
        println!(
            "{set}\treads\t{nm}\ttrue_locus_outweighs_every_other_bin\t{:.4}",
            acc.win[k] as f64 / n
        );
    }
    println!(
        "{set}\treads\tL3+L2x3\ttrue_anchors>=1,2,3,5\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
        acc.any_combo[0] as f64 / n,
        acc.any_combo[1] as f64 / n,
        acc.any_combo[2] as f64 / n,
        acc.any_combo[3] as f64 / n
    );
    println!(
        "{set}\treads\tL3+L2x3\ttrue_locus_outweighs_every_other_bin\t{:.4}",
        acc.win_combo as f64 / n
    );
    for (e, what) in ["internal_exon", "read_end_block"].iter().enumerate() {
        for (b, &(_, _, bn)) in BINS.iter().enumerate() {
            let m = acc.blk_n[e][b].max(1) as f64;
            let mut line = format!("{set}\tblocks\t{what}\t{bn}\tn={}", acc.blk_n[e][b]);
            for (k, kind) in KINDS.iter().enumerate() {
                line += &format!("\t{}={:.4}", kind.name(), acc.blk_hit[k][e][b] as f64 / m);
            }
            for (k, kind) in KINDS.iter().enumerate().skip(1) {
                line += &format!(
                    "\tL3_or_{}={:.4}",
                    kind.name(),
                    acc.blk_hit_l3_or[k][e][b] as f64 / m
                );
            }
            println!("{line}");
        }
    }

    // Short internal exons: is the true exon the only match of its length between an acceptor
    // and a donor dinucleotide inside its intron window? Mismatches allowed: length / 10.
    let mut exons: HashSet<(u32, u32, u32, u32, bool)> = HashSet::new();
    for (name, _) in &reads {
        let tr = &truth[name];
        for i in 1..tr.blocks.len().saturating_sub(1) {
            let (x, y) = tr.blocks[i];
            if y - x < 60 {
                exons.insert((x, y, tr.blocks[i - 1].1, tr.blocks[i + 1].0, tr.gene_minus));
            }
        }
    }
    let exons: Vec<_> = exons.into_iter().collect();
    let t = Instant::now();
    // per bin: n, canonical, unique among canonical pairs, total false hits, window bases
    let res = exons
        .par_iter()
        .fold(
            || [[0u64; 5]; 5],
            |mut r, &(x, y, ws, we, minus)| {
                let k = (y - x) as usize;
                let b = bin_of(k as u32);
                let (acc_m, don_m): (&[u8], &[u8]) = if minus { (b"AC", b"CT") } else { (b"AG", b"GT") };
                let s = &g.seq;
                let ex = &s[x as usize..y as usize];
                let mm = k / 10;
                let canon =
                    &s[x as usize - 2..x as usize] == acc_m && &s[y as usize..y as usize + 2] == don_m;
                let mut hits = 0u64;
                let (lo, hi) = (ws as usize + 2, (we as usize).saturating_sub(k + 2));
                for p in lo..=hi.max(lo) {
                    if p == x as usize || p + k + 2 > s.len() {
                        continue;
                    }
                    if &s[p - 2..p] == acc_m && &s[p + k..p + k + 2] == don_m {
                        let d = s[p..p + k].iter().zip(ex).filter(|(a, b)| a != b).count();
                        hits += (d <= mm) as u64;
                    }
                }
                r[b][0] += 1;
                r[b][1] += canon as u64;
                r[b][2] += (canon && hits == 0) as u64;
                r[b][3] += hits;
                r[b][4] += (we - ws) as u64;
                r
            },
        )
        .reduce(
            || [[0u64; 5]; 5],
            |mut x, y| {
                for i in 0..5 {
                    for j in 0..5 {
                        x[i][j] += y[i][j];
                    }
                }
                x
            },
        );
    let scan = t.elapsed().as_secs_f64();
    for (i, &(_, _, b)) in BINS.iter().enumerate().take(4) {
        let r = res[i];
        let n = r[0].max(1) as f64;
        println!(
            "{set}\tshort_exons\t{b}\tn={}\tcanonical={:.4}\tunique_between_sites={:.4}\tfalse_hits_per_exon={:.3}\tmean_window={:.0}",
            r[0],
            r[1] as f64 / n,
            r[2] as f64 / n,
            r[3] as f64 / n,
            r[4] as f64 / n
        );
    }
    let tot_window: u64 = res.iter().map(|r| r[4]).sum();
    println!(
        "{set}\tshort_exons\tnaive_scan_ns_per_window_base\t{:.2}",
        scan * 1e9 * rayon::current_num_threads() as f64 / tot_window.max(1) as f64
    );
    Ok(())
}
