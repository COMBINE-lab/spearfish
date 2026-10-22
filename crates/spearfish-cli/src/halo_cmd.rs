//! `spearfish halo`: junction halos and isoform resolution per seed scheme (milestone M4,
//! pulled forward to decide the seeding scheme).
//!
//! Seeds are keyed by (label, genomic position of their first transcript base), so the same
//! seed in two isoforms gets the same key, as in the core graph. Three outputs, written to
//! `<prefix>.{halo,iso,genes}.tsv`:
//!
//! - **halo**: for seeds inside a segment shared by >= 2 isoforms, the fraction present in
//!   all of them (isoform-invariant), by distance to the nearest alternative boundary
//!   (the Lemma's halo profile).
//! - **iso**: for ordered pairs (t1, t2) of overlapping isoforms of a gene, restricted to
//!   their common genomic extent, whether t1 has a seed whose key t2 lacks (a
//!   discriminating seed), error-free and after simulated errors in a read from t1. Binned by
//!   how many exonic bases the pair differs by within the overlap.
//! - **genes**: per gene and seed unit, distinct labels against exonic bases and junctions
//!   (Theorem 1 data).

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result, bail};
use rayon::prelude::*;

use crate::annot::{Segments, TxModel, exonic_in, read_gtf, sym_diff};
use crate::rng::SplitMix64;
use crate::seeding::{SeedSet, UnitSeeds, survivors};
use crate::sim::{profile, simulate};

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Transcript FASTA; names are the IDs up to the first '|' or whitespace, as the index
    /// takes them (GENCODE headers: `tx|gene|...`), and must match the GTF's transcript_id.
    transcripts: PathBuf,
    /// Matching GTF (exon records).
    gtf: PathBuf,
    /// Seed scheme spec (repeatable).
    #[arg(short, long = "scheme")]
    schemes: Vec<String>,
    /// File with one seed scheme spec per line.
    #[arg(long)]
    schemes_file: Option<PathBuf>,
    /// Skip LCP levels below this.
    #[arg(long, default_value_t = 3)]
    min_lcp_level: usize,
    /// HPC settings: off, on.
    #[arg(long, value_delimiter = ',', default_value = "off,on")]
    hpc: Vec<String>,
    /// Error profile for the read simulations.
    #[arg(long, default_value = "ont")]
    profile: String,
    /// Error rates for isoform resolution in reads (0 is always included).
    #[arg(long, value_delimiter = ',', default_value = "0.01,0.02,0.06")]
    rates: Vec<f64>,
    /// Position tolerance (raw bases) for seed survival.
    #[arg(long, default_value_t = 3)]
    tol: u32,
    /// At most this many partner isoforms per transcript (sampled deterministically).
    #[arg(long, default_value_t = 20)]
    max_partners: usize,
    #[arg(long, default_value_t = 1)]
    seed: u64,
    /// Worker threads (0 = all cores).
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
    /// Output prefix.
    #[arg(short, long)]
    output: PathBuf,
}

const DIST_EDGES: [u64; 11] = [0, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300];
const DIFF_BINS: [(&str, u64, u64); 6] = [
    ("0_ends_only", 0, 0),
    ("1-10", 1, 10),
    ("11-30", 11, 30),
    ("31-100", 31, 100),
    ("101-300", 101, 300),
    (">300", 301, u64::MAX),
];

fn dist_bin(d: u64) -> usize {
    DIST_EDGES.partition_point(|&e| e <= d) - 1
}

fn diff_bin(d: u64) -> usize {
    DIFF_BINS
        .iter()
        .position(|&(_, lo, hi)| lo <= d && d <= hi)
        .unwrap()
}

struct Dims {
    hpc: usize,
    units: usize,
    rates: usize,
}

#[derive(Clone)]
struct Acc {
    /// [hpc][unit][side][dist bin] -> (occurrences, invariant).
    halo: Vec<(u64, u64)>,
    /// [rate][hpc][unit][diff bin] -> (pairs, resolved, discriminating seeds).
    iso: Vec<(u64, u64, u64)>,
}

impl Acc {
    fn new(d: &Dims) -> Self {
        Acc {
            halo: vec![(0, 0); d.hpc * d.units * 2 * DIST_EDGES.len()],
            iso: vec![(0, 0, 0); d.rates * d.hpc * d.units * DIFF_BINS.len()],
        }
    }

    fn merge(mut self, o: Acc) -> Acc {
        for (a, b) in self.halo.iter_mut().zip(o.halo) {
            a.0 += b.0;
            a.1 += b.1;
        }
        for (a, b) in self.iso.iter_mut().zip(o.iso) {
            a.0 += b.0;
            a.1 += b.1;
            a.2 += b.2;
        }
        self
    }
}

fn halo_idx(d: &Dims, h: usize, u: usize, side: usize, bin: usize) -> usize {
    ((h * d.units + u) * 2 + side) * DIST_EDGES.len() + bin
}

fn iso_idx(d: &Dims, r: usize, h: usize, u: usize, bin: usize) -> usize {
    ((r * d.hpc + h) * d.units + u) * DIFF_BINS.len() + bin
}

/// Per HPC setting, per transcript, per unit: the transcript's seed keys.
type KeySets = Vec<Vec<Vec<HashSet<(u64, u64)>>>>;

struct Tx {
    id: String,
    seq: Vec<u8>,
    model: TxModel,
}

struct Ctx<'a> {
    set: &'a SeedSet,
    hpc: &'a [bool],
    rates: &'a [f64],
    dims: Dims,
    args: &'a Args,
}

/// Analyse one gene; returns its accumulator and its `genes.tsv` rows.
fn gene(cx: &Ctx, gene_id: &str, txs: &[Tx]) -> (Acc, Vec<String>) {
    let d = &cx.dims;
    let mut acc = Acc::new(d);
    let models: Vec<&TxModel> = txs.iter().map(|t| &t.model).collect();
    let segs = Segments::new(&models);
    // seeds[h][t]: per unit, (label, start, end) in transcript coordinates.
    let seeds: Vec<Vec<UnitSeeds>> = cx
        .hpc
        .iter()
        .map(|&h| txs.iter().map(|t| cx.set.seeds(&t.seq, h)).collect())
        .collect();
    // Genomic key of each seed, and per transcript the set of keys.
    let gpos = |t: usize, s: u32| txs[t].model.to_genome(s as u64);
    let keys: KeySets = seeds
        .iter()
        .map(|per_tx| {
            per_tx
                .iter()
                .enumerate()
                .map(|(t, units)| {
                    units
                        .iter()
                        .map(|v| v.iter().map(|&(l, s, _)| (l, gpos(t, s))).collect())
                        .collect()
                })
                .collect()
        })
        .collect();

    // Theorem 1 data.
    let sum_len: u64 = txs.iter().map(|t| t.seq.len() as u64).sum();
    let mut rows = Vec::new();
    for (h, &hpc) in cx.hpc.iter().enumerate() {
        for (u, unit) in cx.set.units.iter().enumerate() {
            let mut labels: Vec<u64> = seeds[h]
                .iter()
                .flat_map(|units| units[u].iter().map(|s| s.0))
                .collect();
            let occ = labels.len();
            labels.sort_unstable();
            labels.dedup();
            rows.push(format!(
                "{gene_id}\t{}\t{}\t{}\t{}\t{sum_len}\t{}\t{}\t{}\t{occ}",
                if hpc { "on" } else { "off" },
                cx.set.name(unit),
                unit.level + 1,
                txs.len(),
                segs.exonic_bases,
                segs.junctions,
                labels.len()
            ));
        }
    }
    if txs.len() < 2 {
        return (acc, rows);
    }

    // Halo profile.
    for h in 0..d.hpc {
        for (t, units) in seeds[h].iter().enumerate() {
            for (u, v) in units.iter().enumerate() {
                for &(l, a, b) in v {
                    let Some(o) = segs.containing(t, a as u64, b as u64) else {
                        continue;
                    };
                    let s = o.seg as usize;
                    let members = &segs.members[s];
                    if members.len() < 2 {
                        continue;
                    }
                    let dl = segs.alt_left[s].then(|| a as u64 - o.tx_start);
                    let dr = segs.alt_right[s].then(|| o.tx_end - b as u64);
                    let (side, dist) = match (dl, dr) {
                        (None, None) => continue,
                        (Some(x), None) => (0, x),
                        (None, Some(y)) => (1, y),
                        (Some(x), Some(y)) => {
                            if x <= y {
                                (0, x)
                            } else {
                                (1, y)
                            }
                        }
                    };
                    let key = (l, gpos(t, a));
                    let inv = members.iter().all(|&m| keys[h][m as usize][u].contains(&key));
                    let e = &mut acc.halo[halo_idx(d, h, u, side, dist_bin(dist))];
                    e.0 += 1;
                    e.1 += inv as u64;
                }
            }
        }
    }

    // Isoform resolution.
    let prof = profile(&cx.args.profile).unwrap();
    let gene_hash = xxhash_rust::xxh3::xxh3_64(gene_id.as_bytes());
    let mut prng = SplitMix64::for_record(cx.args.seed, gene_hash);
    for t1 in 0..txs.len() {
        let (lo1, hi1) = txs[t1].model.span();
        let mut partners: Vec<usize> = (0..txs.len())
            .filter(|&t2| {
                let (lo2, hi2) = txs[t2].model.span();
                t2 != t1 && lo1.max(lo2) < hi1.min(hi2)
            })
            .collect();
        while partners.len() > cx.args.max_partners {
            let i = prng.below(partners.len() as u64) as usize;
            partners.swap_remove(i);
        }
        if partners.is_empty() {
            continue;
        }
        // Survival of t1's seeds in one simulated read per rate (rate 0: all survive).
        let mut rng = SplitMix64::for_record(cx.args.seed, xxhash_rust::xxh3::xxh3_64(txs[t1].id.as_bytes()));
        let reads: Vec<_> = cx
            .rates
            .iter()
            .map(|&r| simulate(&txs[t1].seq, r, &prof, &mut rng))
            .collect();
        // surv[r][h][u][i]
        let surv: Vec<Vec<Vec<Vec<bool>>>> = reads
            .iter()
            .enumerate()
            .map(|(ri, read)| {
                (0..d.hpc)
                    .map(|h| {
                        if cx.rates[ri] == 0.0 {
                            return seeds[h][t1].iter().map(|v| vec![true; v.len()]).collect();
                        }
                        let rs = cx.set.seeds(&read.seq, cx.hpc[h]);
                        (0..d.units)
                            .map(|u| survivors(&seeds[h][t1][u], &rs[u], &read.o2r, cx.args.tol))
                            .collect()
                    })
                    .collect()
            })
            .collect();
        for &t2 in &partners {
            let (lo2, hi2) = txs[t2].model.span();
            let (olo, ohi) = (lo1.max(lo2), hi1.min(hi2));
            let diff = sym_diff(
                &exonic_in(&txs[t1].model, olo, ohi),
                &exonic_in(&txs[t2].model, olo, ohi),
            );
            let bin = diff_bin(diff);
            for h in 0..d.hpc {
                for u in 0..d.units {
                    let disc: Vec<usize> = seeds[h][t1][u]
                        .iter()
                        .enumerate()
                        .filter(|&(_, &(l, a, b))| {
                            let (glo, ghi) = txs[t1].model.genome_span(a as u64, b as u64);
                            glo >= olo && ghi <= ohi && !keys[h][t2][u].contains(&(l, gpos(t1, a)))
                        })
                        .map(|(i, _)| i)
                        .collect();
                    for (r, surv_r) in surv.iter().enumerate() {
                        let e = &mut acc.iso[iso_idx(d, r, h, u, bin)];
                        e.0 += 1;
                        e.1 += disc.iter().any(|&i| surv_r[h][u][i]) as u64;
                        if r == 0 {
                            e.2 += disc.len() as u64;
                        }
                    }
                }
            }
        }
    }
    (acc, rows)
}

pub fn run(args: Args) -> Result<()> {
    let set = SeedSet::load(&args.schemes, args.schemes_file.as_deref(), args.min_lcp_level)?;
    let hpc: Vec<bool> = args
        .hpc
        .iter()
        .map(|h| match h.as_str() {
            "off" => Ok(false),
            "on" => Ok(true),
            _ => bail!("--hpc takes off and/or on, got {h}"),
        })
        .collect::<Result<_>>()?;
    profile(&args.profile).with_context(|| format!("unknown error profile {}", args.profile))?;
    let mut rates = vec![0.0];
    rates.extend(args.rates.iter().copied().filter(|&r| r > 0.0));

    let t0 = Instant::now();
    let models = read_gtf(&args.gtf)?;
    let mut genes: HashMap<String, Vec<Tx>> = HashMap::new();
    let (mut missing, mut mismatched) = (0, 0);
    for (id, _, seq) in spearfish::read_transcripts(&args.transcripts)? {
        let Some(model) = models.get(&id) else {
            missing += 1;
            continue;
        };
        if model.len() != seq.len() as u64 {
            mismatched += 1;
            continue;
        }
        genes.entry(model.gene.clone()).or_default().push(Tx {
            id,
            seq,
            model: model.clone(),
        });
    }
    eprintln!(
        "halo: {} genes, {} transcripts ({missing} not in GTF, {mismatched} length mismatches), loaded in {:.1} s",
        genes.len(),
        genes.values().map(Vec::len).sum::<usize>(),
        t0.elapsed().as_secs_f64()
    );

    let dims = Dims {
        hpc: hpc.len(),
        units: set.units.len(),
        rates: rates.len(),
    };
    let cx = Ctx {
        set: &set,
        hpc: &hpc,
        rates: &rates,
        dims,
        args: &args,
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;
    let mut gene_list: Vec<(&String, &Vec<Tx>)> = genes.iter().collect();
    // Big genes first, for load balance.
    gene_list.sort_by_key(|(_, t)| std::cmp::Reverse(t.len()));
    let t1 = Instant::now();
    let (acc, rows) = pool.install(|| {
        gene_list
            .par_iter()
            .map(|(g, txs)| gene(&cx, g, txs))
            .fold(
                || (Acc::new(&cx.dims), Vec::new()),
                |(a, mut r), (b, rb)| {
                    r.extend(rb);
                    (a.merge(b), r)
                },
            )
            .reduce(
                || (Acc::new(&cx.dims), Vec::new()),
                |(a, mut r), (b, rb)| {
                    r.extend(rb);
                    (a.merge(b), r)
                },
            )
    });
    eprintln!("halo: analysed in {:.1} s", t1.elapsed().as_secs_f64());

    let d = &cx.dims;
    let hname = |h: usize| if hpc[h] { "on" } else { "off" };
    let open = |ext: &str| -> Result<std::io::BufWriter<std::fs::File>> {
        let p = PathBuf::from(format!("{}.{ext}.tsv", args.output.display()));
        Ok(std::io::BufWriter::new(
            std::fs::File::create(&p).with_context(|| format!("creating {}", p.display()))?,
        ))
    };
    let mut f = open("halo")?;
    writeln!(
        f,
        "hpc\tscheme\tlevel\tside\tdist_lo\tdist_hi\toccurrences\tinvariant_frac"
    )?;
    for h in 0..d.hpc {
        for (u, unit) in set.units.iter().enumerate() {
            for side in 0..2 {
                for (b, lo) in DIST_EDGES.iter().enumerate() {
                    let (n, inv) = acc.halo[halo_idx(d, h, u, side, b)];
                    let hi = DIST_EDGES
                        .get(b + 1)
                        .map_or("inf".to_string(), |x| (x - 1).to_string());
                    writeln!(
                        f,
                        "{}\t{}\t{}\t{}\t{}\t{hi}\t{n}\t{:.4}",
                        hname(h),
                        set.name(unit),
                        unit.level + 1,
                        ["left", "right"][side],
                        lo,
                        inv as f64 / n.max(1) as f64
                    )?;
                }
            }
        }
    }
    // Flush each file explicitly: dropping a BufWriter discards write errors.
    f.flush()?;
    let mut f = open("iso")?;
    writeln!(
        f,
        "profile\trate\thpc\tscheme\tlevel\tdiff_bin\tpairs\tresolved_frac\tmean_disc_seeds"
    )?;
    for (r, rate) in rates.iter().enumerate() {
        for h in 0..d.hpc {
            for (u, unit) in set.units.iter().enumerate() {
                for (b, (bname, _, _)) in DIFF_BINS.iter().enumerate() {
                    let (n, res, disc) = acc.iso[iso_idx(d, r, h, u, b)];
                    writeln!(
                        f,
                        "{}\t{rate}\t{}\t{}\t{}\t{bname}\t{n}\t{:.4}\t{:.2}",
                        args.profile,
                        hname(h),
                        set.name(unit),
                        unit.level + 1,
                        res as f64 / n.max(1) as f64,
                        if r == 0 {
                            disc as f64 / n.max(1) as f64
                        } else {
                            f64::NAN
                        }
                    )?;
                }
            }
        }
    }
    f.flush()?;
    let mut f = open("genes")?;
    writeln!(
        f,
        "gene\thpc\tscheme\tlevel\tn_tx\tsum_tx_len\texonic_bases\tjunctions\tdistinct_labels\toccurrences"
    )?;
    for r in rows {
        writeln!(f, "{r}")?;
    }
    f.flush()?;
    Ok(())
}
