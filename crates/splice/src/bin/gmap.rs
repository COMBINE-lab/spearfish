//! Development driver for spliced genome mapping.
//!
//! ```text
//! gmap index <genome.fa> <out.gidx> [--hpc]
//! gmap map <genome.gidx> <reads.fa|fq> [--threads N] [--forward] [--cap N] [--secondary]
//!          [--junctions <BED12 or GTF>] [--sam | --bam [--bam-level N]] [--cigar [--affine]]   (default output: PAF)
//! ```

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use rayon::prelude::*;
use spearfish_splice::bam;
use spearfish_splice::index::{BuildOptions, GenomeIndex};
use spearfish_splice::map::{MapOptions, Mapper, Scratch, Stats, paf_line, sam_header, sam_line};

#[derive(Clone, Copy)]
enum Format {
    Paf,
    Sam,
    Bam(bam::Level),
}

fn main() -> anyhow::Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let flag = |f: &str| a.iter().position(|x| x == f);
    let val = |f: &str| flag(f).map(|i| a[i + 1].clone());
    let threads: usize = val("--threads").map_or(8, |s| s.parse().unwrap());
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?;
    match a.get(1).map(String::as_str) {
        Some("index") => {
            let t = Instant::now();
            let ix = GenomeIndex::build(
                Path::new(&a[2]),
                &BuildOptions {
                    hpc: flag("--hpc").is_some(),
                    ..Default::default()
                },
            )?;
            eprintln!("built in {:.1}s", t.elapsed().as_secs_f64());
            let n = ix.write(Path::new(&a[3]))?;
            eprintln!("wrote {n} bytes in total {:.1}s", t.elapsed().as_secs_f64());
            for (name, size) in ix.sizes() {
                eprintln!("  {name}\t{:.1} MB", size as f64 / 1e6);
            }
        }
        Some("map") => {
            let t = Instant::now();
            let mut ix = GenomeIndex::open(Path::new(&a[2]))?;
            if let Some(bed) = val("--junctions") {
                let j = spearfish_splice::junctions::Junctions::from_file(Path::new(&bed), &ix)?;
                eprintln!("{} known junctions", j.len());
                anyhow::ensure!(
                    !j.is_empty(),
                    "no junction in {bed} lies on a sequence of this index"
                );
                ix.junctions = Some(j);
            }
            let mut opt = MapOptions {
                forward_only: flag("--forward").is_some(),
                pairs: flag("--no-pairs").is_none(),
                ..Default::default()
            };
            if let Some(c) = val("--cap") {
                opt.cap = c.parse()?;
            }
            let secondary = flag("--secondary").is_some();
            // --cigar: base-level alignment of the blocks; --affine: with affine gap costs
            let affine = flag("--affine").is_some();
            let base_level = flag("--cigar").is_some() || affine;
            let mode = if affine {
                spearfish_splice::base::Mode::Affine
            } else {
                spearfish_splice::base::Mode::Edits
            };
            let format = if flag("--bam").is_some() {
                Format::Bam(bam::Level(val("--bam-level").map_or(1, |s| s.parse().unwrap())))
            } else if flag("--sam").is_some() {
                Format::Sam
            } else {
                Format::Paf
            };
            opt.keep_seeds = base_level;
            let ix = std::sync::Arc::new(ix);
            let mapper = Mapper::new(ix.clone(), opt);
            let mut out = BufWriter::with_capacity(1 << 22, std::io::stdout().lock());
            match format {
                Format::Sam => {
                    let mut h = Vec::new();
                    sam_header(&ix, None, &mut h);
                    out.write_all(&h)?;
                }
                Format::Bam(level) => {
                    let (mut h, mut z) = (Vec::new(), Vec::new());
                    bam::header(&ix, None, &mut h);
                    bam::bgzf(&h, level, &mut z);
                    out.write_all(&z)?;
                }
                Format::Paf => {}
            }
            let mut rd = BufReader::with_capacity(1 << 22, File::open(&a[3])?);
            let mut batch: Vec<(String, Vec<u8>)> = Vec::new();
            let mut line = String::new();
            let (mut nreads, mut nmapped) = (0u64, 0u64);
            let mut stats = Stats::default();
            let mut done = false;
            while !done {
                batch.clear();
                while batch.len() < 16384 {
                    line.clear();
                    if rd.read_line(&mut line)? == 0 {
                        done = true;
                        break;
                    }
                    let fastq = line.starts_with('@');
                    let name = line[1..].split_whitespace().next().unwrap_or("").to_string();
                    line.clear();
                    rd.read_line(&mut line)?;
                    let seq = line.trim_end().as_bytes().to_vec();
                    if fastq {
                        line.clear();
                        rd.read_line(&mut line)?;
                        line.clear();
                        rd.read_line(&mut line)?;
                    }
                    batch.push((name, seq));
                }
                // Reads are mapped in groups; each group's records are formatted, and for BAM
                // compressed into complete BGZF blocks, by the thread that mapped them.
                let outs: Vec<(Vec<u8>, u64, Stats)> = batch
                    .par_chunks(256)
                    .map_init(
                        <(Scratch, spearfish_splice::base::Aligner)>::default,
                        |(sc, al), group| {
                            sc.stats = Stats::default();
                            let mut text: Vec<u8> = Vec::new();
                            let mut raw: Vec<u8> = Vec::new();
                            let mut ops: Vec<u32> = Vec::new();
                            let mut mapped = 0u64;
                            for (name, seq) in group {
                                let hits = mapper.map(seq, sc);
                                let mapq = spearfish_splice::map::mapq(&hits);
                                mapped += !hits.is_empty() as u64;
                                for h in hits.iter().take(if secondary { 5 } else { 1 }) {
                                    al.mode = mode;
                                    // (secondary records keep the block-level alignment)
                                    let edits =
                                        (base_level && h.primary).then(|| al.cigar(&ix, seq, h, &mut ops));
                                    let base = edits.map(|e| (&ops[..], e));
                                    match format {
                                        Format::Bam(_) => {
                                            bam::record(&ix, name, seq, h, mapq, &mut ops, edits, &mut raw)
                                        }
                                        Format::Sam => sam_line(&ix, name, seq, h, mapq, base, &mut text),
                                        Format::Paf => {
                                            paf_line(&ix, name, seq.len() as u32, h, mapq, base, &mut text)
                                        }
                                    }
                                }
                            }
                            let bytes = match format {
                                Format::Bam(level) => {
                                    let mut z = Vec::with_capacity(raw.len() / 3);
                                    bam::bgzf(&raw, level, &mut z);
                                    z
                                }
                                _ => text,
                            };
                            (bytes, mapped, sc.stats)
                        },
                    )
                    .collect();
                for (bytes, m, st) in &outs {
                    out.write_all(bytes)?;
                    nmapped += *m;
                    stats.add(st);
                }
                nreads += batch.len() as u64;
            }
            if matches!(format, Format::Bam(_)) {
                out.write_all(&bam::EOF)?;
            }
            out.flush()?;
            eprintln!(
                "mapped {nmapped} of {nreads} reads in {:.1}s",
                t.elapsed().as_secs_f64()
            );
            eprintln!("{stats:?}");
            if base_level {
                use std::sync::atomic::Ordering::Relaxed;
                let (b, g) = (
                    spearfish_splice::base::BASES.load(Relaxed),
                    spearfish_splice::base::AGAIN.load(Relaxed),
                );
                eprintln!("base level: {b} block bases, {g} aligned again with affine gaps");
            }
        }
        _ => anyhow::bail!(
            "usage: gmap index <genome.fa> <out.gidx> | gmap map <genome.gidx> <reads> [--threads N]"
        ),
    }
    Ok(())
}
