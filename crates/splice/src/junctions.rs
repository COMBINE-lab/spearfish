//! Known splice junctions and exons (optional input).
//!
//! Read from a BED12 file of transcript models (as `paftools.js gff2bed` writes) or a GTF. Positions
//! are on the index's concatenated axis. A known splice site costs nothing where an unknown
//! one pays for its motif and a little more for being unknown; known exons and the
//! donor-to-acceptor links are extra candidates where read bases are left unexplained.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};

use crate::index::GenomeIndex;

/// One transcript while reading: sequence offset, reverse strand, exons.
type Model = (u32, bool, Vec<(u32, u32)>);

#[derive(Default)]
pub struct Junctions {
    /// Sorted first-intron-base positions, per transcript strand (0 forward, 1 reverse).
    donors: [Vec<u32>; 2],
    /// Sorted first-exon-base positions after an intron, per transcript strand.
    acceptors: [Vec<u32>; 2],
    /// (left end, right end) of every intron on the genome, sorted, and the same by right end.
    by_left: Vec<(u32, u32)>,
    by_right: Vec<(u32, u32)>,
    /// Exons (start, end), sorted by start.
    exons: Vec<(u32, u32)>,
}

impl Junctions {
    /// Read transcript models from a BED12 file, or from a GTF (by its `.gtf` extension: `exon`
    /// records grouped by `transcript_id`).
    pub fn from_file(path: &Path, ix: &GenomeIndex) -> Result<Junctions> {
        let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let offset: std::collections::HashMap<&str, u32> = ix
            .names()
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), ix.offsets()[i]))
            .collect();
        let gtf = path.extension().is_some_and(|e| e.eq_ignore_ascii_case("gtf"));
        // transcripts as (sequence offset, reverse strand, exons)
        let mut models: Vec<Model> = Vec::new();
        let mut by_id: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
        for line in BufReader::new(f).lines() {
            let line = line?;
            let p: Vec<&str> = line.split('\t').collect();
            if gtf {
                if p.len() < 9 || p[2] != "exon" {
                    continue;
                }
                let Some(&off) = offset.get(p[0]) else { continue };
                let Some(id) = p[8]
                    .split("transcript_id \"")
                    .nth(1)
                    .and_then(|x| x.split('"').next())
                else {
                    continue;
                };
                let (s, e): (u32, u32) = (p[3].parse()?, p[4].parse()?);
                let k = *by_id.entry(id.to_string()).or_insert_with(|| {
                    models.push((off, p[6] == "-", Vec::new()));
                    models.len() - 1
                });
                models[k].2.push((off + s - 1, off + e));
            } else {
                if p.len() < 12 {
                    continue;
                }
                let Some(&off) = offset.get(p[0]) else { continue };
                let start: u32 = p[1].parse()?;
                let sizes = p[10]
                    .split(',')
                    .filter(|s| !s.is_empty())
                    .map(|s| s.parse::<u32>().unwrap_or(0));
                let starts = p[11]
                    .split(',')
                    .filter(|s| !s.is_empty())
                    .map(|s| s.parse::<u32>().unwrap_or(0));
                let exons = sizes
                    .zip(starts)
                    .map(|(z, b)| (off + start + b, off + start + b + z))
                    .collect();
                models.push((off, p[5] == "-", exons));
            }
        }
        let mut j = Junctions::default();
        for (_, minus, mut exons) in models {
            exons.sort_unstable();
            for (k, &(s, e)) in exons.iter().enumerate() {
                j.exons.push((s, e));
                if k > 0 && s > exons[k - 1].1 {
                    let left = exons[k - 1].1;
                    j.by_left.push((left, s));
                    // the cost functions take genome sides (left end, right end of the intron),
                    // whichever of them is the donor on the transcript's strand
                    j.donors[minus as usize].push(left);
                    j.acceptors[minus as usize].push(s);
                }
            }
        }
        for v in j.donors.iter_mut().chain(j.acceptors.iter_mut()) {
            v.sort_unstable();
            v.dedup();
        }
        j.by_left.sort_unstable();
        j.by_left.dedup();
        j.by_right = j.by_left.iter().map(|&(l, r)| (r, l)).collect();
        j.by_right.sort_unstable();
        j.exons.sort_unstable();
        j.exons.dedup();
        Ok(j)
    }

    pub fn len(&self) -> usize {
        self.by_left.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_left.is_empty()
    }

    /// Whether an intron of a transcript on the given strand starts (on the genome) at `p`.
    #[inline]
    pub fn left_site(&self, p: u32, minus: bool) -> bool {
        self.donors[minus as usize].binary_search(&p).is_ok()
    }

    /// Whether an intron of a transcript on the given strand ends (on the genome) at `p`.
    #[inline]
    pub fn right_site(&self, p: u32, minus: bool) -> bool {
        self.acceptors[minus as usize].binary_search(&p).is_ok()
    }

    /// Known exons that lie inside `[lo, hi)`.
    pub fn exons_in(&self, lo: u32, hi: u32) -> &[(u32, u32)] {
        let a = self.exons.partition_point(|e| e.0 < lo);
        let b = a + self.exons[a..].partition_point(|e| e.0 < hi);
        &self.exons[a..b]
    }

    /// Right ends of the known introns whose left end is in `[lo, hi)`: (left, right).
    pub fn from_left(&self, lo: u32, hi: u32) -> &[(u32, u32)] {
        let a = self.by_left.partition_point(|e| e.0 < lo);
        let b = a + self.by_left[a..].partition_point(|e| e.0 < hi);
        &self.by_left[a..b]
    }

    /// Left ends of the known introns whose right end is in `[lo, hi)`: (right, left).
    pub fn from_right(&self, lo: u32, hi: u32) -> &[(u32, u32)] {
        let a = self.by_right.partition_point(|e| e.0 < lo);
        let b = a + self.by_right[a..].partition_point(|e| e.0 < hi);
        &self.by_right[a..b]
    }
}
