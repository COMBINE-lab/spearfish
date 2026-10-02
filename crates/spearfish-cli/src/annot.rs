//! Transcript models from a GTF, and per-gene exon segments.
//!
//! A gene's *segments* are the maximal genomic intervals between consecutive exon boundaries
//! of any of its transcripts, restricted to exonic sequence. Every transcript is then a walk
//! over segments, and two transcripts share a segment exactly when both include that
//! sequence. Segment boundaries where the transcripts containing a segment disagree on the
//! neighbouring segment (alternative splice sites, cassette exons, transcript ends) are
//! *alternative boundaries*: that's where junction halos form.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};

/// One transcript: exons as 0-based half-open genomic intervals in transcript order.
#[derive(Clone, Debug)]
pub struct TxModel {
    pub gene: String,
    pub plus: bool,
    pub exons: Vec<(u64, u64)>,
    /// Transcript offset of each exon's first base.
    offsets: Vec<u64>,
}

impl TxModel {
    pub fn new(gene: String, plus: bool, mut exons: Vec<(u64, u64)>) -> Self {
        exons.sort_unstable();
        if !plus {
            exons.reverse();
        }
        let mut offsets = Vec::with_capacity(exons.len());
        let mut acc = 0;
        for &(s, e) in &exons {
            offsets.push(acc);
            acc += e - s;
        }
        TxModel {
            gene,
            plus,
            exons,
            offsets,
        }
    }

    pub fn len(&self) -> u64 {
        self.exons.iter().map(|(s, e)| e - s).sum()
    }

    /// Genomic extent `[lo, hi)`.
    pub fn span(&self) -> (u64, u64) {
        let lo = self.exons.iter().map(|e| e.0).min().unwrap();
        let hi = self.exons.iter().map(|e| e.1).max().unwrap();
        (lo, hi)
    }

    /// Genomic coordinate of transcript position `p < len`.
    pub fn to_genome(&self, p: u64) -> u64 {
        let i = self.offsets.partition_point(|&o| o <= p) - 1;
        let (s, e) = self.exons[i];
        let d = p - self.offsets[i];
        if self.plus { s + d } else { e - 1 - d }
    }

    /// Genomic extent `[lo, hi)` of the transcript interval `[a, b)` (b > a).
    pub fn genome_span(&self, a: u64, b: u64) -> (u64, u64) {
        let (x, y) = (self.to_genome(a), self.to_genome(b - 1));
        (x.min(y), x.max(y) + 1)
    }
}

fn attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    attrs.split(';').find_map(|kv| {
        let kv = kv.trim();
        let rest = kv.strip_prefix(key)?.strip_prefix(' ')?;
        Some(rest.trim_matches('"'))
    })
}

/// Exon records gathered per transcript: gene, strand (+), exons.
type RawTx = (String, bool, Vec<(u64, u64)>);

/// Transcript models by transcript ID, from the GTF's exon records.
pub fn read_gtf(path: &Path) -> Result<HashMap<String, TxModel>> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut exons: HashMap<String, RawTx> = HashMap::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.splitn(9, '\t').collect();
        if f.len() < 9 || f[2] != "exon" {
            continue;
        }
        let (start, end): (u64, u64) = (f[3].parse()?, f[4].parse()?);
        let plus = match f[6] {
            "+" => true,
            "-" => false,
            s => bail!("unsupported strand {s}"),
        };
        let tx = attr(f[8], "transcript_id").context("exon without transcript_id")?;
        let gene = attr(f[8], "gene_id").context("exon without gene_id")?;
        let e = exons
            .entry(tx.to_string())
            .or_insert_with(|| (gene.to_string(), plus, Vec::new()));
        e.2.push((start - 1, end));
    }
    Ok(exons
        .into_iter()
        .map(|(tx, (gene, plus, ex))| (tx, TxModel::new(gene, plus, ex)))
        .collect())
}

/// A transcript's pass through one segment: transcript interval `[tx_start, tx_end)`.
#[derive(Clone, Copy, Debug)]
pub struct SegOcc {
    pub seg: u32,
    pub tx_start: u64,
    pub tx_end: u64,
}

/// Segment model of one gene.
pub struct Segments {
    /// Per transcript (in the order given), its segments in transcript order.
    pub walks: Vec<Vec<SegOcc>>,
    /// Per segment: the transcripts containing it.
    pub members: Vec<Vec<u32>>,
    /// Per segment: whether its 5' (transcript-order left) / 3' boundary is alternative.
    pub alt_left: Vec<bool>,
    pub alt_right: Vec<bool>,
    /// Distinct exonic bases of the gene.
    pub exonic_bases: u64,
    /// Distinct splice junctions (donor, acceptor) used by any transcript.
    pub junctions: usize,
}

impl Segments {
    pub fn new(txs: &[&TxModel]) -> Self {
        let mut bps: Vec<u64> = txs
            .iter()
            .flat_map(|t| t.exons.iter().flat_map(|&(s, e)| [s, e]))
            .collect();
        bps.sort_unstable();
        bps.dedup();
        // Segment i is [bps[i], bps[i+1]); only segments some exon covers get used.
        let mut walks = Vec::with_capacity(txs.len());
        let mut used = vec![false; bps.len()];
        for t in txs {
            let mut walk = Vec::new();
            let mut off = 0;
            for &(s, e) in &t.exons {
                let (a, b) = (bps.partition_point(|&x| x < s), bps.partition_point(|&x| x < e));
                let segs: Vec<usize> = if t.plus {
                    (a..b).collect()
                } else {
                    (a..b).rev().collect()
                };
                for i in segs {
                    let len = bps[i + 1] - bps[i];
                    used[i] = true;
                    walk.push(SegOcc {
                        seg: i as u32,
                        tx_start: off,
                        tx_end: off + len,
                    });
                    off += len;
                }
            }
            walks.push(walk);
        }
        let n = bps.len();
        let mut members = vec![Vec::new(); n];
        // Neighbour sets in transcript order; u32::MAX marks a transcript end.
        let mut left: Vec<Vec<u32>> = vec![Vec::new(); n];
        let mut right: Vec<Vec<u32>> = vec![Vec::new(); n];
        for (ti, w) in walks.iter().enumerate() {
            for (j, o) in w.iter().enumerate() {
                let s = o.seg as usize;
                members[s].push(ti as u32);
                left[s].push(if j == 0 { u32::MAX } else { w[j - 1].seg });
                right[s].push(w.get(j + 1).map_or(u32::MAX, |x| x.seg));
            }
        }
        let differs = |v: &Vec<u32>| v.iter().any(|&x| x != v[0]);
        let exonic_bases = (0..n.saturating_sub(1))
            .filter(|&i| used[i])
            .map(|i| bps[i + 1] - bps[i])
            .sum();
        let mut junctions: Vec<(u64, u64)> = txs
            .iter()
            .flat_map(|t| {
                t.exons.windows(2).map(|w| {
                    if t.plus {
                        (w[0].1, w[1].0)
                    } else {
                        (w[1].1, w[0].0)
                    }
                })
            })
            .collect();
        junctions.sort_unstable();
        junctions.dedup();
        Segments {
            walks,
            alt_left: left.iter().map(|v| !v.is_empty() && differs(v)).collect(),
            alt_right: right.iter().map(|v| !v.is_empty() && differs(v)).collect(),
            members,
            exonic_bases,
            junctions: junctions.len(),
        }
    }

    /// The segment occurrence of transcript `t` containing transcript interval `[a, b)`, if
    /// the interval lies within one segment.
    pub fn containing(&self, t: usize, a: u64, b: u64) -> Option<SegOcc> {
        let w = &self.walks[t];
        let i = w.partition_point(|o| o.tx_end <= a);
        w.get(i).copied().filter(|o| o.tx_start <= a && b <= o.tx_end)
    }
}

/// Exonic bases of `t` inside `[lo, hi)`, as sorted disjoint intervals.
pub fn exonic_in(t: &TxModel, lo: u64, hi: u64) -> Vec<(u64, u64)> {
    let mut v: Vec<(u64, u64)> = t
        .exons
        .iter()
        .map(|&(s, e)| (s.max(lo), e.min(hi)))
        .filter(|(s, e)| s < e)
        .collect();
    v.sort_unstable();
    v
}

/// Size of the symmetric difference of two sorted disjoint interval sets.
pub fn sym_diff(a: &[(u64, u64)], b: &[(u64, u64)]) -> u64 {
    let total = |v: &[(u64, u64)]| v.iter().map(|(s, e)| e - s).sum::<u64>();
    let (mut i, mut j, mut inter) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        let (s, e) = (a[i].0.max(b[j].0), a[i].1.min(b[j].1));
        if s < e {
            inter += e - s;
        }
        if a[i].1 < b[j].1 { i += 1 } else { j += 1 }
    }
    total(a) + total(b) - 2 * inter
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tx(plus: bool, exons: &[(u64, u64)]) -> TxModel {
        TxModel::new("g".into(), plus, exons.to_vec())
    }

    #[test]
    fn genome_mapping() {
        let p = tx(true, &[(100, 110), (200, 205)]);
        assert_eq!(
            (p.to_genome(0), p.to_genome(9), p.to_genome(10), p.to_genome(14)),
            (100, 109, 200, 204)
        );
        assert_eq!(p.genome_span(8, 12), (108, 202));
        let m = tx(false, &[(100, 110), (200, 205)]);
        // Minus strand: transcript starts at the highest coordinate.
        assert_eq!(
            (m.to_genome(0), m.to_genome(4), m.to_genome(5), m.to_genome(14)),
            (204, 200, 109, 100)
        );
        assert_eq!(m.len(), 15);
    }

    #[test]
    fn segments_of_exon_skip_and_alt_site() {
        // t0: E1 E2 E3; t1: E1 E3 (skip); t2: E1 E2' E3 with E2' starting 5 nt later.
        let t0 = tx(true, &[(0, 10), (20, 30), (40, 50)]);
        let t1 = tx(true, &[(0, 10), (40, 50)]);
        let t2 = tx(true, &[(0, 10), (25, 30), (40, 50)]);
        let s = Segments::new(&[&t0, &t1, &t2]);
        // Segments: [0,10) [20,25) [25,30) [40,50).
        assert_eq!(s.walks[0].len(), 4);
        assert_eq!(s.walks[1].len(), 2);
        assert_eq!(s.walks[2].len(), 3);
        let seg = |t: usize, j: usize| s.walks[t][j].seg as usize;
        // E1: shared by all; its right boundary is alternative (next is [20,25), [40,50) or [25,30)).
        assert_eq!(s.members[seg(0, 0)].len(), 3);
        assert!(s.alt_right[seg(0, 0)] && !s.alt_left[seg(0, 0)]);
        // [25,30): in t0 (after [20,25)) and t2 (after E1): alternative left, shared right.
        assert!(s.alt_left[seg(0, 2)] && !s.alt_right[seg(0, 2)]);
        assert_eq!(s.exonic_bases, 30);
        assert_eq!(s.junctions, 4); // (10,20) (30,40) (10,40) (10,25)
        assert!(s.containing(0, 12, 15).is_some());
        assert!(
            s.containing(0, 12, 16).is_none(),
            "crosses the alt-site boundary at 25"
        );
    }

    #[test]
    fn interval_sym_diff() {
        assert_eq!(sym_diff(&[(0, 10), (20, 30)], &[(0, 10), (25, 30)]), 5);
        assert_eq!(sym_diff(&[(0, 10)], &[]), 10);
    }
}
