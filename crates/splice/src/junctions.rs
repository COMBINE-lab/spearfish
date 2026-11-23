//! Known splice junctions and exons (optional input).
//!
//! Read from transcript models in GTF, GFF3 or BED12 (as `paftools.js gff2bed` writes); the
//! format is told from the content. Positions are on the index's concatenated axis. A known
//! splice site costs nothing where an unknown one pays for its motif and a little more for
//! being unknown; known exons and the donor-to-acceptor links are extra candidates where read
//! bases are left unexplained.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::index::GenomeIndex;

/// One transcript model: the sequence it is on (by index order), whether it is on the reverse
/// strand, and its exons as 0-based half-open intervals on that sequence, in any order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Model {
    pub seq: usize,
    pub reverse: bool,
    pub exons: Vec<(u32, u32)>,
}

/// What a file of transcript models held.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Summary {
    /// `gtf`, `gff3` or `bed12`.
    pub format: &'static str,
    /// Transcript models on sequences of the index.
    pub transcripts: usize,
    /// Distinct introns and exons among them.
    pub introns: usize,
    pub exons: usize,
    /// Records on sequences the index does not hold, and those sequences' names (the first few).
    pub skipped_records: usize,
    pub skipped_sequences: Vec<String>,
}

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

/// The value of `key` among a GTF (`key "value";`) or GFF3 (`key=value;`) record's attributes.
fn attribute<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    for field in attrs.split(';') {
        let field = field.trim();
        if let Some(rest) = field.strip_prefix(key) {
            if let Some(v) = rest.strip_prefix('=') {
                return Some(v.trim());
            }
            if rest.starts_with([' ', '\t']) {
                return Some(rest.trim().trim_matches('"'));
            }
        }
    }
    None
}

impl Junctions {
    /// Read transcript models from a plain-text GTF, GFF3 or BED12 file.
    pub fn from_file(path: &Path, ix: &GenomeIndex) -> Result<Junctions> {
        let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        Ok(Self::from_reader(BufReader::new(f), ix)
            .with_context(|| format!("reading {}", path.display()))?
            .0)
    }

    /// Read transcript models from GTF, GFF3 or BED12 text. GTF and GFF3 `exon` records are
    /// grouped by `transcript_id` or `Parent`; a BED12 line is one transcript. Records on
    /// sequences the index does not hold are counted in the summary and otherwise ignored.
    pub fn from_reader(input: impl BufRead, ix: &GenomeIndex) -> Result<(Junctions, Summary)> {
        let seq_of: HashMap<&str, usize> = ix
            .names()
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), i))
            .collect();
        let mut models: Vec<Model> = Vec::new();
        let mut by_id: HashMap<String, usize> = HashMap::new();
        let mut sum = Summary::default();
        let mut skipped: Vec<String> = Vec::new();
        let mut gff: Option<bool> = None;
        for (n, line) in input.lines().enumerate() {
            let line = line?;
            let n = n + 1;
            if line.starts_with('#') || line.trim().is_empty() {
                gff3_directive(&line, &mut sum);
                continue;
            }
            if line.starts_with("track ") || line.starts_with("browser ") {
                continue;
            }
            let p: Vec<&str> = line.split('\t').collect();
            // GTF and GFF3 have nine columns with a strand in the seventh; BED12 has twelve
            // with numbers in the second and third.
            let is_gff = *gff.get_or_insert_with(|| {
                p.len() == 9 && matches!(p[6], "+" | "-" | "." | "?") && p[3].parse::<u32>().is_ok()
            });
            let num = |s: &str, what: &str| -> Result<u32> {
                s.parse::<u32>()
                    .map_err(|_| anyhow::anyhow!("line {n}: {what} `{s}` is not a number"))
            };
            let mut skip = |name: &str, sum: &mut Summary| {
                sum.skipped_records += 1;
                if skipped.len() < 5 && !skipped.iter().any(|s| s == name) {
                    skipped.push(name.to_string());
                }
            };
            if is_gff {
                if p.len() < 9 {
                    bail!("line {n}: {} columns where GTF and GFF3 have 9", p.len());
                }
                if p[2] != "exon" {
                    continue;
                }
                let Some(&seq) = seq_of.get(p[0]) else {
                    skip(p[0], &mut sum);
                    continue;
                };
                let (id, format) = match attribute(p[8], "transcript_id") {
                    Some(id) => (id, "gtf"),
                    None => match attribute(p[8], "Parent") {
                        Some(id) => (id, "gff3"),
                        None => bail!("line {n}: the exon has neither transcript_id nor Parent"),
                    },
                };
                if sum.format.is_empty() {
                    sum.format = format;
                }
                let (s, e) = (num(p[3], "start")?, num(p[4], "end")?);
                if s == 0 || e < s || e > ix.seq_len(seq) {
                    bail!(
                        "line {n}: exon {s}-{e} is not an interval of {} ({} bases)",
                        p[0],
                        ix.seq_len(seq)
                    );
                }
                // An exon shared by several transcripts lists them all in GFF3.
                for id in id.split(',') {
                    let k = *by_id.entry(format!("{seq}\t{id}")).or_insert_with(|| {
                        models.push(Model {
                            seq,
                            reverse: p[6] == "-",
                            exons: Vec::new(),
                        });
                        models.len() - 1
                    });
                    models[k].exons.push((s - 1, e));
                }
            } else {
                if p.len() < 12 {
                    bail!(
                        "line {n}: {} columns; transcript models are read from GTF, GFF3 or BED12 (12 columns)",
                        p.len()
                    );
                }
                sum.format = "bed12";
                let Some(&seq) = seq_of.get(p[0]) else {
                    skip(p[0], &mut sum);
                    continue;
                };
                let start = num(p[1], "start")?;
                let list = |s: &str, what: &str| -> Result<Vec<u32>> {
                    s.split(',')
                        .filter(|x| !x.is_empty())
                        .map(|x| num(x, what))
                        .collect()
                };
                let (sizes, starts) = (list(p[10], "block size")?, list(p[11], "block start")?);
                if sizes.len() != starts.len() {
                    bail!(
                        "line {n}: {} block sizes and {} block starts",
                        sizes.len(),
                        starts.len()
                    );
                }
                let exons: Vec<(u32, u32)> = sizes
                    .iter()
                    .zip(&starts)
                    .map(|(&z, &b)| (start + b, start + b + z))
                    .collect();
                if exons.iter().any(|&(_, e)| e > ix.seq_len(seq)) {
                    bail!(
                        "line {n}: the transcript ends past the end of {} ({} bases)",
                        p[0],
                        ix.seq_len(seq)
                    );
                }
                models.push(Model {
                    seq,
                    reverse: p[5] == "-",
                    exons,
                });
            }
        }
        sum.skipped_sequences = skipped;
        sum.transcripts = models.len();
        let j = Self::from_models(models, ix);
        sum.introns = j.len();
        sum.exons = j.exons.len();
        Ok((j, sum))
    }

    /// Junctions and exons of transcript models that a caller already holds.
    pub fn from_models(models: impl IntoIterator<Item = Model>, ix: &GenomeIndex) -> Junctions {
        let mut j = Junctions::default();
        for Model {
            seq,
            reverse: minus,
            mut exons,
        } in models
        {
            let off = ix.offsets()[seq];
            exons.sort_unstable();
            for (k, &(s, e)) in exons.iter().enumerate() {
                let (s, e) = (off + s, off + e);
                j.exons.push((s, e));
                if k > 0 && s > off + exons[k - 1].1 {
                    let left = off + exons[k - 1].1;
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
        j
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

    /// The known intron left ends (`right`: right ends) in `[lo, hi]` for a transcript strand.
    pub fn sites_in(&self, lo: u32, hi: u32, minus: bool, right: bool) -> &[u32] {
        let v = if right {
            &self.acceptors[minus as usize]
        } else {
            &self.donors[minus as usize]
        };
        let a = v.partition_point(|&x| x < lo);
        let b = a + v[a..].partition_point(|&x| x <= hi);
        &v[a..b]
    }

    /// Whether `[left, right)` is a known intron.
    pub fn intron(&self, left: u32, right: u32) -> bool {
        self.by_left.binary_search(&(left, right)).is_ok()
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

/// `##gff-version 3` names the format before any record does.
fn gff3_directive(line: &str, sum: &mut Summary) {
    if line.starts_with("##gff-version") && line.split_whitespace().nth(1).is_some_and(|v| v.starts_with('3'))
    {
        sum.format = "gff3";
    }
}
