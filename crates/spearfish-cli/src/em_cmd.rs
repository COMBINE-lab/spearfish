//! `spearfish em`: spearfish's own EM quantification from a mapping TSV (`spearfish map` output, or
//! another mapper's alignments converted to the same columns). Output: transcript, length,
//! estimated reads.
//!
//! Columns are found by the header line's names (`transcript`, `score`, `flags`, and `read`
//! and `qscore` when present), so the 15-column TSV of `spearfish map` with quantification
//! scores and the 14-column TSV without them both read correctly. A file without a header, or
//! whose header does not name those columns, is read by the 14-column positions of
//! [`spearfish::output::TSV_HEADER`]. Hits are weighed by `qscore` when the file has it
//! (`--score-column`).

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::time::Instant;

use crate::em::{EqBuilder, QuantParams, em};
use anyhow::{Context, Result, bail};
use clap::ValueEnum;
use paraseq::Record;

/// Which column weighs a read's hits.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum ScoreColumn {
    /// `qscore` when the header names it, else `score`.
    #[default]
    Auto,
    /// The mapping score.
    Score,
    /// The quantification score (the file must have a `qscore` column).
    Qscore,
}

/// Where the columns `em` reads are on a line (0-based).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Columns {
    read: usize,
    transcript: usize,
    score: usize,
    flags: usize,
    qscore: Option<usize>,
}

/// The 14-column layout of `spearfish map` without quantification scores (and of 0.2.x).
const POSITIONAL: Columns = Columns {
    read: 0,
    transcript: 4,
    score: 5,
    flags: 13,
    qscore: None,
};

impl Columns {
    /// The columns a header line names, or `None` if it does not name `transcript`, `score`
    /// and `flags`. The read name is the `read` column, or the first.
    fn from_header(line: &str) -> Option<Columns> {
        let at = |name: &str| line.split('\t').position(|f| f.trim_end() == name);
        Some(Columns {
            read: at("read").unwrap_or(0),
            transcript: at("transcript")?,
            score: at("score")?,
            flags: at("flags")?,
            qscore: at("qscore"),
        })
    }

    /// The column hits are weighed by, or an error if `col` asks for a missing `qscore`.
    fn weight(&self, col: ScoreColumn) -> Result<usize> {
        Ok(match (col, self.qscore) {
            (ScoreColumn::Score, _) | (ScoreColumn::Auto, None) => self.score,
            (_, Some(q)) => q,
            (ScoreColumn::Qscore, None) => {
                bail!("--score-column qscore: the mapping TSV's header names no `qscore` column")
            }
        })
    }
}

/// What a mapping TSV's first line is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FirstLine {
    /// A header naming the columns.
    Header,
    /// A header that does not name them; the 14-column layout is used.
    UnknownHeader,
    /// A read's line, from a file without a header; the 14-column layout is used.
    Data,
}

/// The columns, and what the first line is. A line that names the columns is a header; one
/// with a number in the 14-column layout's `score` column is a read's line; anything else is a
/// header `em` does not know.
fn layout(first: &str) -> (Columns, FirstLine) {
    if let Some(c) = Columns::from_header(first) {
        return (c, FirstLine::Header);
    }
    let data = first
        .split('\t')
        .nth(POSITIONAL.score)
        .is_some_and(|s| s.parse::<f64>().is_ok());
    (
        POSITIONAL,
        if data {
            FirstLine::Data
        } else {
            FirstLine::UnknownHeader
        },
    )
}

/// One line's read, transcript, weight and flags (`*`, 0 and empty when the line has no such
/// column).
fn fields(line: &str, cols: Columns, weight: usize) -> (&str, &str, f64, &str) {
    let (mut read, mut tx, mut w, mut flags) = ("", "*", "0", "");
    let last = cols.read.max(cols.transcript).max(weight).max(cols.flags);
    for (i, f) in line.split('\t').enumerate().take(last + 1) {
        if i == cols.read {
            read = f;
        }
        if i == cols.transcript {
            tx = f;
        }
        if i == weight {
            w = f;
        }
        if i == cols.flags {
            flags = f;
        }
    }
    (read, tx, w.parse().unwrap_or(0.0), flags)
}

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Transcript FASTA (names and lengths; IDs are cut at the first '|' or whitespace, as
    /// `spearfish index` cuts them).
    transcripts: PathBuf,
    /// Mapping TSV. Columns are found by the header's names (`transcript`, `score`, `flags`,
    /// and `qscore` if present); without a header, or with one that does not name them, by the
    /// 14-column positions of `spearfish map`'s TSV without `qscore`.
    mapping: PathBuf,
    /// Column to weigh hits by: `auto` takes `qscore` when the TSV has it (as quantifiers
    /// should) and `score` otherwise.
    #[arg(long, value_enum, default_value_t = ScoreColumn::Auto)]
    score_column: ScoreColumn,
    /// Score denominator D in exp((s - s_max) / D). This EM's default (2) is what the
    /// mapper's ablations were evaluated with, and also the scale spearfish suggests to
    /// quantifiers such as oarfish (`spearfish::SCORE_DENOM`).
    #[arg(long, default_value_t = QuantParams::default().score_denom)]
    score_denom: f64,
    /// Keep hits scoring at least this fraction of the read's best (in the `--score-column`).
    #[arg(long, default_value_t = QuantParams::default().min_frac)]
    min_frac: f64,
    /// Divide likelihoods by transcript length.
    #[arg(long)]
    length_norm: bool,
    /// Also use reads flagged `low_coverage` (best chain spans < half the read). Hits flagged
    /// `low_score` are never used.
    #[arg(long)]
    keep_low_coverage: bool,
    /// EM: maximum iterations.
    #[arg(long, default_value_t = QuantParams::default().max_iter)]
    max_iter: usize,
    /// EM: stop when the largest relative change (transcripts with >= 0.01 reads) is below
    /// this.
    #[arg(long, default_value_t = QuantParams::default().tol)]
    tol: f64,
    /// Output TSV.
    #[arg(short, long)]
    output: PathBuf,
    /// EM worker threads (0 = all cores).
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
}

pub fn run(args: Args) -> Result<()> {
    let t0 = Instant::now();
    let mut names: Vec<String> = Vec::new();
    let mut lens: Vec<u32> = Vec::new();
    let mut reader = spearfish::open_fastx(Some(&args.transcripts), None)
        .with_context(|| format!("reading {}", args.transcripts.display()))?;
    let mut rset = reader.new_record_set();
    while rset.fill(&mut reader)? {
        for rec in rset.iter() {
            let rec = rec?;
            names.push(spearfish::name_and_gene(&String::from_utf8_lossy(rec.id())).0);
            lens.push(rec.seq().len() as u32);
        }
    }
    let id: HashMap<&str, u32> = names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i as u32))
        .collect();
    let p = QuantParams {
        score_denom: args.score_denom,
        min_frac: args.min_frac,
        length_norm: args.length_norm,
        max_iter: args.max_iter,
        tol: args.tol,
    };
    let mut b = EqBuilder::default();
    let f =
        std::fs::File::open(&args.mapping).with_context(|| format!("opening {}", args.mapping.display()))?;
    let (mut cur, mut hits): (String, Vec<(u32, f64)>) = (String::new(), Vec::new());
    let (mut unknown, mut unmapped, mut skipped) = (0u64, 0u64, 0u64);
    let mut cur_low = false;
    let mut lines = BufReader::new(f).lines().peekable();
    let (cols, first) = match lines.peek() {
        Some(Ok(line)) => layout(line),
        _ => (POSITIONAL, FirstLine::Data),
    };
    if first == FirstLine::UnknownHeader {
        log::warn!(
            "{}: the header names no `transcript`, `score` and `flags` columns; reading spearfish's 14-column layout",
            args.mapping.display()
        );
    }
    if first != FirstLine::Data {
        lines.next().transpose()?;
    }
    let weight = cols.weight(args.score_column)?;
    let weight_name = if Some(weight) == cols.qscore {
        "qscore"
    } else {
        "score"
    };
    for line in lines {
        let line = line?;
        let (read, tx, score, flags) = fields(&line, cols, weight);
        let has = |flag: &str| flags.split(',').any(|f| f == flag);
        // Hits kept under the minimum score (`map --keep-low-score`) are never used.
        let (low, low_score) = (has("low_coverage"), has("low_score"));
        if read != cur {
            if !hits.is_empty() {
                b.add(&hits, &lens, &p);
            } else if cur_low {
                skipped += 1;
            }
            hits.clear();
            cur = read.to_string();
            cur_low = false;
        }
        if tx == "*" {
            unmapped += 1;
            continue;
        }
        if low_score || (low && !args.keep_low_coverage) {
            cur_low = true;
            continue;
        }
        match id.get(tx) {
            Some(&t) => hits.push((t, score)),
            None => unknown += 1,
        }
    }
    if !hits.is_empty() {
        b.add(&hits, &lens, &p);
    } else if cur_low {
        skipped += 1;
    }
    let reads = b.reads;
    let classes = b.finish();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;
    // Created before the EM runs, so a bad path fails early.
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(&args.output).with_context(|| format!("creating {}", args.output.display()))?,
    );
    let (counts, iters) = pool.install(|| em(names.len(), &classes, &p));
    writeln!(out, "transcript\tlength\treads")?;
    for ((n, l), c) in names.iter().zip(&lens).zip(&counts) {
        writeln!(out, "{n}\t{l}\t{c:.4}")?;
    }
    out.flush()?;
    log::info!(
        "quant: {reads} reads used ({unmapped} unmapped, {skipped} low-coverage or low-score skipped, {unknown} hits to unknown transcripts), weights from `{weight_name}`, {} classes, {iters} EM iterations, {:.2} s",
        classes.len(),
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_names_the_columns() {
        // `spearfish map`'s own headers, without and with `qscore`.
        let (cols, first) = layout(spearfish::output::tsv_header(false).trim_end());
        assert_eq!((cols, first), (POSITIONAL, FirstLine::Header));
        let (cols, first) = layout(spearfish::output::tsv_header(true).trim_end());
        assert_eq!(first, FirstLine::Header);
        assert_eq!(
            cols,
            Columns {
                qscore: Some(14),
                ..POSITIONAL
            }
        );
        assert_eq!(cols.weight(ScoreColumn::Auto).unwrap(), 14);
        assert_eq!(cols.weight(ScoreColumn::Qscore).unwrap(), 14);
        assert_eq!(cols.weight(ScoreColumn::Score).unwrap(), 5);
        assert_eq!(POSITIONAL.weight(ScoreColumn::Auto).unwrap(), 5);
        assert!(POSITIONAL.weight(ScoreColumn::Qscore).is_err());
    }

    #[test]
    fn first_line_without_a_known_header() {
        let line = "r1\t1000\t+\t0\tT1\t100.00\t0\t1000\t0\t1000\t10\t500\t0\t-";
        assert_eq!(layout(line), (POSITIONAL, FirstLine::Data));
        assert_eq!(
            layout("read\tlen\ttarget\tas"),
            (POSITIONAL, FirstLine::UnknownHeader)
        );
        assert_eq!(fields(line, POSITIONAL, 5), ("r1", "T1", 100.0, "-"));
        // A short line: no flags.
        assert_eq!(fields("r1\t1\t+\t0\tT1\t7", POSITIONAL, 5), ("r1", "T1", 7.0, ""));
    }
}
