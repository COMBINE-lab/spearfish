//! Per-level core statistics over a FASTA file: occurrences, distinct labels, mean core length
//! and spacing, and the fraction of labels occurring once (a specificity proxy).
//! Usage: fasta_level_stats <fasta> [max_level]
use std::collections::HashMap;

use lcp_core::{LcpParser, Parser};

fn main() {
    let path = std::env::args().nth(1).expect("usage: fasta_level_stats <fasta> [max_level]");
    let k: u8 = std::env::args().nth(2).map(|s| s.parse().unwrap()).unwrap_or(7);
    let text = std::fs::read_to_string(path).unwrap();
    let seqs: Vec<Vec<u8>> = text
        .split('>')
        .skip(1)
        .map(|rec| rec.lines().skip(1).flat_map(|l| l.trim().bytes()).collect())
        .collect();
    let parser = LcpParser::default();
    let mut counts: Vec<HashMap<u64, u32>> = vec![HashMap::new(); k as usize];
    let mut len_sum = vec![0u64; k as usize];
    let mut span_sum = vec![0u64; k as usize];
    let mut gaps = vec![0u64; k as usize];
    let total_bp: u64 = seqs.iter().map(|s| s.len() as u64).sum();
    let t = std::time::Instant::now();
    for s in &seqs {
        let h = parser.parse(s, k);
        for (i, l) in h.levels.iter().enumerate() {
            for c in &l.cores {
                *counts[i].entry(c.label).or_default() += 1;
                len_sum[i] += c.len() as u64;
            }
            if l.cores.len() > 1 {
                span_sum[i] += (l.cores.last().unwrap().start - l.cores[0].start) as u64;
                gaps[i] += l.cores.len() as u64 - 1;
            }
        }
    }
    eprintln!("{} sequences, {:.1} Mbp, parsed in {:.2}s", seqs.len(), total_bp as f64 / 1e6, t.elapsed().as_secs_f64());
    println!("level\toccurrences\tdistinct\tsingleton_frac\tmean_len\tmean_spacing");
    for i in 0..k as usize {
        let occ: u64 = counts[i].values().map(|&c| c as u64).sum();
        let single = counts[i].values().filter(|&&c| c == 1).count();
        println!(
            "{}\t{}\t{}\t{:.3}\t{:.1}\t{:.1}",
            i + 1,
            occ,
            counts[i].len(),
            single as f64 / counts[i].len().max(1) as f64,
            len_sum[i] as f64 / occ.max(1) as f64,
            span_sum[i] as f64 / gaps[i].max(1) as f64
        );
    }
}
