//! Per-level core counts and mean lengths for a random sequence, compared with lcptools,
//! plus single-thread parse throughput. Usage: level_stats [len] [max_level]

extern crate spearfish_lcp as lcp_core;

use lcp_core::{LcpParser, Parser};

fn main() {
    let n: usize = std::env::args()
        .nth(1)
        .map(|s| s.parse().unwrap())
        .unwrap_or(10_000_000);
    let k: u8 = std::env::args().nth(2).map(|s| s.parse().unwrap()).unwrap_or(8);
    let mut x: u64 = 0x1234_5678;
    let seq: Vec<u8> = (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[(x & 3) as usize]
        })
        .collect();
    let parser = LcpParser::default();
    // Parse in 2 kb pieces, like transcripts and reads.
    let t = std::time::Instant::now();
    let mut cores = 0usize;
    for chunk in seq.chunks(2000) {
        cores += parser
            .parse(chunk, k)
            .levels
            .iter()
            .map(|l| l.cores.len())
            .sum::<usize>();
    }
    let dt = t.elapsed().as_secs_f64();
    eprintln!(
        "parse to level {k}: {:.1} Mbp/s single-threaded ({cores} cores)",
        n as f64 / dt / 1e6
    );

    let h = parser.parse(&seq[..1_000_000.min(n)], k);
    let oracle = lcp_oracle_sys::lcptools_boundaries(&seq[..1_000_000.min(n)], k);
    println!("level\tcores\tmean_len\tlcptools_cores");
    for (i, l) in h.levels.iter().enumerate() {
        let c = &l.cores;
        let mean_len = c.iter().map(|c| c.len() as f64).sum::<f64>() / c.len().max(1) as f64;
        println!("{}\t{}\t{:.1}\t{}", i + 1, c.len(), mean_len, oracle[i].len());
    }
}
