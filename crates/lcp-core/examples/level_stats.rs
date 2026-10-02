//! Prints per-level core counts and mean lengths for a random sequence, from lcp-core and
//! (if built with the oracle) lcptools. Quick sanity tool.
use lcp_core::{LcpParser, Parser};

fn main() {
    let n: usize = std::env::args()
        .nth(1)
        .map(|s| s.parse().unwrap())
        .unwrap_or(1_000_000);
    let mut x: u64 = 0x1234_5678;
    let seq: Vec<u8> = (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[(x & 3) as usize]
        })
        .collect();
    let t = std::time::Instant::now();
    let h = LcpParser::default().parse(&seq, 8);
    let dt = t.elapsed().as_secs_f64();
    let oracle = lcp_oracle_sys::lcptools_boundaries(&seq, 8);
    println!("level\tcores\tmean_len\tmean_spacing\tlcptools_cores");
    for (k, l) in h.levels.iter().enumerate() {
        let c = &l.cores;
        let mean_len = c.iter().map(|c| c.len() as f64).sum::<f64>() / c.len().max(1) as f64;
        let spacing = if c.len() > 1 {
            (c.last().unwrap().start - c[0].start) as f64 / (c.len() - 1) as f64
        } else {
            0.0
        };
        println!(
            "{}\t{}\t{:.1}\t{:.1}\t{}",
            k + 1,
            c.len(),
            mean_len,
            spacing,
            oracle[k].len()
        );
    }
    eprintln!("parse to level 8: {:.1} Mbp/s", n as f64 / dt / 1e6);
}
