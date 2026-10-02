//! Single-thread seeding throughput per scheme on random 2 kb sequences.
//! Usage: throughput [spec ...]
use seeds::from_spec;

fn main() {
    let mut specs: Vec<String> = std::env::args().skip(1).collect();
    if specs.is_empty() {
        specs = [
            "lcp:K=3",
            "lcp:K=4",
            "minimizer:k=19,w=23",
            "minimizer:k=15,w=10",
            "closed-syncmer:k=19,s=7,keep=0.54",
            "linked-lcp:b=3,wmin=1,wmax=2,order=2,keep=0.43",
            "randstrobe:k=15,s=5,keep=0.193,wmin=1,wmax=3",
        ]
        .map(String::from)
        .to_vec();
    }
    let mut x: u64 = 0x1234_5678;
    let seq: Vec<u8> = (0..20_000_000)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[(x & 3) as usize]
        })
        .collect();
    for spec in specs {
        let s = from_spec(&spec).unwrap();
        let t = std::time::Instant::now();
        let mut n = 0usize;
        for chunk in seq.chunks(2000) {
            n += s.seeds(chunk).levels.iter().map(|l| l.cores.len()).sum::<usize>();
        }
        let dt = t.elapsed().as_secs_f64();
        println!(
            "{:52} {:7.1} Mbp/s  ({n} seeds)",
            spec,
            seq.len() as f64 / dt / 1e6
        );
    }
}
