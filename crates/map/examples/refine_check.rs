//! The refinement chainer on real data: for each read (both orientations), every candidate's
//! per-transcript chain scores from the refinement chainer (unbounded beam unless given) must
//! equal the oracle's exactly. Also times both chainers and counts identical hits (chain paths
//! can differ only on score ties).
//! Usage: refine_check <index> <reads.fa> [max_reads] [beam]

extern crate spearfish_map as map;

use std::collections::BTreeMap;
use std::sync::Arc;
use std::time::Instant;

use index::Index;
use lcp_core::encode::revcomp;
use map::{Chainer, MapParams, Mapper, OracleChainer, RefineChainer, anchor};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let idx = Arc::new(Index::open(a[1].as_ref()).unwrap());
    let max: usize = a.get(3).map_or(10_000, |s| s.parse().unwrap());
    let beam: Option<usize> = a.get(4).map(|s| s.parse().unwrap());
    let p = MapParams::default();
    let mapper = Mapper::new(idx.clone(), OracleChainer { p: p.chain.clone() }, p.clone()).unwrap();
    let oracle = OracleChainer { p: p.chain.clone() };
    let refine = RefineChainer {
        p: p.chain.clone(),
        beam,
    };
    let text = std::fs::read_to_string(&a[2]).unwrap();
    let (mut to, mut tr) = (0f64, 0f64);
    let (mut cands, mut txs, mut bad_score, mut missing, mut same_hit, mut best_wrong) = (0, 0, 0, 0, 0, 0);
    let mut n = 0;
    let (mut sig_tx, mut sig_uniq, mut place_all, mut place_uniq) = (0usize, 0usize, 0usize, 0usize);
    let mut place_trie = 0usize;
    // Per family-size bucket (transcripts): (candidates, oracle s, refine s).
    let mut buckets = [(0usize, 0f64, 0f64); 8];
    let bucket = |t: usize| match t {
        1 => 0,
        2..=3 => 1,
        4..=7 => 2,
        8..=15 => 3,
        16..=31 => 4,
        32..=63 => 5,
        64..=127 => 6,
        _ => 7,
    };
    for rec in text.split('>').skip(1).take(max) {
        let seq: Vec<u8> = rec.lines().skip(1).flat_map(|l| l.trim().bytes()).collect();
        n += 1;
        for s in [seq.clone(), revcomp(&seq)] {
            let seeds = mapper.seeder.seed(&s);
            let an = anchor::anchors(&idx, &seeds, &p.anchor);
            for mut c in anchor::candidates(&idx, an, &p.anchor) {
                // The candidates the mapper chains: collinear anchors merged, if on.
                if p.merge_collinear {
                    anchor::merge_collinear(&idx, &mut c);
                }
                cands += 1;
                {
                    // Placement lists per transcript, shifted to start at 0: transcripts with
                    // equal lists have identical DPs.
                    let mut by_tx: std::collections::HashMap<u32, Vec<(u32, i64)>> = Default::default();
                    for (ai, an) in c.anchors.iter().enumerate() {
                        let u = idx.node_unitig()[an.node as usize];
                        let off = idx.node_offset()[an.node as usize] as i64;
                        let (ts, ps) = idx.occurrences(u);
                        for (&t, &p) in ts.iter().zip(ps) {
                            by_tx.entry(t).or_default().push((ai as u32, p as i64 + off));
                        }
                    }
                    let mut sigs: std::collections::HashSet<Vec<(u32, i64)>> = Default::default();
                    let mut work_all = 0usize;
                    let mut work_uniq = 0usize;
                    for v in by_tx.values() {
                        let t0 = v[0].1;
                        let sig: Vec<(u32, i64)> = v.iter().map(|&(a, t)| (a, t - t0)).collect();
                        work_all += v.len();
                        if sigs.insert(sig) {
                            work_uniq += v.len();
                        }
                    }
                    // Trie of shifted placement lists: placements in distinct prefixes.
                    let mut lists: Vec<Vec<(u32, i64)>> = by_tx
                        .values()
                        .map(|v| v.iter().map(|&(a, t)| (a, t - v[0].1)).collect())
                        .collect();
                    lists.sort();
                    for (k, l) in lists.iter().enumerate() {
                        let lcp = if k == 0 {
                            0
                        } else {
                            l.iter().zip(&lists[k - 1]).take_while(|(a, b)| a == b).count()
                        };
                        place_trie += l.len() - lcp;
                    }
                    sig_tx += by_tx.len();
                    sig_uniq += sigs.len();
                    place_all += work_all;
                    place_uniq += work_uniq;
                }
                let t0 = Instant::now();
                let o = oracle.chain(&idx, &c);
                let t1 = Instant::now();
                let r = refine.chain(&idx, &c);
                let (dto, dtr) = ((t1 - t0).as_secs_f64(), t1.elapsed().as_secs_f64());
                to += dto;
                tr += dtr;
                let b = &mut buckets[bucket(o.len())];
                b.0 += 1;
                b.1 += dto;
                b.2 += dtr;
                let rm: BTreeMap<u32, &map::TxHit> = r.iter().map(|h| (h.tx, h)).collect();
                let ob = o.iter().map(|h| h.score).fold(f32::NEG_INFINITY, f32::max);
                let rb = r.iter().map(|h| h.score).fold(f32::NEG_INFINITY, f32::max);
                best_wrong += (ob != rb) as usize;
                for h in &o {
                    txs += 1;
                    match rm.get(&h.tx) {
                        None => missing += 1,
                        Some(x) => {
                            if x.score.to_bits() != h.score.to_bits() {
                                bad_score += 1;
                            } else if *x == h {
                                same_hit += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    println!(
        "reads {n} candidates {cands} transcripts {txs} | score mismatches {bad_score} missing {missing} \
         best-score mismatches {best_wrong} | identical hits {:.2}% | oracle {to:.2}s refine {tr:.2}s ({:.2}x)",
        100.0 * same_hit as f64 / txs.max(1) as f64,
        to / tr
    );
    println!(
        "  transcripts {sig_tx} distinct placement lists {sig_uniq} ({:.1}%) | placements {place_all} in distinct lists {place_uniq} ({:.1}%)",
        100.0 * sig_uniq as f64 / sig_tx as f64,
        100.0 * place_uniq as f64 / place_all as f64
    );
    println!(
        "  placements in a trie of shifted lists {place_trie} ({:.1}%)",
        100.0 * place_trie as f64 / place_all as f64
    );
    map::refine::STATS.with(|x| {
        let x = x.borrow();
        if x.0 > 0 {
            println!(
                "  anchors {} candidates/anchor {:.1} states/anchor {:.2} max states {} max transcripts {}",
                x.0,
                x.1 as f64 / x.0 as f64,
                x.2 as f64 / x.0 as f64,
                x.3,
                x.4
            );
        }
    });
    for (i, (c, o, r)) in buckets.iter().enumerate() {
        println!(
            "  family size bucket {i} (<{}): candidates {c} oracle {o:.2}s refine {r:.2}s",
            [2, 4, 8, 16, 32, 64, 128, 99999][i]
        );
    }
}
