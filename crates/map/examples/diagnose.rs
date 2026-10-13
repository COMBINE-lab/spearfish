//! Classify reads whose true transcript is not reported. Reads are named
//! `r<i>|<transcript>|<strand>|...` (as `spearfish simulate` writes them).
//! Usage: diagnose <index> <reads.fa> [max_reads]
use std::collections::HashMap;
use std::sync::Arc;

use index::Index;
use lcp_core::encode::revcomp;
use map::{AnchorParams, MapParams, Mapper, OracleChainer, anchor};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let idx = Arc::new(Index::open(a[1].as_ref()).unwrap());
    let max: usize = a.get(3).map_or(usize::MAX, |s| s.parse().unwrap());
    let p = MapParams::default();
    let mapper = Mapper::new(idx.clone(), OracleChainer { p: p.chain.clone() }, p.clone()).unwrap();
    let variant = |q: MapParams| Mapper::new(idx.clone(), OracleChainer { p: q.chain.clone() }, q).unwrap();
    let no_tb = variant(MapParams {
        tiebreak: None,
        ..p.clone()
    });
    let wide = variant(MapParams {
        delta_rep: 1.0,
        ..p.clone()
    });
    let wide_no_tb = variant(MapParams {
        delta_rep: 1.0,
        tiebreak: None,
        ..p.clone()
    });
    let tx_id: HashMap<&str, u32> = (0..idx.n_transcripts() as u32)
        .map(|t| (idx.tx_name(t), t))
        .collect();
    let text = std::fs::read_to_string(&a[2]).unwrap();
    let mut counts: HashMap<&str, usize> = HashMap::new();
    let mut n = 0;
    let mut gaps: Vec<f32> = Vec::new();
    for rec in text.split('>').skip(1).take(max) {
        let mut lines = rec.lines();
        let name = lines.next().unwrap();
        let seq: Vec<u8> = lines.flat_map(|l| l.trim().bytes()).collect();
        let f: Vec<&str> = name.split('|').collect();
        let truth = tx_id[f[1]];
        let want_rev = f[2] == "-";
        n += 1;
        let r = mapper.map(&seq);
        if r.hits.iter().any(|h| h.tx == truth && h.reverse == want_rev) {
            continue;
        }
        // Re-run the steps on the true strand.
        let strand = if want_rev { revcomp(&seq) } else { seq.clone() };
        let seeds = mapper.seeder.seed(&strand);
        let truth_gene = idx.tx_gene()[truth as usize];
        let in_truth = seeds
            .iter()
            .filter(|s| {
                idx.lookup(s.label)
                    .any(|nd| idx.occurrences(idx.node_unitig()[nd as usize]).0.contains(&truth))
            })
            .count();
        let anchors = anchor::anchors(&idx, &seeds, &mapper.p.anchor);
        let all = anchor::candidates(
            &idx,
            anchors.clone(),
            &AnchorParams {
                f_keep: usize::MAX,
                alpha_fam: 0.0,
                ..mapper.p.anchor.clone()
            },
        );
        let kept = anchor::candidates(&idx, anchors, &mapper.p.anchor);
        let fam_of = |c: &anchor::Candidate| {
            let t0 = idx.occurrences(idx.node_unitig()[c.anchors[0].node as usize]).0[0];
            idx.tx_gene()[t0 as usize]
        };
        let class = if r.hits.is_empty() {
            "unmapped"
        } else if !r.hits.iter().any(|h| h.reverse == want_rev) {
            if in_truth == 0 {
                "wrong strand; no read seed occurs in the truth"
            } else {
                "wrong strand chosen"
            }
        } else if in_truth == 0 {
            "right strand; no read seed occurs in the truth"
        } else if !all.iter().any(|c| fam_of(c) == truth_gene) {
            "right strand; truth gene has no anchors (filtered)"
        } else if !kept.iter().any(|c| fam_of(c) == truth_gene) {
            "right strand; truth gene cut by candidate selection"
        } else {
            "right strand; truth gene kept, transcript lost in chaining/reporting"
        };
        *counts.entry(class).or_default() += 1;
        if class.contains("transcript lost") {
            let rank = |m: &Mapper<OracleChainer>| {
                let r = m.map(&seq);
                let pos = r.hits.iter().position(|h| h.tx == truth);
                (pos, pos.map(|i| (r.hits[i].score, r.hits[0].score, r.hits.len())))
            };
            println!(
                "LOST\t{name}\tno_tb {:?}\twide {:?}\twide_no_tb {:?}",
                rank(&no_tb),
                rank(&wide),
                rank(&wide_no_tb)
            );
        }
        if class.contains("cut by candidate") {
            let rank = all.iter().position(|c| fam_of(c) == truth_gene).unwrap();
            let pass = all
                .iter()
                .filter(|c| c.score >= mapper.p.anchor.alpha_fam * all[0].score)
                .count();
            println!(
                "CUT\t{name}\trank {rank}\tscore {:.1}\tbest {:.1}\tratio {:.2}\tpass_alpha {pass}\tn {}",
                all[rank].score,
                all[0].score,
                all[rank].score / all[0].score,
                all.len()
            );
        }
        if std::env::var_os("DIAG_LIST").is_some() {
            println!("FAIL\t{name}\t{}\t{class}", seq.len());
        }
        if class == "wrong strand chosen" {
            // Best chain per orientation, from the public pipeline pieces.
            let best_on = |s: &[u8]| {
                let seeds = mapper.seeder.seed(s);
                let an = anchor::anchors(&idx, &seeds, &mapper.p.anchor);
                let cands = anchor::candidates(&idx, an, &mapper.p.anchor);
                let mut hits: Vec<map::TxHit> = cands
                    .iter()
                    .flat_map(|c| map::Chainer::chain(&mapper.chainer, &idx, c))
                    .collect();
                hits.sort_by(|a, b| b.score.total_cmp(&a.score));
                let truth_hit = hits
                    .iter()
                    .find(|h| h.tx == truth)
                    .map(|h| (h.score, h.anchors, h.qstart, h.qend));
                (
                    seeds.len(),
                    hits.first().map(|h| {
                        (
                            h.score,
                            h.anchors,
                            h.qstart,
                            h.qend,
                            idx.gene_name(idx.tx_gene()[h.tx as usize]).to_string(),
                        )
                    }),
                    truth_hit,
                )
            };
            let other = revcomp(&strand);
            let (bt, bo) = (best_on(&strand), best_on(&other));
            let gap = bo.1.as_ref().map_or(0.0, |h| h.0) - bt.1.as_ref().map_or(0.0, |h| h.0);
            gaps.push(gap);
            if counts[class] > 15 {
                continue;
            }
            println!(
                "  {name} len {} | TRUE strand: {:?} | OTHER strand: {:?} | truth gene {}",
                seq.len(),
                bt,
                bo,
                idx.gene_name(truth_gene)
            );
        }
    }
    let mut v: Vec<_> = counts.into_iter().collect();
    v.sort_by_key(|a| std::cmp::Reverse(a.1));
    gaps.sort_by(f32::total_cmp);
    let below = |x: f32| gaps.iter().filter(|&&g| g < x).count();
    println!(
        "wrong-strand score gaps (other - true): n {} | <5 {} | <15 {} | <30 {} | <60 {} | median {:?}",
        gaps.len(),
        below(5.0),
        below(15.0),
        below(30.0),
        below(60.0),
        gaps.get(gaps.len() / 2)
    );
    println!("{n} reads");
    for (k, c) in v {
        println!("{c:7}  {:.3}%  {k}", 100.0 * c as f64 / n as f64);
    }
}
