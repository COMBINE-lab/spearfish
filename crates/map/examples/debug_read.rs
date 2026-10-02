//! Print a read's anchors placed in two transcripts and each transcript's oracle chain.
//! Usage: debug_read <index> <reads.fa> <read name> <tx A> <tx B>
use std::sync::Arc;

use index::Index;
use lcp_core::encode::revcomp;
use map::{Chainer, MapParams, Mapper, OracleChainer, anchor};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let idx = Arc::new(Index::open(a[1].as_ref()).unwrap());
    let text = std::fs::read_to_string(&a[2]).unwrap();
    let seq: Vec<u8> = text
        .split('>')
        .find(|r| r.starts_with(&format!("{}\n", a[3])))
        .expect("read not found")
        .lines()
        .skip(1)
        .flat_map(|l| l.trim().bytes())
        .collect();
    let mapper = Mapper::new(
        idx.clone(),
        OracleChainer {
            p: Default::default(),
        },
        MapParams::default(),
    )
    .unwrap();
    let r = mapper.map(&seq);
    let raw = if r.reverse { revcomp(&seq) } else { seq };
    let seeds = mapper.seeder.seed(&raw);
    let anchors = anchor::anchors(&idx, &seeds, &mapper.p.anchor);
    let cands = anchor::candidates(&idx, anchors, &mapper.p.anchor);
    let tx_id = |name: &str| {
        (0..idx.n_transcripts() as u32)
            .find(|&t| idx.tx_name(t) == name)
            .unwrap()
    };
    for name in [&a[4], &a[5]] {
        let t = tx_id(name);
        println!("== {name} (parsed length {})", idx.tx_plen()[t as usize]);
        for c in &cands {
            for (i, an) in c.anchors.iter().enumerate() {
                let u = idx.node_unitig()[an.node as usize];
                let (ts, ps) = idx.occurrences(u);
                let pos: Vec<u32> = ts
                    .iter()
                    .zip(ps)
                    .filter(|&(&x, _)| x == t)
                    .map(|(_, &p)| p + idx.node_offset()[an.node as usize])
                    .collect();
                if !pos.is_empty() {
                    println!(
                        "  a{i:<3} q [{:>5},{:>5})  t {:?}  w {:.1}",
                        an.qs, an.qe, pos, an.weight
                    );
                }
            }
            for h in mapper.chainer.chain(&idx, c).iter().filter(|h| h.tx == t) {
                println!("  chain: {h:?}");
            }
        }
    }
}
