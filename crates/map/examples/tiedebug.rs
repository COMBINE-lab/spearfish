//! For one read, print each reported hit's chain and its whole-read edit distances against
//! the full parsed transcript (infix) and against the transcript's suffix/prefix around the
//! chain, to see where tied isoforms gain or lose edits.
//! Usage: tiedebug <index> <reads.fa|fq> <read name> [max hits]

extern crate spearfish_map as map;

use std::sync::Arc;

use index::Index;
use lcp_core::encode::{codes, revcomp};
use map::{MapParams, Mapper, OracleChainer, align};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let idx = Arc::new(Index::open(a[1].as_ref()).unwrap());
    let max: usize = a.get(4).map_or(12, |s| s.parse().unwrap());
    let text = std::fs::read_to_string(&a[2]).unwrap();
    let raw = find_read(&text, &a[3]).expect("read not found");
    let p = MapParams::default();
    let mapper = Mapper::new(idx.clone(), OracleChainer { p: p.chain.clone() }, p).unwrap();
    let r = mapper.map(&raw);
    let truth = a[3].split('|').nth(1).unwrap_or("");
    for h in r.hits.iter().take(max) {
        let oriented = if h.reverse { revcomp(&raw) } else { raw.clone() };
        let seq = mapper
            .seeder
            .parsed(&oriented)
            .map_or(oriented.clone(), |x| x.seq);
        let codes = codes(&seq);
        let tlen = idx.tx_plen()[h.tx as usize];
        let full = idx.tx_codes(h.tx, 0, tlen);
        let d_full = align::infix_distance(&codes, &full);
        // Read 5' end (before the chain) against the transcript upstream of the chain start,
        // aligned right-to-left (prefix mode on reversed sequences).
        let q5: Vec<u8> = codes[..h.qstart as usize].iter().rev().copied().collect();
        let t5: Vec<u8> = full[..h.tstart as usize].iter().rev().copied().collect();
        let d5 = align::distance(&q5, &t5, align::Mode::Prefix, u32::MAX).unwrap();
        let q3 = &codes[h.qend as usize..];
        let t3 = &full[h.tend as usize..];
        let d3 = align::distance(q3, t3, align::Mode::Prefix, u32::MAX).unwrap();
        let mid = align::distance(
            &codes[h.qstart as usize..h.qend as usize],
            &full[h.tstart as usize..h.tend as usize],
            align::Mode::Global,
            u32::MAX,
        )
        .unwrap();
        let name = idx.tx_name(h.tx);
        println!(
            "{}{name}\t{}\tscore {:.2}\tq [{},{}) of {}\tt [{},{}) of {tlen}\tanchors {}\td_full {d_full}\td5 {d5} (q {} t {})\tmid {mid}\td3 {d3} (q {} t {})",
            if name == truth { "*" } else { " " },
            if h.reverse { '-' } else { '+' },
            h.score,
            h.qstart,
            h.qend,
            codes.len(),
            h.tstart,
            h.tend,
            h.anchors,
            q5.len(),
            t5.len(),
            q3.len(),
            t3.len()
        );
    }
}

/// The sequence of the read whose header's first word is `name`, in FASTA or (four-line)
/// FASTQ text.
fn find_read(text: &str, name: &str) -> Option<Vec<u8>> {
    let mut lines = text.lines();
    while let Some(h) = lines.next() {
        if let Some(id) = h.strip_prefix('@') {
            let seq = lines.next()?;
            lines.nth(1); // `+` and quality lines
            if id.split_whitespace().next() == Some(name) {
                return Some(seq.trim().as_bytes().to_vec());
            }
        } else if let Some(id) = h.strip_prefix('>')
            && id.split_whitespace().next() == Some(name)
        {
            return Some(
                lines
                    .take_while(|l| !l.starts_with('>'))
                    .flat_map(|l| l.trim().bytes())
                    .collect(),
            );
        }
    }
    None
}
