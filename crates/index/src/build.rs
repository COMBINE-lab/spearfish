//! Index construction (design doc, "Construction"), over any base-level seed scheme.
//!
//! 1. Seed every transcript (parallel).
//! 2. Nodes are distinct seed labels; node id = rank of the label in sorted order.
//! 3. Edges are consecutive seeds in a transcript, keyed by (from, to, start offset), so a
//!    unitig's seed offsets are well defined even when seeds overlap or leave gaps
//!    (spec issue 1). Transcript starts and ends are distinct neighbours, so every
//!    transcript enters and leaves unitigs only at their ends.
//! 4. Compact maximal non-branching paths into unitigs (cycles handled).
//! 5. Walk each transcript over unitigs, recording occurrences, and check that the walk
//!    reproduces the transcript's seeds exactly (the round-trip test; the build fails
//!    otherwise).
//! 6. Intern unitig colors; mask unitigs whose color touches more than `f_max` genes; group
//!    transcripts into families (connected through unmasked unitigs).

use std::collections::HashMap;
use std::time::Instant;

use anyhow::{Result, bail};
use rayon::prelude::*;
use seeds::{BaseSeeder, Seed};

use crate::{IndexData, NO_FAMILY};

pub struct BuildParams {
    /// Seed scheme spec; its highest level is the base level.
    pub spec: String,
    pub hpc: bool,
    /// Unitigs whose color touches more genes than this are masked.
    pub f_max: u32,
    /// Build even if a label maps to two different contents (otherwise an error).
    pub allow_collisions: bool,
    pub collapse: Collapse,
}

/// What a graph node is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Collapse {
    /// One node per distinct label, across the whole transcriptome.
    Global,
    /// One node per (label, gene): repeats shared between genes don't merge their graphs.
    Gene,
}

impl Collapse {
    pub fn name(self) -> &'static str {
        match self {
            Collapse::Global => "global",
            Collapse::Gene => "gene",
        }
    }
}

#[derive(Clone, Debug)]
pub struct TxInput {
    pub name: String,
    pub gene: String,
    pub seq: Vec<u8>,
}

#[derive(Debug, Default, Clone)]
pub struct BuildReport {
    pub transcripts: usize,
    pub genes: usize,
    pub unseeded_transcripts: usize,
    pub seed_occurrences: u64,
    pub nodes: usize,
    pub unitigs: usize,
    pub cyclic_unitigs: usize,
    pub unitig_n50_seeds: usize,
    pub occurrences: usize,
    pub colors: usize,
    pub masked_unitigs: usize,
    pub families: usize,
    pub largest_family_tx: usize,
    pub multi_gene_families: usize,
    pub inconsistent_lengths: usize,
    /// Labels seen with two different contents (64-bit collisions or label bugs).
    pub label_collisions: usize,
    /// Why unitigs start: only transcript starts precede the head; several distinct
    /// predecessors; or the unique predecessor branches.
    pub heads_start: usize,
    pub heads_multi_pred: usize,
    pub heads_pred_branches: usize,
    pub seconds: Vec<(&'static str, f64)>,
}

const UNSET: u64 = u64::MAX;
const MULTI: u64 = u64::MAX - 1;
const START: u64 = ((u32::MAX - 2) as u64) << 32;
const END: u64 = ((u32::MAX - 3) as u64) << 32;
const MAX_NODES: usize = (u32::MAX - 4) as usize;

#[inline]
fn enc(node: u32, delta: u32) -> u64 {
    ((node as u64) << 32) | delta as u64
}

#[inline]
fn real(e: u64) -> Option<(u32, u32)> {
    let n = (e >> 32) as u32;
    (e != UNSET && e != MULTI && (n as usize) < MAX_NODES).then_some((n, e as u32))
}

fn note(e: &mut u64, v: u64) {
    if *e == UNSET {
        *e = v;
    } else if *e != v {
        *e = MULTI;
    }
}

struct Uf(Vec<u32>);

impl Uf {
    fn find(&mut self, mut x: u32) -> u32 {
        while self.0[x as usize] != x {
            let p = self.0[x as usize];
            self.0[x as usize] = self.0[p as usize];
            x = p;
        }
        x
    }

    fn union(&mut self, a: u32, b: u32) {
        let (a, b) = (self.find(a), self.find(b));
        if a != b {
            self.0[a.max(b) as usize] = a.min(b);
        }
    }
}

pub fn build(txs: &[TxInput], p: &BuildParams) -> Result<(IndexData, BuildReport)> {
    let mut rep = BuildReport {
        transcripts: txs.len(),
        ..Default::default()
    };
    let seeder = BaseSeeder::new(&p.spec, p.hpc)?;
    let mut t = Instant::now();
    let mut lap = |rep: &mut BuildReport, what: &'static str| {
        rep.seconds.push((what, t.elapsed().as_secs_f64()));
        t = Instant::now();
    };

    // 1. Seeds, parsed lengths, and an independent hash of each seed's content (for the
    // collision check).
    let seeded: Vec<(Vec<Seed>, u32, Vec<u64>)> = txs
        .par_iter()
        .map(|t| {
            let h = seeder.parsed(&t.seq);
            let seq: &[u8] = h.as_ref().map_or(&t.seq, |h| &h.seq);
            let seeds = seeder.seed_parsed(seq);
            let content = seeds
                .iter()
                .map(|s| {
                    xxhash_rust::xxh3::xxh3_64_with_seed(
                        &seq[s.start as usize..s.end as usize],
                        0xc0_11_15_10,
                    )
                })
                .collect();
            (seeds, seq.len() as u32, content)
        })
        .collect();
    rep.seed_occurrences = seeded.iter().map(|s| s.0.len() as u64).sum();
    rep.unseeded_transcripts = seeded.iter().filter(|s| s.0.is_empty()).count();
    lap(&mut rep, "seed");

    // 2. Genes and nodes. A node is a (label, scope) pair, scope = gene or 0; nodes are
    // sorted by label, so each label's nodes are contiguous.
    let mut gene_ids: HashMap<&str, u32> = HashMap::new();
    let mut gene_names: Vec<&str> = Vec::new();
    let tx_gene: Vec<u32> = txs
        .iter()
        .map(|t| {
            *gene_ids.entry(&t.gene).or_insert_with(|| {
                gene_names.push(&t.gene);
                gene_names.len() as u32 - 1
            })
        })
        .collect();
    rep.genes = gene_names.len();
    let scope = |t: usize| {
        if p.collapse == Collapse::Gene {
            tx_gene[t]
        } else {
            0
        }
    };
    let mut keys: Vec<(u64, u32)> = seeded
        .par_iter()
        .enumerate()
        .flat_map_iter(|(t, s)| s.0.iter().map(move |x| (x.label, scope(t))))
        .collect();
    keys.par_sort_unstable();
    keys.dedup();
    if keys.len() > MAX_NODES {
        bail!("too many distinct seeds: {}", keys.len());
    }
    rep.nodes = keys.len();
    let nodes: Vec<Vec<u32>> = seeded
        .par_iter()
        .enumerate()
        .map(|(t, s)| {
            s.0.iter()
                .map(|x| keys.binary_search(&(x.label, scope(t))).unwrap() as u32)
                .collect()
        })
        .collect();
    let labels: Vec<u64> = keys.iter().map(|k| k.0).collect();
    drop(keys);
    lap(&mut rep, "nodes");

    // 3. Unique predecessor/successor (with offset) per node, and seed lengths.
    let n = labels.len();
    let mut pred = vec![UNSET; n];
    let mut succ = vec![UNSET; n];
    let mut node_len = vec![u32::MAX; n];
    let mut node_content = vec![0u64; n];
    for ((s, _, content), v) in seeded.iter().zip(&nodes) {
        for i in 0..v.len() {
            let c = &mut node_content[v[i] as usize];
            if node_len[v[i] as usize] == u32::MAX {
                *c = content[i];
            } else if *c != content[i] {
                rep.label_collisions += 1;
            }
            let len = s[i].end - s[i].start;
            let l = &mut node_len[v[i] as usize];
            if *l == u32::MAX {
                *l = len;
            } else if *l != len {
                rep.inconsistent_lengths += 1;
                *l = (*l).max(len);
            }
            note(
                &mut pred[v[i] as usize],
                if i == 0 {
                    START
                } else {
                    enc(v[i - 1], s[i].start - s[i - 1].start)
                },
            );
            note(
                &mut succ[v[i] as usize],
                if i + 1 == v.len() {
                    END
                } else {
                    enc(v[i + 1], s[i + 1].start - s[i].start)
                },
            );
        }
    }
    drop(node_content);
    if rep.label_collisions > 0 && !p.allow_collisions {
        bail!(
            "{} seed occurrences share a label with different content",
            rep.label_collisions
        );
    }
    lap(&mut rep, "edges");

    // 4. Unitigs. `v` continues `u` iff u's only successor is v and v's only predecessor is
    // u, with the same offset, and u != v.
    let cont = |u: u32| -> Option<(u32, u32)> {
        let (v, d) = real(succ[u as usize])?;
        (v != u && pred[v as usize] == enc(u, d)).then_some((v, d))
    };
    let is_head = |v: u32| match real(pred[v as usize]) {
        Some((u, d)) => !(u != v && succ[u as usize] == enc(v, d)),
        None => true,
    };
    let acyclic = (0..n as u32).filter(|&v| is_head(v)).count();
    for (v, &e) in pred.iter().enumerate() {
        match e {
            START => rep.heads_start += 1,
            MULTI => rep.heads_multi_pred += 1,
            e => {
                if let Some((u, d)) = real(e)
                    && (u as usize == v || succ[u as usize] != enc(v as u32, d))
                {
                    rep.heads_pred_branches += 1;
                }
            }
        }
    }
    let mut node_unitig = vec![u32::MAX; n];
    let mut node_offset = vec![0u32; n];
    let mut unitig_off: Vec<u64> = vec![0];
    let mut unitig_nodes: Vec<u32> = Vec::with_capacity(n);
    let mut unitig_span: Vec<u32> = Vec::new();
    let mut walk = |start: u32, node_unitig: &mut Vec<u32>, cyclic: bool| {
        let u = unitig_span.len() as u32;
        let (mut cur, mut delta) = (start, 0u32);
        loop {
            node_unitig[cur as usize] = u;
            node_offset[cur as usize] = delta;
            unitig_nodes.push(cur);
            match cont(cur) {
                Some((next, d)) if !(cyclic && next == start) => {
                    cur = next;
                    delta += d;
                }
                _ => break,
            }
        }
        unitig_span.push(delta + node_len[cur as usize]);
        unitig_off.push(unitig_nodes.len() as u64);
    };
    for v in 0..n as u32 {
        if is_head(v) {
            walk(v, &mut node_unitig, false);
        }
    }
    for v in 0..n as u32 {
        if node_unitig[v as usize] == u32::MAX {
            walk(v, &mut node_unitig, true);
        }
    }
    rep.unitigs = unitig_span.len();
    rep.cyclic_unitigs = rep.unitigs - acyclic;
    let mut lens: Vec<usize> = unitig_off.windows(2).map(|w| (w[1] - w[0]) as usize).collect();
    lens.sort_unstable_by(|a, b| b.cmp(a));
    let half = lens.iter().sum::<usize>() / 2;
    let mut acc = 0;
    rep.unitig_n50_seeds = lens
        .iter()
        .copied()
        .find(|&l| {
            acc += l;
            acc >= half
        })
        .unwrap_or(0);
    lap(&mut rep, "unitigs");

    // 5. Transcript walks: occurrences, with the round-trip check.
    let walks: Vec<Result<Vec<(u32, u32)>>> = seeded
        .par_iter()
        .zip(&nodes)
        .enumerate()
        .map(|(ti, ((s, _, _), v))| {
            let mut out = Vec::new();
            let mut i = 0;
            while i < v.len() {
                let u = node_unitig[v[i] as usize];
                if node_offset[v[i] as usize] != 0 {
                    bail!(
                        "transcript {} enters unitig {u} mid-way at seed {i}",
                        txs[ti].name
                    );
                }
                let un = &unitig_nodes[unitig_off[u as usize] as usize..unitig_off[u as usize + 1] as usize];
                if i + un.len() > v.len() {
                    bail!("transcript {} ends inside unitig {u}", txs[ti].name);
                }
                for (k, &w) in un.iter().enumerate() {
                    if v[i + k] != w || s[i + k].start - s[i].start != node_offset[w as usize] {
                        bail!(
                            "transcript {} does not follow unitig {u} at seed {}",
                            txs[ti].name,
                            i + k
                        );
                    }
                }
                out.push((u, s[i].start));
                i += un.len();
            }
            Ok(out)
        })
        .collect();
    let mut occ: Vec<(u32, u32, u32)> = Vec::new();
    for (t, w) in walks.into_iter().enumerate() {
        occ.extend(w?.into_iter().map(|(u, pos)| (u, t as u32, pos)));
    }
    occ.par_sort_unstable();
    rep.occurrences = occ.len();
    let nu = rep.unitigs;
    let mut occ_off = vec![0u64; nu + 1];
    for &(u, _, _) in &occ {
        occ_off[u as usize + 1] += 1;
    }
    for i in 0..nu {
        occ_off[i + 1] += occ_off[i];
    }
    let occ_tx: Vec<u32> = occ.iter().map(|o| o.1).collect();
    let occ_pos: Vec<u32> = occ.iter().map(|o| o.2).collect();
    drop(occ);
    lap(&mut rep, "walks");

    // 6. Genes, colors, masks, families.
    let mut color_ids: HashMap<Vec<u32>, u32> = HashMap::new();
    let mut color_off: Vec<u64> = vec![0];
    let mut color_tx: Vec<u32> = Vec::new();
    let mut unitig_color = Vec::with_capacity(nu);
    let mut masked = vec![false; nu];
    let mut genes_buf: Vec<u32> = Vec::new();
    for u in 0..nu {
        let mut set: Vec<u32> = occ_tx[occ_off[u] as usize..occ_off[u + 1] as usize].to_vec();
        set.dedup(); // occurrences are sorted by transcript
        genes_buf.clear();
        genes_buf.extend(set.iter().map(|&t| tx_gene[t as usize]));
        genes_buf.sort_unstable();
        genes_buf.dedup();
        masked[u] = genes_buf.len() as u32 > p.f_max;
        let next = color_ids.len() as u32;
        let id = *color_ids.entry(set).or_insert_with_key(|s| {
            color_tx.extend_from_slice(s);
            color_off.push(color_tx.len() as u64);
            next
        });
        unitig_color.push(id);
    }
    rep.colors = color_ids.len();
    rep.masked_unitigs = masked.iter().filter(|&&m| m).count();
    let mut uf = Uf((0..txs.len() as u32).collect());
    for u in 0..nu {
        if !masked[u] {
            let ts = &occ_tx[occ_off[u] as usize..occ_off[u + 1] as usize];
            for &t in &ts[1..] {
                uf.union(ts[0], t);
            }
        }
    }
    let mut fam_ids: HashMap<u32, u32> = HashMap::new();
    let tx_family: Vec<u32> = (0..txs.len() as u32)
        .map(|t| {
            let r = uf.find(t);
            let next = fam_ids.len() as u32;
            *fam_ids.entry(r).or_insert(next)
        })
        .collect();
    rep.families = fam_ids.len();
    let mut fam_tx = vec![0usize; rep.families];
    let mut fam_genes: Vec<Vec<u32>> = vec![Vec::new(); rep.families];
    for (t, &f) in tx_family.iter().enumerate() {
        fam_tx[f as usize] += 1;
        fam_genes[f as usize].push(tx_gene[t]);
    }
    rep.largest_family_tx = fam_tx.iter().copied().max().unwrap_or(0);
    rep.multi_gene_families = fam_genes
        .iter_mut()
        .map(|g| {
            g.sort_unstable();
            g.dedup();
            g.len()
        })
        .filter(|&l| l > 1)
        .count();
    let unitig_family: Vec<u32> = (0..nu)
        .map(|u| {
            if masked[u] || occ_off[u] == occ_off[u + 1] {
                NO_FAMILY
            } else {
                tx_family[occ_tx[occ_off[u] as usize] as usize]
            }
        })
        .collect();
    let node_flags: Vec<u8> = node_unitig.iter().map(|&u| masked[u as usize] as u8).collect();
    lap(&mut rep, "colors");

    let data = IndexData {
        meta: vec![
            ("spec".into(), seeder.spec().to_string()),
            ("hpc".into(), p.hpc.to_string()),
            ("f_max".into(), p.f_max.to_string()),
            ("collapse".into(), p.collapse.name().to_string()),
            (
                "lcp_parser_version".into(),
                lcp_core::LCP_PARSER_VERSION.to_string(),
            ),
            (
                "lcp_seed".into(),
                format!("{:#x}", lcp_core::LcpConfig::default().seed),
            ),
            (
                "input_xxh3".into(),
                format!("{:016x}", {
                    let mut h = xxhash_rust::xxh3::Xxh3::new();
                    for t in txs {
                        h.update(t.name.as_bytes());
                        h.update(&t.seq);
                    }
                    h.digest()
                }),
            ),
        ],
        labels,
        node_unitig,
        node_offset,
        node_len,
        node_flags,
        unitig_off,
        unitig_nodes,
        unitig_span,
        unitig_color,
        unitig_family,
        occ_off,
        occ_tx,
        occ_pos,
        color_off,
        color_tx,
        tx_gene,
        tx_family,
        tx_len: txs.iter().map(|t| t.seq.len() as u32).collect(),
        tx_plen: seeded.iter().map(|s| s.1).collect(),
        tx_names: txs.iter().map(|t| t.name.clone()).collect(),
        gene_names: gene_names.iter().map(|s| s.to_string()).collect(),
    };
    Ok((data, rep))
}
