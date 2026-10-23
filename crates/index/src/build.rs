//! Index construction, over any base-level seed scheme.
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
    /// Group transcripts into genes from their sequence instead of `TxInput::gene`, for
    /// transcriptomes without a transcript-to-gene annotation (see [`InferGenes`]).
    pub infer_genes: Option<InferGenes>,
}

/// Gene inference: two transcripts share a group if they share a run of `span` consecutive
/// base-level seeds (same labels at the same relative offsets) held by at most `max_tx`
/// transcripts. Isoforms that share an exon share long runs; a single seed shared by chance,
/// or a run inside a repeat held by many transcripts, links nothing.
///
/// A second, attach-only pass then joins each small group (at most `attach_max` transcripts)
/// to a larger group with which it shares runs of `attach_span` seeds, if those runs point to
/// that group unambiguously: short isoforms and fragments often share real exons with their
/// gene but no long run. Large groups never merge in this pass, so it cannot chain groups.
/// `attach_span = 0` skips it.
#[derive(Clone, Copy, Debug)]
pub struct InferGenes {
    pub span: usize,
    pub max_tx: u32,
    pub attach_span: usize,
    pub attach_max: usize,
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
    /// The gene index of each transcript (inferred or given), and the gene names.
    pub tx_gene: Vec<u32>,
    pub gene_names: Vec<String>,
    /// With `infer_genes`: the largest group, in transcripts.
    pub largest_group_tx: usize,
    pub seconds: Vec<(&'static str, f64)>,
}

const UNSET: u64 = u64::MAX;
const MULTI: u64 = u64::MAX - 1;
const START: u64 = ((u32::MAX - 2) as u64) << 32;
const END: u64 = ((u32::MAX - 3) as u64) << 32;
const MAX_NODES: usize = (u32::MAX - 4) as usize;

/// Add one field to an identity hash: its length (little-endian u64), then its bytes.
fn hash_field(h: &mut xxhash_rust::xxh3::Xxh3, bytes: &[u8]) {
    h.update(&(bytes.len() as u64).to_le_bytes());
    h.update(bytes);
}

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

/// Per-transcript results of the seeding pass.
struct Seeded {
    seeds: Vec<Seed>,
    /// Parsed (HPC or raw) length.
    plen: u32,
    /// Independent hash of each seed's content (for the collision check).
    content: Vec<u64>,
    /// 2-bit packed parsed sequence.
    packed: Vec<u64>,
    /// (HPC) The packed raw sequence and its raw-position checkpoints (see [`raw_packing`]).
    raw: Vec<u64>,
    raw_ck: Vec<u32>,
}

/// Packed raw sequence and raw starts of every `RAW_CK`-th parsed base (through the parsed
/// length), if the packed sequence's runs are exactly the parsed runs (non-ACGT bases pack
/// as A and could merge runs).
fn raw_packing(raw: &[u8], h: &lcp_core::hpc::Hpc) -> (Vec<u64>, Vec<u32>) {
    let code = |b: u8| {
        let c = lcp_core::encode::CODE[b as usize];
        if c > 3 { 0 } else { c }
    };
    let plen = h.seq.len() as u32;
    let mut runs = 0u32;
    for (i, &b) in raw.iter().enumerate() {
        if i == 0 || code(b) != code(raw[i - 1]) {
            if runs >= plen || h.to_raw(runs, runs).0 != i as u32 {
                return (Vec::new(), Vec::new());
            }
            runs += 1;
        }
    }
    if runs != plen {
        return (Vec::new(), Vec::new());
    }
    let ck = (0..=plen)
        .step_by(crate::RAW_CK as usize)
        .map(|i| h.to_raw(i, i).0)
        .collect();
    (crate::pack2(raw), ck)
}

/// (run hash, transcript) for every run of `span` consecutive seeds, sorted and deduplicated.
fn seed_runs(seeded: &[Seeded], span: usize) -> Vec<(u64, u32)> {
    let span = span.max(1);
    let mut pairs: Vec<(u64, u32)> = seeded
        .par_iter()
        .enumerate()
        .flat_map_iter(|(t, s)| {
            let seeds = &s.seeds;
            (0..seeds.len().saturating_sub(span - 1)).map(move |i| {
                let mut h = xxhash_rust::xxh3::Xxh3::with_seed(0x6e_6e_65_73);
                for x in &seeds[i..i + span] {
                    h.update(&x.label.to_le_bytes());
                    h.update(&(x.start - seeds[i].start).to_le_bytes());
                }
                (h.digest(), t as u32)
            })
        })
        .collect();
    pairs.par_sort_unstable();
    pairs.dedup();
    pairs
}

/// The attach-only pass of [`InferGenes`]: each group of at most `attach_max` transcripts joins
/// the larger group holding the most of its `attach_span`-seed runs (runs in at most `max_tx`
/// transcripts), if that group holds at least twice as many as any other.
fn attach_small(seeded: &[Seeded], uf: &mut Uf, inf: InferGenes) {
    let n = uf.0.len();
    let root: Vec<u32> = (0..n as u32).map(|t| uf.find(t)).collect();
    let mut size: HashMap<u32, usize> = HashMap::new();
    for &r in &root {
        *size.entry(r).or_default() += 1;
    }
    let small = |r: u32| size[&r] <= inf.attach_max;
    let pairs = seed_runs(seeded, inf.attach_span);
    // (small root, other root) -> runs shared.
    let mut votes: HashMap<(u32, u32), u32> = HashMap::new();
    let mut i = 0;
    while i < pairs.len() {
        let mut j = i + 1;
        while j < pairs.len() && pairs[j].0 == pairs[i].0 {
            j += 1;
        }
        if j - i > 1 && j - i <= inf.max_tx as usize {
            let mut roots: Vec<u32> = pairs[i..j].iter().map(|&(_, t)| root[t as usize]).collect();
            roots.sort_unstable();
            roots.dedup();
            for &a in roots.iter().filter(|&&a| small(a)) {
                for &b in roots.iter().filter(|&&b| b != a && size[&b] > size[&a]) {
                    *votes.entry((a, b)).or_default() += 1;
                }
            }
        }
        i = j;
    }
    let mut best: HashMap<u32, (u32, u32, u32)> = HashMap::new(); // small -> (target, top, second)
    for (&(a, b), &v) in &votes {
        let e = best.entry(a).or_insert((u32::MAX, 0, 0));
        if v > e.1 {
            *e = (b, v, e.1);
        } else if v > e.2 {
            e.2 = v;
        }
    }
    for (a, (b, top, second)) in best {
        if top >= 2 * second.max(1) {
            uf.0[a as usize] = b;
        }
    }
}

/// Union-find over `n` transcripts linking all transcripts that hold a run (`pairs` as
/// [`seed_runs`] gives them) unless more than `max_tx` do. Union by smaller index keeps roots
/// at a group's first transcript.
fn link_runs(n: usize, pairs: &[(u64, u32)], max_tx: u32) -> Uf {
    let mut uf = Uf((0..n as u32).collect());
    let mut i = 0;
    while i < pairs.len() {
        let mut j = i + 1;
        while j < pairs.len() && pairs[j].0 == pairs[i].0 {
            j += 1;
        }
        if j - i > 1 && j - i <= max_tx as usize {
            for &(_, t) in &pairs[i + 1..j] {
                uf.union(pairs[i].1, t);
            }
        }
        i = j;
    }
    uf
}

/// Transcript groups from shared seed runs ([`link_runs`], then [`attach_small`]). Groups are
/// named `group<k>`, numbered in order of their first transcript. Returns (group per
/// transcript, names, largest group size).
fn infer_genes(txs: &[TxInput], seeded: &[Seeded], inf: InferGenes) -> (Vec<u32>, Vec<String>, usize) {
    let n = txs.len();
    let pairs = seed_runs(seeded, inf.span);
    let mut uf = link_runs(n, &pairs, inf.max_tx);
    drop(pairs);
    if inf.attach_span > 0 {
        attach_small(seeded, &mut uf, inf);
    }
    let mut group_of_root: HashMap<u32, u32> = HashMap::new();
    let mut size: Vec<usize> = Vec::new();
    let tx_gene: Vec<u32> = (0..n as u32)
        .map(|t| {
            let r = uf.find(t);
            let g = *group_of_root.entry(r).or_insert_with(|| {
                size.push(0);
                size.len() as u32 - 1
            });
            size[g as usize] += 1;
            g
        })
        .collect();
    let names = (0..size.len()).map(|g| format!("group{g}")).collect();
    let largest = size.iter().copied().max().unwrap_or(0);
    (tx_gene, names, largest)
}

/// Build the index. `txs` is consumed so that the raw sequences can be freed once seeded.
pub fn build(mut txs: Vec<TxInput>, p: &BuildParams) -> Result<(IndexData, BuildReport)> {
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
    let seeded: Vec<Seeded> = txs
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
            let (raw, raw_ck) = h.as_ref().map(|h| raw_packing(&t.seq, h)).unwrap_or_default();
            Seeded {
                seeds,
                plen: seq.len() as u32,
                content,
                packed: crate::pack2(seq),
                raw,
                raw_ck,
            }
        })
        .collect();
    // Past seeding, the raw sequences are needed only for their lengths and the input hash.
    let tx_len: Vec<u32> = txs.iter().map(|t| t.seq.len() as u32).collect();
    // Each field length-prefixed, so that no two different inputs concatenate alike.
    let input_xxh3 = {
        let mut h = xxhash_rust::xxh3::Xxh3::new();
        for t in &txs {
            hash_field(&mut h, t.name.as_bytes());
            hash_field(&mut h, &t.seq);
        }
        h.digest()
    };
    for t in &mut txs {
        t.seq = Vec::new();
    }
    rep.seed_occurrences = seeded.iter().map(|s| s.seeds.len() as u64).sum();
    rep.unseeded_transcripts = seeded.iter().filter(|s| s.seeds.is_empty()).count();
    lap(&mut rep, "seed");

    // 2. Genes and nodes. A node is a (label, scope) pair, scope = gene or 0; nodes are
    // sorted by label, so each label's nodes are contiguous.
    let (tx_gene, gene_names): (Vec<u32>, Vec<String>) = match p.infer_genes {
        Some(inf) => {
            let (tx_gene, names, largest) = infer_genes(&txs, &seeded, inf);
            rep.largest_group_tx = largest;
            (tx_gene, names)
        }
        None => {
            let mut gene_ids: HashMap<&str, u32> = HashMap::new();
            let mut names: Vec<String> = Vec::new();
            let tx_gene = txs
                .iter()
                .map(|t| {
                    *gene_ids.entry(&t.gene).or_insert_with(|| {
                        names.push(t.gene.clone());
                        names.len() as u32 - 1
                    })
                })
                .collect();
            (tx_gene, names)
        }
    };
    rep.genes = gene_names.len();
    // The grouping changes the graph and the mapping, so it is part of the index's identity
    // (`meta`, which the index id hashes) along with the input: each transcript's gene name.
    let genes_xxh3 = {
        let mut h = xxhash_rust::xxh3::Xxh3::new();
        for &g in &tx_gene {
            hash_field(&mut h, gene_names[g as usize].as_bytes());
        }
        h.digest()
    };
    lap(&mut rep, "genes");
    let scope = |t: usize| {
        if p.collapse == Collapse::Gene {
            tx_gene[t]
        } else {
            0
        }
    };
    // Isoforms of a gene share most seeds: dedup the keys of runs of transcripts in gene order
    // before the global sort (fixed-size runs, so one huge gene still splits across threads).
    let mut by_gene: Vec<u32> = (0..txs.len() as u32).collect();
    by_gene.sort_by_key(|&t| tx_gene[t as usize]);
    let mut keys: Vec<(u64, u32)> = by_gene
        .par_chunks(256)
        .flat_map_iter(|ts| {
            let mut k: Vec<(u64, u32)> = ts
                .iter()
                .flat_map(|&t| {
                    let t = t as usize;
                    seeded[t].seeds.iter().map(move |x| (x.label, scope(t)))
                })
                .collect();
            k.sort_unstable();
            k.dedup();
            k
        })
        .collect();
    drop(by_gene);
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
            s.seeds
                .iter()
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
    for (sd, v) in seeded.iter().zip(&nodes) {
        let (s, content) = (&sd.seeds, &sd.content);
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
    let mut heads: Vec<u32> = Vec::new();
    for (v, &e) in pred.iter().enumerate() {
        if is_head(v as u32) {
            heads.push(v as u32);
        }
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
    // The chain of (node, offset) from `start` (a cycle stops before returning to it).
    let chain = |start: u32, cyclic: bool| {
        let (mut cur, mut delta) = (start, 0u32);
        let mut c = vec![(cur, delta)];
        while let Some((next, d)) = cont(cur).filter(|&(next, _)| !(cyclic && next == start)) {
            cur = next;
            delta += d;
            c.push((cur, delta));
        }
        c
    };
    let mut add = |c: &[(u32, u32)], node_unitig: &mut Vec<u32>| {
        let u = unitig_span.len() as u32;
        for &(w, d) in c {
            node_unitig[w as usize] = u;
            node_offset[w as usize] = d;
            unitig_nodes.push(w);
        }
        let &(last, d) = c.last().unwrap();
        unitig_span.push(d + node_len[last as usize]);
        unitig_off.push(unitig_nodes.len() as u64);
    };
    // Chains from heads only read pred and succ, so follow them in parallel (the walks are
    // latency-bound pointer chasing), then number them in head order.
    let chains: Vec<Vec<(u32, u32)>> = heads.par_iter().map(|&h| chain(h, false)).collect();
    for c in &chains {
        add(c, &mut node_unitig);
    }
    drop(chains);
    // What is left are pure cycles: cut each at its first node.
    for v in 0..n as u32 {
        if node_unitig[v as usize] == u32::MAX {
            add(&chain(v, true), &mut node_unitig);
        }
    }
    rep.unitigs = unitig_span.len();
    rep.cyclic_unitigs = rep.unitigs - heads.len();
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
        .map(|(ti, (Seeded { seeds: s, .. }, v))| {
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
            ("input_xxh3".into(), format!("{input_xxh3:016x}")),
            ("genes_xxh3".into(), format!("{genes_xxh3:016x}")),
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
        tx_len,
        tx_plen: seeded.iter().map(|s| s.plen).collect(),
        // The per-transcript packings, each starting on a word boundary.
        tx_seq: crate::concat(seeded.iter().map(|s| &s.packed[..])),
        tx_raw: if p.hpc {
            crate::concat(seeded.iter().map(|s| &s.raw[..]))
        } else {
            Default::default()
        },
        tx_raw_ck: if p.hpc {
            crate::concat(seeded.iter().map(|s| &s.raw_ck[..]))
        } else {
            Default::default()
        },
        tx_names: txs.into_iter().map(|t| t.name).collect(),
        gene_names: gene_names.clone(),
        manifest: None,
    };
    rep.tx_gene = data.tx_gene.clone();
    rep.gene_names = gene_names;
    Ok((data, rep))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every transcript holding a run joins one group, however the groups it links were rooted.
    #[test]
    fn runs_link_all_their_groups() {
        // Runs 1 and 2 make groups {0, 10} and {1, 11}; run 3, held by 5, 10 and 11, joins
        // both with {5}.
        let pairs = [(1, 0), (1, 10), (2, 1), (2, 11), (3, 5), (3, 10), (3, 11)];
        let mut uf = link_runs(12, &pairs, 64);
        let roots: Vec<u32> = [0, 1, 5, 10, 11].iter().map(|&t| uf.find(t)).collect();
        assert_eq!(roots, [0; 5]);
        // Runs held by more than max_tx transcripts link nothing.
        let mut uf = link_runs(12, &pairs, 2);
        assert_ne!(uf.find(5), uf.find(10));
        assert_eq!(uf.find(11), 1);
    }
}
