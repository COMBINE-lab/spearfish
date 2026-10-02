//! The spearfish index: a colored, compacted graph over base-level seeds.
//!
//! Nodes are distinct seed labels; unitigs are maximal non-branching chains of seeds with
//! fixed relative offsets; each transcript is a walk over whole unitigs, stored as
//! occurrences `(transcript, position of the unitig's first seed)`. A seed in unitig `u` at
//! offset `δ` therefore sits at `pos + δ` in every transcript `t` with an occurrence
//! `(t, pos)` of `u`, which is all the chainer needs. Coordinates are in the parsed
//! (homopolymer-compressed, if HPC is on) transcript. See `build` for construction.

mod build;
pub mod store;

use std::path::Path;

use anyhow::{Context, Result};

pub use build::{BuildParams, BuildReport, Collapse, TxInput, build};
use store::{Store, Writer};

/// `unitig_family` value for masked or unused unitigs.
pub const NO_FAMILY: u32 = u32::MAX;
/// `node_flags` bit: the node's unitig is masked (color touches too many genes).
pub const FLAG_MASKED: u8 = 1;

/// Owned index contents, as produced by `build`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexData {
    pub meta: Vec<(String, String)>,
    pub labels: Vec<u64>,
    pub node_unitig: Vec<u32>,
    pub node_offset: Vec<u32>,
    pub node_len: Vec<u32>,
    pub node_flags: Vec<u8>,
    pub unitig_off: Vec<u64>,
    pub unitig_nodes: Vec<u32>,
    pub unitig_span: Vec<u32>,
    pub unitig_color: Vec<u32>,
    pub unitig_family: Vec<u32>,
    pub occ_off: Vec<u64>,
    pub occ_tx: Vec<u32>,
    pub occ_pos: Vec<u32>,
    pub color_off: Vec<u64>,
    pub color_tx: Vec<u32>,
    pub tx_gene: Vec<u32>,
    pub tx_family: Vec<u32>,
    pub tx_len: Vec<u32>,
    pub tx_plen: Vec<u32>,
    /// Parsed transcript sequences, 2-bit packed (32 bases per word, base i at bits
    /// 2(i mod 32)); each transcript starts on a word boundary: (word offsets, words).
    pub tx_seq: (Vec<u64>, Vec<u64>),
    pub tx_names: Vec<String>,
    pub gene_names: Vec<String>,
}

/// 2-bit pack an ASCII sequence (non-ACGT as A), 32 bases per word, base i at bits
/// 2(i mod 32).
pub(crate) fn pack2(seq: &[u8]) -> Vec<u64> {
    let mut words = vec![0u64; seq.len().div_ceil(32)];
    for (i, &b) in seq.iter().enumerate() {
        let c = lcp_core::encode::CODE[b as usize];
        let c = if c > 3 { 0 } else { c as u64 };
        words[i / 32] |= c << (2 * (i % 32));
    }
    words
}

fn pool(strings: &[String]) -> (Vec<u64>, Vec<u8>) {
    let mut off = vec![0u64];
    let mut bytes = Vec::new();
    for s in strings {
        bytes.extend_from_slice(s.as_bytes());
        off.push(bytes.len() as u64);
    }
    (off, bytes)
}

impl IndexData {
    fn writer(&self) -> Writer {
        let mut w = Writer::default();
        let meta: String = self.meta.iter().map(|(k, v)| format!("{k}={v}\n")).collect();
        w.add("meta", meta.as_bytes());
        w.add("labels", &self.labels);
        w.add("node_unitig", &self.node_unitig);
        w.add("node_offset", &self.node_offset);
        w.add("node_len", &self.node_len);
        w.add("node_flags", &self.node_flags);
        w.add("unitig_off", &self.unitig_off);
        w.add("unitig_nodes", &self.unitig_nodes);
        w.add("unitig_span", &self.unitig_span);
        w.add("unitig_color", &self.unitig_color);
        w.add("unitig_family", &self.unitig_family);
        w.add("occ_off", &self.occ_off);
        w.add("occ_tx", &self.occ_tx);
        w.add("occ_pos", &self.occ_pos);
        w.add("color_off", &self.color_off);
        w.add("color_tx", &self.color_tx);
        w.add("tx_gene", &self.tx_gene);
        w.add("tx_family", &self.tx_family);
        w.add("tx_len", &self.tx_len);
        w.add("tx_plen", &self.tx_plen);
        w.add("tx_seq_off", &self.tx_seq.0);
        w.add("tx_seq", &self.tx_seq.1);
        let (o, b) = pool(&self.tx_names);
        w.add("tx_name_off", &o);
        w.add("tx_names", &b);
        let (o, b) = pool(&self.gene_names);
        w.add("gene_name_off", &o);
        w.add("gene_names", &b);
        w
    }

    /// Write to `path`; returns bytes written.
    pub fn write(&self, path: &Path) -> Result<u64> {
        self.writer().write(path)
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.writer().to_bytes()
    }
}

/// A loaded (memory-mapped) index.
pub struct Index {
    store: Store,
    meta: Vec<(String, String)>,
}

macro_rules! section {
    ($name:ident, $t:ty) => {
        pub fn $name(&self) -> &[$t] {
            self.store.get::<$t>(stringify!($name)).unwrap()
        }
    };
}

impl Index {
    pub fn open(path: &Path) -> Result<Self> {
        Self::from_store(Store::open(path)?)
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        Self::from_store(Store::from_bytes(bytes)?)
    }

    fn from_store(store: Store) -> Result<Self> {
        let text = std::str::from_utf8(store.get::<u8>("meta")?)?;
        let meta = text
            .lines()
            .filter_map(|l| l.split_once('='))
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect();
        let idx = Index { store, meta };
        // Touch every section once so a malformed index fails here, not mid-mapping.
        for name in [
            "labels",
            "node_unitig",
            "node_offset",
            "node_len",
            "node_flags",
            "unitig_off",
            "unitig_nodes",
            "unitig_span",
            "unitig_color",
            "unitig_family",
            "occ_off",
            "occ_tx",
            "occ_pos",
            "color_off",
            "color_tx",
            "tx_gene",
            "tx_family",
            "tx_len",
            "tx_plen",
            "tx_seq_off",
            "tx_seq",
            "tx_name_off",
            "tx_names",
            "gene_name_off",
            "gene_names",
        ] {
            idx.store
                .sizes()
                .iter()
                .find(|(n, _)| n == name)
                .with_context(|| format!("index lacks {name}"))?;
        }
        Ok(idx)
    }

    pub fn meta(&self, key: &str) -> Option<&str> {
        self.meta.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
    }

    section!(labels, u64);
    section!(node_unitig, u32);
    section!(node_offset, u32);
    section!(node_len, u32);
    section!(node_flags, u8);
    section!(unitig_off, u64);
    section!(unitig_nodes, u32);
    section!(unitig_span, u32);
    section!(unitig_color, u32);
    section!(unitig_family, u32);
    section!(occ_off, u64);
    section!(occ_tx, u32);
    section!(occ_pos, u32);
    section!(color_off, u64);
    section!(color_tx, u32);
    section!(tx_gene, u32);
    section!(tx_family, u32);
    section!(tx_len, u32);
    section!(tx_plen, u32);
    section!(tx_seq_off, u64);
    section!(tx_seq, u64);

    pub fn n_nodes(&self) -> usize {
        self.labels().len()
    }

    pub fn n_unitigs(&self) -> usize {
        self.unitig_span().len()
    }

    pub fn n_transcripts(&self) -> usize {
        self.tx_gene().len()
    }

    /// Nodes carrying a seed label (empty if the index lacks it). With global collapsing at
    /// most one; with per-gene collapsing one per gene containing the label.
    #[inline]
    pub fn lookup(&self, label: u64) -> std::ops::Range<u32> {
        let l = self.labels();
        let lo = l.partition_point(|&x| x < label);
        let hi = lo + l[lo..].partition_point(|&x| x == label);
        lo as u32..hi as u32
    }

    /// 2-bit codes (A=0, C=1, G=2, T=3) of the parsed transcript `t` over `[lo, hi)`.
    pub fn tx_codes(&self, t: u32, lo: u32, hi: u32) -> Vec<u8> {
        let words = &self.tx_seq()[self.tx_seq_off()[t as usize] as usize..];
        (lo..hi)
            .map(|i| ((words[(i / 32) as usize] >> (2 * (i % 32))) & 3) as u8)
            .collect()
    }

    /// Nodes of unitig `u`, in order.
    pub fn nodes_of(&self, u: u32) -> &[u32] {
        let off = self.unitig_off();
        &self.unitig_nodes()[off[u as usize] as usize..off[u as usize + 1] as usize]
    }

    /// Occurrences of unitig `u`: transcripts (sorted) and the position of its first seed.
    pub fn occurrences(&self, u: u32) -> (&[u32], &[u32]) {
        let off = self.occ_off();
        let r = off[u as usize] as usize..off[u as usize + 1] as usize;
        (&self.occ_tx()[r.clone()], &self.occ_pos()[r])
    }

    /// Transcripts of color class `c` (sorted).
    pub fn color(&self, c: u32) -> &[u32] {
        let off = self.color_off();
        &self.color_tx()[off[c as usize] as usize..off[c as usize + 1] as usize]
    }

    fn pooled(&self, off: &str, bytes: &str, i: usize) -> &str {
        let o = self.store.get::<u64>(off).unwrap();
        let b = self.store.get::<u8>(bytes).unwrap();
        std::str::from_utf8(&b[o[i] as usize..o[i + 1] as usize]).unwrap()
    }

    pub fn tx_name(&self, t: u32) -> &str {
        self.pooled("tx_name_off", "tx_names", t as usize)
    }

    pub fn gene_name(&self, g: u32) -> &str {
        self.pooled("gene_name_off", "gene_names", g as usize)
    }

    /// Section sizes in bytes, by name.
    pub fn sizes(&self) -> Vec<(String, usize)> {
        self.store.sizes()
    }

    /// Every transcript's seeds, reconstructed from its unitig walk as (node, position),
    /// sorted by position. Used to check the index against fresh seeding.
    pub fn transcript_seeds(&self) -> Vec<Vec<(u32, u32)>> {
        let mut out = vec![Vec::new(); self.n_transcripts()];
        for u in 0..self.n_unitigs() as u32 {
            let (ts, ps) = self.occurrences(u);
            for (&t, &p) in ts.iter().zip(ps) {
                out[t as usize].extend(
                    self.nodes_of(u)
                        .iter()
                        .map(|&v| (v, p + self.node_offset()[v as usize])),
                );
            }
        }
        for v in &mut out {
            v.sort_unstable_by_key(|&(n, p)| (p, n));
        }
        out
    }
}
