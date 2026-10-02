//! Seed schemes behind one trait, so the index and chainer are seed-agnostic and the
//! LCP-versus-minimizer/syncmer/strobemer ablations are configuration changes.
//!
//! A scheme returns a [`Hierarchy`]: single-level schemes (minimizers, syncmers, strobemers)
//! fill level 1 only. All schemes are forward-strand (non-canonical), like LCP, and skip
//! non-ACGT characters. Labels are 64-bit content hashes. Schemes are built from spec strings
//! such as `minimizer:k=15,w=10` (see [`from_spec`]).

mod kmer;
mod linked;
mod spec;

use lcp_core::{Hierarchy, LcpParser, Parser};

pub use kmer::{Minimizer, Randstrobe, Syncmer, SyncmerKind};
pub use linked::LinkedCores;
pub use spec::{SpecError, from_spec};

pub trait SeedScheme: Send + Sync {
    fn seeds(&self, seq: &[u8]) -> Hierarchy;
    /// Canonical spec string; `from_spec(&s.name())` rebuilds the scheme.
    fn name(&self) -> String;
}

/// LCP cores at levels `1..=max_level`.
pub struct LcpSeeds {
    pub parser: LcpParser,
    pub max_level: u8,
}

impl SeedScheme for LcpSeeds {
    fn seeds(&self, seq: &[u8]) -> Hierarchy {
        self.parser.parse(seq, self.max_level)
    }

    fn name(&self) -> String {
        format!("lcp:K={}", self.max_level)
    }
}
