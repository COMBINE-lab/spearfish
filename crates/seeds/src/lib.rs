//! Seed schemes behind one trait, so the index and chainer are seed-agnostic and the
//! LCP-versus-minimizer/syncmer/strobemer ablations are configuration changes.
//!
//! A scheme returns a [`Hierarchy`]: single-level schemes (minimizers, syncmers) fill
//! level 1 only.

use lcp_core::{Hierarchy, LcpParser, Parser};

pub trait SeedScheme: Send + Sync {
    fn seeds(&self, seq: &[u8]) -> Hierarchy;
    /// Short name recorded in index metadata.
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
