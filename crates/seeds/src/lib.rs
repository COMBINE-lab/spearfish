//! Seed schemes behind one trait, so the index and chainer are seed-agnostic and the
//! LCP-versus-minimizer/syncmer/strobemer ablations are configuration changes.
//!
//! A scheme returns a [`Hierarchy`]: single-level schemes (minimizers, syncmers, strobemers)
//! fill level 1 only. All schemes are forward-strand (non-canonical), like LCP, and skip
//! non-ACGT characters. Labels are 64-bit content hashes. Schemes are built from spec strings
//! such as `minimizer:k=15,w=10` (see [`from_spec`]).

mod base;
mod kmer;
mod linked;
mod spec;

use lcp_core::{Hierarchy, LcpParser, Parser};

pub use base::{BaseSeeder, Seed};
pub use kmer::{Minimizer, Randstrobe, Syncmer, SyncmerKind};
pub use linked::LinkedCores;
pub use spec::{SpecError, from_spec};

pub trait SeedScheme: Send + Sync {
    fn seeds(&self, seq: &[u8]) -> Hierarchy;
    /// The highest level's seeds, appended to `out` as (label, start, end) in start order.
    fn top_into(&self, seq: &[u8], out: &mut Vec<(u64, u32, u32)>) {
        let h = self.seeds(seq);
        if let Some(l) = h.levels.last() {
            out.extend(l.cores.iter().map(|c| (c.label, c.start, c.end)));
        }
    }
    /// Canonical spec string; `from_spec(&s.name())` rebuilds the scheme.
    fn name(&self) -> String;
}

/// LCP cores at levels `1..=max_level`.
pub struct LcpSeeds {
    pub parser: LcpParser,
    pub max_level: u8,
}

thread_local! {
    static PARSE: std::cell::RefCell<lcp_core::ParseScratch> = std::cell::RefCell::new(Default::default());
}

impl SeedScheme for LcpSeeds {
    fn seeds(&self, seq: &[u8]) -> Hierarchy {
        self.parser.parse(seq, self.max_level)
    }

    fn top_into(&self, seq: &[u8], out: &mut Vec<(u64, u32, u32)>) {
        PARSE.with(|sc| {
            self.parser
                .top_cores(seq, self.max_level, &mut sc.borrow_mut(), out)
        });
    }

    fn name(&self) -> String {
        format!("lcp:K={}", self.max_level)
    }
}
