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
    /// The highest level's seeds, appended to `out` as (label, start, end) in start order, at
    /// most one per start.
    fn top_into(&self, seq: &[u8], out: &mut Vec<(u64, u32, u32)>) {
        let h = self.seeds(seq);
        if let Some(l) = h.levels.last() {
            out.extend(l.cores.iter().map(|c| (c.label, c.start, c.end)));
        }
    }
    /// [`Self::top_into`] of a sequence given as 2-bit codes (A=0, C=1, G=2, T=3).
    fn top_codes_into(&self, codes: &[u8], out: &mut Vec<(u64, u32, u32)>) {
        let seq: Vec<u8> = codes.iter().map(|&c| b"ACGT"[c as usize]).collect();
        self.top_into(&seq, out)
    }
    /// The (label, start) of [`Self::top_codes_into`]'s seeds, sorted by (start, end, label)
    /// and without duplicates (as [`BaseSeeder`] returns them), appended to `out`. Schemes that
    /// can skip computing the ends override it.
    fn top_starts_codes_into(&self, codes: &[u8], out: &mut Vec<(u64, u32)>) {
        let mut v = Vec::new();
        self.top_codes_into(codes, &mut v);
        v.sort_unstable_by_key(|&(l, s, e)| (s, e, l));
        v.dedup();
        out.extend(v.iter().map(|&(l, s, _)| (l, s)));
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

    fn top_codes_into(&self, codes: &[u8], out: &mut Vec<(u64, u32, u32)>) {
        PARSE.with(|sc| {
            self.parser
                .top_cores_codes(codes, self.max_level, &mut sc.borrow_mut(), out)
        });
    }

    fn top_starts_codes_into(&self, codes: &[u8], out: &mut Vec<(u64, u32)>) {
        // Cores come in strictly increasing start order: nothing to sort or remove.
        PARSE.with(|sc| {
            self.parser
                .top_starts_codes(codes, self.max_level, &mut sc.borrow_mut(), out)
        });
    }

    fn name(&self) -> String {
        format!("lcp:K={}", self.max_level)
    }
}
