//! Locally consistent parsing (LCP) behind an abstract [`Parser`] interface.
//!
//! Everything downstream (index, seeding, chaining) is written against [`Parser`] and
//! [`Hierarchy`], so the LCP variant can be swapped. See the design doc, "LCP preliminaries".

mod bits;
pub mod encode;
mod lcp;
pub mod locality;

pub use lcp::{LcpConfig, LcpParser};

/// One core of a parse: a labelled interval of the input.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Core {
    /// 64-bit content label. Equal labels at the same level imply equal content
    /// (up to hash collisions, which the index build counts).
    pub label: u64,
    /// Content span `[start, end)` in input (nucleotide) coordinates.
    pub start: u32,
    pub end: u32,
    /// Children: for level 1 the base range (equal to the span); for level k > 1 an
    /// index range `[child_lo, child_hi)` into the level k-1 cores. Children of
    /// adjacent cores may overlap.
    pub child_lo: u32,
    pub child_hi: u32,
}

impl Core {
    #[inline]
    pub fn len(&self) -> u32 {
        self.end - self.start
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }
}

/// The cores of one level, in increasing start order.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Level {
    pub cores: Vec<Core>,
    /// Start of the context that deterministic coin tossing (DCT) made each core depend on
    /// (`ctx_start[i] <= cores[i].start`). This is the start lcptools reports, and a lower
    /// bound for the left context radius. Parallel to `cores`.
    pub ctx_start: Vec<u32>,
}

/// A parse at levels `1..=levels.len()`. `levels[0]` is level 1.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Hierarchy {
    /// Length of the parsed input.
    pub len: u32,
    pub levels: Vec<Level>,
}

impl Hierarchy {
    /// Cores at level `k` (1-based). Empty if `k` is out of range.
    pub fn level(&self, k: u8) -> &[Core] {
        match (k as usize).checked_sub(1).and_then(|i| self.levels.get(i)) {
            Some(l) => &l.cores,
            None => &[],
        }
    }

    pub fn num_levels(&self) -> u8 {
        self.levels.len() as u8
    }
}

/// Parser identity, recorded in index metadata so index and reads are guaranteed to use
/// bit-identical parses.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParserMeta {
    pub name: &'static str,
    pub version: u32,
    pub seed: u64,
    pub dct_iters: u8,
}

pub trait Parser: Send + Sync {
    /// Parse `seq` (ASCII nucleotides) up to `max_level`. The result always has exactly
    /// `max_level` levels, some possibly empty. Non-ACGT characters split the input: no core
    /// spans one.
    fn parse(&self, seq: &[u8], max_level: u8) -> Hierarchy;

    fn meta(&self) -> ParserMeta;
}
