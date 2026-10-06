//! The base-level seeder shared by index construction and read mapping, so both seed
//! through bit-identical code: a scheme's highest reported level, optionally after
//! homopolymer compression. Coordinates are in the parsed (HPC, if on) sequence.

use lcp_core::hpc::Hpc;

use crate::{SeedScheme, SpecError, from_spec};

/// One seed: content label and span `[start, end)` in parsed coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Seed {
    pub label: u64,
    pub start: u32,
    pub end: u32,
}

pub struct BaseSeeder {
    scheme: Box<dyn SeedScheme>,
    spec: String,
    pub hpc: bool,
}

impl BaseSeeder {
    pub fn new(spec: &str, hpc: bool) -> Result<Self, SpecError> {
        let scheme = from_spec(spec)?;
        Ok(BaseSeeder {
            spec: scheme.name(),
            scheme,
            hpc,
        })
    }

    /// Canonical scheme spec.
    pub fn spec(&self) -> &str {
        &self.spec
    }

    /// The sequence that gets seeded: `raw`, or its homopolymer compression.
    pub fn parsed(&self, raw: &[u8]) -> Option<Hpc> {
        self.hpc.then(|| Hpc::new(raw))
    }

    /// Seeds of `raw`, sorted by (start, end), in parsed coordinates.
    pub fn seed(&self, raw: &[u8]) -> Vec<Seed> {
        let h = self.parsed(raw);
        let seq: &[u8] = h.as_ref().map_or(raw, |h| &h.seq);
        self.seed_parsed(seq)
    }

    /// Seeds of an already parsed (compressed, if HPC is on) sequence.
    pub fn seed_parsed(&self, seq: &[u8]) -> Vec<Seed> {
        let mut v = Vec::new();
        self.seed_parsed_into(seq, &mut v, &mut Vec::new());
        v
    }

    /// `seed_parsed` into `out` (cleared first), with `tmp` as reusable scratch.
    pub fn seed_parsed_into(&self, seq: &[u8], out: &mut Vec<Seed>, tmp: &mut Vec<(u64, u32, u32)>) {
        tmp.clear();
        self.scheme.top_into(seq, tmp);
        out.clear();
        out.extend(tmp.iter().map(|&(label, start, end)| Seed { label, start, end }));
        // Parsers emit seeds nearly in start order; sort only if needed, with the stable sort,
        // which merges the presorted runs (equal keys are equal seeds, so the order is the
        // same as any sort's).
        if !out.is_sorted_by_key(|s| (s.start, s.end, s.label)) {
            out.sort_by_key(|s| (s.start, s.end, s.label));
        }
        out.dedup();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn top_level_and_hpc() {
        let raw = b"ACGTTTTGCAAAGGCTTACGATCGATCGGGATCGATCGTAGCTAGCATCGACTGACTAGCCCGATTACG".repeat(20);
        let s = BaseSeeder::new("lcp:K=3", false).unwrap();
        let seeds = s.seed(&raw);
        assert!(!seeds.is_empty());
        assert!(
            seeds
                .windows(2)
                .all(|w| (w[0].start, w[0].end) <= (w[1].start, w[1].end))
        );
        let h = BaseSeeder::new("lcp:K=3", true).unwrap();
        let hs = h.seed(&raw);
        let plen = Hpc::new(&raw).seq.len() as u32;
        assert!(hs.iter().all(|x| x.end <= plen));
        assert_eq!(s.spec(), "lcp:K=3");
    }
}
