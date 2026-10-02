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
        let hier = self.scheme.seeds(seq);
        let top = hier.levels.last().map_or(&[][..], |l| &l.cores[..]);
        let mut v: Vec<Seed> = top
            .iter()
            .map(|c| Seed {
                label: c.label,
                start: c.start,
                end: c.end,
            })
            .collect();
        v.sort_unstable_by_key(|s| (s.start, s.end, s.label));
        v.dedup();
        v
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
