//! A set of seed schemes evaluated side by side, shared by the analysis subcommands.
//! Each scheme contributes one *unit* per level it reports (LCP: one per level).

use std::path::Path;

use anyhow::{Context, Result, bail};
use lcp_core::hpc::Hpc;
use seeds::{SeedScheme, from_spec};

#[derive(Clone, Debug)]
pub struct Unit {
    pub scheme: usize,
    pub level: usize,
}

pub struct SeedSet {
    pub schemes: Vec<Box<dyn SeedScheme>>,
    pub units: Vec<Unit>,
}

/// Seeds `(label, start, end)` in raw coordinates, per unit.
pub type UnitSeeds = Vec<Vec<(u64, u32, u32)>>;

impl SeedSet {
    /// Schemes from `specs` plus one spec per line of `file` (`#` comments allowed). LCP
    /// levels below `min_lcp_level` are skipped.
    pub fn load(specs: &[String], file: Option<&Path>, min_lcp_level: usize) -> Result<Self> {
        let mut specs = specs.to_vec();
        if let Some(f) = file {
            let text = std::fs::read_to_string(f).with_context(|| format!("reading {}", f.display()))?;
            specs.extend(
                text.lines()
                    .map(|l| l.split('#').next().unwrap().trim())
                    .filter(|l| !l.is_empty())
                    .map(String::from),
            );
        }
        if specs.is_empty() {
            bail!("no seed schemes given (use --scheme or --schemes-file)");
        }
        let schemes: Vec<Box<dyn SeedScheme>> =
            specs.iter().map(|s| from_spec(s)).collect::<Result<_, _>>()?;
        let mut units = Vec::new();
        for (i, s) in schemes.iter().enumerate() {
            let is_lcp = s.name().starts_with("lcp:");
            // A probe sequence reveals how many levels the scheme reports.
            for level in 0..s.seeds(b"ACGTACGTAC").levels.len() {
                if !is_lcp || level + 1 >= min_lcp_level {
                    units.push(Unit { scheme: i, level });
                }
            }
        }
        Ok(SeedSet { schemes, units })
    }

    pub fn name(&self, u: &Unit) -> String {
        self.schemes[u.scheme].name()
    }

    pub fn seeds(&self, raw: &[u8], hpc: bool) -> UnitSeeds {
        let h = hpc.then(|| Hpc::new(raw));
        let seq: &[u8] = h.as_ref().map_or(raw, |h| &h.seq);
        let hier: Vec<_> = self.schemes.iter().map(|s| s.seeds(seq)).collect();
        self.units
            .iter()
            .map(|u| {
                hier[u.scheme].levels[u.level]
                    .cores
                    .iter()
                    .map(|c| {
                        let (s, e) = h.as_ref().map_or((c.start, c.end), |h| h.to_raw(c.start, c.end));
                        (c.label, s, e)
                    })
                    .collect()
            })
            .collect()
    }
}

/// Is there an entry `(label, p)` in `sorted` with `|p - pos| <= tol`?
pub fn has_near(sorted: &[(u64, u32)], label: u64, pos: u32, tol: u32) -> bool {
    let i = sorted.partition_point(|&(l, p)| (l, p) < (label, pos.saturating_sub(tol)));
    sorted.get(i).is_some_and(|&(l, p)| l == label && p <= pos + tol)
}

/// For each transcript seed `(label, start, _)`: does the read have the same label within
/// `tol` of the mapped start `o2r[start]`?
pub fn survivors(orig: &[(u64, u32, u32)], read: &[(u64, u32, u32)], o2r: &[u32], tol: u32) -> Vec<bool> {
    let mut sorted: Vec<(u64, u32)> = read.iter().map(|&(l, s, _)| (l, s)).collect();
    sorted.sort_unstable();
    orig.iter()
        .map(|&(l, s, _)| has_near(&sorted, l, o2r[s as usize], tol))
        .collect()
}
