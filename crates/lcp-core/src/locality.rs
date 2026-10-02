//! Empirical context radii (design doc, Property 2).
//!
//! A level-k core at `[s, e)` depends only on `S[s - λ_k, e + ρ_k)`. We measure λ_k and ρ_k
//! by applying single substitutions and recording how far from the edit the changed cores
//! lie.

use std::collections::HashSet;

use crate::Hierarchy;

/// Largest observed left and right context distances, per level (index 0 = level 1).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Radii {
    pub lambda: Vec<u32>,
    pub rho: Vec<u32>,
}

impl Radii {
    pub fn new(levels: usize) -> Self {
        Radii {
            lambda: vec![0; levels],
            rho: vec![0; levels],
        }
    }

    pub fn merge(&mut self, other: &Radii) {
        for (a, b) in self.lambda.iter_mut().zip(&other.lambda) {
            *a = (*a).max(*b);
        }
        for (a, b) in self.rho.iter_mut().zip(&other.rho) {
            *a = (*a).max(*b);
        }
    }
}

/// Radii implied by comparing the parse of a sequence with the parse after a substitution
/// at `p` (both parses of equal-length sequences). A core counts as changed if its
/// `(label, start, end)` occurs in one parse but not the other. For a changed core at
/// `[s, e)`, an edit at `p < s` witnesses `λ >= s - p` and an edit at `p >= e` witnesses
/// `ρ >= p - e + 1`.
pub fn substitution_radii(orig: &Hierarchy, mutated: &Hierarchy, p: u32) -> Radii {
    let levels = orig.levels.len().min(mutated.levels.len());
    let mut r = Radii::new(levels);
    for k in 0..levels {
        let key = |c: &crate::Core| (c.label, c.start, c.end);
        let a: HashSet<_> = orig.levels[k].cores.iter().map(key).collect();
        let b: HashSet<_> = mutated.levels[k].cores.iter().map(key).collect();
        for &(_, s, e) in a.symmetric_difference(&b) {
            if p < s {
                r.lambda[k] = r.lambda[k].max(s - p);
            } else if p >= e {
                r.rho[k] = r.rho[k].max(p - e + 1);
            }
        }
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{LcpParser, Parser};

    fn xorshift(x: &mut u64) -> u64 {
        *x ^= *x << 13;
        *x ^= *x >> 7;
        *x ^= *x << 17;
        *x
    }

    /// Level-1 cores are decided by rules that look one symbol back and three ahead, so a
    /// substitution can only change level-1 cores within a few bases.
    #[test]
    fn level1_radii_are_small() {
        let p = LcpParser::default();
        let mut x = 42u64;
        let mut total = Radii::new(4);
        for _ in 0..200 {
            let seq: Vec<u8> = (0..2000)
                .map(|_| b"ACGT"[(xorshift(&mut x) & 3) as usize])
                .collect();
            let h = p.parse(&seq, 4);
            for _ in 0..20 {
                let pos = (xorshift(&mut x) % seq.len() as u64) as usize;
                let mut m = seq.clone();
                m[pos] = b"ACGT"[((crate::encode::CODE[seq[pos] as usize] as u64 + 1 + xorshift(&mut x) % 3)
                    % 4) as usize];
                let hm = p.parse(&m, 4);
                total.merge(&substitution_radii(&h, &hm, pos as u32));
            }
        }
        eprintln!("observed radii: {total:?}");
        assert!(total.lambda[0] <= 4 && total.rho[0] <= 4, "{total:?}");
        // Radii grow with level.
        for k in 1..4 {
            assert!(
                total.lambda[k] + total.rho[k] >= total.lambda[k - 1] + total.rho[k - 1],
                "{total:?}"
            );
        }
    }
}
