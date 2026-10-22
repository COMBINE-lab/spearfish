//! Empirical context radii.
//!
//! A level-k core at `[s, e)` depends only on `S[s - λ_k, e + ρ_k)`. We measure λ_k and ρ_k
//! by applying single edits (substitution, insertion, deletion) and recording how far from
//! the edit the changed cores lie.

use std::collections::HashSet;

use crate::{Core, Hierarchy};

/// A single edit, in coordinates of the original sequence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Edit {
    /// Replace the base at `p`.
    Sub(u32),
    /// Insert a base before original position `p` (it becomes position `p` of the mutated
    /// sequence).
    Ins(u32),
    /// Delete the base at `p`.
    Del(u32),
}

impl Edit {
    /// The original interval the edit touches: one base for a substitution or deletion, an
    /// empty interval at `p` for an insertion.
    pub fn span(self) -> (u32, u32) {
        match self {
            Edit::Sub(p) | Edit::Del(p) => (p, p + 1),
            Edit::Ins(p) => (p, p),
        }
    }

    /// Apply the edit; `base` is the substituted or inserted base.
    pub fn apply(self, seq: &[u8], base: u8) -> Vec<u8> {
        let mut m = seq.to_vec();
        match self {
            Edit::Sub(p) => m[p as usize] = base,
            Edit::Ins(p) => m.insert(p as usize, base),
            Edit::Del(p) => {
                m.remove(p as usize);
            }
        }
        m
    }

    /// Map a core start in the mutated sequence to original coordinates.
    fn start_to_orig(self, x: u32) -> u32 {
        match self {
            Edit::Sub(_) => x,
            Edit::Ins(p) => {
                if x <= p {
                    x
                } else {
                    x - 1
                }
            }
            Edit::Del(p) => {
                if x < p {
                    x
                } else {
                    x + 1
                }
            }
        }
    }

    /// Map a core end (exclusive) in the mutated sequence to original coordinates.
    fn end_to_orig(self, x: u32) -> u32 {
        match self {
            Edit::Sub(_) => x,
            Edit::Ins(p) => {
                if x <= p {
                    x
                } else {
                    x - 1
                }
            }
            Edit::Del(p) => {
                if x <= p {
                    x
                } else {
                    x + 1
                }
            }
        }
    }

    /// The interval the edit occupies in the mutated sequence: the new base for a
    /// substitution or insertion, an empty interval for a deletion.
    fn mutated_span(self) -> (u32, u32) {
        match self {
            Edit::Sub(p) | Edit::Ins(p) => (p, p + 1),
            Edit::Del(p) => (p, p),
        }
    }
}

/// Where a changed core lies relative to the edit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Placement {
    /// The core overlaps the edit.
    Overlap,
    /// The core starts `d >= 1` positions right of the edit, so `λ >= d` (1 = adjacent).
    Right(u32),
    /// The core ends `d >= 1` positions left of the edit, so `ρ >= d` (1 = adjacent).
    Left(u32),
}

/// Placement of core `[s, e)` relative to an edit occupying `[a, b)` in the same coordinates.
fn place(s: u32, e: u32, (a, b): (u32, u32)) -> Placement {
    if b <= s {
        Placement::Right(s - b + 1)
    } else if a >= e {
        Placement::Left(a - e + 1)
    } else {
        Placement::Overlap
    }
}

/// Compare the parses before and after `edit` and call `visit(level_index, placement)` for
/// every changed core. A core is changed if its `(label, start, end)`, with mutated
/// coordinates mapped back to the original, occurs in only one of the two parses. Each
/// changed core is placed in its own parse's coordinates, so a mutated core that starts on
/// an inserted base counts as overlapping the edit. Distances agree across the two
/// coordinate systems.
pub fn edit_effects(
    orig: &Hierarchy,
    mutated: &Hierarchy,
    edit: Edit,
    mut visit: impl FnMut(usize, Placement),
) {
    let levels = orig.levels.len().min(mutated.levels.len());
    let key = |c: &Core| (c.label, edit.start_to_orig(c.start), edit.end_to_orig(c.end));
    for k in 0..levels {
        let before: HashSet<(u64, u32, u32)> = orig.levels[k]
            .cores
            .iter()
            .map(|c| (c.label, c.start, c.end))
            .collect();
        let after: HashSet<(u64, u32, u32)> = mutated.levels[k].cores.iter().map(key).collect();
        for c in &orig.levels[k].cores {
            if !after.contains(&(c.label, c.start, c.end)) {
                visit(k, place(c.start, c.end, edit.span()));
            }
        }
        for c in &mutated.levels[k].cores {
            if !before.contains(&key(c)) {
                visit(k, place(c.start, c.end, edit.mutated_span()));
            }
        }
    }
}

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

    /// Maximum radii witnessed by one edit.
    pub fn from_edit(orig: &Hierarchy, mutated: &Hierarchy, edit: Edit) -> Self {
        let mut r = Radii::new(orig.levels.len().min(mutated.levels.len()));
        edit_effects(orig, mutated, edit, |k, p| match p {
            Placement::Right(d) => r.lambda[k] = r.lambda[k].max(d),
            Placement::Left(d) => r.rho[k] = r.rho[k].max(d),
            Placement::Overlap => {}
        });
        r
    }
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

    fn random_seq(x: &mut u64, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[(xorshift(x) & 3) as usize]).collect()
    }

    #[test]
    fn coordinate_mapping() {
        // Insertion before 5: mutated starts > 5 shift back by one; a core starting at the
        // inserted base starts at 5 in original coordinates.
        let e = Edit::Ins(5);
        assert_eq!(
            (e.start_to_orig(3), e.start_to_orig(5), e.start_to_orig(6)),
            (3, 5, 5)
        );
        assert_eq!((e.end_to_orig(5), e.end_to_orig(6), e.end_to_orig(9)), (5, 5, 8));
        // Deletion of 5: mutated position 5 is original 6.
        let d = Edit::Del(5);
        assert_eq!((d.start_to_orig(4), d.start_to_orig(5)), (4, 6));
        assert_eq!((d.end_to_orig(5), d.end_to_orig(6)), (5, 7));
        assert_eq!(d.apply(b"ACGTACGT", b'N'), b"ACGTAGT");
        assert_eq!(e.apply(b"ACGTACGT", b'N'), b"ACGTANCGT");
    }

    #[test]
    fn unchanged_far_from_edit() {
        // An edit in the middle of a long sequence leaves cores far away untouched, so all
        // changed cores lie within a bounded distance, for every edit type.
        let p = LcpParser::default();
        let mut x = 7u64;
        for _ in 0..50 {
            let seq = random_seq(&mut x, 3000);
            let h = p.parse(&seq, 4);
            for edit in [Edit::Sub(1500), Edit::Ins(1500), Edit::Del(1500)] {
                let base = if seq[1500] == b'A' { b'C' } else { b'A' };
                let hm = p.parse(&edit.apply(&seq, base), 4);
                let r = Radii::from_edit(&h, &hm, edit);
                assert!(
                    r.lambda.iter().chain(&r.rho).all(|&d| d < 1000),
                    "{edit:?}: {r:?}"
                );
            }
        }
    }

    /// Level-1 cores are decided by rules that look one symbol back and three ahead, so an
    /// edit can only change level-1 cores within a few bases.
    #[test]
    fn level1_radii_are_small() {
        let p = LcpParser::default();
        let mut x = 42u64;
        let mut total = Radii::new(4);
        for _ in 0..200 {
            let seq = random_seq(&mut x, 2000);
            let h = p.parse(&seq, 4);
            for i in 0..30 {
                let pos = (xorshift(&mut x) % seq.len() as u64) as u32;
                let base = b"ACGT"[(xorshift(&mut x) & 3) as usize];
                let edit = match i % 3 {
                    0 if base != seq[pos as usize] => Edit::Sub(pos),
                    1 => Edit::Ins(pos),
                    _ => Edit::Del(pos),
                };
                let hm = p.parse(&edit.apply(&seq, base), 4);
                total.merge(&Radii::from_edit(&h, &hm, edit));
            }
        }
        eprintln!("observed radii: {total:?}");
        assert!(total.lambda[0] <= 4 && total.rho[0] <= 4, "{total:?}");
        for k in 1..4 {
            assert!(
                total.lambda[k] + total.rho[k] >= total.lambda[k - 1] + total.rho[k - 1],
                "{total:?}"
            );
        }
    }
}
