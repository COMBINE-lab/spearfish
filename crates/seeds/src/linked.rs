//! Linked LCP cores: long, specific seeds made by chaining nearby base-level cores, in the
//! style of randstrobes. A contiguous upper-level core is destroyed by an error anywhere in
//! its span; a linked seed only by an error in (or near) one of its member cores, since the
//! gaps between members aren't part of the seed.
//!
//! Each base core `j` whose label hash is in the lowest `keep` fraction starts a seed. Its
//! next member is the candidate among those `w_min..=w_max` positions downstream minimising
//! popcount(acc ^ label), where `acc` mixes the members chosen so far; `order` members in all.
//! Candidates are all base cores (`sampled = false`), or only the sampled ones
//! (`sampled = true`, as randstrobes link subsampled syncmers): then a core created by an
//! error competes for partner only if it is itself sampled. Window positions count
//! candidates. Label: hash of the member labels. Span: first member's start to last
//! member's end.

use lcp_core::{Core, Hierarchy, LcpParser, Level, Parser};
use xxhash_rust::xxh3::xxh3_64_with_seed;

use crate::SeedScheme;

#[derive(Clone, Debug)]
pub struct LinkedCores {
    pub parser: LcpParser,
    /// Level of the member cores.
    pub base: u8,
    pub w_min: usize,
    pub w_max: usize,
    /// Number of members (2 = pairs).
    pub order: usize,
    /// Fraction of base cores that start a seed.
    pub keep: f64,
    /// Link only among sampled cores.
    pub sampled: bool,
}

const LINK_SEED: u64 = 0x11_4ced_c0de;

impl LinkedCores {
    /// Linked seeds over base-level cores `all`. Member indices (`child_lo`, `child_hi - 1`)
    /// refer to `all`.
    pub fn link(&self, all: &[Core], out: &mut Vec<Core>) {
        let threshold = if self.keep >= 1.0 {
            u64::MAX
        } else {
            (self.keep * u64::MAX as f64) as u64
        };
        // Candidate pool and each candidate's index in `all`.
        let idx: Vec<usize> = if self.sampled {
            (0..all.len()).filter(|&i| all[i].label <= threshold).collect()
        } else {
            (0..all.len()).collect()
        };
        let cores: Vec<Core> = idx.iter().map(|&i| all[i]).collect();
        let mut members: Vec<usize> = Vec::with_capacity(self.order);
        let mut buf: Vec<u8> = Vec::with_capacity(8 * self.order);
        'seeds: for j in 0..cores.len() {
            if cores[j].label > threshold {
                continue;
            }
            members.clear();
            members.push(j);
            let mut acc = cores[j].label;
            let mut cur = j;
            for _ in 1..self.order {
                let lo = cur + self.w_min;
                if lo >= cores.len() {
                    continue 'seeds;
                }
                let hi = (cur + self.w_max).min(cores.len() - 1);
                let m = (lo..=hi)
                    .min_by_key(|&m| ((acc ^ cores[m].label).count_ones(), m))
                    .unwrap();
                members.push(m);
                acc = acc.rotate_left(17) ^ cores[m].label;
                cur = m;
            }
            buf.clear();
            for &m in &members {
                buf.extend_from_slice(&cores[m].label.to_le_bytes());
            }
            let first = &cores[members[0]];
            let end = members.iter().map(|&m| cores[m].end).max().unwrap();
            out.push(Core {
                label: xxh3_64_with_seed(&buf, LINK_SEED ^ self.base as u64),
                start: first.start,
                end,
                child_lo: idx[members[0]] as u32,
                child_hi: idx[*members.last().unwrap()] as u32 + 1,
            });
        }
    }
}

impl SeedScheme for LinkedCores {
    fn seeds(&self, seq: &[u8]) -> Hierarchy {
        let h = self.parser.parse(seq, self.base);
        let mut cores = Vec::new();
        self.link(h.level(self.base), &mut cores);
        let ctx_start = cores.iter().map(|c| c.start).collect();
        Hierarchy {
            len: h.len,
            levels: vec![Level { cores, ctx_start }],
        }
    }

    fn name(&self) -> String {
        let keep = if self.keep < 1.0 {
            format!(",keep={}", self.keep)
        } else {
            String::new()
        };
        let pool = if self.sampled { ",pool=sampled" } else { "" };
        format!(
            "linked-lcp:b={},wmin={},wmax={},order={}{keep}{pool}",
            self.base, self.w_min, self.w_max, self.order
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x & 3) as usize]
            })
            .collect()
    }

    fn scheme(keep: f64) -> LinkedCores {
        LinkedCores {
            parser: LcpParser::default(),
            base: 3,
            w_min: 1,
            w_max: 4,
            order: 2,
            keep,
            sampled: false,
        }
    }

    #[test]
    fn sampled_pool_links_sampled_cores() {
        let seq = random(100_000, 14);
        let cores = LcpParser::default().parse(&seq, 3).level(3).to_vec();
        let sch = LinkedCores {
            sampled: true,
            ..scheme(0.4)
        };
        let threshold = (0.4 * u64::MAX as f64) as u64;
        let seeds = sch.seeds(&seq);
        assert!(!seeds.level(1).is_empty());
        for s in seeds.level(1) {
            let (a, b) = (s.child_lo as usize, s.child_hi as usize - 1);
            assert!(cores[a].label <= threshold && cores[b].label <= threshold);
            assert_eq!(s.start, cores[a].start);
        }
    }

    #[test]
    fn density_scales_with_keep() {
        let seq = random(300_000, 11);
        let base = LcpParser::default().parse(&seq, 3).level(3).len() as f64;
        let all = scheme(1.0).seeds(&seq).level(1).len() as f64;
        let some = scheme(0.4).seeds(&seq).level(1).len() as f64;
        assert!(all / base > 0.99);
        assert!((some / base - 0.4).abs() < 0.02, "{}", some / base);
    }

    #[test]
    fn members_and_span() {
        let seq = random(20_000, 12);
        let cores = LcpParser::default().parse(&seq, 3).level(3).to_vec();
        for s in scheme(1.0).seeds(&seq).level(1) {
            let (a, b) = (s.child_lo as usize, s.child_hi as usize - 1);
            assert!((1..=4).contains(&(b - a)));
            assert_eq!(s.start, cores[a].start);
            assert!(s.end >= cores[b].end);
        }
    }

    #[test]
    fn gap_errors_do_not_break_seeds() {
        // Substitute a base in the middle of the gap between a seed's two members. A
        // contiguous core over the same span always breaks; a linked seed often survives,
        // because the gap isn't part of it.
        let seq = random(50_000, 13);
        let sch = scheme(1.0);
        let cores = LcpParser::default().parse(&seq, 3).level(3).to_vec();
        let seeds = sch.seeds(&seq).level(1).to_vec();
        let (mut tried, mut survived) = (0, 0);
        for s in seeds.iter().step_by(7) {
            let (a, b) = (&cores[s.child_lo as usize], &cores[s.child_hi as usize - 1]);
            if b.start < a.end + 20 {
                continue;
            }
            let p = ((a.end + b.start) / 2) as usize;
            let mut m = seq.clone();
            m[p] = if m[p] == b'A' { b'C' } else { b'A' };
            tried += 1;
            if sch
                .seeds(&m)
                .level(1)
                .iter()
                .any(|c| (c.label, c.start, c.end) == (s.label, s.start, s.end))
            {
                survived += 1;
            }
        }
        assert!(tried > 50, "{tried}");
        assert!(survived as f64 / tried as f64 > 0.3, "{survived}/{tried}");
    }
}
