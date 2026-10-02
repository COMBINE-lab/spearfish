//! Small deterministic RNG for simulations: SplitMix64, seeded per record so results don't
//! depend on thread scheduling.

#[derive(Clone, Debug)]
pub struct SplitMix64(u64);

impl SplitMix64 {
    /// Independent stream for record `index` under a global `seed`.
    pub fn for_record(seed: u64, index: u64) -> Self {
        let mut r = SplitMix64(seed ^ index.wrapping_mul(0xd1b5_4a32_d192_ed03));
        r.next_u64();
        r
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// Uniform in `0..n` (n > 0), by multiply-shift.
    pub fn below(&mut self, n: u64) -> u64 {
        ((self.next_u64() as u128 * n as u128) >> 64) as u64
    }

    /// Uniform in `[0, 1)`.
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A random base other than `b` (case-insensitive), as uppercase ASCII.
    pub fn other_base(&mut self, b: u8) -> u8 {
        let alts: Vec<u8> = b"ACGT"
            .iter()
            .copied()
            .filter(|&c| c != b.to_ascii_uppercase())
            .collect();
        alts[self.below(alts.len() as u64) as usize]
    }

    pub fn base(&mut self) -> u8 {
        b"ACGT"[self.below(4) as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_and_in_range() {
        let mut a = SplitMix64::for_record(1, 7);
        let mut b = SplitMix64::for_record(1, 7);
        for _ in 0..1000 {
            let (x, y) = (a.below(10), b.below(10));
            assert_eq!(x, y);
            assert!(x < 10);
            assert_ne!(a.other_base(b'a'), b'A');
            b.other_base(b'a');
        }
    }
}
