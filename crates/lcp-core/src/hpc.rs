//! Homopolymer compression (HPC): collapse each run of one base to a single base.
//!
//! Applied identically to transcripts and reads before parsing, so all internal coordinates
//! are HPC coordinates. [`Hpc::to_raw`] maps spans back to raw coordinates for output.
//! Runs are compared case-insensitively; the first base of each run is kept as is.

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Hpc {
    /// The compressed sequence.
    pub seq: Vec<u8>,
    /// `raw_start[i]` is the raw position where compressed base `i`'s run starts;
    /// `raw_start[seq.len()]` is the raw length.
    raw_start: Vec<u32>,
}

impl Hpc {
    pub fn new(raw: &[u8]) -> Self {
        // Branchless: every base is written to the next free slot, which is kept only when the
        // base starts a run (runs compare case-insensitively; a run keeps its first base).
        let n = raw.len();
        let mut seq = vec![0u8; n];
        let mut raw_start = vec![0u32; n + 1];
        let mut k = 0usize;
        let mut prev = u16::MAX; // matches no base
        for (i, &b) in raw.iter().enumerate() {
            let u = b.to_ascii_uppercase() as u16;
            // SAFETY: k <= i (k starts at 0 and grows by at most 1 per base), i < n = seq.len(),
            // and raw_start.len() = n + 1.
            unsafe {
                *seq.get_unchecked_mut(k) = b;
                *raw_start.get_unchecked_mut(k) = i as u32;
            }
            k += (u != prev) as usize;
            prev = u;
        }
        seq.truncate(k);
        raw_start[k] = n as u32;
        raw_start.truncate(k + 1);
        Hpc { seq, raw_start }
    }

    /// `Hpc::new(&revcomp(raw))`, from this compression of `raw`, provided `raw` is all
    /// A/C/G/T (either case): runs map to runs in reverse order, each starting with the
    /// complement (upper case, as `revcomp` writes it) of its base. Other symbols all become
    /// `N` under `revcomp` and can merge runs, so callers must check.
    pub fn revcomp_acgt(&self) -> Self {
        let n = self.raw_len();
        let seq = self
            .seq
            .iter()
            .rev()
            .map(|&b| b"TGCA"[(crate::encode::CODE[b as usize] & 3) as usize])
            .collect();
        let raw_start = self.raw_start.iter().rev().map(|&r| n - r).collect();
        Hpc { seq, raw_start }
    }

    pub fn raw_len(&self) -> u32 {
        *self.raw_start.last().unwrap()
    }

    /// Raw span `[s, e)` covered by the compressed span `[start, end)`: from the start of
    /// the first run to the end of the last.
    #[inline]
    pub fn to_raw(&self, start: u32, end: u32) -> (u32, u32) {
        (self.raw_start[start as usize], self.raw_start[end as usize])
    }

    /// Length of the run that compressed base `i` stands for.
    #[inline]
    pub fn run_len(&self, i: u32) -> u32 {
        self.raw_start[i as usize + 1] - self.raw_start[i as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        /// The mirrored compression equals compressing the reverse complement.
        #[test]
        fn revcomp_mirrors(raw in "[ACGTacgt]{0,200}") {
            let raw = raw.as_bytes();
            let a = Hpc::new(raw).revcomp_acgt();
            let b = Hpc::new(&crate::encode::revcomp(raw));
            prop_assert_eq!(&a.seq, &b.seq);
            prop_assert_eq!(&a.raw_start, &b.raw_start);
        }
    }

    #[test]
    fn compresses_runs() {
        let h = Hpc::new(b"AAACGGGTtA");
        assert_eq!(h.seq, b"ACGTA");
        assert_eq!(h.raw_len(), 10);
        assert_eq!(h.to_raw(1, 3), (3, 7)); // "CGGG"
        assert_eq!(h.to_raw(3, 5), (7, 10)); // "TtA"
        assert_eq!(h.run_len(3), 2);
    }

    #[test]
    fn empty() {
        let h = Hpc::new(b"");
        assert!(h.seq.is_empty());
        assert_eq!(h.raw_len(), 0);
        assert_eq!(h.to_raw(0, 0), (0, 0));
    }

    proptest! {
        /// Expanding each compressed base by its run length gives back the raw sequence
        /// (case-insensitively), and no two adjacent compressed bases are equal.
        #[test]
        fn roundtrip(raw in prop::collection::vec(prop::sample::select(b"AACGTTN".to_vec()), 0..500)) {
            let h = Hpc::new(&raw);
            let mut back = Vec::new();
            for (i, &b) in h.seq.iter().enumerate() {
                back.extend(std::iter::repeat_n(b, h.run_len(i as u32) as usize));
            }
            prop_assert_eq!(back, raw);
            prop_assert!(h.seq.windows(2).all(|w| w[0] != w[1]));
        }
    }
}
