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
        let mut seq = Vec::with_capacity(raw.len());
        let mut raw_start = Vec::with_capacity(raw.len() + 1);
        let mut prev: Option<u8> = None;
        for (i, &b) in raw.iter().enumerate() {
            let u = b.to_ascii_uppercase();
            if prev != Some(u) {
                seq.push(b);
                raw_start.push(i as u32);
                prev = Some(u);
            }
        }
        raw_start.push(raw.len() as u32);
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
