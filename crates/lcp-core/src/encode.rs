//! 2-bit nucleotide encoding and helpers.
//!
//! The code order A=0, C=1, G=2, T=3 matches lcptools' default alphabet, which matters
//! because the parsing rules compare codes (local minima/maxima).

/// Code for each ASCII byte, or `INVALID`.
pub const INVALID: u8 = 0xff;

const fn build_table() -> [u8; 256] {
    let mut t = [INVALID; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    // RNA reads (direct RNA basecalls use U).
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;
    t
}

pub static CODE: [u8; 256] = build_table();

/// Reverse complement of an ASCII nucleotide sequence. Non-ACGT bytes map to `N`.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match CODE[b as usize] {
            0 => b'T',
            1 => b'G',
            2 => b'C',
            3 => b'A',
            _ => b'N',
        })
        .collect()
}

/// Maximal runs of ACGT characters, as `(offset, codes)`.
pub fn acgt_segments(seq: &[u8]) -> impl Iterator<Item = (usize, Vec<u8>)> + '_ {
    let mut i = 0;
    std::iter::from_fn(move || {
        while i < seq.len() && CODE[seq[i] as usize] == INVALID {
            i += 1;
        }
        if i >= seq.len() {
            return None;
        }
        let start = i;
        let mut codes = Vec::new();
        while i < seq.len() {
            let c = CODE[seq[i] as usize];
            if c == INVALID {
                break;
            }
            codes.push(c);
            i += 1;
        }
        Some((start, codes))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn segments_split_on_non_acgt() {
        let segs: Vec<_> = acgt_segments(b"NNACgtNAC-T").collect();
        assert_eq!(segs, vec![(2, vec![0, 1, 2, 3]), (7, vec![0, 1]), (10, vec![3])]);
    }

    #[test]
    fn revcomp_basic() {
        assert_eq!(revcomp(b"AACGTN"), b"NACGTT");
    }
}
