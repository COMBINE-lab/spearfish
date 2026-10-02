//! Variable-length bit strings for deterministic coin tossing.
//!
//! Each item is a bit string stored LSB-first: bit 0 is the low bit of the *last* symbol,
//! as in lcptools, so "first difference from the right" is a trailing-zero count.

#[derive(Default)]
pub(crate) struct BitArena {
    words: Vec<u64>,
    /// Word offset of each item; `off[i]..off[i + 1]` are its words.
    off: Vec<usize>,
    len: Vec<u32>,
    cur_len: u32,
    /// Bits of the current item not yet flushed to `words` (the low `acc_bits` bits).
    acc: u64,
    acc_bits: u32,
}

#[derive(Clone, Copy)]
pub(crate) struct BitView<'a> {
    words: &'a [u64],
    len: u32,
}

impl BitView<'_> {
    #[inline]
    fn word(&self, i: usize) -> u64 {
        self.words.get(i).copied().unwrap_or(0)
    }

    #[inline]
    fn bit(&self, p: u32) -> u64 {
        if p >= self.len {
            0
        } else {
            (self.word((p / 64) as usize) >> (p % 64)) & 1
        }
    }
}

impl BitArena {
    pub(crate) fn with_capacity(items: usize) -> Self {
        let mut off = Vec::with_capacity(items + 1);
        off.push(0);
        BitArena {
            words: Vec::with_capacity(items),
            off,
            len: Vec::with_capacity(items),
            cur_len: 0,
            acc: 0,
            acc_bits: 0,
        }
    }

    /// Empty the arena, keeping its allocations.
    pub(crate) fn clear(&mut self) {
        self.words.clear();
        self.off.clear();
        self.off.push(0);
        self.len.clear();
        self.cur_len = 0;
        self.acc = 0;
        self.acc_bits = 0;
    }

    pub(crate) fn len(&self) -> usize {
        self.len.len()
    }

    /// Append `width` low bits of `v` above the bits already pushed for the current item.
    #[inline]
    pub(crate) fn push_bits(&mut self, v: u64, width: u32) {
        debug_assert!(width <= 64 && (width == 64 || v >> width == 0));
        if width == 0 {
            return;
        }
        let free = 64 - self.acc_bits;
        if width < free {
            self.acc |= v << self.acc_bits;
            self.acc_bits += width;
        } else {
            // Fill the accumulator, flush it, keep the rest.
            self.words.push(self.acc | (v << self.acc_bits));
            self.acc = if free < 64 { v >> free } else { 0 };
            self.acc_bits = width - free;
            if self.acc_bits == 0 {
                self.acc = 0;
            }
        }
        self.cur_len += width;
    }

    /// Close the current item.
    #[inline]
    pub(crate) fn finish_item(&mut self) {
        if self.acc_bits > 0 {
            self.words.push(self.acc);
        }
        self.acc = 0;
        self.acc_bits = 0;
        self.off.push(self.words.len());
        self.len.push(self.cur_len);
        self.cur_len = 0;
    }

    #[inline]
    pub(crate) fn get(&self, i: usize) -> BitView<'_> {
        BitView {
            words: &self.words[self.off[i]..self.off[i + 1]],
            len: self.len[i],
        }
    }
}

/// Width of a DCT value when it is concatenated into the next level's bit strings
/// (lcptools: bit length, at least 2).
#[inline]
pub(crate) fn value_width(v: u64) -> u32 {
    (64 - v.leading_zeros()).max(2)
}

/// One DCT step: compress `right` relative to its left neighbour. Returns `2p + b`, where
/// `p` is the index of the first differing bit from the right (capped at the shorter length)
/// and `b` is `right`'s bit at `p` (0 past its end).
#[inline]
pub(crate) fn dct_value(left: BitView<'_>, right: BitView<'_>) -> u64 {
    let m = left.len.min(right.len);
    let mut p = m;
    for w in 0..m.div_ceil(64) as usize {
        let x = left.word(w) ^ right.word(w);
        if x != 0 {
            p = (64 * w as u32 + x.trailing_zeros()).min(m);
            break;
        }
    }
    2 * p as u64 + right.bit(p)
}

/// [`dct_value`] for two single-word bit strings.
#[inline]
pub(crate) fn dct_value_words(left: u64, left_len: u32, right: u64, right_len: u32) -> u64 {
    let w = [left];
    let r = [right];
    dct_value(
        BitView {
            words: &w,
            len: left_len,
        },
        BitView {
            words: &r,
            len: right_len,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn arena(items: &[&[(u64, u32)]]) -> BitArena {
        let mut a = BitArena::with_capacity(items.len());
        for it in items {
            for &(v, w) in it.iter() {
                a.push_bits(v, w);
            }
            a.finish_item();
        }
        a
    }

    #[test]
    fn push_across_word_boundary() {
        let a = arena(&[&[(0x3, 2); 40]]);
        let v = a.get(0);
        assert_eq!(v.len, 80);
        assert_eq!(v.words, &[u64::MAX, 0xffff]);
    }

    #[test]
    fn dct_first_difference() {
        // left = ...0b1010, right = ...0b0110: lowest differing bit is bit 2, right has 1 there.
        let a = arena(&[&[(0b1010, 4)], &[(0b0110, 4)]]);
        assert_eq!(dct_value(a.get(0), a.get(1)), 2 * 2 + 1);
        assert_eq!(dct_value_words(0b1010, 4, 0b0110, 4), 5);
    }

    #[test]
    fn dct_suffix_match_uses_next_bit_of_right() {
        // right = 0b1_01 (3 bits), left = 0b01 (2 bits): all 2 shared bits match, p = 2,
        // right's bit 2 is 1.
        assert_eq!(dct_value_words(0b01, 2, 0b101, 3), 5);
        // right shorter: the bit past its end is 0.
        assert_eq!(dct_value_words(0b101, 3, 0b01, 2), 4);
    }

    #[test]
    fn dct_difference_in_second_word() {
        let a = arena(&[&[(0, 64), (0b1, 2)], &[(0, 64), (0b0, 2)]]);
        assert_eq!(dct_value(a.get(0), a.get(1)), 2 * 64);
    }
}
