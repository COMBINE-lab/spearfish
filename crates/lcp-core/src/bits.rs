//! Variable-length bit strings for deterministic coin tossing.
//!
//! Each item is a bit string stored LSB-first: bit 0 is the low bit of the *last* symbol,
//! as in lcptools, so "first difference from the right" is a trailing-zero count.

/// Items are almost always at most 64 bits (level-1 cores under 32 bases; a few small DCT
/// values per core above), so each item's low word and length sit in flat arrays and only
/// longer items keep their further words aside.
#[derive(Default)]
pub(crate) struct BitArena {
    /// Low 64 bits of each item (0 if empty).
    head: Vec<u64>,
    len: Vec<u32>,
    /// Words past the first of items longer than 64 bits, and per such item (item, index of
    /// its second word in `tail`), in item order.
    tail: Vec<u64>,
    long: Vec<(u32, u32)>,
    cur_len: u32,
    /// The current item's first word once flushed.
    cur_head: Option<u64>,
    /// Bits of the current item not yet flushed (the low `acc_bits` bits).
    acc: u64,
    acc_bits: u32,
}

impl BitArena {
    pub(crate) fn with_capacity(items: usize) -> Self {
        BitArena {
            head: Vec::with_capacity(items),
            len: Vec::with_capacity(items),
            ..Default::default()
        }
    }

    /// Empty the arena, keeping its allocations.
    pub(crate) fn clear(&mut self) {
        self.head.clear();
        self.len.clear();
        self.tail.clear();
        self.long.clear();
        self.cur_len = 0;
        self.cur_head = None;
        self.acc = 0;
        self.acc_bits = 0;
    }

    pub(crate) fn len(&self) -> usize {
        self.len.len()
    }

    #[inline]
    fn flush(&mut self, w: u64) {
        match self.cur_head {
            None => {
                self.cur_head = Some(w);
                self.long.push((self.len.len() as u32, self.tail.len() as u32));
            }
            Some(_) => self.tail.push(w),
        }
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
            self.flush(self.acc | (v << self.acc_bits));
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
        match self.cur_head.take() {
            None => self.head.push(self.acc),
            Some(h) => {
                self.head.push(h);
                if self.acc_bits > 0 {
                    self.tail.push(self.acc);
                }
                // A flushed item of exactly 64 bits has no tail: not long after all.
                if self.cur_len <= 64 {
                    self.long.pop();
                }
            }
        }
        self.acc = 0;
        self.acc_bits = 0;
        self.len.push(self.cur_len);
        self.cur_len = 0;
    }

    /// A whole item of `width <= 64` bits (`v`'s low bits), with no item open.
    #[inline]
    pub(crate) fn push_item(&mut self, v: u64, width: u32) {
        debug_assert!(width <= 64 && self.cur_len == 0 && self.cur_head.is_none());
        self.head.push(v);
        self.len.push(width);
    }

    /// A whole item of DCT values, the last value lowest, each [`value_width`] bits wide.
    #[inline]
    pub(crate) fn push_values(&mut self, vals: &[u64]) {
        let (mut acc, mut nb) = (0u64, 0u32);
        for &v in vals.iter().rev() {
            let w = value_width(v);
            if nb + w > 64 {
                for &v in vals.iter().rev() {
                    self.push_bits(v, value_width(v));
                }
                self.finish_item();
                return;
            }
            acc |= v << nb;
            nb += w;
        }
        self.push_item(acc, nb);
    }

    /// Word `w` of item `i` (0 past its end).
    fn word(&self, i: usize, w: usize) -> u64 {
        if w == 0 {
            return self.head[i];
        }
        if (self.len[i] as usize) <= 64 * w {
            return 0;
        }
        let k = self.long.partition_point(|&(it, _)| (it as usize) < i);
        self.tail[self.long[k].1 as usize + w - 1]
    }

    /// One DCT step for item `i` (`i >= 1`) against item `i - 1`: [`dct_value`] of the two.
    #[cfg(test)]
    pub(crate) fn dct(&self, i: usize) -> u64 {
        self.dct_with(i, self.head[i - 1], self.head[i], self.len[i - 1], self.len[i])
    }

    /// The DCT step of every item `i >= 1` against item `i - 1`, into `out[i - 1]`
    /// (`out.len() == self.len() - 1`).
    pub(crate) fn dct_all(&self, out: &mut [u64]) {
        let pairs = self.head.windows(2).zip(self.len.windows(2));
        for (i, (v, (h, l))) in out.iter_mut().zip(pairs).enumerate() {
            *v = self.dct_with(i + 1, h[0], h[1], l[0], l[1]);
        }
    }

    /// The DCT step for item `i` given its and item `i - 1`'s low words and lengths.
    #[inline(always)]
    fn dct_with(&self, i: usize, lh: u64, rh: u64, ll: u32, rl: u32) -> u64 {
        let m = ll.min(rl);
        let x = lh ^ rh;
        let p = if x != 0 {
            x.trailing_zeros().min(m)
        } else if m <= 64 {
            m
        } else {
            return self.dct_slow(i);
        };
        if p >= 64 {
            return self.dct_slow(i);
        }
        2 * p as u64 + (((rh >> p) & 1) & (p < rl) as u64)
    }

    #[cold]
    fn dct_slow(&self, i: usize) -> u64 {
        let m = self.len[i - 1].min(self.len[i]);
        let mut p = m;
        for w in 0..m.div_ceil(64) as usize {
            let x = self.word(i - 1, w) ^ self.word(i, w);
            if x != 0 {
                p = (64 * w as u32 + x.trailing_zeros()).min(m);
                break;
            }
        }
        let b = if p < self.len[i] {
            (self.word(i, (p / 64) as usize) >> (p % 64)) & 1
        } else {
            0
        };
        2 * p as u64 + b
    }
}

/// Width of a DCT value when it is concatenated into the next level's bit strings
/// (lcptools: bit length, at least 2).
#[inline]
pub(crate) fn value_width(v: u64) -> u32 {
    (64 - v.leading_zeros()).max(2)
}

/// One DCT step on single-word bit strings: compress `right` relative to its left neighbour.
/// Returns `2p + b`, where `p` is the index of the first differing bit from the right (capped
/// at the shorter length) and `b` is `right`'s bit at `p` (0 past its end).
#[inline]
pub(crate) fn dct_value_words(left: u64, left_len: u32, right: u64, right_len: u32) -> u64 {
    let m = left_len.min(right_len);
    let x = left ^ right;
    let p = if x != 0 { x.trailing_zeros().min(m) } else { m };
    let b = if p < right_len && p < 64 {
        (right >> p) & 1
    } else {
        0
    };
    2 * p as u64 + b
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
        assert_eq!(a.len[0], 80);
        assert_eq!((a.word(0, 0), a.word(0, 1), a.word(0, 2)), (u64::MAX, 0xffff, 0));
    }

    #[test]
    fn dct_first_difference() {
        // left = ...0b1010, right = ...0b0110: lowest differing bit is bit 2, right has 1 there.
        let a = arena(&[&[(0b1010, 4)], &[(0b0110, 4)]]);
        assert_eq!(a.dct(1), 2 * 2 + 1);
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
        assert_eq!(a.dct(1), 2 * 64);
        // Long items among short ones keep their own words.
        let a = arena(&[
            &[(1, 2)],
            &[(0, 64), (0b10, 2)],
            &[(0, 64)],
            &[(0, 64), (0b11, 2)],
        ]);
        // Bit 0 differs (1 vs 0), and right's bit 0 is 0.
        assert_eq!(a.dct(1), 0);
        assert_eq!(a.dct(2), 2 * 64);
        // All 64 shared bits agree; right's bit 64 is 1.
        assert_eq!(a.dct(3), 2 * 64 + 1);
        assert_eq!(a.word(3, 1), 0b11);
        assert_eq!(a.word(1, 1), 0b10);
    }

    /// The flat layout gives what the bits spell out, item by item.
    #[test]
    fn dct_matches_bit_strings() {
        let mut x = 7u64;
        let mut rnd = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        for _ in 0..2000 {
            let items: Vec<Vec<(u64, u32)>> = (0..6)
                .map(|_| {
                    let k = (rnd() % 30) as usize;
                    // Few distinct values, so neighbours often share long suffixes.
                    (0..k)
                        .map(|_| (rnd() % 2, 1 + (rnd() % 2) as u32))
                        .map(|(v, w)| (v & ((1 << w) - 1), w))
                        .collect()
                })
                .collect();
            let refs: Vec<&[(u64, u32)]> = items.iter().map(|v| v.as_slice()).collect();
            let a = arena(&refs);
            let bits: Vec<Vec<u8>> = items
                .iter()
                .map(|it| {
                    it.iter()
                        .flat_map(|&(v, w)| (0..w).map(move |b| ((v >> b) & 1) as u8))
                        .collect()
                })
                .collect();
            let mut all = vec![0; items.len() - 1];
            a.dct_all(&mut all);
            for i in 1..items.len() {
                let (l, r) = (&bits[i - 1], &bits[i]);
                let m = l.len().min(r.len());
                let p = (0..m).find(|&k| l[k] != r[k]).unwrap_or(m);
                let want = 2 * p as u64 + r.get(p).copied().unwrap_or(0) as u64;
                assert_eq!(a.dct(i), want, "{:?}", items);
                assert_eq!(all[i - 1], want, "{:?}", items);
            }
        }
    }
}
