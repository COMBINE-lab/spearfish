//! Edit distance of a pattern against the best-matching substring of a text (pattern fully
//! aligned, text ends free), with Myers' bit-parallel algorithm in multi-word blocks (as in
//! edlib). Inputs are 2-bit codes (0..4). Cost O(⌈m/64⌉ · n) word operations.

/// Minimum edit distance between `pattern` and any substring of `text`.
pub fn infix_distance(pattern: &[u8], text: &[u8]) -> u32 {
    let m = pattern.len();
    if m == 0 {
        return 0;
    }
    let blocks = m.div_ceil(64);
    // Peq[c][b]: bit i set where pattern[64b + i] == c.
    let mut peq = vec![[0u64; 4]; blocks];
    for (i, &c) in pattern.iter().enumerate() {
        peq[i / 64][(c & 3) as usize] |= 1 << (i % 64);
    }
    let last_bit = (m - 1) % 64;
    let mut pv = vec![u64::MAX; blocks];
    let mut mv = vec![0u64; blocks];
    // D[m][0] = m; D[0][j] = 0 (free start in the text).
    let mut score = m as i64;
    let mut best = score;
    for &c in text {
        let c = (c & 3) as usize;
        let mut hin: i64 = 0;
        for b in 0..blocks {
            let mut eq = peq[b][c];
            let (p, mm) = (pv[b], mv[b]);
            let xv = eq | mm;
            if hin < 0 {
                eq |= 1;
            }
            let xh = ((eq & p).wrapping_add(p) ^ p) | eq;
            let mut ph = mm | !(xh | p);
            let mut mh = p & xh;
            let bit = if b + 1 == blocks { last_bit } else { 63 };
            let hout = ((ph >> bit) & 1) as i64 - ((mh >> bit) & 1) as i64;
            ph <<= 1;
            mh <<= 1;
            if hin < 0 {
                mh |= 1;
            } else if hin > 0 {
                ph |= 1;
            }
            pv[b] = mh | !(xv | ph);
            mv[b] = ph & xv;
            hin = hout;
        }
        score += hin;
        best = best.min(score);
    }
    best as u32
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn naive(p: &[u8], t: &[u8]) -> u32 {
        let (m, n) = (p.len(), t.len());
        let mut prev: Vec<u32> = (0..=m as u32).collect();
        let mut best = prev[m];
        for j in 1..=n {
            let mut cur = vec![0u32; m + 1]; // D[0][j] = 0
            for i in 1..=m {
                let sub = prev[i - 1] + (p[i - 1] != t[j - 1]) as u32;
                cur[i] = sub.min(prev[i] + 1).min(cur[i - 1] + 1);
            }
            best = best.min(cur[m]);
            prev = cur;
        }
        best
    }

    #[test]
    fn examples() {
        assert_eq!(infix_distance(&[0, 1, 2, 3], &[3, 3, 0, 1, 2, 3, 3]), 0);
        assert_eq!(infix_distance(&[0, 1, 2, 3], &[0, 1, 3]), 1);
        assert_eq!(infix_distance(&[0, 1, 2, 3], &[]), 4);
        assert_eq!(infix_distance(&[], &[0, 1]), 0);
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(500))]

        #[test]
        fn matches_naive_dp(
            p in prop::collection::vec(0u8..4, 0..300),
            t in prop::collection::vec(0u8..4, 0..300),
        ) {
            prop_assert_eq!(infix_distance(&p, &t), naive(&p, &t));
        }

        /// Related sequences (the realistic case): a mutated copy of a substring.
        #[test]
        fn matches_naive_on_related(
            t in prop::collection::vec(0u8..4, 100..400),
            edits in prop::collection::vec((0usize..400, 0u8..4, 0u8..3), 0..20),
            a in 0usize..50, len in 60usize..300,
        ) {
            let end = (a + len).min(t.len());
            let mut p: Vec<u8> = t[a.min(end)..end].to_vec();
            for (pos, base, kind) in edits {
                if p.is_empty() { break; }
                let i = pos % p.len();
                match kind { 0 => p[i] = base, 1 => p.insert(i, base), _ => { p.remove(i); } }
            }
            prop_assert_eq!(infix_distance(&p, &t), naive(&p, &t));
        }
    }
}
