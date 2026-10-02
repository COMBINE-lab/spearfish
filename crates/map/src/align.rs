//! Edit distance with Myers' bit-parallel algorithm in multi-word blocks and Ukkonen's
//! cutoff (as in edlib). Inputs are 2-bit codes (0..4). Three modes:
//!
//! - `Infix`: pattern fully aligned, both text ends free (best substring of the text).
//! - `Prefix`: pattern fully aligned starting at text position 0, text end free.
//! - `Global`: pattern and text fully aligned.
//!
//! With a threshold `k`, only blocks of pattern rows that can still score ≤ k are computed
//! (the band), and the result is `None` if the distance exceeds `k`. Cost
//! O(⌈min(m, k)/64⌉ · n) word operations in practice, O(⌈m/64⌉ · n) at worst.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Infix,
    Prefix,
    Global,
}

const W: usize = 64;

/// One Myers block step; returns the horizontal delta out of row `bit`.
#[inline(always)]
fn block(pv: &mut u64, mv: &mut u64, mut eq: u64, hin: i32, bit: usize) -> i32 {
    let (p, m) = (*pv, *mv);
    let xv = eq | m;
    if hin < 0 {
        eq |= 1;
    }
    let xh = ((eq & p).wrapping_add(p) ^ p) | eq;
    let mut ph = m | !(xh | p);
    let mut mh = p & xh;
    let hout = ((ph >> bit) & 1) as i32 - ((mh >> bit) & 1) as i32;
    ph <<= 1;
    mh <<= 1;
    if hin < 0 {
        mh |= 1;
    } else if hin > 0 {
        ph |= 1;
    }
    *pv = mh | !(xv | ph);
    *mv = ph & xv;
    hout
}

/// Edit distance in `mode`, or `None` if it exceeds `k`.
pub fn distance(pattern: &[u8], text: &[u8], mode: Mode, k: u32) -> Option<u32> {
    let m = pattern.len();
    let n = text.len();
    if m == 0 {
        let d = if mode == Mode::Global { n as u32 } else { 0 };
        return (d <= k).then_some(d);
    }
    let k = k.min((m + n) as u32) as i32;
    let blocks = m.div_ceil(W);
    let last_bit = (m - 1) % W;
    let rows = |b: usize| if b + 1 == blocks { last_bit + 1 } else { W };
    let bit = |b: usize| if b + 1 == blocks { last_bit } else { W - 1 };
    let mut peq = vec![[0u64; 4]; blocks];
    for (i, &c) in pattern.iter().enumerate() {
        peq[i / W][(c & 3) as usize] |= 1 << (i % W);
    }
    // Horizontal delta entering row 0: D[0][j] is 0 (free text start) or j.
    let top: i32 = if mode == Mode::Infix { 0 } else { 1 };
    let mut pv = vec![u64::MAX; blocks];
    let mut mv = vec![0u64; blocks];
    // score[b] = D[bottom row of block b][current column]; column 0 is D[i][0] = i.
    let mut score: Vec<i32> = (0..blocks).map(|b| (b * W + rows(b)) as i32).collect();
    // Active blocks 0..=last: rows that can still be <= k (D[i][0] = i).
    let mut last: isize = (blocks.min((k as usize + 1).div_ceil(W)) as isize) - 1;
    let mut best: i32 = if last as usize + 1 == blocks && score[blocks - 1] <= k {
        score[blocks - 1]
    } else {
        i32::MAX
    };
    for (j, &c) in text.iter().enumerate() {
        let c = (c & 3) as usize;
        let mut hin = top;
        for b in 0..=last as usize {
            hin = block(&mut pv[b], &mut mv[b], peq[b][c], hin, bit(b));
            score[b] += hin;
        }
        let lb = last as usize;
        // Ukkonen: add the next block if its top can still be <= k, else drop exhausted
        // blocks from the bottom.
        if lb + 1 < blocks && score[lb] - hin <= k && ((peq[lb + 1][c] & 1) != 0 || hin < 0) {
            let nb = lb + 1;
            pv[nb] = u64::MAX;
            mv[nb] = 0;
            let h = block(&mut pv[nb], &mut mv[nb], peq[nb][c], hin, bit(nb));
            score[nb] = score[lb] - hin + rows(nb) as i32 + h;
            last += 1;
        } else {
            while last >= 0 && score[last as usize] >= k + W as i32 {
                last -= 1;
            }
            if last < 0 {
                // No row can be <= k from here on (only possible in the anchored modes):
                // Global needs the last column, so it fails; the others keep their best.
                return if mode == Mode::Global {
                    None
                } else {
                    (best <= k).then_some(best as u32)
                };
            }
        }
        if last as usize + 1 == blocks {
            let s = score[blocks - 1];
            match mode {
                Mode::Infix | Mode::Prefix => best = best.min(s),
                Mode::Global => {
                    if j + 1 == n {
                        best = s;
                    }
                }
            }
        } else if mode == Mode::Global && j + 1 == n {
            best = i32::MAX;
        }
    }
    if mode == Mode::Global && n == 0 {
        best = m as i32;
    }
    (best <= k).then_some(best as u32)
}

/// Minimum edit distance between `pattern` and any substring of `text` (no cutoff).
pub fn infix_distance(pattern: &[u8], text: &[u8]) -> u32 {
    distance(pattern, text, Mode::Infix, u32::MAX).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn naive(p: &[u8], t: &[u8], mode: Mode) -> u32 {
        let (m, n) = (p.len(), t.len());
        let mut prev: Vec<u32> = (0..=m as u32).collect();
        let mut best = if mode == Mode::Global && n > 0 {
            u32::MAX
        } else {
            prev[m]
        };
        for j in 1..=n {
            let mut cur = vec![0u32; m + 1];
            cur[0] = if mode == Mode::Infix { 0 } else { j as u32 };
            for i in 1..=m {
                let sub = prev[i - 1] + (p[i - 1] != t[j - 1]) as u32;
                cur[i] = sub.min(prev[i] + 1).min(cur[i - 1] + 1);
            }
            match mode {
                Mode::Global => {
                    if j == n {
                        best = cur[m]
                    }
                }
                _ => best = best.min(cur[m]),
            }
            prev = cur;
        }
        best
    }

    fn check(p: &[u8], t: &[u8], k: u32) -> Result<(), TestCaseError> {
        for mode in [Mode::Infix, Mode::Prefix, Mode::Global] {
            let want = naive(p, t, mode);
            let got = distance(p, t, mode, k);
            prop_assert_eq!(got, (want <= k).then_some(want), "mode {:?} k {}", mode, k);
        }
        Ok(())
    }

    #[test]
    fn examples() {
        assert_eq!(infix_distance(&[0, 1, 2, 3], &[3, 3, 0, 1, 2, 3, 3]), 0);
        assert_eq!(infix_distance(&[0, 1, 2, 3], &[0, 1, 3]), 1);
        assert_eq!(infix_distance(&[0, 1, 2, 3], &[]), 4);
        assert_eq!(
            distance(&[0, 1, 2, 3], &[3, 3, 0, 1, 2, 3, 3], Mode::Prefix, 10),
            Some(2)
        );
        assert_eq!(
            distance(&[0, 1, 2, 3], &[3, 3, 0, 1, 2, 3, 3], Mode::Global, 10),
            Some(3)
        );
        assert_eq!(
            distance(&[0, 1, 2, 3], &[3, 3, 0, 1, 2, 3, 3], Mode::Global, 2),
            None
        );
        assert_eq!(distance(&[], &[0, 1], Mode::Global, 5), Some(2));
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(600))]

        #[test]
        fn random_vs_naive(
            p in prop::collection::vec(0u8..4, 0..300),
            t in prop::collection::vec(0u8..4, 0..300),
            k in prop_oneof![Just(u32::MAX), 0u32..400],
        ) {
            check(&p, &t, k)?;
        }

        /// Related sequences (the realistic case): a mutated copy of a substring, with
        /// thresholds around the true distance.
        #[test]
        fn related_vs_naive(
            t in prop::collection::vec(0u8..4, 100..500),
            edits in prop::collection::vec((0usize..500, 0u8..4, 0u8..3), 0..30),
            a in 0usize..50, len in 60usize..450, slack in 0u32..20,
        ) {
            let end = (a + len).min(t.len());
            let mut p: Vec<u8> = t[a.min(end)..end].to_vec();
            for (pos, base, kind) in edits {
                if p.is_empty() { break; }
                let i = pos % p.len();
                match kind { 0 => p[i] = base, 1 => p.insert(i, base), _ => { p.remove(i); } }
            }
            let d = naive(&p, &t, Mode::Infix);
            for k in [d.saturating_sub(1), d, d + slack, u32::MAX] {
                check(&p, &t, k)?;
            }
        }
    }
}
