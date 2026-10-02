//! A simple long-read error simulator with an exact original-to-read coordinate map.
//!
//! Each base is deleted, substituted or kept, then an insertion may follow it. Indel
//! probabilities are raised inside homopolymer runs (length >= 3) by `hp_bias`, and scaled
//! down elsewhere so the mean per-base rate stays as requested. Insertions inside runs
//! extend the run, as ONT and HiFi errors typically do.

use crate::rng::SplitMix64;

#[derive(Clone, Copy, Debug)]
pub struct ErrorProfile {
    pub name: &'static str,
    /// Fractions of errors that are substitutions, insertions and deletions (sum to 1).
    pub sub: f64,
    pub ins: f64,
    pub del: f64,
    /// Relative indel rate inside homopolymer runs.
    pub hp_bias: f64,
}

pub const PROFILES: [ErrorProfile; 3] = [
    ErrorProfile {
        name: "uniform",
        sub: 1.0 / 3.0,
        ins: 1.0 / 3.0,
        del: 1.0 / 3.0,
        hp_bias: 1.0,
    },
    // ONT cDNA-like: deletions dominate indels; indels concentrate in homopolymers.
    ErrorProfile {
        name: "ont",
        sub: 0.40,
        ins: 0.25,
        del: 0.35,
        hp_bias: 4.0,
    },
    // HiFi-like: few substitutions; most errors are homopolymer-length indels.
    ErrorProfile {
        name: "hifi",
        sub: 0.20,
        ins: 0.40,
        del: 0.40,
        hp_bias: 10.0,
    },
];

pub fn profile(name: &str) -> Option<ErrorProfile> {
    PROFILES.iter().copied().find(|p| p.name == name)
}

pub struct SimRead {
    pub seq: Vec<u8>,
    /// `o2r[i]`: read position of original base `i`, or of the next surviving base if `i`
    /// was deleted; `o2r[len]` is the read length.
    pub o2r: Vec<u32>,
    pub edits: u32,
}

pub fn simulate(orig: &[u8], rate: f64, p: &ErrorProfile, rng: &mut SplitMix64) -> SimRead {
    let n = orig.len();
    // Homopolymer membership: base i lies in a run of length >= 3.
    let mut hp = vec![false; n];
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && orig[j].eq_ignore_ascii_case(&orig[i]) {
            j += 1;
        }
        if j - i >= 3 {
            hp[i..j].fill(true);
        }
        i = j;
    }
    let f_hp = hp.iter().filter(|&&h| h).count() as f64 / n.max(1) as f64;
    let scale = 1.0 / (1.0 - f_hp + p.hp_bias * f_hp);
    let (p_sub, p_ins, p_del) = (rate * p.sub, rate * p.ins, rate * p.del);

    let mut seq = Vec::with_capacity(n + n / 8);
    let mut o2r = Vec::with_capacity(n + 1);
    let mut edits = 0;
    for i in 0..n {
        o2r.push(seq.len() as u32);
        let bias = scale * if hp[i] { p.hp_bias } else { 1.0 };
        let r = rng.unit();
        if r < p_del * bias {
            edits += 1;
        } else if r < p_del * bias + p_sub {
            seq.push(rng.other_base(orig[i]));
            edits += 1;
        } else {
            seq.push(orig[i]);
        }
        if rng.unit() < p_ins * bias {
            seq.push(if hp[i] { orig[i] } else { rng.base() });
            edits += 1;
        }
    }
    o2r.push(seq.len() as u32);
    SimRead { seq, o2r, edits }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rate_and_map() {
        let mut rng = SplitMix64::for_record(5, 0);
        let orig: Vec<u8> = (0..200_000).map(|_| rng.base()).collect();
        for prof in PROFILES {
            for rate in [0.0, 0.05, 0.15] {
                let r = simulate(&orig, rate, &prof, &mut rng);
                let realized = r.edits as f64 / orig.len() as f64;
                assert!((realized - rate).abs() < 0.01, "{} {rate}: {realized}", prof.name);
                assert_eq!(r.o2r.len(), orig.len() + 1);
                assert_eq!(*r.o2r.last().unwrap() as usize, r.seq.len());
                assert!(r.o2r.windows(2).all(|w| w[0] <= w[1]));
                if rate == 0.0 {
                    assert_eq!(r.seq, orig);
                }
            }
        }
    }
}
