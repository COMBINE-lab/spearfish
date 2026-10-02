//! A mergeable integer histogram with an overflow bucket.

#[derive(Clone, Debug)]
pub struct Hist {
    /// `counts[v]` for `v < cap`; `counts[cap]` collects everything `>= cap`.
    counts: Vec<u64>,
    sum: u128,
    n: u64,
    max: u64,
}

impl Hist {
    pub fn new(cap: usize) -> Self {
        Hist {
            counts: vec![0; cap + 1],
            sum: 0,
            n: 0,
            max: 0,
        }
    }

    #[inline]
    pub fn add(&mut self, v: u64) {
        let cap = self.counts.len() - 1;
        self.counts[(v as usize).min(cap)] += 1;
        self.sum += v as u128;
        self.n += 1;
        self.max = self.max.max(v);
    }

    pub fn merge(&mut self, o: &Hist) {
        assert_eq!(self.counts.len(), o.counts.len());
        for (a, b) in self.counts.iter_mut().zip(&o.counts) {
            *a += b;
        }
        self.sum += o.sum;
        self.n += o.n;
        self.max = self.max.max(o.max);
    }

    pub fn n(&self) -> u64 {
        self.n
    }

    pub fn mean(&self) -> f64 {
        if self.n == 0 {
            0.0
        } else {
            self.sum as f64 / self.n as f64
        }
    }

    pub fn max(&self) -> u64 {
        self.max
    }

    /// Smallest value `v` with at least a fraction `q` of observations `<= v`. Values in the
    /// overflow bucket report as the cap.
    pub fn quantile(&self, q: f64) -> u64 {
        if self.n == 0 {
            return 0;
        }
        let target = ((q * self.n as f64).ceil() as u64).max(1);
        let mut acc = 0;
        for (v, &c) in self.counts.iter().enumerate() {
            acc += c;
            if acc >= target {
                return v as u64;
            }
        }
        (self.counts.len() - 1) as u64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantiles_and_merge() {
        let mut a = Hist::new(100);
        for v in 1..=10 {
            a.add(v);
        }
        let mut b = Hist::new(100);
        b.add(1000);
        a.merge(&b);
        assert_eq!(a.n(), 11);
        assert_eq!(a.quantile(0.5), 6);
        assert_eq!(a.quantile(1.0), 100);
        assert_eq!(a.max(), 1000);
        assert!((a.mean() - 1055.0 / 11.0).abs() < 1e-9);
    }
}
