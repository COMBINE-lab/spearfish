//! SIMD instruction-set levels: what the CPU supports ([`detected`]) and what the kernels run
//! ([`active`]), shared by every crate with hand-written vector kernels.
//!
//! Kernels dispatch on [`active`] once per call (a parse, a DP, a batch), into an entry
//! function compiled for that level as a whole, so nothing is checked inside hot loops. All
//! levels give bit-identical results; `SPEARFISH_SIMD=scalar|avx2|avx512` lowers the level
//! (for testing and benchmarking), capped at what the CPU supports.

use std::sync::OnceLock;

/// An instruction-set level a kernel can be compiled for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Level {
    /// The target's baseline (SSE2 on x86-64, NEON on aarch64), autovectorized at most.
    Scalar,
    /// x86-64 AVX2 with FMA, BMI1/2, POPCNT and LZCNT (Haswell, Zen 1 and later).
    Avx2,
    /// [`Level::Avx2`] plus AVX-512 F/BW/VL (Skylake-SP, Ice Lake, Zen 4 and later).
    Avx512,
}

impl Level {
    /// Every level, lowest first.
    pub const ALL: [Level; 3] = [Level::Scalar, Level::Avx2, Level::Avx512];

    /// Lower-case name, as `SPEARFISH_SIMD` takes it.
    pub fn name(&self) -> &'static str {
        match self {
            Level::Scalar => "scalar",
            Level::Avx2 => "avx2",
            Level::Avx512 => "avx512",
        }
    }

    /// The level named `s` (case-insensitive), as [`Level::name`] spells it.
    pub fn from_name(s: &str) -> Option<Level> {
        Level::ALL
            .into_iter()
            .find(|l| l.name().eq_ignore_ascii_case(s.trim()))
    }

    /// Whether a CPU at this level can run kernels compiled for `other`.
    pub fn includes(self, other: Level) -> bool {
        matches!(
            (self, other),
            (_, Level::Scalar) | (Level::Avx2 | Level::Avx512, Level::Avx2) | (Level::Avx512, Level::Avx512)
        )
    }
}

impl std::fmt::Display for Level {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

/// The highest level this CPU supports (checked once).
pub fn detected() -> Level {
    static L: OnceLock<Level> = OnceLock::new();
    *L.get_or_init(detect)
}

#[cfg(target_arch = "x86_64")]
fn detect() -> Level {
    use std::arch::is_x86_feature_detected as has;
    let avx2 = has!("avx2") && has!("fma") && has!("bmi1") && has!("bmi2") && has!("popcnt") && has!("lzcnt");
    if avx2 && has!("avx512f") && has!("avx512bw") && has!("avx512vl") {
        Level::Avx512
    } else if avx2 {
        Level::Avx2
    } else {
        Level::Scalar
    }
}

#[cfg(not(target_arch = "x86_64"))]
fn detect() -> Level {
    Level::Scalar
}

/// Every level this CPU can run, lowest first (what tests and benchmarks iterate over).
pub fn supported() -> Vec<Level> {
    let d = detected();
    Level::ALL.into_iter().filter(|&l| d.includes(l)).collect()
}

/// The level the kernels run: [`detected`], lowered by `SPEARFISH_SIMD` if set. Resolved once
/// per process (a warning, if any, is printed then).
pub fn active() -> Level {
    static L: OnceLock<Level> = OnceLock::new();
    *L.get_or_init(|| {
        let var = std::env::var("SPEARFISH_SIMD").ok();
        let (l, warning) = resolve(var.as_deref(), detected());
        if let Some(w) = warning {
            eprintln!("spearfish: {w}");
        }
        l
    })
}

/// The level for `SPEARFISH_SIMD` = `var` on a CPU at `detected`, and a warning to print if
/// the request can't be honoured (more than the CPU supports, or not a level).
pub fn resolve(var: Option<&str>, detected: Level) -> (Level, Option<String>) {
    let Some(v) = var.filter(|v| !v.trim().is_empty()) else {
        return (detected, None);
    };
    match Level::from_name(v) {
        Some(l) if detected.includes(l) => (l, None),
        Some(l) => (
            detected,
            Some(format!(
                "SPEARFISH_SIMD={} requested, but this CPU supports only {}; using {}",
                l.name(),
                detected.name(),
                detected.name()
            )),
        ),
        None => (
            detected,
            Some(format!(
                "ignoring SPEARFISH_SIMD={v:?} (expected scalar, avx2 or avx512); using {}",
                detected.name()
            )),
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_round_trip() {
        for l in Level::ALL {
            assert_eq!(Level::from_name(l.name()), Some(l));
            assert_eq!(Level::from_name(&l.name().to_uppercase()), Some(l));
        }
        assert_eq!(Level::from_name("sse"), None);
    }

    #[test]
    fn resolve_caps_and_rejects() {
        use Level::*;
        assert_eq!(resolve(None, Avx512), (Avx512, None));
        assert_eq!(resolve(Some(""), Avx2), (Avx2, None));
        assert_eq!(resolve(Some("scalar"), Avx512), (Scalar, None));
        assert_eq!(resolve(Some("AVX2"), Avx512), (Avx2, None));
        assert_eq!(resolve(Some("avx512"), Avx512), (Avx512, None));
        // More than the CPU has: a warning, and the detected level.
        let (l, w) = resolve(Some("avx512"), Avx2);
        assert_eq!(l, Avx2);
        assert!(w.unwrap().contains("supports only avx2"));
        let (l, w) = resolve(Some("avx2"), Scalar);
        assert_eq!(l, Scalar);
        assert!(w.is_some());
        // Not a level: a warning, never a panic.
        let (l, w) = resolve(Some("fast"), Avx512);
        assert_eq!(l, Avx512);
        assert!(w.unwrap().contains("ignoring"));
    }

    #[test]
    fn levels_are_ordered_by_inclusion() {
        let s = supported();
        assert_eq!(s[0], Level::Scalar);
        assert_eq!(*s.last().unwrap(), detected());
        assert!(detected().includes(active()));
        for a in Level::ALL {
            assert!(a.includes(a) && a.includes(Level::Scalar));
        }
        assert!(!Level::Avx2.includes(Level::Avx512));
    }
}
