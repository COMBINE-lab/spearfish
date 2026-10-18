//! The dispatcher honours `SPEARFISH_SIMD`: each case re-runs this test binary with the
//! variable set, and the child checks the level `active()` resolved (it is fixed per process).

use std::process::Command;

use lcp_core::simd::{self, Level};

/// Runs in the child (with `SPEARFISH_SIMD_EXPECT` set): `active()` is the expected level.
#[test]
fn child_reports_active() {
    if let Ok(want) = std::env::var("SPEARFISH_SIMD_EXPECT") {
        assert_eq!(simd::active().name(), want);
    }
}

/// `active()` in a child process with `SPEARFISH_SIMD` = `var`, checked against `want`; returns
/// the child's stderr.
fn child(var: &str, want: Level) -> String {
    let out = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "child_reports_active",
            "--nocapture",
            "--test-threads=1",
        ])
        .env("SPEARFISH_SIMD", var)
        .env("SPEARFISH_SIMD_EXPECT", want.name())
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&out.stderr).into_owned();
    assert!(
        out.status.success(),
        "SPEARFISH_SIMD={var} should give {want}:\n{}{err}",
        String::from_utf8_lossy(&out.stdout)
    );
    err
}

#[test]
fn env_selects_each_supported_level() {
    for l in simd::supported() {
        let err = child(l.name(), l);
        assert!(!err.contains("spearfish:"), "unexpected warning for {l}: {err}");
        child(&l.name().to_uppercase(), l);
    }
}

#[test]
fn env_above_detected_warns_and_caps() {
    let d = simd::detected();
    if let Some(&above) = Level::ALL.iter().find(|&&l| !d.includes(l)) {
        let err = child(above.name(), d);
        assert!(err.contains("supports only"), "{err}");
    }
}

#[test]
fn env_invalid_warns_and_falls_back() {
    let err = child("turbo", simd::detected());
    assert!(err.contains("ignoring SPEARFISH_SIMD"), "{err}");
    // Empty means unset.
    let err = child("", simd::detected());
    assert!(!err.contains("spearfish:"), "{err}");
}
