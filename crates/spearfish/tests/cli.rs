//! End-to-end checks of the `spearfish` binary on the toy fixture.

use std::path::PathBuf;
use std::process::Command;

fn toy() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/toy/transcripts.fa")
}

fn lcp_stats(extra: &[&str]) -> Vec<Vec<String>> {
    let out = Command::new(env!("CARGO_BIN_EXE_spearfish"))
        .arg("lcp-stats")
        .arg(toy())
        .args(["-K", "4", "-t", "2"])
        .args(extra)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .map(|l| l.split('\t').map(String::from).collect())
        .collect()
}

#[test]
fn lcp_stats_on_toy() {
    for extra in [&[][..], &["--hpc"][..]] {
        let rows = lcp_stats(extra);
        assert_eq!(rows[0][0], "level");
        assert_eq!(rows.len(), 5, "header plus 4 levels");
        let col = |name: &str| rows[0].iter().position(|c| c == name).unwrap();
        let occ: Vec<u64> = rows[1..]
            .iter()
            .map(|r| r[col("occurrences")].parse().unwrap())
            .collect();
        let len: Vec<f64> = rows[1..]
            .iter()
            .map(|r| r[col("len_mean")].parse().unwrap())
            .collect();
        // Fewer, longer cores at each level.
        assert!(occ.windows(2).all(|w| w[1] < w[0]), "{occ:?}");
        assert!(len.windows(2).all(|w| w[1] > w[0]), "{len:?}");
    }
}

#[test]
fn hpc_lengths_are_raw() {
    // With HPC, core lengths are mapped back to raw coordinates, so they are at least as
    // long as the compressed cores; the level-1 mean must exceed its uncompressed value.
    let plain = lcp_stats(&[]);
    let hpc = lcp_stats(&["--hpc"]);
    let c = plain[0].iter().position(|c| c == "len_mean").unwrap();
    let a: f64 = plain[1][c].parse().unwrap();
    let b: f64 = hpc[1][c].parse().unwrap();
    assert!(b > a, "plain {a}, hpc {b}");
}

#[test]
fn locality_on_toy() {
    let out = Command::new(env!("CARGO_BIN_EXE_spearfish"))
        .arg("locality")
        .arg(toy())
        .args(["-K", "3", "-t", "2", "--edits", "5"])
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let text = String::from_utf8(out.stdout).unwrap();
    let rows: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows.len(), 1 + 3 * 3, "header plus 3 edit types x 3 levels");
    let col = |name: &str| rows[0].iter().position(|c| *c == name).unwrap();
    for r in &rows[1..] {
        // 6 toy transcripts x 5 edits.
        assert_eq!(r[col("edits")], "30");
        let changed: f64 = r[col("changed_mean")].parse().unwrap();
        assert!(changed > 0.0 && changed < 20.0, "{r:?}");
    }
}
