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
