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

#[test]
fn survival_on_toy() {
    let out = Command::new(env!("CARGO_BIN_EXE_spearfish"))
        .arg("survival")
        .arg(toy())
        .args([
            "-s",
            "lcp:K=3",
            "-s",
            "minimizer:k=15,w=10",
            "--rates",
            "0,0.05",
            "--profiles",
            "ont",
        ])
        .args(["--hpc", "off,on", "-t", "2"])
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let text = String::from_utf8(out.stdout).unwrap();
    let rows: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
    // 2 rates x 2 HPC settings x (3 LCP levels + 1 minimizer set).
    assert_eq!(rows.len(), 1 + 2 * 2 * 4);
    let col = |name: &str| rows[0].iter().position(|c| *c == name).unwrap();
    for r in &rows[1..] {
        let s: f64 = r[col("survival")].parse().unwrap();
        if r[col("rate")] == "0" {
            // No errors: every seed survives.
            assert_eq!(s, 1.0, "{r:?}");
        } else {
            assert!(s < 1.0, "{r:?}");
        }
    }
}

#[test]
fn halo_on_toy() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("halo_toy");
    std::fs::create_dir_all(&dir).unwrap();
    let prefix = dir.join("toy");
    let out = Command::new(env!("CARGO_BIN_EXE_spearfish"))
        .arg("halo")
        .arg(toy())
        .arg(toy().with_file_name("annotation.gtf"))
        .args([
            "-s",
            "lcp:K=3",
            "-s",
            "closed-syncmer:k=19,s=7",
            "--hpc",
            "off",
            "--rates",
            "0.02",
            "-t",
            "2",
            "-o",
        ])
        .arg(&prefix)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let read = |ext: &str| std::fs::read_to_string(format!("{}.{ext}.tsv", prefix.display())).unwrap();
    // Toy genes: GENEA (exon skip, 120 nt) and GENEB (alt 3' site, 36 nt) give 2 ordered pairs
    // each, in the 101-300 and 31-100 bins; both are resolved error-free by every seed type.
    let iso = read("iso");
    let resolved: Vec<&str> = iso
        .lines()
        .filter(|l| l.contains("\t0\toff\t") && (l.contains("\t31-100\t") || l.contains("\t101-300\t")))
        .collect();
    assert_eq!(resolved.len(), 2 * 2, "{iso}");
    for l in resolved {
        let f: Vec<&str> = l.split('\t').collect();
        assert_eq!((f[6], f[7]), ("2", "1.0000"), "{l}");
    }
    // genes.tsv: GENEA has 1270 distinct exonic bases and 4 junctions.
    assert!(
        read("genes")
            .lines()
            .any(|l| l.starts_with("TOYG0000001\toff\tlcp:K=3\t3\t2\t2420\t1270\t4\t"))
    );
    assert!(read("halo").lines().count() > 1);
}
