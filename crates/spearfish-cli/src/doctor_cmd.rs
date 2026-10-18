//! `spearfish doctor`: report the build and the machine, run a self-test, and optionally check
//! an index.
//!
//! The self-test runs at the SIMD level in use (`SPEARFISH_SIMD=scalar|avx2|avx512` lowers
//! it). It builds an index of the toy transcripts in memory and maps reads made from
//! them (forward, reverse-complemented, and with substitutions) with every preset. The
//! output's digest must match the reference: the SIMD kernels this machine selects compute
//! exactly what the scalar code does, so any difference is a bug.

use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Result, bail};
use spearfish::provenance::BuildInfo;
use spearfish::{Index, IndexOptions, MapSettings, Preset};

/// The toy transcripts (`fixtures/toy/transcripts.fa`).
const TOY: &str = include_str!("../../../fixtures/toy/transcripts.fa");

/// Digest of the self-test's output, per preset (see [`self_test`]).
pub const EXPECTED: [(Preset, &str); 3] = [
    (Preset::Ont, "5612907c8a21abe5"),
    (Preset::Drna, "4328f84a7c92686b"),
    (Preset::Hifi, "61ac6d16cfab9738"),
];

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Also check that this index loads and suits this build.
    #[arg(long)]
    index: Option<PathBuf>,
    /// Print the report as JSON.
    #[arg(long)]
    json: bool,
}

fn toy_transcripts() -> Vec<(String, String, Vec<u8>)> {
    let mut txs: Vec<(String, String, Vec<u8>)> = Vec::new();
    for line in TOY.lines() {
        if let Some(id) = line.strip_prefix('>') {
            let (name, gene) = spearfish::name_and_gene(id);
            txs.push((name, gene, Vec::new()));
        } else if let Some(t) = txs.last_mut() {
            t.2.extend_from_slice(line.trim().as_bytes());
        }
    }
    txs
}

/// Reads from each transcript: as is, reverse-complemented, and with a substitution every
/// 29 bases (deterministic).
fn toy_reads(txs: &[(String, String, Vec<u8>)]) -> Vec<Vec<u8>> {
    let mut reads = Vec::new();
    for (_, _, seq) in txs {
        reads.push(seq.clone());
        reads.push(lcp_core::encode::revcomp(seq));
        let mut m = seq.clone();
        for (i, b) in m.iter_mut().enumerate().skip(13).step_by(29) {
            *b = match (*b, i % 3) {
                (b'A', 0) | (b'G', 1) => b'C',
                (b'A', _) => b'G',
                (b'C', _) => b'T',
                (b'G', _) => b'A',
                _ => b'A',
            };
        }
        reads.push(m);
    }
    reads
}

/// Map the toy reads with `preset`; returns the xxh3-64 of every hit's fields.
pub fn self_test(preset: Preset) -> Result<String> {
    let txs = toy_transcripts();
    let reads = toy_reads(&txs);
    let (index, _) = Index::build(txs, &IndexOptions::default())?;
    let mapper = MapSettings::builder(preset).mapper(Arc::new(index))?;
    let mut h = xxhash_rust::xxh3::Xxh3::new();
    for r in &reads {
        let m = mapper.map(r);
        h.update(format!("{} {} {}\n", m.len, m.low_coverage, m.rescued).as_bytes());
        for x in &m.hits {
            h.update(
                format!(
                    "{} {} {:.4} {} {} {} {} {} {} {:.4} {}\n",
                    x.tx,
                    x.reverse,
                    x.score,
                    x.query_start,
                    x.query_end,
                    x.target_start,
                    x.target_end,
                    x.anchors,
                    x.anchored,
                    x.max_discrepancy,
                    x.extended
                )
                .as_bytes(),
            );
        }
    }
    Ok(format!("{:016x}", h.digest()))
}

pub fn run(args: Args) -> Result<()> {
    let build = BuildInfo::current();
    let mut report = serde_json::json!({
        "spearfish": build,
        "os": std::env::consts::OS,
        "arch": std::env::consts::ARCH,
        "cpus_available": crate::runinfo::available_cpus(),
        "simd": spearfish::provenance::simd_level(),
        "simd_detected": spearfish::provenance::simd_detected(),
        "score_denom": spearfish::SCORE_DENOM,
    });
    let mut ok = true;
    let mut tests = serde_json::Map::new();
    for (p, want) in EXPECTED {
        let got = self_test(p)?;
        let pass = got == want;
        ok &= pass;
        tests.insert(
            p.name().into(),
            serde_json::json!({"digest": got, "expected": want, "pass": pass}),
        );
    }
    report["self_test"] = tests.into();
    if let Some(path) = &args.index {
        let r = (|| -> Result<serde_json::Value> {
            let idx = Index::open(path)?;
            idx.raw().check_compatible()?;
            Ok(serde_json::json!({
                "path": path.display().to_string(),
                "index_id": idx.id(),
                "transcripts": idx.n_transcripts(),
                "raw_sequence": idx.has_raw_sequence(),
                "built_by": idx.manifest().map(|m| m.built_by.short()),
                "compatible": true,
            }))
        })();
        report["index"] = match r {
            Ok(v) => v,
            Err(e) => {
                ok = false;
                serde_json::json!({"path": path.display().to_string(), "compatible": false, "error": e.to_string()})
            }
        };
    }
    if args.json {
        println!("{}", serde_json::to_string_pretty(&report)?);
    } else {
        println!("spearfish {}", build.short());
        println!("target      {} ({})", build.target, build.profile);
        println!("cpus        {}", report["cpus_available"]);
        println!(
            "simd        {} (CPU supports {})",
            report["simd"].as_str().unwrap_or(""),
            report["simd_detected"].as_str().unwrap_or("")
        );
        for (p, t) in report["self_test"].as_object().unwrap() {
            let pass = t["pass"].as_bool() == Some(true);
            println!(
                "self-test   {p:5} {}",
                if pass {
                    "ok".to_string()
                } else {
                    format!(
                        "FAILED (digest {}, expected {})",
                        t["digest"].as_str().unwrap_or(""),
                        t["expected"].as_str().unwrap_or("")
                    )
                }
            );
        }
        if let Some(i) = report.get("index") {
            match i["compatible"].as_bool() {
                Some(true) => println!(
                    "index       ok: {} transcripts, id {}, built by {}",
                    i["transcripts"],
                    i["index_id"].as_str().unwrap_or(""),
                    i["built_by"]
                        .as_str()
                        .unwrap_or("an older spearfish (no manifest)")
                ),
                _ => println!("index       FAILED: {}", i["error"].as_str().unwrap_or("")),
            }
        }
    }
    if !ok {
        bail!("doctor found problems (see above)");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The reference digests: update them (from `spearfish doctor --json`) only together
    /// with a deliberate change to mapping output.
    #[test]
    fn self_test_matches_reference() {
        for (p, want) in EXPECTED {
            assert_eq!(self_test(p).unwrap(), want, "{p:?}");
        }
    }
}
