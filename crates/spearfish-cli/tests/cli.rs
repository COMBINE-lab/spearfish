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

#[test]
fn index_and_stats_on_toy() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    let idx = dir.join("toy_cli.sfi");
    let out = Command::new(env!("CARGO_BIN_EXE_spearfish"))
        .arg("index")
        .arg(toy())
        .args(["-s", "lcp:K=3", "-t", "2", "-o"])
        .arg(&idx)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let out = Command::new(env!("CARGO_BIN_EXE_spearfish"))
        .arg("index-stats")
        .arg(&idx)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let text = String::from_utf8(out.stdout).unwrap();
    let get = |k: &str| {
        text.lines()
            .find_map(|l| l.strip_prefix(&format!("{k}\t")))
            .unwrap()
            .to_string()
    };
    assert_eq!(get("spec"), "lcp:K=3");
    assert_eq!(get("hpc"), "true");
    assert_eq!(get("collapse"), "gene");
    assert_eq!(get("transcripts"), "6");
    assert_eq!(get("compatible"), "yes");
    assert!(get("built_by").starts_with("spearfish "));

    // The manifest is stored in the index and written next to it, and names its inputs.
    let side: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("toy_cli.sfi.json")).unwrap()).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_spearfish"))
        .args(["inspect", "--json", "--sections"])
        .arg(&idx)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let j: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(j["manifest"], side);
    assert_eq!(j["index_id"], side["index_id"]);
    assert_eq!(side["input"]["transcripts"], 6);
    assert_eq!(side["input"]["genes_from"], "fasta");
    assert_eq!(side["options"]["scheme"], "lcp:K=3");
    assert!(
        j["sections"]
            .as_array()
            .unwrap()
            .iter()
            .any(|s| s["name"] == "manifest")
    );
}

/// Commands around mapping: doctor, completions, man pages, parameter printing.
#[test]
fn utility_commands() {
    let bin = env!("CARGO_BIN_EXE_spearfish");
    let ok = |args: &[&str]| {
        let o = Command::new(bin).args(args).output().unwrap();
        assert!(
            o.status.success(),
            "{args:?}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        String::from_utf8(o.stdout).unwrap()
    };
    let d: serde_json::Value = serde_json::from_str(&ok(&["doctor", "--json"])).unwrap();
    for p in ["ont", "drna", "hifi", "ont-fast", "drna-fast", "hifi-fast"] {
        assert_eq!(d["self_test"][p]["pass"], true, "{d}");
    }
    assert!(ok(&["completions", "bash"]).contains("spearfish"));
    let man = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("man");
    ok(&["man", "-q", "-o", man.to_str().unwrap()]);
    for page in ["spearfish.1", "spearfish-map.1", "spearfish-index.1"] {
        assert!(man.join(page).exists(), "{page}");
    }
    let p: serde_json::Value = serde_json::from_str(&ok(&[
        "map",
        "--print-params",
        "--preset",
        "hifi",
        "--extend-ends",
    ]))
    .unwrap();
    assert_eq!(p["preset"], "hifi");
    assert!(p["params"]["ends"].is_object());
    assert!(p["params"]["rescue"].is_null());
    let v = ok(&["--version"]);
    assert!(v.starts_with("spearfish "), "{v}");
    // Two outputs can't share stdout.
    let o = Command::new(bin)
        .args(["map", "x.sfi", "-", "-o", "-", "--paf", "-"])
        .output()
        .unwrap();
    assert!(!o.status.success());
    // RAD can't go to stdout (checked before the index is opened, so x.sfi need not exist).
    let o = Command::new(bin)
        .args(["map", "x.sfi", "-", "--rad", "-"])
        .output()
        .unwrap();
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("RAD output needs a file"));
}

#[test]
fn simulate_and_map_on_toy() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    let (idx, reads, out_tsv) = (
        dir.join("toy_map.sfi"),
        dir.join("toy_reads.fa"),
        dir.join("toy_map.tsv"),
    );
    let run = |args: &[&std::ffi::OsStr]| {
        let o = Command::new(env!("CARGO_BIN_EXE_spearfish"))
            .args(args)
            .output()
            .unwrap();
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    };
    run(&[
        "index".as_ref(),
        toy().as_os_str(),
        "-o".as_ref(),
        idx.as_os_str(),
    ]);
    run(&[
        "simulate".as_ref(),
        toy().as_os_str(),
        "-n".as_ref(),
        "200".as_ref(),
        "--rate".as_ref(),
        "0.03".as_ref(),
        "-o".as_ref(),
        reads.as_os_str(),
    ]);
    let out_rad = dir.join("toy_map.rad");
    run(&[
        "map".as_ref(),
        idx.as_os_str(),
        reads.as_os_str(),
        "-t".as_ref(),
        "2".as_ref(),
        "-o".as_ref(),
        out_tsv.as_os_str(),
        "--rad".as_ref(),
        out_rad.as_os_str(),
        "--paf".as_ref(),
        dir.join("toy_map.paf").as_os_str(),
    ]);
    // Every read maps, its source transcript is among its top-scoring hits (truncated reads
    // can tie between isoforms), and each read's lines are contiguous (reads may appear in
    // any order).
    let text = std::fs::read_to_string(&out_tsv).unwrap();
    let mut reads: Vec<(String, Vec<(f64, String)>)> = Vec::new();
    for l in text.lines().skip(1) {
        let f: Vec<&str> = l.split('\t').collect();
        if reads.last().is_none_or(|r| r.0 != f[0]) {
            reads.push((f[0].to_string(), Vec::new()));
        }
        reads
            .last_mut()
            .unwrap()
            .1
            .push((f[5].parse().unwrap(), f[4].to_string()));
    }
    assert_eq!(
        reads.len(),
        200,
        "each read must appear as exactly one contiguous block"
    );
    let names: std::collections::HashSet<&str> = reads.iter().map(|r| r.0.as_str()).collect();
    assert_eq!(names.len(), 200);
    let mut correct = 0;
    for (name, hits) in &reads {
        let parts: Vec<&str> = name.split('|').collect();
        let best = hits.iter().map(|h| h.0).fold(f64::MIN, f64::max);
        correct += hits.iter().any(|h| h.0 == best && h.1 == parts[1]) as usize;
    }
    assert!(correct >= 195, "{correct}/200 top hits correct");

    // The RAD file holds the same hits as the TSV (it has no read names, so compare the
    // multisets of hits), and its header the targets and counts.
    let key = |t: &str, strand: &str, score: f32, ts: u32, te: u32| {
        format!("{t}\t{strand}\t{score:.2}\t{ts}\t{te}")
    };
    let mut tsv_hits: Vec<String> = text
        .lines()
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let n = |i: usize| f[i].parse::<u32>().unwrap();
            key(f[4], f[2], f[5].parse().unwrap(), n(8), n(9))
        })
        .collect();
    let (h, mut r) = long_read_rad::open(&out_rad, std::num::NonZeroUsize::new(2).unwrap()).unwrap();
    assert_eq!(h.meta.mapper, "spearfish");
    assert_eq!(h.meta.score_denom, spearfish::SCORE_DENOM);
    assert_eq!(h.counts.reads, 200);
    assert_eq!(h.counts.mapped, 200);
    assert_eq!(h.ref_names.len(), h.ref_lengths.len());
    // The toy headers give genes, which group its isoforms.
    let info = long_read_rad::IndexInfo::from_file_tags(&r.file_tag_map).unwrap();
    assert_eq!(info.genes_from, "fasta");
    assert_eq!(info.num_transcripts, h.ref_names.len() as u64);
    assert!(
        info.num_genes > 0 && info.num_genes < info.num_transcripts,
        "{info:?}"
    );
    assert!(info.gene_groups);
    let rad_hits = std::sync::Mutex::new(Vec::new());
    r.process_parallel(std::num::NonZeroUsize::new(2).unwrap(), |mc| {
        for c in mc.iter() {
            for rec in &c.reads {
                for a in rec.alns() {
                    let strand = if a.fw { "+" } else { "-" };
                    let t = &h.ref_names[a.tid as usize];
                    rad_hits
                        .lock()
                        .unwrap()
                        .push(key(t, strand, a.score, a.ref_start, a.ref_end));
                }
            }
        }
    })
    .unwrap();
    let mut rad_hits = rad_hits.into_inner().unwrap();
    tsv_hits.sort();
    rad_hits.sort();
    assert_eq!(rad_hits, tsv_hits);

    // PAF: one line per TSV hit (all reads map), with the same spans.
    let paf = std::fs::read_to_string(dir.join("toy_map.paf")).unwrap();
    let mut paf_hits: Vec<String> = paf
        .lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert!(f.len() >= 12, "{l}");
            let score: f32 = f
                .iter()
                .find_map(|x| x.strip_prefix("sc:f:"))
                .unwrap()
                .parse()
                .unwrap();
            key(f[5], f[4], score, f[7].parse().unwrap(), f[8].parse().unwrap())
        })
        .collect();
    paf_hits.sort();
    assert_eq!(paf_hits, tsv_hits);

    // The run record: the index, the settings, every parameter and the counts.
    let meta: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("toy_map.tsv.meta.json")).unwrap()).unwrap();
    assert_eq!(meta["kind"], "spearfish-map-run");
    assert_eq!(meta["stats"]["reads"], 200);
    assert_eq!(meta["stats"]["mapped"], 200);
    assert_eq!(meta["settings"]["preset"], "ont");
    assert_eq!(meta["parameters"]["preset"], "ont");
    assert_eq!(meta["index"]["index_id"].as_str().unwrap(), h.meta.index_id);
    assert!(meta["index"]["manifest"]["built_by"]["version"].is_string());
    assert_eq!(meta["index"]["genes"], info.num_genes);
    assert_eq!(meta["index"]["gene_groups"], true);
    assert_eq!(meta["index"]["manifest"]["input"]["genes_from"], "fasta");
    // The toy headers give genes: no warnings.
    assert!(meta.get("warnings").is_none(), "{}", meta["warnings"]);

    // Output to a device: no run record beside it; the first regular output file gets it.
    #[cfg(unix)]
    {
        let null = std::ffi::OsStr::new("/dev/null");
        let reads_fa = dir.join("toy_reads.fa");
        run(&[
            "map".as_ref(),
            idx.as_os_str(),
            reads_fa.as_os_str(),
            "-o".as_ref(),
            null,
        ]);
        let paf = dir.join("toy_null.paf");
        run(&[
            "map".as_ref(),
            idx.as_os_str(),
            reads_fa.as_os_str(),
            "-o".as_ref(),
            null,
            "--paf".as_ref(),
            paf.as_os_str(),
        ]);
        assert!(dir.join("toy_null.paf.meta.json").is_file());
    }
}

/// `em` skips reads flagged `low_coverage` whatever other flags the line carries.
/// The simulator's harder reads: fragments, random sequence at the ends, reads of nothing but
/// random sequence, and unexpressed transcripts. The same seed gives the same reads; random
/// reads stay unmapped and the rest still map to their source.
#[test]
fn bump_alloc_is_on_by_default_and_takes_an_optional_size() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    let (idx, reads) = (dir.join("toy_bump.sfi"), dir.join("toy_bump.fa"));
    let meta = dir.join("toy_bump.meta.json");
    let run = |args: &[&std::ffi::OsStr]| {
        Command::new(env!("CARGO_BIN_EXE_spearfish"))
            .args(args)
            .output()
            .unwrap()
    };
    let o = run(&[
        "index".as_ref(),
        toy().as_os_str(),
        "-o".as_ref(),
        idx.as_os_str(),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let o = run(&[
        "simulate".as_ref(),
        toy().as_os_str(),
        "-n".as_ref(),
        "100".as_ref(),
        "-o".as_ref(),
        reads.as_os_str(),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    // The options stand before the index, so a value `--bump-alloc` took without `=` would be
    // the index.
    let map = |flags: &[&str]| {
        let mut args: Vec<&std::ffi::OsStr> = vec!["map".as_ref()];
        args.extend(flags.iter().map(std::ffi::OsStr::new));
        args.extend([
            idx.as_os_str(),
            reads.as_os_str(),
            "-t".as_ref(),
            "2".as_ref(),
            "-v".as_ref(),
            "--meta".as_ref(),
            meta.as_os_str(),
        ]);
        let o = run(&args);
        let record: serde_json::Value = if o.status.success() {
            serde_json::from_str(&std::fs::read_to_string(&meta).unwrap()).unwrap()
        } else {
            serde_json::Value::Null
        };
        (o, record)
    };
    let lines = |o: &std::process::Output| {
        let mut l: Vec<String> = String::from_utf8_lossy(&o.stdout)
            .lines()
            .map(String::from)
            .collect();
        l.sort();
        l
    };
    let (off, _) = map(&["--no-bump-alloc"]);
    assert!(off.status.success());
    assert!(lines(&off).len() > 100);
    // Where this target has the arenas, the limit in force: the default is on, and says
    // nothing where the arenas are missing; a size asked for says so there.
    let arenas = cfg!(all(
        target_pointer_width = "64",
        any(target_os = "linux", target_os = "macos")
    ));
    for (flags, limit, asked) in [
        (&[][..], Some(8192), false),
        (&["--bump-alloc"][..], Some(8192), true),
        (&["--bump-alloc=8KB"][..], Some(8192), true),
        (&["--bump-alloc=1KB"][..], Some(1024), true),
        (&["--bump-alloc=500"][..], Some(512), true),
        (&["--bump-alloc=32kb"][..], Some(32768), true),
        (&["--no-bump-alloc", "--bump-alloc=1KB"][..], Some(1024), true),
        (&["--bump-alloc=1KB", "--no-bump-alloc"][..], None, false),
        (&["--no-bump-alloc"][..], None, false),
        (&["--bump-alloc=off"][..], None, false),
        (&["--bump-alloc=none"][..], None, false),
        (&["--bump-alloc=0"][..], None, false),
    ] {
        let (o, record) = map(flags);
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(o.status.success(), "{flags:?}: {err}");
        assert_eq!(lines(&o), lines(&off), "{flags:?}");
        let recorded = &record["stats"]["bump_alloc"];
        match limit {
            Some(bytes) if arenas => {
                assert!(
                    err.contains(&format!("bump-alloc: blocks of up to {bytes} bytes")),
                    "{flags:?}: {err}"
                );
                assert_eq!(recorded["max_block"], bytes, "{flags:?}: {record}");
                assert!(
                    recorded["arenas_peak"].as_u64().unwrap() >= 1,
                    "{flags:?}: {record}"
                );
            }
            Some(_) if asked => {
                assert!(err.contains("--bump-alloc is not available"), "{flags:?}: {err}");
                assert!(recorded.is_null(), "{flags:?}: {record}");
            }
            _ => {
                assert!(!err.contains("bump-alloc"), "{flags:?}: {err}");
                assert!(recorded.is_null(), "{flags:?}: {record}");
            }
        }
    }
    for (flag, why) in [
        ("--bump-alloc=64KB", "at most 32 KB"),
        ("--bump-alloc=8MB", "expected a size"),
        ("--bump-alloc=on", "expected a size"),
        ("--bump-alloc=", "expected a size"),
        ("--no-bump-alloc=1KB", "unexpected value"),
    ] {
        let (o, _) = map(&[flag]);
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(!o.status.success(), "{flag}");
        assert!(err.contains(why), "{flag}: {err}");
    }
}

#[test]
fn simulated_random_reads_stay_unmapped() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    let run = |args: &[&std::ffi::OsStr]| {
        let o = Command::new(env!("CARGO_BIN_EXE_spearfish"))
            .args(args)
            .output()
            .unwrap();
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    };
    let (toy, idx) = (toy(), dir.join("hard_map.sfi"));
    run(&["index".as_ref(), toy.as_os_str(), "-o".as_ref(), idx.as_os_str()]);
    let simulate = |reads: &std::path::Path, truth: &std::path::Path| {
        let mut args: Vec<&std::ffi::OsStr> = vec!["simulate".as_ref(), toy.as_os_str()];
        args.extend(
            [
                "-n",
                "400",
                "--rate",
                "0.05",
                "--rate-sd",
                "0.02",
                "--expressed",
                "0.5",
                "--fragment",
                "0.5",
                "--min-fragment",
                "300",
                "--junk",
                "40",
                "--unmappable",
                "0.25",
                "--seed",
                "3",
                "-o",
            ]
            .map(std::ffi::OsStr::new),
        );
        args.extend([reads.as_os_str(), "--truth".as_ref(), truth.as_os_str()]);
        run(&args);
    };
    let (reads, truth) = (dir.join("hard_reads.fa"), dir.join("hard_truth.tsv"));
    let (reads2, truth2) = (dir.join("hard_reads2.fa"), dir.join("hard_truth2.tsv"));
    simulate(&reads, &truth);
    simulate(&reads2, &truth2);
    let fa = std::fs::read_to_string(&reads).unwrap();
    assert_eq!(fa, std::fs::read_to_string(&reads2).unwrap());
    let truth = std::fs::read_to_string(&truth).unwrap();
    assert_eq!(truth, std::fs::read_to_string(&truth2).unwrap());

    // Names: `r<i>|<transcript>|<strand>|<start>|<end>|<head>|<tail>`, or `r<i>|*|*|0|0` for
    // a random read; the truth counts the others, on expressed transcripts only.
    let names: Vec<Vec<&str>> = fa
        .lines()
        .filter_map(|l| l.strip_prefix('>'))
        .map(|n| n.split('|').collect())
        .collect();
    assert_eq!(names.len(), 400);
    let random = names.iter().filter(|n| n[1] == "*").count();
    assert!((50..150).contains(&random), "{random} random reads");
    assert!(names.iter().all(|n| n.len() == if n[1] == "*" { 5 } else { 7 }));
    let rows: Vec<Vec<&str>> = truth.lines().skip(1).map(|l| l.split('\t').collect()).collect();
    let count = |r: &Vec<&str>| r[2].parse::<usize>().unwrap();
    assert_eq!(rows.iter().map(count).sum::<usize>(), 400 - random);
    let unexpressed = rows
        .iter()
        .filter(|r| r[3].parse::<f64>().unwrap() == 0.0)
        .count();
    assert!(unexpressed > 0 && unexpressed < rows.len());
    assert!(
        rows.iter()
            .all(|r| count(r) == 0 || r[3].parse::<f64>().unwrap() > 0.0)
    );

    let out = dir.join("hard_map.tsv");
    run(&[
        "map".as_ref(),
        idx.as_os_str(),
        reads.as_os_str(),
        "-t".as_ref(),
        "2".as_ref(),
        "-q".as_ref(),
        "-o".as_ref(),
        out.as_os_str(),
    ]);
    let (mut real, mut correct) = (0, 0);
    for l in std::fs::read_to_string(&out).unwrap().lines().skip(1) {
        let f: Vec<&str> = l.split('\t').collect();
        let name: Vec<&str> = f[0].split('|').collect();
        if name[1] == "*" {
            assert!(f[4] == "*" && f[13].starts_with("unmapped"), "{l}");
        } else if f[3] == "0" {
            real += 1;
            correct += (f[4] != "*" && !f[13].contains("low_score")) as usize;
        }
    }
    assert_eq!(real, 400 - random);
    assert!(correct * 100 >= real * 95, "{correct}/{real} reads mapped");
}

#[test]
fn em_skips_low_coverage_reads_with_other_flags() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    let (map, out) = (dir.join("flags_map.tsv"), dir.join("flags_quant.tsv"));
    let head = "read\tread_len\tstrand\trank\ttranscript\tscore\tqstart\tqend\ttstart\ttend\tanchors\tanchored\tmax_disc\tflags\n";
    let line = |r: &str, t: &str, flags: &str| {
        format!("{r}\t1000\t+\t0\t{t}\t100.00\t0\t1000\t0\t1000\t10\t500\t0\t{flags}\n")
    };
    let tsv = [
        head.to_string(),
        line("r1", "TOYT0000001", "-"),
        line("r2", "TOYT0000002", "low_coverage,anchor_span"),
        line("r3", "TOYT0000002", "rescued,low_coverage"),
        line("r4", "TOYT0000003", "anchor_span"),
        // A kept chance-level hit (`--keep-low-score`) never counts.
        line("r5", "TOYT0000004", "low_score"),
    ]
    .concat();
    std::fs::write(&map, tsv).unwrap();
    let o = Command::new(env!("CARGO_BIN_EXE_spearfish"))
        .args([
            "em".as_ref(),
            toy().as_os_str(),
            map.as_os_str(),
            "-o".as_ref(),
            out.as_os_str(),
            "-q".as_ref(),
        ])
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    // `-q` silences the summary.
    assert!(o.stderr.is_empty(), "{}", String::from_utf8_lossy(&o.stderr));
    let counts: std::collections::HashMap<String, f64> = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_string(), f[2].parse().unwrap())
        })
        .collect();
    assert_eq!(counts["TOYT0000001"], 1.0);
    assert_eq!(counts["TOYT0000002"], 0.0);
    assert_eq!(counts["TOYT0000003"], 1.0);
    assert_eq!(counts["TOYT0000004"], 0.0);
}

/// Run `em` on the toy transcripts and the mapping TSV `tsv` (written to `name`), with `extra`
/// arguments: the exit status, stderr, and the output file's text.
fn run_em(name: &str, tsv: &str, extra: &[&str]) -> (bool, String, String) {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("em_columns");
    std::fs::create_dir_all(&dir).unwrap();
    let (map, out) = (
        dir.join(format!("{name}.tsv")),
        dir.join(format!("{name}.quant.tsv")),
    );
    std::fs::write(&map, tsv).unwrap();
    let _ = std::fs::remove_file(&out);
    let o = Command::new(env!("CARGO_BIN_EXE_spearfish"))
        .args([
            "em".as_ref(),
            toy().as_os_str(),
            map.as_os_str(),
            "-o".as_ref(),
            out.as_os_str(),
        ])
        .args(extra)
        .output()
        .unwrap();
    (
        o.status.success(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
        std::fs::read_to_string(&out).unwrap_or_default(),
    )
}

/// `em`'s output as transcript → estimated reads.
fn em_counts(quant: &str) -> std::collections::HashMap<String, f64> {
    quant
        .lines()
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_string(), f[2].parse().unwrap())
        })
        .collect()
}

const EM_HEAD_14: &str = "read\tread_len\tstrand\trank\ttranscript\tscore\tqstart\tqend\ttstart\ttend\tanchors\tanchored\tmax_disc\tflags\n";
const EM_HEAD_15: &str = "read\tread_len\tstrand\trank\ttranscript\tscore\tqstart\tqend\ttstart\ttend\tanchors\tanchored\tmax_disc\tflags\tqscore\n";

/// A TSV line (14 columns, or 15 with `qscore`).
fn em_line(read: &str, tx: &str, score: &str, flags: &str, qscore: Option<&str>) -> String {
    let mut l = format!("{read}\t1000\t+\t0\t{tx}\t{score}\t0\t1000\t0\t1000\t10\t500\t0\t{flags}");
    if let Some(q) = qscore {
        l.push('\t');
        l.push_str(q);
    }
    l.push('\n');
    l
}

/// `em` finds the flags in a 14-column (0.2.x) and a 15-column (`qscore`) TSV alike: hits
/// flagged `low_score` are never used, and `low_coverage` ones only with
/// `--keep-low-coverage`. (0.3.0 and 0.3.1 read the last column as the flags, `qscore` in a
/// 15-column file, and so used them all.)
#[test]
fn em_reads_flags_in_14_and_15_column_tsvs() {
    let reads = [
        ("r1", "TOYT0000001", "-"),
        ("r2", "TOYT0000002", "low_coverage"),
        ("r3", "TOYT0000003", "rescued,low_coverage"),
        ("r4", "TOYT0000004", "low_score"),
        ("r5", "TOYT0000005", "anchor_span"),
    ];
    let tsv = |q: bool| -> String {
        let head = if q { EM_HEAD_15 } else { EM_HEAD_14 };
        let mut t = head.to_string();
        for (r, tx, fl) in reads {
            t.push_str(&em_line(r, tx, "100.00", fl, q.then_some("100.00")));
        }
        // An unmapped read's line.
        t.push_str(&em_line("r6", "*", "0", "unmapped", q.then_some("0")));
        t
    };
    for q in [false, true] {
        let name = if q { "flags15" } else { "flags14" };
        let (ok, log, quant) = run_em(name, &tsv(q), &[]);
        assert!(ok, "{log}");
        assert!(
            log.contains("2 reads used (1 unmapped, 3 low-coverage or low-score skipped"),
            "{log}"
        );
        let c = em_counts(&quant);
        assert_eq!(
            [1, 2, 3, 4, 5, 6].map(|i| c[&format!("TOYT000000{i}")]),
            [1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            "{name}"
        );
        let (ok, log, quant) = run_em(name, &tsv(q), &["--keep-low-coverage"]);
        assert!(ok, "{log}");
        let c = em_counts(&quant);
        assert_eq!(
            [1, 2, 3, 4, 5, 6].map(|i| c[&format!("TOYT000000{i}")]),
            [1.0, 1.0, 1.0, 0.0, 1.0, 0.0],
            "{name} --keep-low-coverage"
        );
    }
}

/// `em` weighs hits by `qscore` when the TSV has it (`--score-column auto`), by `score` with
/// `--score-column score`, and refuses `--score-column qscore` without the column. Columns are
/// found by name, in any order.
#[test]
fn em_weighs_hits_by_qscore() {
    // By score, TOYT0000001 is best and TOYT0000002 within `--min-frac`; by qscore (the end fit
    // moved TOYT0000002 up), TOYT0000002 is best by 5 points.
    let tsv = [
        EM_HEAD_15.to_string(),
        em_line("r1", "TOYT0000001", "100.00", "-", Some("100.00")),
        em_line("r1", "TOYT0000002", "98.00", "-", Some("105.00")),
    ]
    .concat();
    let (ok, log, quant) = run_em("qscore_auto", &tsv, &[]);
    assert!(ok && log.contains("weights from `qscore`"), "{log}");
    let c = em_counts(&quant);
    assert!(c["TOYT0000002"] > 0.99 && c["TOYT0000001"] < 0.01, "{quant}");
    let (ok, log, quant) = run_em("qscore_score", &tsv, &["--score-column", "score"]);
    assert!(ok && log.contains("weights from `score`"), "{log}");
    let c = em_counts(&quant);
    assert!(c["TOYT0000001"] > 0.99 && c["TOYT0000002"] < 0.01, "{quant}");

    // r1's two hits tie by score and are 2 points apart by qscore; r2 and r3 make the two
    // transcripts otherwise equally abundant. By score r1 splits 0.5 / 0.5; by qscore it is
    // exp(2 / 2) times as likely on TOYT0000002, which gets e / (1 + e) of it, about 0.73.
    let rows = [
        ("r1", "TOYT0000001", "100.00", "100.00"),
        ("r1", "TOYT0000002", "100.00", "102.00"),
        ("r2", "TOYT0000001", "90.00", "90.00"),
        ("r3", "TOYT0000002", "90.00", "90.00"),
    ];
    let tie: String = std::iter::once(EM_HEAD_15.to_string())
        .chain(rows.iter().map(|&(r, t, s, q)| em_line(r, t, s, "-", Some(q))))
        .collect();
    let s = em_counts(&run_em("tie_score", &tie, &["--score-column", "score"]).2);
    assert_eq!((s["TOYT0000001"], s["TOYT0000002"]), (1.5, 1.5));
    let (ok, log, quant) = run_em("tie_qscore", &tie, &["--score-column", "qscore"]);
    assert!(ok, "{log}");
    let q = em_counts(&quant);
    assert!(q["TOYT0000002"] > 1.6 && q["TOYT0000001"] < 1.4, "{q:?}");
    assert_eq!(em_counts(&run_em("tie_auto", &tie, &[]).2), q);
    // The same columns in another order, with no others: found by name.
    let moved: String = std::iter::once("qscore\tread\tflags\ttranscript\tscore\n".to_string())
        .chain(
            rows.iter()
                .map(|&(r, t, s, q)| format!("{q}\t{r}\t-\t{t}\t{s}\n")),
        )
        .collect();
    let (ok, log, quant) = run_em("tie_reordered", &moved, &[]);
    assert!(ok, "{log}");
    assert_eq!(em_counts(&quant), q);

    // No `qscore` column: auto weighs by score; asking for qscore is an error.
    let tsv14 = [
        EM_HEAD_14.to_string(),
        em_line("r1", "TOYT0000001", "100.00", "-", None),
    ]
    .concat();
    let (ok, log, _) = run_em("noq_auto", &tsv14, &[]);
    assert!(ok && log.contains("weights from `score`"), "{log}");
    let (ok, log, _) = run_em("noq_qscore", &tsv14, &["--score-column", "qscore"]);
    assert!(!ok && log.contains("no `qscore` column"), "{log}");
}

/// A 14-column TSV gives exactly 0.3.1's output (checked against the 0.3.1 binary), with a
/// header and without one. Without a header the first line is a read's; with a header `em`
/// does not know, the 14-column layout is read and a warning logged.
#[test]
fn em_on_14_columns_matches_0_3_1() {
    let rows = [
        ("r1", "TOYT0000001", "100.00", "-"),
        ("r1", "TOYT0000002", "99.00", "-"),
        ("r1", "TOYT0000003", "96.50", "-"),
        ("r2", "TOYT0000002", "80.00", "-"),
        ("r2", "TOYT0000001", "79.00", "-"),
        ("r3", "TOYT0000004", "120.00", "anchor_span"),
        ("r3", "TOYT0000005", "117.50", "anchor_span"),
        ("r4", "TOYT0000003", "60.00", "rescued"),
        ("r5", "TOYT0000006", "50.00", "low_coverage"),
        ("r5", "TOYT0000005", "49.00", "low_coverage"),
        ("r6", "NOT_A_TX", "70.00", "-"),
        ("r7", "*", "0", "unmapped"),
        ("r8", "TOYT0000001", "90.00", "-"),
        ("r8", "TOYT0000006", "90.00", "-"),
        ("r9", "TOYT0000004", "30.00", "low_score"),
    ];
    let body: String = rows
        .iter()
        .map(|&(r, t, s, f)| em_line(r, t, s, f, None))
        .collect();
    // `spearfish em` 0.3.1 on the same file.
    let v031 = "transcript\tlength\treads\nTOYT0000001\t1270\t2.9411\nTOYT0000002\t1150\t0.0000\nTOYT0000003\t1100\t1.0589\nTOYT0000004\t1064\t1.0000\nTOYT0000005\t1050\t0.0000\nTOYT0000006\t1050\t0.0000\n";
    let v031_keep = "transcript\tlength\treads\nTOYT0000001\t1270\t2.5923\nTOYT0000002\t1150\t0.0000\nTOYT0000003\t1100\t1.0667\nTOYT0000004\t1064\t1.0000\nTOYT0000005\t1050\t0.0000\nTOYT0000006\t1050\t1.3409\n";
    for (name, tsv) in [
        ("golden14", format!("{EM_HEAD_14}{body}")),
        ("golden14_noheader", body.clone()),
        ("golden14_otherheader", format!("a\tb\tc\n{body}")),
    ] {
        let (ok, log, quant) = run_em(name, &tsv, &[]);
        assert!(ok, "{log}");
        assert_eq!(quant, v031, "{name}");
        assert_eq!(
            log.contains("header names no"),
            name == "golden14_otherheader",
            "{log}"
        );
        let (ok, log, quant) = run_em(name, &tsv, &["--keep-low-coverage"]);
        assert!(ok, "{log}");
        assert_eq!(quant, v031_keep, "{name} --keep-low-coverage");
    }
}

/// The toy transcripts with Ensembl-style headers (`>TOYT0000001 cdna gene:TOYG0000001`): the
/// ID ends at the first whitespace.
fn described_toy(path: &std::path::Path) {
    let fa: String = std::fs::read_to_string(toy())
        .unwrap()
        .lines()
        .map(|l| match l.strip_prefix('>') {
            Some(id) => {
                let f: Vec<&str> = id.split('|').collect();
                format!(">{} cdna gene:{}\n", f[0], f[1])
            }
            None => format!("{l}\n"),
        })
        .collect();
    std::fs::write(path, fa).unwrap();
}

/// simulate, em and halo name transcripts as the index does when headers carry descriptions.
#[test]
fn description_headers_name_transcripts_like_the_index() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("described");
    std::fs::create_dir_all(&dir).unwrap();
    let fa = dir.join("transcripts.fa");
    described_toy(&fa);
    let run = |args: &[&std::ffi::OsStr]| {
        let o = Command::new(env!("CARGO_BIN_EXE_spearfish"))
            .args(args)
            .output()
            .unwrap();
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        String::from_utf8(o.stderr).unwrap()
    };
    let toy_names: std::collections::HashSet<String> = (1..=6).map(|i| format!("TOYT000000{i}")).collect();

    // simulate: read names are r<i>|<transcript>|<strand>|<start>|<end>, and the truth names
    // are the index's names.
    let (reads, truth) = (dir.join("reads.fa"), dir.join("truth.tsv"));
    run(&[
        "simulate".as_ref(),
        fa.as_os_str(),
        "-n".as_ref(),
        "20".as_ref(),
        "-o".as_ref(),
        reads.as_os_str(),
        "--truth".as_ref(),
        truth.as_os_str(),
    ]);
    for l in std::fs::read_to_string(&reads).unwrap().lines() {
        if let Some(id) = l.strip_prefix('>') {
            let f: Vec<&str> = id.split('|').collect();
            assert!(f.len() == 5 && toy_names.contains(f[1]), "{id}");
        }
    }
    let truth_names: std::collections::HashSet<String> = std::fs::read_to_string(&truth)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();
    assert_eq!(truth_names, toy_names);

    // em: hits name the index's transcripts, which must be found.
    let (map, quant) = (dir.join("map.tsv"), dir.join("quant.tsv"));
    let line = |r: &str, t: &str| format!("{r}\t1000\t+\t0\t{t}\t100.00\t0\t1000\t0\t1000\t10\t500\t0\t-\n");
    std::fs::write(
        &map,
        [
            "read\tread_len\tstrand\trank\ttranscript\tscore\tqstart\tqend\ttstart\ttend\tanchors\tanchored\tmax_disc\tflags\n".to_string(),
            line("r1", "TOYT0000001"),
            line("r2", "TOYT0000003"),
        ]
        .concat(),
    )
    .unwrap();
    let log = run(&[
        "em".as_ref(),
        fa.as_os_str(),
        map.as_os_str(),
        "-o".as_ref(),
        quant.as_os_str(),
    ]);
    assert!(
        log.contains("2 reads used") && log.contains(" 0 hits to unknown"),
        "{log}"
    );
    let counts: std::collections::HashMap<String, f64> = std::fs::read_to_string(&quant)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_string(), f[2].parse().unwrap())
        })
        .collect();
    assert_eq!(
        counts.keys().cloned().collect::<std::collections::HashSet<_>>(),
        toy_names
    );
    assert_eq!((counts["TOYT0000001"], counts["TOYT0000003"]), (1.0, 1.0));

    // halo: every transcript is found in the GTF.
    let log = run(&[
        "halo".as_ref(),
        fa.as_os_str(),
        toy().with_file_name("annotation.gtf").as_os_str(),
        "-s".as_ref(),
        "lcp:K=3".as_ref(),
        "--hpc".as_ref(),
        "off".as_ref(),
        "-t".as_ref(),
        "2".as_ref(),
        "-o".as_ref(),
        dir.join("halo").as_os_str(),
    ]);
    assert!(log.contains("6 transcripts (0 not in GTF"), "{log}");
}

/// Headers that are not UTF-8 are read lossily, in transcripts and in reads.
#[test]
fn non_utf8_headers() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("non_utf8");
    std::fs::create_dir_all(&dir).unwrap();
    let toy_fa = std::fs::read(toy()).unwrap();
    // A transcript header with byte 0xE9 (Latin-1 'é').
    let mut fa = b">TOYT0000000_\xe9|TOYG0000000\n".to_vec();
    let first_seq: Vec<u8> = toy_fa
        .split(|&b| b == b'\n')
        .skip(1)
        .take_while(|l| !l.starts_with(b">"))
        .flat_map(|l| l.iter().copied())
        .collect();
    fa.extend_from_slice(&first_seq);
    fa.push(b'\n');
    fa.extend_from_slice(&toy_fa);
    let (tx, idx) = (dir.join("transcripts.fa"), dir.join("toy.sfi"));
    std::fs::write(&tx, &fa).unwrap();
    let run = |args: &[&std::ffi::OsStr]| {
        let o = Command::new(env!("CARGO_BIN_EXE_spearfish"))
            .args(args)
            .output()
            .unwrap();
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    };
    run(&["index".as_ref(), tx.as_os_str(), "-o".as_ref(), idx.as_os_str()]);
    let mut reads = b">read_\xe9 sample=caf\xe9\n".to_vec();
    reads.extend_from_slice(&first_seq);
    reads.push(b'\n');
    let (rd, tsv, paf) = (dir.join("reads.fa"), dir.join("out.tsv"), dir.join("out.paf"));
    std::fs::write(&rd, reads).unwrap();
    run(&[
        "map".as_ref(),
        idx.as_os_str(),
        rd.as_os_str(),
        "-o".as_ref(),
        tsv.as_os_str(),
        "--paf".as_ref(),
        paf.as_os_str(),
    ]);
    let text = std::fs::read_to_string(&tsv).unwrap();
    let rows: Vec<Vec<&str>> = text.lines().skip(1).map(|l| l.split('\t').collect()).collect();
    assert!(
        !rows.is_empty() && rows.iter().all(|r| r[0] == "read_\u{fffd}"),
        "{text}"
    );
    assert!(rows.iter().any(|r| r[4] == "TOYT0000000_\u{fffd}"), "{text}");
    assert!(
        std::fs::read_to_string(&paf)
            .unwrap()
            .starts_with("read_\u{fffd}\t")
    );
}

/// Paths that are not valid Unicode work as arguments, and are recorded lossily.
/// Linux only: macOS file systems refuse names that are not UTF-8, so the case can't arise there.
#[cfg(target_os = "linux")]
#[test]
fn non_unicode_path_arguments() {
    use std::os::unix::ffi::OsStrExt;
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("non_unicode");
    std::fs::create_dir_all(&dir).unwrap();
    let path = |name: &[u8]| dir.join(std::ffi::OsStr::from_bytes(name));
    let (idx, tsv, rad) = (path(b"idx_\xe9.sfi"), path(b"o_\xe9.tsv"), path(b"o_\xe9.rad"));
    let run = |args: &[&std::ffi::OsStr]| {
        let o = Command::new(env!("CARGO_BIN_EXE_spearfish"))
            .args(args)
            .output()
            .unwrap();
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    };
    run(&[
        "index".as_ref(),
        toy().as_os_str(),
        "-o".as_ref(),
        idx.as_os_str(),
    ]);
    run(&[
        "map".as_ref(),
        idx.as_os_str(),
        toy().as_os_str(),
        "-o".as_ref(),
        tsv.as_os_str(),
        "--rad".as_ref(),
        rad.as_os_str(),
    ]);
    // The run records hold the arguments with U+FFFD for the invalid byte.
    let has = |record: &[u8], arg: &str| {
        let j: serde_json::Value = serde_json::from_slice(&std::fs::read(path(record)).unwrap()).unwrap();
        let args = j["command_line"].as_array().unwrap();
        assert!(
            args.iter().any(|a| a.as_str().unwrap().ends_with(arg)),
            "{args:?}"
        );
    };
    has(b"idx_\xe9.sfi.json", "idx_\u{fffd}.sfi");
    has(b"o_\xe9.tsv.meta.json", "o_\u{fffd}.tsv");
    let h = long_read_rad::read_header(&rad).unwrap();
    assert!(h.meta.cmdline.contains("o_\u{fffd}.rad"), "{}", h.meta.cmdline);
}

/// A failed write of the last buffered output is an error, not a silent truncation.
#[cfg(target_os = "linux")]
#[test]
fn final_flush_errors_are_reported() {
    use std::ffi::OsStr;
    let (tx, full) = (toy(), OsStr::new("/dev/full"));
    let idx = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("full_t2g.sfi");
    let s = |x: &'static str| OsStr::new(x);
    for args in [
        vec![s("lcp-stats"), tx.as_os_str(), s("-K"), s("4"), s("-o"), full],
        vec![s("simulate"), tx.as_os_str(), s("-n"), s("5"), s("-o"), full],
        vec![
            s("index"),
            tx.as_os_str(),
            s("-o"),
            idx.as_os_str(),
            s("--write-t2g"),
            full,
        ],
    ] {
        let o = Command::new(env!("CARGO_BIN_EXE_spearfish"))
            .args(&args)
            .output()
            .unwrap();
        assert!(!o.status.success(), "{args:?}");
        assert!(
            String::from_utf8_lossy(&o.stderr).contains("No space left"),
            "{args:?}"
        );
    }
}

/// An index whose every transcript is its own gene (headers without gene IDs, no --t2g or
/// --infer-genes): `index` and `map` both warn, under `-q` too, `map` once per run, and the run
/// record lists the warning. With genes, `-q` stays silent. The RAD header and the run record
/// carry the gene labels' provenance in both cases.
#[test]
fn warns_on_index_without_gene_groups() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("geneless");
    std::fs::create_dir_all(&dir).unwrap();
    let bin = env!("CARGO_BIN_EXE_spearfish");
    let run = |args: &[&std::ffi::OsStr]| {
        let o = Command::new(bin).args(args).output().unwrap();
        assert!(
            o.status.success(),
            "{args:?}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        String::from_utf8(o.stderr).unwrap()
    };
    let os = |s: &'static str| std::ffi::OsStr::new(s);
    // The toy transcripts, headers cut to the transcript name.
    let plain = dir.join("plain.fa");
    let fa: String = std::fs::read_to_string(toy())
        .unwrap()
        .lines()
        .map(|l| match l.strip_prefix('>') {
            Some(h) => format!(">{}\n", h.split('|').next().unwrap()),
            None => format!("{l}\n"),
        })
        .collect();
    std::fs::write(&plain, fa).unwrap();
    let reads = dir.join("reads.fa");
    run(&[
        os("simulate"),
        toy().as_os_str(),
        os("-n"),
        os("20"),
        os("-q"),
        os("-o"),
        reads.as_os_str(),
    ]);

    let warning = "has no gene groups: every transcript is its own gene";
    let idx = dir.join("plain.sfi");
    for quiet in [false, true] {
        let mut args = vec![os("index"), plain.as_os_str(), os("-o"), idx.as_os_str()];
        if quiet {
            args.push(os("-q"));
        }
        let err = run(&args);
        assert_eq!(err.matches(warning).count(), 1, "quiet {quiet}: {err}");
        assert!(
            err.contains("--t2g FILE") && err.contains("--infer-genes"),
            "{err}"
        );
        assert_eq!(quiet, err.lines().count() == 1, "{err}");
    }

    let out = dir.join("plain.tsv");
    let rad = dir.join("plain.rad");
    let err = run(&[
        os("map"),
        idx.as_os_str(),
        reads.as_os_str(),
        os("-q"),
        os("-t"),
        os("2"),
        os("-o"),
        out.as_os_str(),
        os("--rad"),
        rad.as_os_str(),
    ]);
    assert_eq!(err.matches(warning).count(), 1, "{err}");
    assert!(err.contains("rebuild it with --t2g FILE"), "{err}");
    let meta: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("plain.tsv.meta.json")).unwrap()).unwrap();
    let w = meta["warnings"].as_array().unwrap();
    assert_eq!(w.len(), 1, "{w:?}");
    assert_eq!(w[0]["code"], "no_gene_groups");
    assert!(w[0]["message"].as_str().unwrap().contains(warning));
    // The RAD header and the run record say the index has no gene groups, and why.
    let n_tx = long_read_rad::read_header(&rad).unwrap().ref_names.len() as u64;
    assert_eq!(
        long_read_rad::read_index_info(&rad).unwrap(),
        Some(long_read_rad::IndexInfo::new("none", n_tx, n_tx, false))
    );
    assert_eq!(meta["index"]["genes"], n_tx);
    assert_eq!(meta["index"]["gene_groups"], false);

    // Inferred genes group the toy's isoforms: no warning.
    let inferred = dir.join("inferred.sfi");
    let err = run(&[
        os("index"),
        plain.as_os_str(),
        os("--infer-genes"),
        os("-q"),
        os("-o"),
        inferred.as_os_str(),
    ]);
    assert!(err.is_empty(), "{err}");
    let err = run(&[
        os("map"),
        inferred.as_os_str(),
        reads.as_os_str(),
        os("-q"),
        os("-o"),
        out.as_os_str(),
        os("--rad"),
        rad.as_os_str(),
    ]);
    assert!(err.is_empty(), "{err}");
    let meta: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("plain.tsv.meta.json")).unwrap()).unwrap();
    assert!(meta.get("warnings").is_none());
    let info = long_read_rad::read_index_info(&rad).unwrap().unwrap();
    assert_eq!(info.genes_from, "inferred");
    assert_eq!(info.num_transcripts, n_tx);
    assert!(info.num_genes < n_tx, "{info:?}");
    assert!(info.gene_groups);
    assert_eq!(meta["index"]["genes"], info.num_genes);
    assert_eq!(meta["index"]["gene_groups"], true);
}

fn toy_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/toy")
        .join(name)
}

fn spearfish(args: &[&std::ffi::OsStr]) -> (bool, String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_spearfish"))
        .args(args)
        .output()
        .unwrap();
    (
        out.status.success(),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

/// `index --genome` builds a genome index, and `map-spliced` maps to it: the toy transcripts
/// come back as their annotated exons, in PAF, SAM and BAM, and the run is recorded.
#[test]
fn genome_index_and_spliced_map_on_toy() {
    use std::ffi::OsStr;
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("genome_cli");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let (idx, genome, gtf, reads) = (
        dir.join("toy.sfg"),
        toy_file("genome.fa"),
        toy_file("annotation.gtf"),
        toy(),
    );
    let o = |s: &'static str| OsStr::new(s);
    let (ok, _, err) = spearfish(&[
        o("index"),
        o("--genome"),
        genome.as_os_str(),
        o("-o"),
        idx.as_os_str(),
    ]);
    assert!(ok, "{err}");
    assert!(err.contains("genome") && err.contains("sequences 1"), "{err}");
    assert!(dir.join("toy.sfg.json").exists());

    let (ok, out, err) = spearfish(&[o("inspect"), idx.as_os_str()]);
    assert!(ok, "{err}");
    assert!(
        out.contains("kind\tgenome") && out.contains("bases\t18090"),
        "{out}"
    );
    let (ok, out, err) = spearfish(&[o("inspect"), o("--json"), idx.as_os_str()]);
    assert!(ok, "{err}");
    let j: serde_json::Value = serde_json::from_str(&out).unwrap();
    assert_eq!(j["kind"], "genome");
    assert_eq!(j["manifest"]["kind"], "spearfish-genome-index");
    let (ok, out, err) = spearfish(&[o("doctor"), o("--index"), idx.as_os_str()]);
    assert!(ok, "{err}");
    assert!(out.contains("ok: genome, 1 sequences"), "{out}");

    // PAF on stdout by default; with the annotation, exact reads give the annotated exons.
    let map = |extra: &[&OsStr]| {
        let mut a = vec![o("map-spliced"), idx.as_os_str(), reads.as_os_str()];
        a.extend_from_slice(extra);
        spearfish(&a)
    };
    let (ok, paf, err) = map(&[o("-a"), gtf.as_os_str(), o("-t"), o("2")]);
    assert!(ok, "{err}");
    assert!(err.contains("annotation: 6 transcripts, 11 introns"), "{err}");
    let cg: std::collections::HashMap<&str, &str> = paf
        .lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!((f[4], f[5], f[11]), ("+", "toy", "60"), "{l}");
            (f[0].split('|').next().unwrap(), l.rsplit("cg:Z:").next().unwrap())
        })
        .collect();
    assert_eq!(cg.len(), 6);
    assert_eq!(cg["TOYT0000001"], "400M1176N250M1114N120M1321N500M");
    assert_eq!(cg["TOYT0000002"], "400M1176N250M2555N500M");
    assert_eq!(cg["TOYT0000004"], "350M793N264M1194N450M");

    // SAM and BAM files, with a base-level alignment, and a run record beside the first.
    let (sam, bam) = (dir.join("out.sam"), dir.join("out.bam"));
    let (ok, out, err) = map(&[
        o("--annotation"),
        gtf.as_os_str(),
        o("--preset"),
        o("drna"),
        o("--sam"),
        sam.as_os_str(),
        o("--bam"),
        bam.as_os_str(),
        o("--base-level"),
        o("affine"),
        o("--secondary"),
        o("2"),
    ]);
    assert!(ok && out.is_empty(), "{err}");
    let sam_text = std::fs::read_to_string(&sam).unwrap();
    assert!(
        sam_text.starts_with("@HD\tVN:1.6\tSO:unsorted\tGO:query\n@SQ\tSN:toy\tLN:18090\n@PG\tID:spearfish")
    );
    let line = sam_text.lines().find(|l| l.starts_with("TOYT0000001")).unwrap();
    let rec: Vec<&str> = line.split('\t').collect();
    assert_eq!(
        (rec[1], rec[2], rec[3], rec[5]),
        ("0", "toy", "858", "400M1176N250M1114N120M1321N500M")
    );
    assert!(rec.contains(&"NM:i:0") && rec.contains(&"XS:A:+"), "{rec:?}");
    let bam_bytes = std::fs::read(&bam).unwrap();
    assert_eq!(&bam_bytes[..4], &[0x1f, 0x8b, 8, 4]);
    // (the empty BGZF block that ends a BAM file)
    assert_eq!(&bam_bytes[bam_bytes.len() - 28..][..4], &[0x1f, 0x8b, 8, 4]);
    let meta: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("out.sam.meta.json")).unwrap()).unwrap();
    assert_eq!(meta["kind"], "spearfish-genome-map-run");
    assert_eq!(meta["index"]["kind"], "genome");
    assert_eq!(meta["annotation"]["transcripts"], 6);
    assert_eq!(meta["preset"], "drna");
    assert_eq!(meta["parameters"]["sense_reads"], true);
    assert_eq!(meta["parameters"]["max_hits"], 3);
    assert_eq!(
        (meta["stats"]["reads"].as_u64(), meta["stats"]["spliced"].as_u64()),
        (Some(6), Some(6))
    );

    // GFF3 and BED12 give the same alignments as the GTF.
    let (gff, bed) = (dir.join("toy.gff3"), dir.join("toy.bed"));
    let (mut g, mut b) = (String::from("##gff-version 3\n"), String::new());
    let mut tx: std::collections::BTreeMap<String, Vec<(u32, u32)>> = Default::default();
    for l in std::fs::read_to_string(&gtf).unwrap().lines() {
        let p: Vec<&str> = l.split('\t').collect();
        if p.len() == 9 && p[2] == "exon" {
            let id = p[8]
                .split("transcript_id \"")
                .nth(1)
                .unwrap()
                .split('"')
                .next()
                .unwrap();
            g += &format!("toy\tx\texon\t{}\t{}\t.\t+\t.\tID=e;Parent={id}\n", p[3], p[4]);
            tx.entry(id.into())
                .or_default()
                .push((p[3].parse::<u32>().unwrap() - 1, p[4].parse().unwrap()));
        }
    }
    for (id, ex) in &tx {
        let (s, e) = (ex[0].0, ex[ex.len() - 1].1);
        let sizes: Vec<String> = ex.iter().map(|x| (x.1 - x.0).to_string()).collect();
        let starts: Vec<String> = ex.iter().map(|x| (x.0 - s).to_string()).collect();
        b += &format!(
            "toy\t{s}\t{e}\t{id}\t0\t+\t{s}\t{e}\t0\t{}\t{},\t{},\n",
            ex.len(),
            sizes.join(","),
            starts.join(",")
        );
    }
    std::fs::write(&gff, g).unwrap();
    std::fs::write(&bed, b).unwrap();
    let sorted = |s: &str| {
        let mut v: Vec<String> = s.lines().map(String::from).collect();
        v.sort();
        v
    };
    for (file, format) in [(&gff, "gff3"), (&bed, "bed12")] {
        let (ok, out, err) = map(&[o("-a"), file.as_os_str()]);
        assert!(ok, "{err}");
        assert!(err.contains(&format!("({format})")), "{err}");
        assert_eq!(sorted(&out), sorted(&paf), "{format}");
    }

    // Each mapping command refuses the other kind of index and names the right command.
    let (ok, _, err) = spearfish(&[o("map"), idx.as_os_str(), reads.as_os_str()]);
    assert!(
        !ok && err.contains("is a genome index") && err.contains("map-spliced"),
        "{err}"
    );
    let tidx = dir.join("toy.sfi");
    let (ok, _, err) = spearfish(&[o("index"), reads.as_os_str(), o("-o"), tidx.as_os_str()]);
    assert!(ok, "{err}");
    let (ok, _, err) = spearfish(&[o("map-spliced"), tidx.as_os_str(), reads.as_os_str()]);
    assert!(
        !ok && err.contains("is a transcriptome index") && err.contains("spearfish map"),
        "{err}"
    );
    let (ok, _, err) = spearfish(&[
        o("index"),
        o("--genome"),
        genome.as_os_str(),
        o("-o"),
        o("x"),
        o("--infer-genes"),
    ]);
    assert!(!ok && err.contains("cannot be used with"), "{err}");
    // An annotation of another assembly stops the run.
    let other = dir.join("other.gtf");
    std::fs::write(&other, "chr1\tx\texon\t10\t50\t.\t+\t.\ttranscript_id \"t1\";\n").unwrap();
    let (ok, _, err) = map(&[o("-a"), other.as_os_str()]);
    assert!(
        !ok && err.contains("no transcript is on a sequence of the index"),
        "{err}"
    );
}
