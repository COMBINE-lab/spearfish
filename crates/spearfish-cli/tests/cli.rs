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
    for p in ["ont", "drna", "hifi"] {
        assert_eq!(d["self_test"][p]["pass"], true, "{d}");
    }
    assert!(ok(&["completions", "bash"]).contains("spearfish"));
    let man = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("man");
    ok(&["man", "-q", "-o", man.to_str().unwrap()]);
    for page in [
        "spearfish.1",
        "spearfish-map.1",
        "spearfish-index.1",
        "spearfish-quant.1",
    ] {
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
    assert_eq!(h.meta.score_denom, 1.0);
    assert_eq!(h.counts.reads, 200);
    assert_eq!(h.counts.mapped, 200);
    assert_eq!(h.ref_names.len(), h.ref_lengths.len());
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
