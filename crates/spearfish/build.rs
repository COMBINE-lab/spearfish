//! Record the git revision the library is built from (`SPEARFISH_GIT`), for provenance.
//! Packagers building outside a git checkout can set `SPEARFISH_GIT` themselves.

use std::process::Command;

fn git(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    out.status
        .success()
        .then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
        .filter(|s| !s.is_empty())
}

fn main() {
    println!("cargo:rerun-if-env-changed=SPEARFISH_GIT");
    let rev = std::env::var("SPEARFISH_GIT")
        .ok()
        .or_else(|| git(&["describe", "--always", "--dirty", "--abbrev=12"]))
        .unwrap_or_else(|| "unknown".into());
    println!("cargo:rustc-env=SPEARFISH_GIT={rev}");
    // Rebuild when HEAD moves or the working tree's index changes (commits, checkouts).
    // (Only existing files: Cargo reruns the script on every build for a missing one.)
    let watch = |p: String| {
        if std::path::Path::new(&p).exists() {
            println!("cargo:rerun-if-changed={p}");
        }
    };
    if let Some(dir) = git(&["rev-parse", "--absolute-git-dir"]) {
        watch(format!("{dir}/HEAD"));
        watch(format!("{dir}/index"));
        if let (Some(r), Some(common)) = (
            git(&["symbolic-ref", "-q", "HEAD"]),
            git(&["rev-parse", "--path-format=absolute", "--git-common-dir"]),
        ) {
            watch(format!("{common}/{r}"));
            watch(format!("{common}/packed-refs"));
        }
    }
    println!(
        "cargo:rustc-env=SPEARFISH_TARGET={}",
        std::env::var("TARGET").unwrap_or_default()
    );
    println!(
        "cargo:rustc-env=SPEARFISH_PROFILE={}",
        std::env::var("PROFILE").unwrap_or_default()
    );
}
