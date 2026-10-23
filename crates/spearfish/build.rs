//! Record the git revision the library is built from (`SPEARFISH_GIT`), for provenance: as
//! set in the environment (packagers building outside a git checkout set it), else the commit
//! a published crate was packaged from (`.cargo_vcs_info.json`, which `cargo package` writes),
//! else `git describe` of the checkout, else `unknown`.

use std::process::Command;

fn git(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    out.status
        .success()
        .then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
        .filter(|s| !s.is_empty())
}

/// The commit a published crate was packaged from, from the `.cargo_vcs_info.json` that
/// `cargo package` writes next to its manifest (`{"git": {"sha1": "...", "dirty": true}, ...}`),
/// abbreviated as `git describe` abbreviates it.
fn packaged_commit(crate_dir: &str) -> Option<String> {
    let path = format!("{crate_dir}/.cargo_vcs_info.json");
    let info = std::fs::read_to_string(&path).ok()?;
    println!("cargo:rerun-if-changed={path}");
    let rest = &info[info.find("\"sha1\"")? + 6..];
    let rest = &rest[rest.find('"')? + 1..];
    let sha = &rest[..rest.find('"')?];
    if sha.len() < 12 || !sha.chars().all(|c| c.is_ascii_hexdigit()) {
        return None;
    }
    let dirty = info.contains("\"dirty\": true") || info.contains("\"dirty\":true");
    Some(format!("{}{}", &sha[..12], if dirty { "-dirty" } else { "" }))
}

fn main() {
    println!("cargo:rerun-if-env-changed=SPEARFISH_GIT");
    let crate_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_default();
    // A packaged crate's own record comes before git: its source may sit inside some other
    // repository (a vendored dependency), whose revision is not the crate's.
    let rev = std::env::var("SPEARFISH_GIT")
        .ok()
        .or_else(|| packaged_commit(&crate_dir))
        .or_else(|| git(&["describe", "--always", "--dirty", "--abbrev=12"]))
        .unwrap_or_else(|| "unknown".into());
    println!("cargo:rustc-env=SPEARFISH_GIT={rev}");
    // Rebuild when HEAD moves or the working tree's index changes (commits, checkouts), and
    // when the workspace's sources change, which may make the tree dirty without touching
    // the git index. (Only existing paths: Cargo reruns the script on every build for a
    // missing one.)
    let watch = |p: String| {
        if std::path::Path::new(&p).exists() {
            println!("cargo:rerun-if-changed={p}");
        }
    };
    for p in ["crates", "Cargo.toml", "Cargo.lock"] {
        watch(format!("{crate_dir}/../../{p}"));
    }
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
