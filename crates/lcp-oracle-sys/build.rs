use std::path::PathBuf;

fn main() {
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    // Same configuration as `make LABEL=64 POS=32 DCT=1`.
    let cfg = std::fs::read_to_string("lcptools/config.h.in")
        .unwrap()
        .replace("@LABEL@", "64")
        .replace("@POS@", "32")
        .replace("@DCT@", "1");
    std::fs::write(out.join("config.h"), cfg).unwrap();

    cc::Build::new()
        .files(["lcptools/encoding.c", "lcptools/core.c", "lcptools/lps.c"])
        .include(&out)
        .include("lcptools")
        .opt_level(2)
        .warnings(false)
        .compile("lcptools");
    println!("cargo:rustc-link-lib=m");
    println!("cargo:rerun-if-changed=lcptools");
}
