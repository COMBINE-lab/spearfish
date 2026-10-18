---
title: Installation
description: Build spearfish from source, and check the build.
---

spearfish is written in Rust. Until the first release, it is installed from source.

## Requirements

- Rust 1.89 or newer, installed with [rustup](https://rustup.rs).
- Linux or macOS, on x86-64 or arm64 (Apple Silicon). The alignment and parsing kernels pick
  AVX-512, AVX2 or scalar code at run time, and every choice gives exactly the same output.

## Build

```sh
git clone https://github.com/COMBINE-lab/spearfish.git
cd spearfish
cargo install --locked --path crates/spearfish-cli
```

This puts the `spearfish` binary in `~/.cargo/bin`. To build without installing:

```sh
cargo build --release
./target/release/spearfish --help
```

## Check the build

```sh
spearfish doctor
```

`doctor` prints the version and git revision, the target, the CPUs available and the SIMD
kernels in use, then maps reads from a built-in toy transcriptome with every preset and checks
the result against a reference. With `--index transcripts.sfi` it also checks that an index
loads and suits this build. `--json` prints the same report as JSON, to attach to a bug report.

The full test suite is `cargo test --workspace --release`.

## Shell completions and man pages

```sh
spearfish completions bash > ~/.local/share/bash-completion/completions/spearfish
spearfish completions zsh  > ~/.zfunc/_spearfish        # also fish, elvish, powershell
spearfish man -o ~/.local/share/man/man1                 # spearfish.1, spearfish-map.1, ...
```

## For quantification

`spearfish quant` runs [oarfish](https://github.com/COMBINE-lab/oarfish), which must be a build
with spearfish support (its `--help` lists `spearfish` as a `--mapper`). spearfish looks for it
on the `PATH`, in `$OARFISH`, or where `--oarfish` points.
