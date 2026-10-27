---
title: Installation
description: Install spearfish, and check the build.
---

spearfish installs from bioconda, or with Cargo from crates.io or from a checkout.

## Requirements

- Linux or macOS, on x86-64 or arm64 (Apple Silicon). The alignment and parsing kernels pick
  AVX-512, AVX2 or scalar code at run time on x86-64, and NEON or scalar code on arm64, and
  every choice gives exactly the same output.

## Install with conda

```sh
conda install -c bioconda -c conda-forge spearfish
```

[bioconda](https://bioconda.github.io/recipes/spearfish/README.html) has builds for Linux and macOS,
on x86-64 and arm64. They are generic builds: the SIMD kernels are still chosen at run time, so a
conda build is as fast as one compiled for your CPU.

## Install with Cargo

With Rust 1.89 or newer, installed with [rustup](https://rustup.rs):

```sh
cargo install --locked spearfish-cli
```

This puts the `spearfish` binary in `~/.cargo/bin`. From a checkout instead:

```sh
git clone https://github.com/COMBINE-lab/spearfish.git
cd spearfish
cargo install --locked --path crates/spearfish-cli
```

To build without installing:

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
the result against a reference. `SPEARFISH_SIMD=scalar|avx2|avx512|neon` makes spearfish use at
most that level, for example to compare speeds. With `--index transcripts.sfi` it also checks that an index
loads and suits this build. `--json` prints the same report as JSON, to attach to a bug report.

The full test suite is `cargo test --workspace --release`.

## Shell completions and man pages

```sh
spearfish completions bash > ~/.local/share/bash-completion/completions/spearfish
spearfish completions zsh  > ~/.zfunc/_spearfish        # also fish, elvish, powershell
spearfish man -o ~/.local/share/man/man1                 # spearfish.1, spearfish-map.1, ...
```

## For quantification

Abundances come from [oarfish](https://github.com/COMBINE-lab/oarfish), installed separately:
version 0.11.0 or newer has spearfish support (its `--help` lists `spearfish` as a `--mapper`).

```sh
conda install -c bioconda -c conda-forge oarfish
```

See [Quantification](/spearfish/guide/quantification/).
