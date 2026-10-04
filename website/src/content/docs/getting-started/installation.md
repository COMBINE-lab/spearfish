---
title: Installation
description: Build spearfish from source.
---

spearfish is written in Rust. Until the first release, it is installed from source.

## Requirements

- A recent stable Rust toolchain (1.88 or newer), installed with [rustup](https://rustup.rs).
- spearfish is developed and tested on Linux x86-64. Its alignment kernels pick SIMD
  instructions (AVX-512, AVX2, or NEON on ARM) at run time; other platforms should build but are
  untested.

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
cargo test --workspace --release
```
