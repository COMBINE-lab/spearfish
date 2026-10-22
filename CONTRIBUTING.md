# Contributing

## Build and check

Before every push, run what CI runs:

```sh
cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
python3 fixtures/toy/make_toy.py && git diff --exit-code fixtures/toy
```

CI also runs these on arm64 Linux and Apple Silicon, runs `spearfish doctor`, and checks that
the workspace builds with the minimum Rust version (`rust-version` in `Cargo.toml`).

## Profiling

Release builds carry no debug info. For profilers, build with line tables, and with frame
pointers for call graphs (`perf record --call-graph fp`):

```sh
RUSTFLAGS="-C force-frame-pointers=yes" cargo build --profile profiling -p spearfish-cli
```

The binary is `target/profiling/spearfish`, optimized as a release build is.

## Output must not change by accident

Changes that are not meant to change mappings (speedups, refactors, interface work) must leave
the output byte-identical. Checks, from cheapest:

- `cargo test` includes `spearfish doctor`'s self-test, which maps reads from the toy
  transcriptome with every preset and compares a digest of the hits with a reference
  (`EXPECTED` in `crates/spearfish-cli/src/doctor_cmd.rs`). Update the reference only together
  with a deliberate change to mapping, and say so in the commit message.
- On real data, compare sorted TSV output of the old and new builds, and
  `spearfish inspect --sections` of indexes they build (the per-section checksums must agree).

The SIMD kernels must agree exactly with the scalar code; each kernel's property test checks
every level the machine supports against a scalar reference. CI runs the tests on each
runner's own level and again with `SPEARFISH_SIMD` forcing scalar and avx2 (x86-64) or scalar
(Apple Silicon). Hosted x86-64 runners may lack AVX-512, so CI may not run the AVX-512
kernels at all (each run's summary shows the level it used): run `cargo test` on an AVX-512
machine after changing them. The levels live in
`lcp_core::simd`; adding one (e.g. NEON) means a `Level` variant, its detection, and an arm in
each kernel's dispatch.

## Adding a mapping setting

1. Add the parameter to its struct in `crates/map` (with its default).
2. Add an optional field to `MapSettings` (`crates/spearfish/src/settings.rs`), apply it in
   `resolve()`, and add its name to `MapSettings::FIELDS` (a test fails otherwise).
3. Add a flag in `crates/spearfish-cli/src/map_args.rs`, with `[default: ...]` in its help;
   tests check that every setting has a flag and that the documented defaults are the real
   ones. Experimental settings get `hide = true`.
4. Note it in `CHANGELOG.md`, and in the website's mapping guide if users should know it.

## Index format

The index is a flat container of named sections (`crates/index/src/store.rs`). Adding a
section that older readers can ignore needs no version change; changing how an existing
section is laid out or interpreted needs `FORMAT_VERSION` bumped. Changing how the LCP parser
cuts sequences needs `LCP_PARSER_VERSION` bumped, so that indexes built before are rejected
instead of silently mismatching.

## Developer switches

Hidden commands (`em`, `simulate`, `lcp-stats`, `locality`, `survival`, `halo`) and hidden
mapping flags exist for method development; `spearfish <command> --help` does not list them,
but they are in `crates/spearfish-cli/src`. Each mapping flag sets the `MapSettings` field of
its name (`--no-x` turns `x` off; `--beam` goes with `--chainer`), documented there, and
`spearfish map --print-params` shows the values a run would use.

Renamed flags keep their old names as hidden aliases: `--delta-rep`, `--f-keep`, `--f-tie`,
`--f-max`, `--h-max`, `--g-max`, `--alpha`, `--window`, `--beta`, `--no-merge`, `--kappa`,
`--lambda`, `--affine-within`, `--unit-tiebreak`, `--whole-read-tiebreak`, `--orient-edit`,
`--flank-min`, `--no-tie-clip`, `--tie-band`, `--tie-pad` and `--tie-pad-div` for `map` (a
test in `map_args.rs` maps each to its current name), `--f-max` (now `--mask-genes`) for
`index`, and `index-stats` for `inspect`. Benchmark scripts still use `--no-merge`,
`--unit-tiebreak`, `--whole-read-tiebreak`, `--no-tie-clip` and `--orient-edit`.

Environment variables:

| Variable | Effect |
| --- | --- |
| `SPEARFISH_SIMD=scalar\|avx2\|avx512\|neon` | Use at most this SIMD level (to test or time a lower one; a level the CPU lacks, such as an x86-64 one on arm64, falls back with a warning). |
| `SPEARFISH_STATS=1` | Count internal events while mapping; printed at the end and stored in the run record. |
| `REFINE_STATS=1` | Statistics of the experimental refine chainer. |
| `SPEARFISH_GIT` | At build time: the revision to record when building outside a git checkout. |

## Commits

The repository is public: commit messages describe changes and relative effects ("fewer
instructions per read"), not unpublished benchmark datasets or numbers.
