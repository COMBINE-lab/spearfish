# Contributing

## Build and check

Before every push, run what CI runs:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
python3 fixtures/toy/make_toy.py && git diff --exit-code fixtures/toy
```

CI also runs these on arm64 Linux and Apple Silicon, runs `spearfish doctor`, and checks that
the workspace builds with the minimum Rust version (`rust-version` in `Cargo.toml`).

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
every level the machine supports against a scalar reference.

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
but they are in `crates/spearfish-cli/src`. Environment variables:

| Variable | Effect |
| --- | --- |
| `SPEARFISH_STATS=1` | Count internal events while mapping; printed at the end and stored in the run record. |
| `REFINE_STATS=1` | Statistics of the experimental refine chainer. |
| `SPEARFISH_GIT` | At build time: the revision to record when building outside a git checkout. |

## Commits

The repository is public: commit messages describe changes and relative effects ("fewer
instructions per read"), not unpublished benchmark datasets or numbers.
