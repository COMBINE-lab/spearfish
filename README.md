# spearfish

Fast long-read RNA-seq mapping to a transcriptome with locally consistent parsing (LCP).

For each read (ONT cDNA, ONT direct RNA, PacBio Kinnex/Iso-Seq), spearfish reports a scored
set of compatible transcripts for downstream EM quantification, without base-level
alignment. LCP cores are collapsed into a transcriptome *core graph* that factors out
splicing and paralog redundancy. A single multi-transcript chaining DP with lazy color
refinement computes exact per-transcript chain scores while sharing work across isoforms.

**Status:** early development. See `PLAN.md` for the schedule and `DECISIONS.md` for
design decisions and the evidence behind them.

## Layout

| Path | Contents |
| --- | --- |
| `crates/lcp-core` | LCP parser (lcptools-compatible boundaries, 64-bit content labels, hierarchy) |
| `crates/lcp-oracle-sys` | test-only build of [lcptools](https://github.com/BilkentCompGen/lcptools) (BSD-3-Clause) |
| `crates/seeds` | seed-scheme trait: LCP, and minimizers, syncmers and strobemers for ablations |
| `crates/index`, `map`, `quant` | core-graph index, read mapping and chaining, EM (in progress) |
| `crates/spearfish` | the `spearfish` binary |
| `fixtures/toy` | toy transcriptome with an exon skip, an alternative 3′ site and a paralog pair |
| `bench/` | data acquisition and benchmark workflow |

## Building

```sh
cargo build --release
cargo test --workspace
```

## License

BSD-3-Clause.
