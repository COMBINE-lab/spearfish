<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="website/src/assets/spearfish-logo-dark.svg">
    <img src="website/src/assets/spearfish-logo.svg" alt="spearfish" width="420">
  </picture>
</p>

# spearfish

Fast long-read RNA-seq mapping to a transcriptome with locally consistent parsing (LCP).

For each read (ONT cDNA, ONT direct RNA, PacBio Kinnex/Iso-Seq), spearfish performs *delta
mapping*: it reports every compatible transcript, the read's position on it, and its score
difference from the best candidate, for downstream EM quantification. LCP cores are collapsed
into a transcriptome *core graph* that factors out splicing and paralog redundancy, chaining
computes exact per-transcript chain scores while sharing work across isoforms, and base-level
alignment is restricted to the read segments where tied isoforms differ.

Documentation: the Astro/Starlight site in `website/` (`cd website && npm install && npm run dev`).

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
