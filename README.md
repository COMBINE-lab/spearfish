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

Documentation: https://combine-lab.github.io/spearfish/ (source in `website/`; preview locally
with `cd website && npm install && npm run dev`).

**Status:** pre-release. Interfaces may still change; see [CHANGELOG.md](CHANGELOG.md).

## Install

Rust 1.89 or newer ([rustup](https://rustup.rs)):

```sh
cargo install --locked --path crates/spearfish-cli
spearfish doctor          # build, CPU, SIMD kernels, and a mapping self-test
```

spearfish runs on Linux and macOS, x86-64 and arm64. Its alignment and parsing kernels pick
AVX-512, AVX2 or scalar code at run time on x86-64, and NEON or scalar code on arm64; every
choice gives the same output.

## Use

```sh
spearfish index transcripts.fa.gz -o transcripts.sfi             # also writes transcripts.sfi.json
spearfish map transcripts.sfi reads.fq.gz --preset ont -o reads.tsv   # and/or --paf, --rad
spearfish quant transcripts.sfi reads.fq.gz --preset ont -o sample     # abundances via oarfish
```

| Command | Does |
| --- | --- |
| `index` | Build the index from a transcript FASTA (gene IDs from headers, `--t2g`, or `--infer-genes`). |
| `map` | Map reads: TSV (`-o`), PAF (`--paf`) and/or RAD (`--rad`) output, plus a run record. |
| `quant` | Map to RAD and quantify with [oarfish](https://github.com/COMBINE-lab/oarfish). |
| `inspect` | Show what an index holds, what built it, and whether this build can use it. |
| `doctor` | Report the build and machine and run a self-test; `--index` also checks an index. |
| `completions`, `man` | Shell completions and man pages. |

Presets (`--preset ont|drna|hifi`) set every default; any setting can be changed with its
flag, and `spearfish map --print-params` prints the parameters a run would use. `-h` lists
the common flags, `--help` all of them.

Every index records how it was built (`<index>.json`, and `spearfish inspect`), and every
mapping run writes `<output>.meta.json`: the build, command line, index, settings, every
parameter used, and run statistics.

## As a library

```rust
use std::sync::Arc;
use spearfish::{Index, MapSettings, Preset};

let index = Arc::new(Index::open("transcripts.sfi".as_ref())?);
let mapper = MapSettings::builder(Preset::Hifi).extend_ends(true).mapper(index.clone())?;
for hit in mapper.map(b"ACGT...").hits {
    println!("{} {} {}", index.tx_name(hit.tx), hit.score, hit.reverse);
}
```

## Layout

| Path | Contents |
| --- | --- |
| `crates/lcp-core` | LCP parser (lcptools-compatible boundaries, 64-bit content labels, hierarchy) |
| `crates/lcp-oracle-sys` | test-only build of [lcptools](https://github.com/BilkentCompGen/lcptools) (BSD-3-Clause) |
| `crates/seeds` | seed-scheme trait: LCP, and minimizers, syncmers and strobemers for ablations |
| `crates/index` | core-graph index and its on-disk container |
| `crates/map` | read mapping: anchoring, chaining, tie-break, rescue, span extension |
| `crates/quant` | spearfish's own EM (a baseline for method development) |
| `crates/rad` | `bulk_long_read` RAD records, shared with oarfish |
| `crates/spearfish` | the library: index, settings, mapper, provenance |
| `crates/spearfish-cli` | the `spearfish` binary |
| `fixtures/toy` | toy transcriptome with an exon skip, an alternative 3′ site and a paralog pair |

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

BSD-3-Clause.
