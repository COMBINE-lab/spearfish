<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="website/src/assets/spearfish-logo-dark.svg">
    <img src="website/src/assets/spearfish-logo.svg" alt="spearfish" width="420">
  </picture>
</p>

# spearfish

Fast long-read RNA-seq mapping with locally consistent parsing (LCP): to a transcriptome, or
spliced to a genome.

For each read (ONT cDNA, ONT direct RNA, PacBio Kinnex/Iso-Seq), spearfish performs *delta
mapping*: it reports every compatible transcript, the read's position on it, and its score
difference from the best candidate, for downstream EM quantification. LCP cores are collapsed
into a transcriptome *core graph* that factors out splicing and paralog redundancy, chaining
computes exact per-transcript chain scores while sharing work across isoforms, and base-level
alignment is restricted to the read segments where tied isoforms differ.

Documentation: https://combine-lab.github.io/spearfish/ (source in `website/`; preview locally
with `cd website && npm install && npm run dev`).

**Status:** 0.3.0. Each hit now carries a quantification score next to its mapping score, and
the library and RAD records gained fields for it: see [CHANGELOG.md](CHANGELOG.md). Interfaces may
change before 1.0.

## Install

From [bioconda](https://bioconda.github.io/recipes/spearfish/README.html) (Linux and macOS,
x86-64 and arm64):

```sh
conda install -c bioconda -c conda-forge spearfish
```

Or with Rust 1.91 or newer ([rustup](https://rustup.rs)), from crates.io or from a checkout:

```sh
cargo install --locked spearfish-cli              # or: cargo install --locked --path crates/spearfish-cli
```

Then check the build:

```sh
spearfish doctor          # build, CPU, SIMD kernels, and a mapping self-test
```

spearfish runs on Linux and macOS, x86-64 and arm64. Its alignment and parsing kernels pick
AVX-512, AVX2 or scalar code at run time on x86-64, and NEON or scalar code on arm64; every
choice gives the same output.

## Use

```sh
spearfish index transcripts.fa.gz -o transcripts.sfi             # also writes transcripts.sfi.json
spearfish map transcripts.sfi reads.fq.gz --preset ont -o reads.tsv   # and/or --paf, --rad
```

| Command | Does |
| --- | --- |
| `index` | Build the index from a transcript FASTA (gene IDs from headers, `--t2g`, or `--infer-genes`), or from a genome FASTA with `--genome`. |
| `map` | Map reads to transcripts: TSV (`-o`), PAF (`--paf`) and/or RAD (`--rad`) output, plus a run record. |
| `map-spliced` | Map reads, spliced, to a genome index: PAF (`--paf`), SAM (`--sam`) and/or BAM (`--bam`) output, plus a run record. |
| `inspect` | Show what an index holds, what built it, and whether this build can use it. |
| `doctor` | Report the build and machine and run a self-test; `--index` also checks an index. |
| `completions`, `man` | Shell completions and man pages. |

Presets (`--preset ont|drna|hifi`) set every default (`ont-fast`, `drna-fast` and `hifi-fast`
map as 0.1.0 did: faster, less accurate); any setting can be changed with its flag, and
`spearfish map --print-params` prints the parameters a run would use. `-h` lists the common
flags, `--help` all of them.

Every index records how it was built (`<index>.json`, and `spearfish inspect`), and every
mapping run writes `<output>.meta.json`: the build, command line, index, settings, every
parameter used, and run statistics.

## Spliced mapping to a genome

Reads can also be mapped to a genome, with introns. `index --genome` builds the genome index and
`map-spliced` maps to it.

```sh
spearfish index --genome genome.fa.gz -o genome.sfg
spearfish map-spliced genome.sfg reads.fq.gz --preset ont --bam reads.bam   # or --sam, --paf (default: PAF on stdout)
spearfish map-spliced genome.sfg reads.fq.gz -a annotation.gtf.gz --bam reads.bam   # known transcripts are optional
```

The index needs no annotation and stores none. Known transcripts (`-a`: GTF, GFF3 or BED12; or a BED list of introns, as minimap2 reads and writes)
are given when mapping; junctions outside them are still found. By default an alignment is
its exon blocks, with introns as `N` in the CIGAR, which is what quantification and transcript
discovery read; `--base-level edits|affine` adds a base-level alignment of each block. See the
[guide](https://combine-lab.github.io/spearfish/guide/genome/).

## Quantification

spearfish maps; [oarfish](https://github.com/COMBINE-lab/oarfish) turns its mappings into
transcript abundances. With an oarfish that has spearfish support (0.11.0 or newer, e.g.
`conda install -c bioconda -c conda-forge oarfish`; its `--help` lists `spearfish` as a `--mapper`), either let oarfish map with spearfish in process:

```sh
oarfish --reads reads.fq.gz --annotated transcripts.fa.gz --seq-tech ont-cdna --mapper spearfish -o sample
```

or map once to a RAD file and quantify that (to keep the mappings, or quantify them several ways):

```sh
spearfish map transcripts.sfi reads.fq.gz --preset ont --rad sample.rad
oarfish -a sample.rad -o sample
```

oarfish 0.12 maps in process with spearfish 0.2. It also reads RAD files written by 0.3.0, but
weighs alignments by the mapping score; using the quantification score (`qscore`) needs an oarfish
release that embeds spearfish 0.3.

## As a library

`spearfish = "0.3"` in `Cargo.toml`:

```rust
use std::sync::Arc;
use spearfish::{Index, MapSettings, Preset};

let index = Arc::new(Index::open("transcripts.sfi".as_ref())?);
let mapper = MapSettings::builder(Preset::Hifi).extend_ends(true).mapper(index.clone())?;
for hit in mapper.map(b"ACGT...").hits {
    println!("{} {} {}", index.tx_name(hit.tx), hit.score, hit.reverse);
}
```

Spliced mapping to a genome is the `genome` module:

```rust
use spearfish::{Preset, genome::{GenomeIndex, SpliceOptions, SplicedMapper}};

let mut index = GenomeIndex::open("genome.sfg".as_ref())?;
index.load_annotation("annotation.gtf".as_ref())?;      // optional
let mapper = SplicedMapper::new(index, &SpliceOptions::preset(Preset::Ont));
if let Some(best) = mapper.map(b"ACGT...").alignments.first() {
    for exon in best.blocks() {
        println!("{} {}-{}", mapper.index().seq_name(best.seq), exon.target_start, exon.target_end);
    }
}
```

## Layout

| Path (crate) | Contents |
| --- | --- |
| `crates/lcp-core` (`spearfish-lcp`) | LCP parser (lcptools-compatible boundaries, 64-bit content labels, hierarchy) |
| `crates/lcp-oracle-sys` | test-only build of [lcptools](https://github.com/BilkentCompGen/lcptools) (BSD-3-Clause; not published) |
| `crates/seeds` (`spearfish-seeds`) | seed-scheme trait: LCP, and minimizers, syncmers and strobemers for ablations |
| `crates/index` (`spearfish-index`) | core-graph index and its on-disk container |
| `crates/map` (`spearfish-map`) | read mapping: anchoring, chaining, tie-break, rescue, span extension |
| `crates/splice` (`spearfish-splice`) | spliced mapping to a genome: genome index, junction placement, PAF/SAM/BAM records |
| `crates/rad` (`long-read-rad`) | `bulk_long_read` RAD records, shared with oarfish |
| `crates/spearfish` (`spearfish`) | the library: index, settings, mapper, output formats, provenance |
| `crates/spearfish-cli` (`spearfish-cli`) | the `spearfish` binary |
| `crates/nudgy` (`nudgy`) | the allocator behind `spearfish map`'s arenas (`--bump-alloc`): per-thread bump arenas in front of any global allocator (independent of spearfish, with its own version) |
| `crates/py-spearfish` | the Python package (`spearfish` on PyPI; built with maturin) |
| `fixtures/toy` | toy transcriptome with an exon skip, an alternative 3′ site and a paralog pair, and its genome and annotation |

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

BSD-3-Clause.
