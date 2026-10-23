# Changelog

## Unreleased

### Command line

- Every mapping setting is an optional flag applied over the preset: a flag left out keeps
  the preset's value. On/off settings use `--x` / `--no-x` (`--extend-ends`,
  `--no-extend-ends`, `--rescue`, `--no-rescue`, `--tail-trim`, `--no-tail-trim`,
  `--tie-group`, `--no-tie-group`, `--clip-ends`, `--no-clip-ends`, `--raw-ends`,
  `--no-raw-ends`, `--end-diagonal`, `--no-end-diagonal`). The forms `--extend-ends true|false`
  and `--raw-ends true|false` are gone; `spearfish index --hpc false` is now `--no-hpc`.
- Flags renamed for consistency; the old names still work as hidden aliases: `--delta-rep` →
  `--report-within`, `--f-keep`/`--f-tie`/`--f-max` → `--keep-genes`/`--keep-genes-tie`/
  `--max-genes`, `--h-max` → `--max-seed-genes`, `--alpha`/`--window` → `--chain-alpha`/
  `--chain-window`, `--kappa`/`--lambda`/`--affine-within` → `--tie-kappa`/`--tie-lambda`/
  `--tie-affine-within`, `--tie-band`/`--tie-pad`/`--tie-pad-div` → `--band`/`--pad`/
  `--pad-div`, `--no-tie-clip` → `--no-clip-ends`; `spearfish index --f-max` →
  `--mask-genes`. `--rescue-full` is now `--rescue-window full`.
- `--xdrop`, `--end-bonus`, `--band` and `--pad` apply to span extension as well as to the
  tie-break (their defaults were already equal).
- Settings that had no flag now have one (`--tie-margin-rel`, `--tie-cap`, `--tail-run`,
  `--min-coverage`, `--min-gene-frac`, `--seed-span-cap`, `--chain-max-gap`,
  `--chain-max-overlap`, `--rescue-max-tx-len`, `--rescue-other-always`), and the preset-only
  ones can be turned both ways (e.g. `--rescue` for hifi, `--no-tie-group`).
- `-h` shows the common flags, `--help` every flag and the settings that differ between
  presets. `spearfish map --print-params` prints the parameters a run would use.
- `spearfish map --paf FILE`: PAF output.
- `spearfish quant` now maps reads to RAD and quantifies them with oarfish. spearfish's own EM
  over a mapping TSV is `spearfish em` (hidden; a baseline for method development).
- New: `spearfish inspect` (`index-stats` is an alias), `spearfish doctor`,
  `spearfish completions`, `spearfish man`. Global `-q` / `-v`.
- The research commands (`simulate`, `lcp-stats`, `locality`, `survival`, `halo`) and the
  experimental mapping flags are hidden from help; they still work.
- Errors from a reader closing the pipe (`| head`) no longer fail the run.
- Reads and transcripts can be plain or gzip- (BGZF included), bzip2- or zstd-compressed.
  xz is no longer read: it needed a C library (liblzma) for a format reads rarely come in, and
  an xz file now fails with a message saying so.

### Provenance

- Indexes carry a manifest (build, command line, options, input checksums, where genes came
  from, construction statistics), also written next to the index as `<index>.json`. Indexes
  without one still load.
- Mapping checks that the index was built with this build's LCP parser version and seed, and
  asks for a rebuild otherwise (these were recorded but not checked).
- Every mapping run writes a run record, `<first output>.meta.json` (or `--meta FILE`): build,
  command line, index and manifest, settings changed, every parameter used, read counts,
  throughput, CPU time and peak memory.
- RAD files record the parameters as JSON, the spearfish version with its git revision, and an
  `index_id` tag.
- `--version` shows the git revision, target and SIMD level.
- The index id now covers the gene grouping: indexes built from the same transcripts with
  different genes (headers, `--t2g` or `--infer-genes`) used to share an id. The input hash
  also length-prefixes names and sequences, so moving bases between a name and its sequence
  changes it. Existing indexes keep their recorded id; a rebuild gets a new one.

### Library

- `MapSettings`: a preset plus the settings changed, with a builder
  (`MapSettings::builder(Preset::Hifi).extend_ends(true).mapper(index)`) and `resolve()` to
  the complete `MapOptions`. `IndexOptions` has a builder and now covers gene inference.
- `SCORE_DENOM` (1) is the one score scale; `MapOptions` and `Mapper` report it. The docs that
  said D = 2 were out of date.
- `MapOptions` records its preset and serializes; `Mapper::options()`, `Index::manifest()`,
  `Index::id()`. `Strands`, `TieBreak`, `MapParams` and the other parameter types are
  re-exported. The existing API (`MapOptions::preset`, `Mapper::new`, `Index::build`) is
  unchanged.
- Removed: `MapOptions::tiebreak_mut` and `chain_mut` (use the public `params` fields), and
  the `spearfish::quant` re-export (the EM crate serves `spearfish em` only).

### Performance

- AVX2 kernels for the affine alignment rows, the batched edit distances and the LCP scan:
  CPUs with AVX2 but not AVX-512 map 16 to 34% faster, close to AVX-512 speed.
- SIMD dispatch happens once per DP pass or parse call instead of per row or per 64-base
  block, and the whole LCP parse is compiled per level; a few percent fewer instructions at
  every level.
- NEON kernels on aarch64 (Apple Silicon): the affine alignment rows, the batched edit
  distances (four texts in lockstep), the LCP scan of the bases, and the chaining predecessor
  scores with their maximum.
- `SPEARFISH_SIMD` lowers the SIMD level used; `spearfish doctor` reports it.
- `spearfish index` builds about a third faster with about a third less peak memory
  (GENCODE with spike-ins, 32 threads): node keys are deduplicated within genes before the
  global sort, unitig chains are followed in parallel, and the transcript sequences are freed
  once seeded. Building an index in memory (`Index::build`, as oarfish does from a FASTA) lays
  the image out once instead of copying it twice. Index sections are byte-identical.
- Rescue does less work per candidate window: genes come from the strand's anchor selection,
  window seeds probe a table of the read's seed labels instead of being sorted and merged,
  points are sorted as packed keys, and windows are parsed straight from their 2-bit codes,
  for labels and starts only. 11-14% fewer instructions on ONT cDNA and simulated 90%-accuracy
  ONT reads, 6% on direct RNA, 2% on simulated 95%-accuracy ONT reads and under 1% where little
  is rescued (HiFi, PacBio), with identical output.
- The tie-break's Global and Prefix edit distances (direct RNA's segment costs) run in the
  batched SIMD kernels, as infix ones already did, and the affine DP rows write their own band
  sentinel and run whole AVX-512 vectors: 17% fewer instructions on direct RNA and 0.5-3% on
  the other read types, with identical output.
- Chaining does less work per transcript: the DP scores predecessor windows in whole vector
  blocks and finds the best one in the same pass, a transcript whose placements match the
  previous transcript's reuses its chain, and placements are stored relative to their
  transcript's first. 6-13% fewer instructions across read types, with identical output.
- Less work per read outside the alignment and chaining kernels: the read is compressed once
  (the library reuses the mapper's homopolymer runs for query spans), TSV and PAF lines are
  written without `core::fmt`, packed bases are decoded 8 at a time, LCP core ends and the
  first DCT step skip redundant work, and one orientation's hits are sorted by cached keys.
  9-12% fewer instructions (ONT cDNA, direct RNA, PacBio and simulated reads) and 4-8% fewer
  cycles single-threaded, with identical output.

### Fixes

- Near-tied read orientations and tied hits without a shared anchor are compared by a clipped
  alignment that scored the read between the outermost anchors in compressed bases and the read
  ends in raw ones, so whichever hit left more of the read to its ends gained. With raw read
  ends the whole read is now scored in raw bases. Reads from transcripts whose sequence lies
  inside another gene's (pseudogenes in a host's intron, often antisense) picked the wrong gene
  half as often: on the simulated sets, top-1 99.16% -> 99.37% (HiFi) and 98.51% -> 98.67%
  (ONT), reads with a best hit in the wrong gene 0.30% -> 0.09% and 0.40% -> 0.23%.
- `--orient-tail` (default 30 points): of two near-tied read orientations, the one whose raw
  read ends in a poly(A) tail is preferred. Where sequence is shared antisense (SIRV antisense
  isoforms, nested pseudogenes) both orientations align equally well, and the tail is the
  evidence of strand: on ONT PCR-cDNA SIRV reads, top hits on the strand opposite the read's
  tail fell from 6,144 to 2,159. Mapping cost is unchanged (instructions within ±0.5%).
- `--orient-margin-rel` (default 0.35): both read orientations go to the tie-break when the
  weaker one's chain scores within 35% of the stronger one's. Before, only `--tie-margin`
  points of slack applied, so an orientation that happened to seed fewer level-3 cores (LCP
  parses the two strands differently) was dropped before the alignment could compare them.
  Only the orientation test widens; ties within one orientation are unchanged.
- `--rescue-genes-len` (default 350 compressed bases): rescue keeps candidate genes in
  proportion to `350 / read length`, between 1 and 4 times `--max-genes`. A short read
  (pseudogenes and small RNAs are often under 500 nt) seeds few level-3 cores, and its true
  gene was often not among the first `--max-genes` candidates. `0` turns this off.
- Together, on the simulated ONT 95% set, top-1 98.67% -> 98.72% (pseudogene 98.17% -> 98.42%,
  other biotypes 97.62% -> 97.99%), reads with a best hit in the wrong gene 0.23% -> 0.19%; on
  ONT 90%, 95.95% -> 96.05% (pseudogene 94.64% -> 95.07%, other 91.99% -> 93.05%); on HiFi,
  99.37% -> 99.38% (pseudogene 99.17% -> 99.27%). Real-data spike-in accuracy is unchanged
  within ±0.001 Spearman. Mapping cost rises 0.4-1.8% in instructions (largest on ONT cDNA)
  and at most 1.5% in wall time.

- `--orient-edit` and `--unit-tiebreak` (experimental): near-tied read orientations were
  compared by edit distance over the leading orientation's anchored span but over the other
  orientation's whole read, adapters and primers included, so the chain-score leader almost
  always stayed. Both are now scored over the whole read: on the simulated ONT 90% set, top-1
  95.75% -> 96.07% (`--orient-edit`) and 95.74% -> 96.05% (`--unit-tiebreak`). The presets
  don't use this comparison.
- A poly(A) tail after an A-run inside the read (an A-rich 3' UTR, a genomic A-stretch) was
  missed: only the first A-run was looked at, and with more than 80 read bases after it the
  read had no tail. The search now goes on past such runs, for tail trimming, `--tail-bonus`
  and `--orient-tail`. 0.1-0.35% of real reads change mapping; simulated and spike-in accuracy
  are unchanged.
- On indexes built with `--no-hpc`, poly(A) tails were never found, so `--tail-trim`,
  `--tail-bonus` and `--orient-tail` did nothing. On such an index, simulated ONT reads with
  tails and adapters now reach 99.01% top-1 instead of 98.35%, and SIRV E2 Spearman on ONT
  PCR-cDNA rises from 0.706 to 0.758. HPC indexes are unaffected.
- With raw read ends (`ont`), a read end's band was laid out for a window computed from the
  raw read end, while the candidate transcript windows are sized in compressed bases and are
  often longer in raw ones: isoforms that agree near the read end got different bands and
  could not share DP rows. Every candidate now gets the band of the longest window. About 0.1%
  of ONT reads change mapping, accuracy is unchanged, and ONT mapping takes 0.4-2% fewer
  instructions.

- `spearfish map -o /dev/null` (or any output that is a device or pipe) no longer fails at the
  end trying to write `/dev/null.meta.json`: the run record goes next to the first output that
  is a regular file, and is skipped if there is none and `--meta` is not given.

- aarch64 (Apple Silicon, ARM Linux): the affine alignment rows were not computed at all,
  because a configuration attribute removed the scalar code along with the x86 path.

- Settings that have no effect are now also reported for `--no-clip-ends` with the tie-break
  off, and for `--end-bonus`, `--xdrop`, `--band`, `--pad` and `--pad-div` with both the
  tie-break and span extension off (e.g. `--preset hifi --no-tiebreak`).
- Errors about input and output files name the file (a missing FASTA for `spearfish index`
  or `spearfish em`, an output that can't be created). `spearfish index --write-t2g` and
  `--t2g`, and `spearfish em -o`, fail on a bad path before the index build or the EM rather
  than after it.
- `spearfish quant` keeps the RAD file when oarfish fails and names it, so `oarfish -a` can be
  rerun without mapping again; with `-v` the error no longer points to a log file that was not
  written (oarfish logs to stderr then).
- A damaged raw-sequence section made an index load as one without raw sequence (which
  changes mapping output), and sections of mismatched sizes loaded and then panicked while
  mapping. Loading now fails with the section's name.
- Index writing staged the file at the output's name with its extension replaced by `.tmp`, so
  builds to `x.lcp` and `x.mm` at the same time shared one temporary file, and an unrelated
  `x.tmp` was moved away. It now stages to `<output>.partial-<pid>` and removes it on failure.
- Reads that spell T both as T and as U got wrong query coordinates. The mapper reads U as T,
  so a T next to a U is one homopolymer run, but the query span was converted back to raw
  bases with the runs of the read as given. Basecallers write one spelling per read, so
  mapping real reads is unchanged.

### Other

- gzip input decodes with zlib-rs.
- Index writing streams sections to disk instead of building the file in memory.
- The workspace declares its minimum Rust version, 1.89.
- CI runs on x86-64 Linux, arm64 Linux and Apple Silicon, and checks the minimum Rust version.
- Release builds (`cargo install`) no longer embed debug info, so the binary is about 8 times
  smaller; backtraces lose file and line numbers (panic messages keep them). The `profiling`
  profile builds with line tables (see CONTRIBUTING.md).
