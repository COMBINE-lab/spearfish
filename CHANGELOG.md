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

### Other

- gzip input decodes with zlib-rs.
- Index writing streams sections to disk instead of building the file in memory.
- The workspace declares its minimum Rust version, 1.89.
- CI runs on x86-64 Linux, arm64 Linux and Apple Silicon, and checks the minimum Rust version.
