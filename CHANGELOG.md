# Changelog

## 0.3.0 (unreleased)

- **A second score per hit, the quantification score (`qscore`), with an end fit for reads
  without a poly(A) tail** (`--end-fit`, 5 points; on in every preset except the `-fast` ones).
  The poly(A) tail is spearfish's evidence for where a read ends and, between near-tied
  orientations, for its strand. Without it (poly(A)-trimmed reads, or reads that lost their
  tail), isoforms that differ only in their 3' end and antisense transcripts overlapping on
  opposite strands tie, and the quantifier tends to give such reads to whichever transcript has
  more unique reads. The quantification score is the mapping score plus `--end-fit` points for
  the candidates whose annotated 3' end the read reaches, when the read has no tail and the
  candidates' ends lie farther apart than `--end-fit-tol` (4 compressed bases). The mapping
  score is unchanged and still decides the ranks, the best set, the primary hit and the
  orientation, so mapping output is the same as in 0.2.1; only the weights a quantifier gives a
  read's hits change. Where a read ends is measured from the alignment of its right end, or,
  in the `drna` preset (whose spans are not extended), from each candidate's chain
  (`--end-fit-chain`; `--no-end-fit-chain` aligns instead, at a higher mapping cost; span
  extension implies alignment). On simulated ONT cDNA and real spike-in data with
  poly(A)-trimmed reads, much of the accuracy lost to trimming comes back; on reads with tails
  results are about the same. Mapping costs well under 1% more. Output: a 15th TSV column
  `qscore`, a PAF tag `qs:f`, an optional RAD alignment tag `qscore` (file tag `quant_scores`),
  `Hit::quant_score` in the library and `Hit.quant_score` / `Mapper.tsv_header` in Python; all
  are left out when the end fit is off (`--end-fit 0`, the `-fast` presets), so output is then
  as before. See the new docs page "Mapping and quantification scores".

- **Read orientation: antisense pairs that are one sequence under the read are compared, not
  tied** (`--orient-same-exact`, `--orient-same-end-cut`, `--orient-tail-after`; on in every
  preset except the `-fast` ones). When a read's two orientations have best hits that are one
  sequence on opposite strands, as with a short antisense isoform contained in a longer one,
  0.2.1 reported the two orientations tied and let only the poly(A) tail separate them. The
  read often holds more: bases one transcript has of its own (a short isoform's first bases),
  and a tail the detector missed. Such a pair is now compared by exact local alignments of the
  whole read; the shared sequence scores the same in both, so the difference is what the read
  holds beyond it. The alignments are skipped where the two transcripts are one sequence under
  the whole read as well, a read end earns the end bonus there only where both transcripts
  reach it (as `--end-bonus-cut` does within an orientation), and the pair stands at the better
  of the two chain scores. In every orientation comparison, a poly(A) tail followed by more
  than 80 read bases now counts too, up to 120, where it starts at the hit transcript's 3′ end:
  cDNA reads that kept their primer and sequencing adapter after the tail have that much. On
  spike-in isoforms antisense to a longer one, many more reads go to the short isoform.
  By default (`--orient-quant`) the longer tail limit and the exact comparison act through the
  quantification score only: the orientation, the hits and their mapping scores are as in
  0.2.1, and each hit's `qscore` moves by what the two would have moved its score. So the
  mapping output, best sets included, is unchanged, and a quantifier that reads `qscore`
  weighs the read by the new evidence; `--no-orient-quant` applies them to the mapping score.
  With them on, a `qscore` is written even with `--end-fit 0`. Mapping costs slightly more.
  `--no-orient-same-exact --orient-tail-after=80` turns the change off.

## 0.2.1 (2026-11-12)

- **Tie-break: the end bonus no longer rewards bases past another candidate's transcript end**
  (`--end-bonus-cut`, on in every preset except the `-fast` ones). A read end earns
  `--end-bonus` points when its alignment reaches the read's end, and the bonus is larger than
  a mismatch costs. A read that ends at its transcript's end followed by one or two unaligned
  bases (adapter, primer, basecalling junk) could reach its end only in a longer isoform,
  through those bases aligned as mismatches into the isoform's continuation, and the longer
  isoform won by a point or more where the read fits both equally. When some candidates'
  windows stop at their transcript's end and others go on, the bonus now counts only for
  reaching the read's end within the shortest such window; past it, a longer candidate is
  credited for the bases it aligns and nothing more. Such reads now tie between the two
  isoforms, as with minimap2, and leave the choice to the quantifier. On simulated ONT cDNA
  reads at several error rates with sparse expression, fewer reads are placed on unexpressed
  isoforms and abundance estimates follow the truth more closely; on spike-in isoforms that are
  a prefix of a longer one, reads move to the shorter isoform. Mapping costs the same.
  `--no-end-bonus-cut` maps as 0.2.0 did.
- **`spearfish map` takes its small allocations from per-thread arenas by default** (the
  `--bump-alloc` of 0.2.0, with its 8 KB limit). Mapping uses less CPU time, with the same
  output, at the cost of a higher peak memory: each thread's arena holds at most 4 MB (at most
  1 GB in all), and where transparent huge pages are always on, the general allocator's peak
  can rise by a few megabytes a thread more (`MIMALLOC_ALLOW_THP=0` leaves only the arenas'
  share). `--no-bump-alloc` turns them off, as `--bump-alloc=off` or `none` still does; of
  the two options, the last one given counts. `--bump-alloc=<size>` sets the limit as before,
  and still warns where the arenas are missing; by default, mapping there goes without them
  and says nothing. The run record says what the arenas did (`stats.bump_alloc`: the limit,
  the most threads that had an arena at once, and how often a thread found none free or its
  own full), and leaves the field out when they were off. Only the `spearfish` binary
  changes: the library never turns the arenas on, the Python package does not use them, and
  a program that embeds the mapper keeps its own choice.
- **Mapping asks the caches for more of the index ahead of use.** The raw-position lookups of
  the tie-break and of the output, the comparison of tied candidates' blocks, a label's
  directory bucket beyond its first line, and the occurrence lists of unitigs shared by many
  isoforms were read one cache miss after another; they are now prefetched in batches. The
  output is the same and mapping takes less time, most on high-accuracy reads. This applies
  to the library and the Python package as well as the binary.
- **Tie-break: the search for a candidate's own A-run at a poly(A) read end carries its count
  from one start position to the next** instead of counting each window again. Same result,
  slightly less work on cDNA and HiFi reads.

## 0.2.0 (2026-11-11)

- **A 0.2, not a 0.1.1: mapping results and the library change.** The default presets map
  differently from 0.1.0 (the entries below), and code written against the 0.1.0 library may
  not compile: the mapping presets gained variants (`Preset::OntFast`, `DrnaFast`,
  `HifiFast`), and `MappedRead`, `IndexOptions` and the option structs of `spearfish-map` and
  `spearfish-index` gained fields. A program that asks for `spearfish = "0.1"` therefore
  keeps 0.1.0's behaviour until it asks for 0.2. So that later additions are not breaking,
  `Preset`, `Hit`, `MappedRead` and `IndexOptions` are now `#[non_exhaustive]`, as
  `MapSettings` already was: outside the crate, match a `Preset` with a wildcard arm, read
  the fields of a `Hit` or `MappedRead` without building one by a literal, and build
  `IndexOptions` with `IndexOptions::builder()` or by changing fields of
  `IndexOptions::default()`. The structs of the lower-level crates (`spearfish-map`,
  `spearfish-index`) stay plain, and may gain fields in any 0.x release; `MapSettings` is the
  stable way to set them.
- **Reads with only a chance match are unmapped.** Every read with a chain used to be
  reported, so reads from outside the transcriptome (unannotated loci, contaminants, random
  sequence) were placed on whichever transcript shared a seed with them. `spearfish map` now
  reports a read unmapped, flagged `low_score`, when its best chain scores under
  `--min-score` + `--min-score-ln-len` · ln(read length / 1000) (36 and 4 by default, for
  every preset). A best hit that only the rescue found, and that scores near the threshold or
  spans little of the read, must also align around its chain (`--min-score-verify`): the
  rescue's seeds are dense enough to chain past the threshold by chance on low-complexity
  reads. On simulated reads this removes nearly every placement of a random read and most
  wrong-gene placements, for under 0.1 points of correctly placed reads. `--keep-low-score`
  keeps the hits (flagged; RAD `read_flags` bit 2), and `--min-score 0` restores the old
  behaviour. The library has `MappedRead::low_score` and the settings `min_score`,
  `min_score_ln_len`, `min_score_anchors`, `min_score_verify` and `keep_low_score`; the run
  record counts `low_score` reads. Callers no longer need a score filter of their own.
- **Doubtful reads get rescue-level work (`ont` and `drna` presets).** A read at 10% error or
  more often has too few base-level seeds on its true transcript: the chain covers a small
  part of the read, a paralog or pseudogene with one more surviving seed wins, or the wrong
  orientation does. Reads whose best chain does not cover them well, or that have a rival on
  the other orientation, are now (`--rescue-settled` sets who is exempt)
  - chained again at the rescue level in their best transcripts (`--rescue-rechain`);
  - seeded against the whole transcriptome, by a table of every transcript's runs of three
    consecutive rescue seeds, which names transcripts the base-level seeds never reached
    (`--rescue-global`, `--rescue-global-cands`, `--rescue-extra-cands`);
  - chained on the other orientation in the transcripts that share seeds, reverse-complemented,
    with the best hit (`--rescue-antisense`), so that an antisense transcript of the same
    locus competes on equal terms.

  On simulated ONT reads, correctly placed reads rise by about 0.4 points at 5% error and 1.5
  points at 10%, and transcript abundances from oarfish follow the truth more closely at
  every error rate from 3% to 14%. The cost: ONT cDNA reads take 1.1 to 1.15 times the CPU
  time they took (direct RNA reads within 5% of it), and mapping holds about 1.5 GB more of
  the index in memory for a human transcriptome: the table, the transcripts' rescue seeds
  (`--rescue-store`) and their antisense mates, all read in place from the index file. The
  `hifi` preset, which has no rescue, costs no more than it did.
- **The index holds the rescue tables.** `spearfish index` writes the table, the
  transcripts' rescue seeds and each transcript's antisense mates (the transcripts that hold
  stretches of its reverse complement) into the index, as sections that 0.1.0 skips. Mapping
  reads them in place, so it starts as fast as it did and keeps no copy of them of its own.
  They are a little over half of the file (2.8 GB for a human transcriptome) and add two to
  three seconds to the build on 32 threads. `spearfish index --no-rescue-tables`
  (`rescue_tables` in the library and in Python) leaves them out; with such an index, or one
  built by 0.1.0, the `ont` and `drna` presets make the table and the seeds at every start,
  in about half a second on 32 threads, and a transcript's mates when a read first needs
  them, and the peak is about 1.4 GB higher. Nothing else in the index changed, and its
  format version is the same.
- **`ont-fast`, `drna-fast` and `hifi-fast` presets** map as 0.1.0 did: the same hits and
  scores, every read with a chain reported; of the mapping changes listed here, only the
  score scale applies to them. On noisy reads they are faster and less accurate than the
  presets they are named for, and they never read the rescue tables.
- Tie-break (all presets unless noted):
  - Read orientations that nearly tie are compared by exact local alignment around each
    chain (`--orient-exact`; off for `hifi`), and an orientation's poly(A) tail counts only
    for the hits it was trimmed for.
  - With raw read ends, the stretches between shared anchors are also aligned in raw bases,
    so that homopolymer-length differences between isoforms and paralogs count
    (`--raw-mid`; off for `hifi`).
  - A read's poly(A) tail is trimmed only past the A's its transcript has there
    (`--tail-templated`): transcripts that end in a templated A-run no longer lose them.
  - Tied hits without a shared anchor are grouped by the anchors they share (`--tie-group`,
    now on).
  - Transcripts within `--report-within-abs` score points of the best (30) are reported even
    when `--report-within` of a low best score would drop them.
- Candidate genes: a seed that overlaps the gene's previous seed on the read now counts for the
  part that extends past it, where it counted for nothing. Reads with many errors, whose few
  seeds on the true gene are often consecutive overlapping ones, keep that gene as a candidate
  more often.
- Rescue: the other orientation is also re-seeded when its base-level chain is more than a
  chance match by itself (`--min-score`), not only when it scores within
  `--rescue-other-margin` of the chosen orientation's. When neither chain covers the read, the
  chosen orientation's lead says little.
- Faster, with the same hits and scores (all presets, the `-fast` ones included):
  - the tie-break's banded affine alignment runs on 16-bit cells when a pass's scores are sure
    to fit them (AVX-512: 32 cells per register instead of 16, with narrow bands held in
    registers from row to row);
  - a read's whole-transcriptome votes for an orientation are kept, and reused when a later
    step asks for no more candidates than the first did;
  - the index's run-start select is branch-free, and output lines format their integers
    without a call per number.
- `SCORE_DENOM`, the scale RAD output suggests for turning score differences into
  probabilities (`score_denom`), is 2, where it was 1: a hit `d` below its read's best counts
  `exp(-d / 2)` times as much. 1 was chosen on simulated sets in which every transcript is
  expressed. When most are not, as in real samples, it hands too many reads to unexpressed
  relatives of the expressed transcripts, and a larger scale gives more accurate abundances.
  2 takes most of that gain for a small loss when every transcript is expressed. 3, oarfish's
  default for minimap2 alignments, gains a little more on such sets, but loses more when
  every transcript is expressed and takes the reads of spike-in isoforms that are nearly
  contained in another isoform.
- `spearfish simulate` can make harder read sets: `--expressed` (fraction of transcripts with
  any abundance), `--fragment` and `--min-fragment` (internal fragments), `--rate-sd` (per-read
  error rates), `--junk` (random sequence at the read ends) and `--unmappable` (reads of random
  sequence, named `r<i>|*|*|0|0`). Without them it writes the reads it wrote before.
- The `low_coverage` flag's documentation said oarfish skips such reads by default; it does
  not, and most of them are correctly placed.
- **An optional allocator for mapping (`spearfish map --bump-alloc`, the new crate
  [`nudgy`](https://crates.io/crates/nudgy)).** Mapping a read makes a few hundred small
  allocations and frees nearly all of them before the next read. With `--bump-alloc`, blocks
  of up to 8 KB come from memory that belongs to the allocating thread (64 KB chunks, each
  used again from its start once its blocks are freed), and everything else from mimalloc as
  before. Mapping takes about 3% less CPU time, for 3 to 8 MB of memory a thread; the output
  is the same. It is off by default, and there on 64-bit Linux and macOS. The allocator is a
  crate of its own, `nudgy` (`crates/nudgy`, published with its own version), and does not
  depend on spearfish: `nudgy::Bump<A>` wraps any global allocator, a program that embeds
  the mapper declares it as its global allocator and calls `nudgy::enable()`, and until that
  call every allocation goes to the wrapped allocator. `nudgy::enable_with` takes another
  limit than 8 KB (up to 32 KB), chosen at run time, and so does the option:
  `--bump-alloc=1KB` (`--bump-alloc=off` or `none` is no arenas). `nudgy::stats()` says how
  many arenas were in use and how often a thread went without one (`spearfish map -v` prints
  it). The crate carries its checks: loom models of what its threads share, a test for Miri,
  stress tests and a benchmark of the allocation paths.
- **`hifi` preset: the poly(A) tail orders the two orientations** (`--tail-first`, on for
  `hifi`). A read that ends in a tail of 20 bases or more in exactly one orientation is mapped
  in that orientation first, and in the other only when the first one's best hit spans less
  than 80% of the read (`--tail-first-span`, `--tail-first-run`). On HiFi libraries that keep
  their tails this is about 20% fewer instructions per mapped read; about 0.2% of reads report
  a different (nearly always smaller) hit set and 0.01% a different top hit, with mapping
  accuracy on simulated reads and spike-in quantification unchanged. Reads without such a tail
  are mapped as before. `--no-tail-first` restores the previous behaviour; the other presets
  are unchanged, and `--tail-first` turns it on for them.
- `spearfish map` says when the reads look oriented: with both orientations mapped and at
  least 99% of the mapped reads on one strand, the summary names the `--strands` value that
  maps them in one orientation, and the run record counts the mapped reads whose best hit is
  forward (`forward`). `--strands fwd` or `rev` takes precedence over `--tail-first`.
- `spearfish --version` and every output's provenance name the commit a crates.io build was
  packaged from (`.cargo_vcs_info.json`), where they said `unknown`.
- `bump_and_publish.sh` resumes an interrupted release: run it again from the same checkout,
  and it continues from the unpushed release commit. While publishing, it waits and retries
  when crates.io is slow to index a new crate or rate-limits new crate names.

## 0.1.0 (2026-10-23)

The first release.

- **Crates.** The library is `spearfish`, and the `spearfish` binary is `spearfish-cli`
  (`cargo install --locked spearfish-cli`). They build on `spearfish-lcp`, `spearfish-seeds`,
  `spearfish-index`, `spearfish-map` and `long-read-rad`, which are published with them; the
  first four were renamed from `lcp-core`, `seeds`, `index` and `map`, which other crates own on
  crates.io (the code still imports them by their old names).
- **Python.** The `spearfish` package (`crates/py-spearfish`).
- **Quantification is oarfish's.** oarfish either maps with spearfish in process
  (`oarfish --mapper spearfish`) or quantifies a RAD file from `spearfish map --rad`
  (`oarfish -a`).
- **Removed before release: `spearfish quant`.** It mapped to a temporary RAD file and ran an
  external `oarfish -a`: a third route to the same estimates, which left the mapping settings
  oarfish's models need to the user. It was last present in commit
  [`2bf22ad`](https://github.com/COMBINE-lab/spearfish/tree/2bf22ad) (`crates/spearfish-cli/src/quant_cmd.rs`),
  from which it can be restored. spearfish's own EM stays as the hidden `spearfish em`, a
  baseline for method development; its crate (`quant`) is now a module of the command-line
  crate.

## Development before 0.1.0

The changes below were made during development, before any release, and are kept for the record.

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

- `spearfish::output`: the TSV and PAF writers (`TSV_HEADER`, `tsv_lines`, `paf_lines`,
  `read_flags`) and, with the new default `rad` feature, RAD output for oarfish (`rad_writer`,
  `rad_meta`, `rad_record`; the format crate is re-exported as `spearfish::rad`), moved from
  the CLI so that other front ends write exactly what `spearfish map` writes.
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
- `Index::build_with`: an in-memory build that records its input in the manifest, as
  `Index::build_to_file_with` does.

### Python

- Python bindings (`crates/py-spearfish`, built with maturin; one abi3 wheel for CPython
  3.9+):
  - `Index.load` / `Index.build`, the latter from a FASTA, a dict, or tuples, with `t2g` or
    `infer_genes`.
  - `Mapper(index, preset, **settings)`: every `MapSettings` field as a keyword.
  - Mapping one read (`map`) or many in parallel with the GIL released (`map_batch`,
    `map_records`, and `map_file`, which reads the file in Rust).
  - Results as objects, NumPy columns, or a pandas data frame. Output matches
    `spearfish map`'s.
  - `MappingBatch.to_tsv`/`to_paf`/`write_tsv`/`write_paf` and `RadWriter` write exactly what
    `spearfish map` writes, through `spearfish::output`.
  - `map_file` and `map_records` map the next chunk while Python works on the current one;
    `map_file` also reads on a thread of its own.
  - Release workflow (`publish-python.yml`): wheels for Linux x86-64/arm64 and macOS, tested
    and then published to PyPI on a `py-spearfish-v*` tag (TestPyPI by hand).

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
- Compressed input decompresses on a thread of its own, ahead of the reader, instead of
  inside paraseq's reader lock, where it held up the mapping threads. zstd input now maps
  about as fast as plain input; gzip and bzip2 are limited by their single-threaded decoding.
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
