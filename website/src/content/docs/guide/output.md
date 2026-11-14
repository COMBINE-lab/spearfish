---
title: Output format
description: spearfish map's TSV, PAF and RAD outputs, and the run record next to them.
---

`spearfish map` writes a tab-separated file with a header line and one line per reported
(read, transcript) pair, best first. A read's lines are contiguous; reads appear in the order
they finish, not input order. An unmapped read gets one line with transcript `*` and the flag
`unmapped` (`unmapped,low_score` if a chain was found and rejected).

| Column | Meaning |
| --- | --- |
| `read` | Read name. |
| `read_len` | Read length in bases. |
| `strand` | `+` if the read maps as given, `-` if its reverse complement does. Per line: a read whose orientations nearly tie (antisense transcripts) is reported on both. |
| `rank` | Position of this transcript in the read's list, 0 = best. |
| `transcript` | Transcript name (the FASTA ID up to the first `|`). |
| `score` | Chain score after the base-level tie-break. Higher is better; differences between a read's transcripts are on minimap2's alignment-score scale. |
| `qstart`, `qend` | Span on the read, in raw read coordinates on its forward strand. |
| `tstart`, `tend` | Span on the transcript, in raw transcript coordinates. |
| `anchors` | Seeds in the read's chain on this transcript. |
| `anchored` | Read bases covered by those seeds (homopolymer-compressed). |
| `max_disc` | Largest gap discrepancy between consecutive anchors in the chain. |
| `flags` | Comma-separated flags, or `-`. |

## Flags

| Flag | Meaning |
| --- | --- |
| `rescued` | The hits come from the low-coverage rescue (lower-level LCP seeds). |
| `low_coverage` | The best chain's anchors span less than half the read (`--min-coverage`), measured before spans are extended to the read ends. Most such reads are correctly placed (short matches inside longer reads, or reads with many errors), so the flag alone is no reason to drop a read. |
| `low_score` | The best chain is no better than a chance match (`--min-score`): the read is reported unmapped. With `--keep-low-score` its hits are written with this flag. |
| `anchor_span` | The span runs from the first to the last anchor instead of to the read's ends (extension off, or the hit scores well below the best). |
| `unmapped` | No transcript was found, or the best one was rejected (`low_score`). |

### Mapped or not

A read that matches nothing in the transcriptome (a read from an unannotated locus, a
contaminant, an adapter chimera, low-quality sequence) still shares a seed or two with some
transcript by chance. spearfish reports a read unmapped when its best chain scores under
`--min-score` + `--min-score-ln-len` · ln(read length / 1000): 36 for a 1,000-base read, 27 for
100 bases, 45 for 10,000. The threshold grows with the logarithm of the read's length because
the best chance chain does, as the largest of the chance matches a longer read has more room
for. Every preset uses the same rule, except the `-fast` ones, which report every read with a
chain.

The rescue's seeds are denser, and on low-complexity and foreign reads their chance matches
chain past the threshold more often. A best hit that only the rescue found, and that scores
within 30 points of the threshold or spans less than half the read, must therefore also align
around its chain: its clipped alignment score has to reach 40 plus one point per compressed
read base it spans, or the read is unmapped (`--min-score-verify`).

The defaults were set on a human transcriptome with the default index. Random sequence of 300
to 8,000 bases all but never passes them (no read in 80,000 with the `ont` preset); reversed
(not complemented) real ONT cDNA reads, which keep their repeats and low-complexity runs, pass
about 3 times in 10,000. Nearly every read with a true placement of 150 bases or more passes
at ONT error rates; the reads that fail are mostly shorter than that, and a few are short
reads with many errors, which an aligner can still place. Chain scores depend on the index
scheme, so another scheme needs its own threshold. `--min-score 0` turns the rejection off,
and the run record counts the reads it removed (`low_score`).

## PAF output

With `--paf FILE`, the same hits are written as [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md)
lines, best first per read, as minimap2 writes them without base-level alignment (`-c`).
Unmapped reads are left out.

| Column | Value |
| --- | --- |
| 1–4 | Read name, length, `qstart`, `qend`. |
| 5 | Strand. |
| 6–9 | Transcript name, length, `tstart`, `tend`. |
| 10 | Read bases covered by the hit's seeds (homopolymer-compressed), in place of exact matches. |
| 11 | The longer of the read and transcript spans. |
| 12 | 255: mapping quality is not computed (tied isoforms are reported as hits instead). |

Tags: `tp:A:P` on the best hit and `tp:A:S` on the others, `AS:i` the score rounded, `sc:f`
the score, `cm:i` the anchors in the chain, `rk:i` the rank, and `fl:Z` the TSV flags when there
are any.

## RAD output

With `--rad mappings.rad`, `spearfish map` also writes its mappings in the binary
[RAD format](https://combine-lab.github.io/libradicl/) used by alevin-fry and piscem. The file
holds one record per mapped read with all of its reported transcripts, in chunks that are each
compressed (`--rad-codec`, LZ4 by default). Reads can be decoded in parallel, chunk by chunk.
On 2.2 million PacBio reads (42.2 million hits) the TSV takes 4.6 GB, the RAD file 0.66 GB
with LZ4 and 0.47 GB with zstd, and writing it does not slow mapping.

The record layout (`known_rad_type` `bulk_long_read`) does not depend on spearfish, so other
long-read mappers can write it too:

| Level | Tag | Type | Meaning |
| --- | --- | --- | --- |
| read | `read_len` | u32 | Read length in bases. |
| read | `read_flags` | u8 | Bit 0: `low_coverage`; bit 1: `rescued`; bit 2: `low_score` (only with `--keep-low-score`; such a read's alignments should not be counted). |
| alignment | `ref_ori` | u32 | Transcript id in the low 31 bits; bit 31 set if the read maps as given (`+`). |
| alignment | `ref_start`, `ref_end` | u32 | Span on the transcript (`tstart`, `tend`). |
| alignment | `query_start`, `query_end` | u32 | Span on the read (`qstart`, `qend`). |
| alignment | `score` | f32 | The score column. |
| alignment | `aln_flags` | u8 | Bit 0: `anchor_span`. |

File-level tags give the transcript lengths (`ref_lengths`), the mapper, its version and git
revision, command line, the settings and every parameter used (`params`, JSON), the preset, the
sequencing technology, the score scale (`score_denom`, 2 for spearfish), whether spans were
extended to the read ends (`ends_extended`), the orientations mapped, the index's identifier
(`index_id`), and the numbers of reads, mapped reads and bases. Unmapped reads are counted but get no
record, and records carry no read names.

## Run records

Every `spearfish map` run writes a JSON record next to its first output file
(`reads.tsv.meta.json`; `--meta FILE` to choose; none if `--meta` is not given and no output is a
regular file, as when everything goes to standard output or `/dev/null`). It holds:

- `spearfish`: version, git revision, target and build profile; `simd`: the kernels used;
- `command_line`, `working_directory`, start and finish times, threads and CPUs available;
- `index`: its path, identifier and manifest (what built it, from which transcripts);
- `reads` and `outputs`;
- `settings`: the preset and the settings given; `parameters`: every parameter used;
  `score_denom`;
- `warnings`, when there are any: each with a `code` and the `message` shown. Today the one
  warning is `no_gene_groups`, for an index whose every transcript is its own gene
  ([Gene membership](/spearfish/guide/indexing/#gene-membership));
- `stats`: reads, mapped, unmapped, multi-mapped, rescued, low-coverage and low-score reads,
  mapped reads whose best hit is forward (`forward`), hits, bases, mean read length, hits per mapped read, index load and mapping time, reads and megabases per
  second, CPU time and peak memory, and, when the per-thread arenas were on, what they did
  (`bump_alloc`: the largest block they served, the most threads that had an arena at once,
  and how often a thread found none free or its own full).

## Using the scores

Only differences between a read's scores carry meaning: they are the "delta" of
[delta mapping](/spearfish/concepts/delta-mapping/). An EM turns them into conditional
probabilities, P(read | transcript) ∝ `exp((s - s_max) / D)`. spearfish suggests D = 2
(`spearfish::SCORE_DENOM`), which the RAD file records and oarfish reads there; for minimap2's
alignment scores oarfish uses 3.
