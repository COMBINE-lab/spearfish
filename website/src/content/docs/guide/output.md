---
title: Output format
description: spearfish map's TSV, PAF and RAD outputs, and the run record next to them.
---

`spearfish map` writes a tab-separated file with a header line and one line per reported
(read, transcript) pair, best first. A read's lines are contiguous; reads appear in the order
they finish, not input order. An unmapped read gets one line with transcript `*` and the flag
`unmapped`.

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
| `low_coverage` | The best chain covers less than half the read (`--min-coverage`). oarfish skips such reads by default. |
| `anchor_span` | The span runs from the first to the last anchor instead of to the read's ends (extension off, or the hit scores well below the best). |
| `unmapped` | No transcript was found. |

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
| read | `read_flags` | u8 | Bit 0: `low_coverage`; bit 1: `rescued`. |
| alignment | `ref_ori` | u32 | Transcript id in the low 31 bits; bit 31 set if the read maps as given (`+`). |
| alignment | `ref_start`, `ref_end` | u32 | Span on the transcript (`tstart`, `tend`). |
| alignment | `query_start`, `query_end` | u32 | Span on the read (`qstart`, `qend`). |
| alignment | `score` | f32 | The score column. |
| alignment | `aln_flags` | u8 | Bit 0: `anchor_span`. |

File-level tags give the transcript lengths (`ref_lengths`), the mapper, its version and git
revision, command line, the settings and every parameter used (`params`, JSON), the preset, the
sequencing technology, the score scale (`score_denom`, 1 for spearfish), whether spans were
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
- `stats`: reads, mapped, unmapped, multi-mapped, rescued and low-coverage reads, hits, bases,
  mean read length, hits per mapped read, index load and mapping time, reads and megabases per
  second, CPU time and peak memory.

## Using the scores

Only differences between a read's scores carry meaning: they are the "delta" of
[delta mapping](/spearfish/concepts/delta-mapping/). An EM turns them into conditional
probabilities, P(read | transcript) ∝ `exp((s - s_max) / D)`. oarfish uses D = 1 for
spearfish's scores, which the RAD file records; for minimap2's alignment scores it uses 3.
