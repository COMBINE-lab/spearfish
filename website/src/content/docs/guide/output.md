---
title: Output format
description: The columns of spearfish map's TSV output.
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
| `low_coverage` | The best chain covers less than half the read. `spearfish quant` skips such reads by default, as oarfish does. |
| `anchor_span` | The span runs from the first to the last anchor instead of to the read's ends (extension off, or the hit scores well below the best). |
| `unmapped` | No transcript was found. |

## Using the scores

Only differences between a read's scores carry meaning: they are the "delta" of
[delta mapping](/spearfish/concepts/delta-mapping/). An EM turns them into conditional
probabilities, for example P(read | transcript) ∝ `exp((s - s_max) / D)`, which is what
`spearfish quant` does with D = 2.
