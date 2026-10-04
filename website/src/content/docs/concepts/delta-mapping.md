---
title: Delta mapping
description: What delta mapping keeps, what it gives up, and why that is what quantification needs.
---

Fast transcriptome mapping methods are defined by what they give up relative to full alignment.

| | Candidate transcripts | Position on each | Score | Alignment path |
| --- | --- | --- | --- | --- |
| Alignment | yes | yes | absolute | yes |
| Quasi-mapping | yes | yes | no | no |
| Pseudoalignment | yes | no | no | no |
| **Delta mapping** | yes | yes | **relative to the best candidate** | no |

**Delta mapping (Δ-mapping)** keeps every transcript a read is compatible with, the read's
position on each, and each transcript's score difference from the best candidate. It gives up
only the part of each score that every candidate shares.

## Why the shared part can go

Suppose a read matches two isoforms that differ by one skipped exon. Along the exons both
isoforms contain, any sequencing error costs both of them the same. Aligning those stretches
changes both scores by the same amount and leaves their difference alone. Only the region around
the skipped exon can change which isoform the read prefers.

So a delta mapper never aligns sequence common to all candidates. It aligns only where they
differ, and it computes everything else once per distinct piece of sequence instead of once per
transcript.

## Why the difference is enough

Quantifiers use scores only through differences. An EM over transcripts weights a read's
candidates by a likelihood such as `exp((s - s_max) / D)`, where `s_max` is the read's best
score: adding a constant to every score of a read changes nothing. The deltas are exactly what
the EM consumes.

## What you give up

- **Alignment paths.** There are no CIGAR strings, so variant calling or visualizing alignments
  needs a full aligner.
- **Absolute scores.** You cannot compare a read's score with another read's, only a read's
  candidates with each other.

Positions are kept, which is what lets positional models, such as models of 5′ truncation, work
from delta-mapped reads.

## How spearfish computes it

[Inside spearfish](/spearfish/concepts/inside/) describes the three pieces: seeds that are the
same in every isoform that shares the sequence, chaining shared across isoforms, and base-level
alignment restricted to where tied isoforms differ.
