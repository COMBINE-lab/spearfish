---
title: Introduction
description: What spearfish computes, and when to use it.
---

spearfish maps long RNA-seq reads (Oxford Nanopore cDNA and direct RNA, PacBio HiFi) to an
annotated transcriptome. For each read it reports every transcript the read is compatible with,
where the read lies on each one, and a score; these scored sets are the input that
transcript-quantification methods resolve by expectation-maximization (EM).

## Why a transcriptome needs its own mapper

A transcriptome is an unusual reference. Human GENCODE v50 has 670,670 transcripts, but they are
assembled from a shared pool of exons, and a seed from a typical read occurs in about 24 candidate
transcripts. A genome mapper treats each transcript as a separate target, as it would treat
chromosomes, so it finds, chains and aligns the same exon once for every isoform that contains it.
Faster genome mappers make each of those alignments cheaper; spearfish avoids repeating them.

## Delta mapping

spearfish performs **delta mapping (Δ-mapping)**. For each read, it computes:

- every transcript the read is compatible with,
- the read's position on each,
- each transcript's score difference from the best candidate.

It gives up only the part of each score that every candidate shares. Sequence common to all
candidates adds the same amount to every score, so it never needs to be aligned. Compare
pseudoalignment, which gives up positions and scores, and quasi-mapping, which gives up scores.
[Delta mapping](/spearfish/concepts/delta-mapping/) explains the idea in more detail.

## When to use spearfish

- **Transcript quantification from long reads.** This is the application spearfish is built
  for. Its output feeds an EM directly.
- **Fast compatibility sets with positions.** If you need which transcripts a read could come
  from and where, but not base-level alignments, spearfish gives you that quickly.

spearfish does not produce base-level alignments (CIGAR strings), and it maps to a given
transcriptome, so it does not discover novel isoforms. For those, use a spliced aligner such as
minimap2 against the genome.
