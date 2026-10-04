---
title: Inside spearfish
description: The core graph, shared chaining and the base-level tie-break, and the intuition behind each.
---

spearfish computes delta mappings in three stages. Each one exploits the same fact: the isoforms
of a gene are built from shared exons, so most of the work for one isoform is also the work for
the others.

## Seeds that are the same wherever their context is the same

spearfish seeds with **locally consistent parsing (LCP)**. LCP cuts a sequence into *cores*
whose boundaries and labels depend only on a small neighborhood. A core in the middle of an exon
therefore comes out identical, at the same offset, in every isoform that contains the exon.

Collapsing identical cores within a gene gives the **core graph**. It stores each exon interior
once and adds a constant number of nodes per splice junction, however many isoforms use the exon.
The cores that differ between isoforms sit in a narrow *halo* around each junction, and those are
exactly the seeds that tell isoforms apart.

Reads are homopolymer-compressed before seeding, because long-read errors concentrate in
homopolymer runs.

## Chaining once per distinct placement

A read's seed hits become isoform evidence through **colinear chaining** in transcript
coordinates: an exon skip shows up as a gap that is consistent in the isoforms that skip the exon
and inconsistent in those that include it.

The chaining dynamic program for a transcript depends on the transcript only through where it
places the read's anchors, so transcripts that place them the same way get the same result.
spearfish chains every candidate transcript exactly, but shares the computation among
transcripts whose placements agree, keeping them in a prefix trie and splitting them only where
they diverge.

## Base-level evidence only where isoforms disagree

When transcripts tie on chain score, spearfish breaks the tie with affine-gap alignment, but only
of the read segments between anchors shared by all tied transcripts, and only where their
sequences differ. Read ends are aligned as extensions, so adapters and other untrimmed sequence
are clipped instead of forced into the alignment, and they are aligned in uncompressed sequence,
where template-switch Gs and isoform ends a few bases apart can be told apart.

## Rescue for reads that seed poorly

Reads whose best chain covers less than half of the read are re-seeded with lower-level LCP cores,
which survive more errors, but only against the transcripts of candidate genes, which keeps them
specific.

## More detail

The manuscript describing spearfish, with proofs that the shared chaining is exact, is in
preparation.
