---
title: Building an index
description: Index a transcriptome with spearfish index.
---

```sh
spearfish index transcripts.fa.gz -o transcripts.sfi -t 16
```

The input is a transcript FASTA, plain or compressed. spearfish reads gene membership from
GENCODE-style headers (`transcript|gene|...`); with other headers each transcript is its own
gene. Gene membership matters: seeds are collapsed within a gene, which is where isoforms share
their exons.

For human GENCODE v50 (670,670 transcripts) the index is 1.3 GB and builds in about 6 seconds on
32 threads. Map against an index built from the same transcript FASTA you quantify against.

## What the index holds

- The **core graph**: the transcriptome's seeds, collapsed so that each exon interior is stored
  once per gene instead of once per isoform, and joined into unitigs.
- For every unitig, the transcripts that contain it and where.
- The homopolymer-compressed sequence of every transcript, and its raw sequence, used to align
  read ends where isoforms differ.

[Inside spearfish](/spearfish/concepts/inside/) explains why this structure is small.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `-t, --threads` | all cores | Worker threads. |
| `-s, --scheme` | `lcp:K=3` | Seed scheme. LCP level-3 cores are the default; minimizers (`minimizer:k=19,w=23`) give similar accuracy. |
| `--hpc` | `true` | Homopolymer-compress transcripts, and reads at mapping time. |
| `--collapse` | `gene` | Collapse seeds per (label, gene) or globally (`global`, larger and slower). |
| `--f-max` | `20` | Mask unitigs shared by more than this many genes (repeats). |

`spearfish index-stats` reports on a built index.
