---
title: Building an index
description: Index a transcriptome with spearfish index.
---

```sh
spearfish index transcripts.fa.gz -o transcripts.sfi -t 16
```

The input is a transcript FASTA, plain or compressed.

## Gene membership

spearfish collapses seeds within a gene, which is where isoforms share their exons, and it picks
candidate genes for each read before scoring their transcripts. So it needs to know which
transcripts belong together. There are three ways to tell it:

- **GENCODE-style headers** (`transcript|gene|...`): used automatically.
- **A transcript-to-gene table**, `--t2g t2g.tsv`, with one `transcript gene` pair per line
  (whitespace-separated; further columns are ignored). Use this for Ensembl, RefSeq or any other
  annotation.
- **Inferred from sequence**, `--infer-genes`, for transcriptomes without an annotation (other
  species, de novo assemblies). Transcripts that share a run of 12 consecutive seeds (about 300
  bases) are grouped, ignoring runs found in more than 64 transcripts, which are repeats; then
  short transcripts left on their own join the one group they share shorter runs (about 100
  bases) with. On human GENCODE these groups give the same mapping accuracy as the annotated
  genes. They merge gene families whose copies are near-identical, which reads cannot tell apart
  either.

Without any of these, every transcript is its own gene and `spearfish index` prints a warning.
Isoforms are then not collapsed: the GENCODE index grows from 1.3 to 3.1 GB, and fewer reads are
assigned to the right transcript (top-1 95.9% falls to 92.5% on simulated ONT reads), because the
per-read candidate set holds a few transcripts instead of a few genes. `--write-t2g out.tsv`
writes the grouping that was used.

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
| `--t2g` | | Transcript-to-gene table (see above). |
| `--infer-genes [N]` | off (N = 64) | Infer genes from shared seed runs, ignoring runs in more than N transcripts. |
| `--infer-span` | `12` | Consecutive seeds a shared run must span. |
| `--infer-attach-span` | `4` | Shorter runs by which small groups attach to a larger one (0 = off). |
| `--infer-attach-max` | `3` | Largest group, in transcripts, that may attach. |
| `--write-t2g` | | Write the transcript-to-gene grouping used. |

`spearfish index-stats` reports on a built index.
