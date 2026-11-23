---
title: Spliced mapping to a genome
description: Build a genome index and map long RNA-seq reads to it with introns, with or without an annotation.
---

spearfish maps long RNA-seq reads to a transcriptome, and it can also map them, spliced, to a
genome. Use the genome when you want to find transcripts that are not annotated, or when the
tool that reads the alignments works in genome coordinates.

```sh
spearfish index --genome genome.fa.gz -o genome.sfg
spearfish map-spliced genome.sfg reads.fastq.gz --preset ont -t 16 --bam reads.bam
```

`index --genome` says that the FASTA is a genome, and `map-spliced` maps reads to the index it
builds. Mapping to a transcriptome stays with `index` and `map`.

## What is computed

A downstream tool reads three things from a spliced alignment: the locus, the strand, and the
coordinates of the exon blocks. spearfish computes those and, by default, nothing else. It
anchors exons with sparse seeds, then settles each junction from the few read bases beside
it. It looks for an exon between two anchored ones when read bases are left over, and extends
the read ends. No base is aligned across an intron, and by default no base is aligned inside
an exon either. That is why the mapping is fast.

If you need a base-level alignment, ask for it with `--base-level`:

| `--base-level` | Each exon block is reported as | Use it for |
| --- | --- | --- |
| `blocks` (default) | one `M` run, with one `I` or `D` if the read and genome lengths differ | quantification, transcript discovery, anything that reads exon structure |
| `edits` | the alignment with the fewest edits | counting mismatches, `NM` |
| `affine` | an alignment with affine gap costs, as minimap2 places indels | variant calling and other tools that read indels |

## Building the index

```sh
spearfish index --genome genome.fa.gz -o genome.sfg -t 16
```

The FASTA can be plain or compressed (gzip, bzip2, zstd). Sequence names are the headers up to
the first white space. The index needs no annotation and stores none, so one index serves
every annotation and every run. A human genome index is a few gigabytes and takes a few
minutes to build. `spearfish inspect genome.sfg` shows what a file holds and what built it.

A genome index and a transcriptome index are different files. `map` and `map-spliced` each
refuse the other's index and name the command to use.

## Using an annotation

Known transcripts are optional and are given when mapping:

```sh
spearfish map-spliced genome.sfg reads.fastq.gz -a annotation.gtf.gz --bam reads.bam
```

`-a` (`--annotation`) takes GTF, GFF3 or BED12, plain or compressed; the format is told from
the content. From the file spearfish takes the exons of each transcript, and uses them in two
ways. A known splice site costs nothing where an unknown one pays for its motif, so a read
that fits a known junction is placed on it. And where read bases are left unexplained, known
exons and junctions are tried first. Junctions that are not in the file are still found, at
the cost their motif gives them.

Sequence names in the annotation must be those of the genome (`chr1` is not `1`). Records on
other sequences are counted and reported; if no transcript is on a sequence of the index, the
run stops with an error.

## Presets

| Preset | Reads | Effect for a genome |
| --- | --- | --- |
| `ont` (default) | ONT cDNA | both orientations; the transcript strand comes from the splice sites |
| `drna` | ONT direct RNA | reads are the transcript's sense strand, so splice sites are sought on one strand |
| `hifi` | PacBio HiFi, Kinnex, Iso-Seq | as `ont` |

## Output

Any combination of the three; PAF goes to standard output when none is named. Reads are
written as they finish, not in input order, and the records of a read are adjacent.

| Option | Output |
| --- | --- |
| `--paf FILE` | PAF with a `cg:Z` CIGAR in which introns are `N` |
| `--sam FILE` | SAM |
| `--bam FILE` | BAM, unsorted (`--bam-level 0` to `9`, default 1) |

The primary record of a read carries its sequence and, for FASTQ input, its base qualities;
secondary records carry neither. SAM and BAM also get a record (flag 4) for each read that
does not map, so the file holds every read; PAF has mapped reads only. `--secondary N` adds up
to N secondary alignments per read.

Qualities nearly double the size of a BAM file. A tool that does not read them, such as a
quantifier, loses nothing with `--no-qualities`, and `--no-unmapped` leaves out the reads
that did not map.

| Tag | Meaning |
| --- | --- |
| `AS:i` | score; comparable between the alignments of one read |
| `cm:i` | seed matches in the chain |
| `ts:A` (PAF), `XS:A` (SAM, BAM) | transcript strand from the splice sites: `+`, `-`, or `.` in PAF when unknown |
| `NM:i` | edits, with `--base-level edits` or `affine` |
| `tp:A` | PAF: `P` primary, `S` secondary |
| `jo:Z` | PAF: one letter per junction for how it was found, for diagnosis |

Mapping quality runs from 0 to 60. It is 60 for a read with one candidate locus and falls as
the second-best locus approaches the best.

Each run writes a record next to its first output file (`reads.bam.meta.json`, or `--meta`):
the build, the command line, the index, the annotation, every parameter, and the counts of
reads, mapped reads and spliced reads.

## Quantification from genome alignments

[oarfish](https://github.com/COMBINE-lab/oarfish) projects genome alignments onto annotated
transcripts and estimates their abundances:

```sh
spearfish map-spliced genome.sfg reads.fastq.gz -a annotation.gtf --secondary 4 \
    --no-qualities --no-unmapped --bam reads.bam
oarfish --genome-alignments reads.bam --annotation annotation.gtf --genome-fasta genome.fa -o sample
```

If all you want is the abundance of annotated transcripts, mapping to the transcriptome
([Quantification](/spearfish/guide/quantification/)) is the shorter route.

## Limits

- Long reads only: ONT cDNA, ONT direct RNA, PacBio.
- Genomes up to 4.29 Gbp in total.
- A read is placed at one locus with introns up to `--max-intron` (200,000 bases by default).
  Chimeric reads, fusion transcripts and circular RNA are not reported as such.
