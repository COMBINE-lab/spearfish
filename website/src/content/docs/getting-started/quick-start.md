---
title: Quick start
description: Index a transcriptome, map reads and quantify, on the toy example and on real data.
---

## On the toy example

The repository includes a toy transcriptome with an exon skip, an alternative 3′ splice site
and a pair of paralogs. Simulate reads from it (with `simulate`, a developer command not
listed in `--help`), then index and map:

```sh
spearfish simulate fixtures/toy/transcripts.fa -n 1000 -o reads.fa --truth truth.tsv
spearfish index fixtures/toy/transcripts.fa -o toy.sfi
spearfish map toy.sfi reads.fa -o map.tsv
```

The simulated read names carry their origin, and `truth.tsv` the true counts.

`map.tsv` has one line per reported transcript for each read:

```
read                    read_len strand rank transcript  score  qstart qend tstart tend anchors anchored max_disc flags
r0|TOYT0000003|-|91|1100     997 -      0    TOYT0000003 344.00 0      997  91     1100 11      371      4        -
r0|TOYT0000003|-|91|1100     997 -      1    TOYT0000004 253.50 11     868  221    1053 10      321      31       anchor_span
```

Next to it, `map.tsv.meta.json` records the run: the spearfish build, the command line, the
index, every parameter used and the read counts ([Output](/spearfish/guide/output/#run-records)).

## On real data

Index a transcriptome once, then map each sample with the preset for its read type:

```sh
spearfish index gencode.v50.transcripts.fa.gz -o gencode_v50.sfi -t 16
spearfish map gencode_v50.sfi sample.fastq.gz --preset ont -t 16 -o sample.tsv
```

Use `--preset drna` for ONT direct RNA and `--preset hifi` for PacBio HiFi
([Mapping reads](/spearfish/guide/mapping/)). Add `--paf sample.paf` for PAF output, or
`--rad sample.rad` for oarfish.

For transcript abundances, quantify the RAD file with
[oarfish](https://github.com/COMBINE-lab/oarfish), or let oarfish map with spearfish itself:

```sh
spearfish map gencode_v50.sfi sample.fastq.gz --preset ont -t 16 --rad sample.rad
oarfish -a sample.rad -j 16 -o sample
```

This writes oarfish's estimates to `sample.quant` ([Quantification](/spearfish/guide/quantification/)).
