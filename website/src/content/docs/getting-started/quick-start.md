---
title: Quick start
description: Index a transcriptome, map reads and quantify, on the toy example and on real data.
---

## On the toy example

The repository includes a toy transcriptome with an exon skip, an alternative 3′ splice site
and a pair of paralogs. Simulate reads from it, then index, map and quantify:

```sh
spearfish index fixtures/toy/transcripts.fa -o toy.sfi
spearfish simulate fixtures/toy/transcripts.fa -n 1000 -o reads.fa --truth truth.tsv
spearfish map toy.sfi reads.fa -o map.tsv
spearfish quant fixtures/toy/transcripts.fa map.tsv -o quant.tsv
```

`map.tsv` has one line per reported transcript for each read:

```
read                    read_len strand rank transcript  score  qstart qend tstart tend anchors anchored max_disc flags
r0|TOYT0000003|-|91|1100     997 -      0    TOYT0000003 344.00 0      997  91     1100 11      371      4        -
r0|TOYT0000003|-|91|1100     997 -      1    TOYT0000004 253.50 11     868  221    1053 10      321      31       anchor_span
```

and `quant.tsv` has the estimated read count of each transcript. The simulated read names carry
their origin, and `truth.tsv` the true counts, so you can compare.

## On real data

Index a transcriptome once, then map each sample with the preset for its read type:

```sh
spearfish index gencode.v50.transcripts.fa.gz -o gencode_v50.sfi -t 16
spearfish map gencode_v50.sfi sample.fastq.gz --preset ont -t 16 -o sample.map.tsv
spearfish quant gencode.v50.transcripts.fa.gz sample.map.tsv -t 16 -o sample.quant.tsv
```

Use `--preset drna` for ONT direct RNA and `--preset hifi` for PacBio HiFi
([Mapping reads](/spearfish/guide/mapping/)).
