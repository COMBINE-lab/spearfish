---
title: Quantification
description: Estimate transcript abundances from spearfish's mappings with oarfish.
---

spearfish reports scored transcript sets; an EM turns them into abundances.
[oarfish](https://github.com/COMBINE-lab/oarfish) does this, with coverage and positional models.
It needs a build with spearfish support, which oarfish 0.11.0 and newer have (`conda install -c bioconda
-c conda-forge oarfish`): its `--help` lists `spearfish` as a `--mapper`.

## In one step: oarfish maps with spearfish

```sh
oarfish --reads reads.fastq.gz --annotated transcripts.fa.gz --seq-tech ont-cdna \
    --mapper spearfish -j 16 -o sample
```

oarfish builds a spearfish index from the transcripts in process (or opens one written by
`spearfish index`, given with `--index transcripts.sfi`), maps the reads with the preset that
`--seq-tech` names (`ont-cdna`, `ont-drna` or `pac-bio-hifi`), and quantifies them, writing
`sample.quant` and `sample.meta_info.json`. It sets the mapping options its models need itself:
with `--model-truncation`, for instance, read spans are extended to the read ends. This is the
recommended route.

## In two steps: map to RAD, then quantify

```sh
spearfish map transcripts.sfi reads.fastq.gz --preset ont -t 16 --rad sample.rad
oarfish -a sample.rad -o sample -j 16
```

The RAD file holds every read's scored hits, and oarfish takes the transcripts, the score
scale, the preset and whether spans reach the read ends from it, so no other options are
needed. Map this way to keep the mappings: to quantify them several ways without mapping
again, or to map and quantify on different machines. Mapping settings are yours to choose
here, so map with `--extend-ends` for oarfish's truncation model (`--model-truncation`); the
`ont` preset already extends them.
