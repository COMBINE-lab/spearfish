---
title: Quantification
description: Estimate transcript abundances from spearfish's mappings with oarfish.
---

spearfish reports scored transcript sets; an EM turns them into abundances.
[oarfish](https://github.com/COMBINE-lab/oarfish) does this, with coverage and positional models.

## spearfish quant

```sh
spearfish quant transcripts.sfi reads.fastq.gz --preset ont -t 16 -o sample
```

maps the reads to a RAD file and runs `oarfish -a` on it. oarfish takes the transcripts, the
score scale, the preset and whether spans reach the read ends from the file, so no other options
are needed. Its outputs are written under the prefix (`sample.quant`, `sample.meta_info.json`,
...), the mapping's run record as `sample.map.meta.json`, and oarfish's log as
`sample.oarfish.log` (`-v` shows it on the terminal instead). The RAD file is removed afterwards
unless `--keep-rad` is given.

Mapping settings work as in `spearfish map`. Arguments after `--` go to oarfish:

```sh
spearfish quant transcripts.sfi reads.fastq.gz -o sample --extend-ends -- --model-truncation
```

oarfish's truncation model uses where reads start and end on transcripts, so map with
`--extend-ends` for it (the `ont` preset already does).

`spearfish quant` needs an oarfish build with spearfish support; give its path with
`--oarfish` or `$OARFISH` if it is not on the `PATH`.

## Step by step

The same, in two steps:

```sh
spearfish map transcripts.sfi reads.fastq.gz --preset ont -t 16 --rad sample.rad
oarfish -a sample.rad -o sample -j 16
```

oarfish can also map reads itself with spearfish (`oarfish --mapper spearfish`), building or
opening a spearfish index in process.
