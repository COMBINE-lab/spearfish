---
title: Quantification
description: Estimate transcript abundances from spearfish's output.
---

spearfish reports scored transcript sets; an EM turns them into abundances. There are two ways
to run one.

## spearfish quant

```sh
spearfish quant transcripts.fa.gz reads.map.tsv -t 16 -o reads.quant.tsv
```

The transcript FASTA supplies names and lengths. The output lists each transcript's length and
estimated number of reads. The EM groups reads with the same scored transcript set into
equivalence classes and runs in parallel over them.

| Option | Default | Meaning |
| --- | --- | --- |
| `--score-denom` | `2` | D in `exp((s - s_max) / D)`. Use 5, oarfish's default, for minimap2 alignment scores. |
| `--min-frac` | `0.95` | Keep hits scoring at least this fraction of the read's best. |
| `--keep-low-coverage` | off | Also use reads flagged `low_coverage`. |
| `--length-norm` | off | Divide likelihoods by transcript length. Long-read protocols sequence about one read per molecule, so the default is off. |

## With oarfish

[oarfish](https://github.com/COMBINE-lab/oarfish) is a long-read transcript quantifier with
coverage and positional models. Two ways to use it with spearfish are in development and not
yet part of a released oarfish: a mapping backend that maps reads in process and passes the
scored sets to its EM, and reading a RAD file from `spearfish map --rad` with
`oarfish -a mappings.rad`. The RAD file's tags tell oarfish the transcripts and the score scale,
so no further options are needed. For oarfish's truncation model, map with
`--extend-ends true` so that spans reach the read ends.
