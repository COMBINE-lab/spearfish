---
title: Mapping reads
description: Map long reads with spearfish map, and choose the preset for your read type.
---

```sh
spearfish map transcripts.sfi reads.fastq.gz --preset ont -t 16 -o reads.map.tsv
```

Reads can be FASTA or FASTQ, plain or compressed; `-` reads standard input. Direct RNA reads
that spell uracil as `U` are handled directly. Output goes to standard output unless `-o` is
given (see [Output format](/spearfish/guide/output/)).

## Presets

Choose the preset that matches how the reads were made. It sets the base-level scoring and the
orientations mapped.

| Preset | Reads | What it sets |
| --- | --- | --- |
| `ont` (default) | ONT cDNA (PCR-cDNA and direct cDNA) | Both orientations; minimap2 `map-ont` scores; read ends aligned in uncompressed sequence; spans extended to the read ends. |
| `drna` | ONT direct RNA | Sense orientation only; read ends in compressed sequence (no template-switch Gs to correct); spans end at the outermost anchors. |
| `hifi` | PacBio HiFi (Kinnex, Iso-Seq) | minimap2 `map-hifi` scores; tie groups compared by clipped scores; no low-coverage rescue, since accurate reads rarely need it. |

## Options you may want to change

| Option | Default | Meaning |
| --- | --- | --- |
| `-t, --threads` | all cores | Worker threads. |
| `--rad` | | Also write the mappings as a RAD file (see [Output format](/spearfish/guide/output/#rad-output)). Without `-o`, only the RAD file is written. |
| `--rad-codec` | `lz4` | Per-chunk compression of the RAD file: `lz4`, `zstd` or `none`. |
| `--delta-rep` | `0.3` | Report transcripts scoring within this fraction of the read's best. Quantifiers usually keep a tighter set (95% of the best), so the default is generous. |
| `--extend-ends` | preset | Extend each reported span past its outermost anchors to where the read's alignment ends. Needed by models that use read positions, such as a truncation model. |
| `--extend-within` | `0.05` | Extend only hits within this fraction of the best score. |
| `--strands` | preset | Override the orientations mapped: `both`, `fwd` or `rev`. |
| `--no-rescue` | | Turn off re-seeding of reads whose best chain covers less than half the read. |

`spearfish map --help` lists the remaining options, which tune chaining and the tie-break and are
mainly for experiments.
