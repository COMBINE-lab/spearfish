---
title: Mapping reads
description: Map long reads with spearfish map, choose the preset for your read type, and change settings.
---

```sh
spearfish map transcripts.sfi reads.fastq.gz --preset ont -t 16 -o reads.tsv
```

Reads can be FASTA or FASTQ, plain or compressed (gzip, bzip2, zstd); `-` reads standard
input. Direct RNA reads that spell uracil as `U` are handled directly.

## Outputs

Any combination of:

| Option | Output |
| --- | --- |
| `-o FILE` | TSV, one line per reported transcript ([Output](/spearfish/guide/output/)). The default, to standard output, when no output is given. |
| `--paf FILE` | PAF, the same hits ([PAF](/spearfish/guide/output/#paf-output)). |
| `--rad FILE` | RAD records for oarfish ([RAD](/spearfish/guide/output/#rad-output)). |

`-` sends the TSV or the PAF to standard output. Each run also writes a run record,
`<first output file>.meta.json`, or where `--meta` says ([Run records](/spearfish/guide/output/#run-records)).

## Presets

Choose the preset that matches how the reads were made. It sets the defaults of every other
setting.

| Preset | Reads | What it sets |
| --- | --- | --- |
| `ont` (default) | ONT cDNA (PCR-cDNA and direct cDNA) | Both orientations; minimap2 `map-ont` scores; read ends aligned in uncompressed sequence; spans extended to the read ends; low-coverage rescue. |
| `drna` | ONT direct RNA | Sense orientation only; read ends in compressed sequence (no template-switch Gs to correct); shorter poly(A) tails recognized (10 bases); spans end at the outermost anchors. |
| `hifi` | PacBio HiFi (Kinnex, Iso-Seq) | Both orientations, the one ending in a poly(A) tail first; minimap2 `map-hifi` scores; tied hits without a shared anchor grouped by the anchors they share; read-end bands kept on the anchor diagonal; no rescue, since accurate reads rarely need it; spans end at the outermost anchors. |

With every preset, a poly(A) tail at a read's 3′ end is found (a run of A's with at most one
other base) and trimmed before the read end is aligned, and isoforms whose annotated 3′ end lies
where the tail starts gain a few points (`--tail-bonus`). When a read's two orientations align
almost equally well (sequence shared antisense, as with antisense isoforms or a pseudogene inside
another gene), the orientation that ends in the poly(A) tail is preferred (`--orient-tail`).
LCP parsing is not strand-symmetric, so a read's two orientations draw different seeds; both are
compared whenever the other orientation's best chain scores within 35% of the chosen one's
(`--orient-margin-rel`).

The `hifi` preset also uses the tail to save work (`--tail-first`). A read that ends in a
poly(A) tail of 20 bases or more in one orientation, and not in the other, is mapped in that
orientation first, and the other orientation is mapped only if the first one's best hit spans
less than 80% of the read. Most reads of a library that keeps its tails then need one
orientation, which is about a fifth less mapping work, and well under 1% of reads report a
different set of transcripts. Reads without such a tail are mapped in both orientations as
before, so a library whose tails were trimmed loses nothing and gains nothing; if its reads were
also oriented (Iso-Seq FLNC reads, for example), `--strands fwd` maps the one orientation
directly. `--no-tail-first` always maps both.

## Changing settings

Every setting has a flag, and a flag you leave out keeps the preset's value; a preset's choices
are never overwritten by a generic default. On/off settings come in pairs, `--x` and `--no-x`
(`--rescue` turns the rescue on for HiFi reads, `--no-extend-ends` turns span extension off for
ONT reads). To see what a run will use:

```sh
spearfish map --preset hifi --extend-ends --print-params
```

prints every parameter as JSON. `spearfish map -h` lists the common settings; `--help` lists
all of them, with their defaults and a table of the settings that differ between presets.

Settings you may want to change:

| Option | Default | Meaning |
| --- | --- | --- |
| `-t, --threads` | all CPUs | Worker threads. |
| `--report-within` | `0.3` | Report transcripts scoring within this fraction of the read's best. Quantifiers usually keep a tighter set (95% of the best), so the default is generous. |
| `--extend-ends`, `--no-extend-ends` | preset | Extend each reported span past its outermost anchors to where the read's alignment ends. Needed by models that use read positions, such as oarfish's truncation model. |
| `--extend-within` | `0.05` | Extend only hits within this fraction of the best score. |
| `--strands` | preset | Orientations mapped: `both`, `fwd` or `rev`. |
| `--tail-first`, `--no-tail-first` | preset | Map the orientation that ends in a poly(A) tail first, and the other only when the first one's best hit spans less than `--tail-first-span` (0.8) of the read. |
| `--rescue`, `--no-rescue` | preset | Re-seed reads whose best chain covers less than half the read at a lower LCP level. Reads shorter than `--rescue-genes-len` (350 compressed bases) re-seed proportionally more candidate genes. |
| `--tail-trim`, `--no-tail-trim` | on | Trim poly(A) tails before aligning read ends. |
| `--tail-run` | preset | Shortest poly(A) tail recognized. |
| `--clip-ends`, `--no-clip-ends` | on | Clip read ends that don't align (adapters, primers). |

The remaining settings tune gene selection, chaining, the tie-break and the rescue; they are
listed under their headings in `--help`.

## Messages

Progress and a summary go to standard error. `-q` shows errors only, `-v` adds detail. On a
terminal, a progress line shows the reads mapped so far.
