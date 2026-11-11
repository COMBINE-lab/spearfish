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
| `ont` (default) | ONT cDNA (PCR-cDNA and direct cDNA) | Both orientations; minimap2 `map-ont` scores; read ends aligned in uncompressed sequence; spans extended to the read ends; low-coverage rescue and rescue-level work for doubtful reads (see [Doubtful reads](#doubtful-reads)). |
| `drna` | ONT direct RNA | Sense orientation only; read ends in compressed sequence (no template-switch Gs to correct); shorter poly(A) tails recognized (10 bases); spans end at the outermost anchors. |
| `hifi` | PacBio HiFi (Kinnex, Iso-Seq) | Both orientations, the one ending in a poly(A) tail first; minimap2 `map-hifi` scores; read-end bands kept on the anchor diagonal; no rescue and no raw-base or exact-alignment comparisons of near ties (`--raw-mid`, `--orient-exact`), since accurate reads rarely need them; spans end at the outermost anchors. |
| `ont-fast`, `drna-fast`, `hifi-fast` | As the preset each is named for | That preset as spearfish 0.1.0 mapped it: every read with a chain is reported ([no read is rejected](/spearfish/guide/output/#mapped-or-not)), doubtful reads get no extra work, and near ties are compared as they were then. ONT cDNA reads map in about a fifth to a quarter less CPU time and 1.5 to 2 GB less memory, less accurately when reads have many errors; for direct RNA and HiFi reads the time is about the same. |

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

## Doubtful reads

Most reads are settled by their base-level chain: it covers the read well, and the other
orientation has nothing close. The rest are mostly reads with many errors (10% and more), on
which few base-level seeds survive, so that the chain covers a small part of the read and a
paralog, a pseudogene or the other orientation can win by one seed. With the `ont` and `drna`
presets these reads get more work at the rescue level, whose seeds are shorter and denser:

- their best transcripts are chained again with rescue seeds (`--rescue-rechain`);
- they are seeded against the whole transcriptome, through a table of every transcript's runs
  of three consecutive rescue seeds, which names transcripts the base-level seeds never
  reached (`--rescue-global`);
- the other orientation is chained in the transcripts that share seeds, reverse-complemented,
  with the best hit (`--rescue-antisense`).

`--rescue-settled` sets how well a chain must cover its read for the read to skip all this.
The table, every transcript's rescue seeds (`--rescue-store`) and the transcripts' antisense
mates are read in place from the index, where `spearfish index` keeps them: about 1.5 GB of the
file for a human transcriptome, which is the memory these presets take over the `-fast` ones.
It is the file's own pages, which processes mapping against the same index share. With an
index that lacks them (built with `--no-rescue-tables`, or by spearfish 0.1.0) the table and
the seeds are made when mapping starts, which takes about half a second on 32 threads, and a
transcript's mates when a read first needs them; the peak is then about 1.4 GB higher.
`--no-rescue-global --no-rescue-store` gives most of the memory back, at some cost in
accuracy on high-error reads; `--no-rescue` turns the rescue off altogether, tables included.

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
| `--min-score` | `36` | Report a read unmapped when its best chain scores under this for a 1,000-base read: no better than a chance match (see [Mapped or not](/spearfish/guide/output/#mapped-or-not)). `0` maps every read with a chain. |
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

Where Linux has transparent huge pages always on, the allocator's memory uses them, and each
2 MB stretch that mapping touches is resident whole. If memory is tight, `MIMALLOC_ALLOW_THP=0`
in the environment turns that off, which lowers the peak by several hundred megabytes.

## Messages

Progress and a summary go to standard error. `-q` shows errors only, `-v` adds detail. On a
terminal, a progress line shows the reads mapped so far.
