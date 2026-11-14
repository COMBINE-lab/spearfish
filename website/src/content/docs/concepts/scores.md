---
title: Mapping and quantification scores
description: The two scores spearfish reports for each hit, what each decides, and how the end fit that separates them is computed.
---

Every hit spearfish reports carries two scores.

## The mapping score (`score`)

The mapping score is spearfish's alignment score after the tie-break: the chain score, with
near-tied candidates re-ranked by base-level alignment of the read where they differ. It
includes the read-end evidence that a poly(A) tail gives:

- **3' end**: a candidate whose annotated 3' end lies where the read's poly(A) tail starts gains
  `--tail-bonus` points;
- **strand**: when both read orientations align about equally well, the orientation whose read
  ends in a poly(A) tail gains `--orient-tail` points.

The mapping score alone decides everything about the mapping itself:

- which hit is the primary (rank 0) and the order of the hits;
- the **best set**: the hits at the best mapping score;
- the read's orientation;
- which hits are reported at all.

## The quantification score (`qscore`)

The quantification score is the mapping score plus an **end-fit adjustment** and, between a
read's two orientations, the **orientation evidence** described below. It is the score a
quantifier should weigh hits by: oarfish turns the difference between a hit's score and the
read's best into a weight, `exp(-delta / D)`, and it uses the quantification score for that.

For most hits the two scores are equal. The end fit changes them only when all of the
following hold:

1. **No poly(A) tail.** The read has no detectable poly(A) tail in either orientation (for
   example, because the reads were poly(A)-trimmed, or the tail was lost). When there is a tail,
   its evidence is already part of the mapping score, and the quantification score equals it.
2. **Competing candidates.** At least two hits lie within the tie-break margin of the best
   mapping score. These can be isoforms on the same strand, or hits of the two read
   orientations (antisense transcripts).
3. **Distinguishable ends.** Some candidates' annotated 3' ends lie where the read ends, and
   others' lie farther away than the measurement error (`--end-fit-tol`).

The candidates whose annotated 3' end the read reaches then gain `--end-fit` points in their
quantification score. The others keep their mapping score.

## How the end fit is measured

For each candidate, spearfish uses the alignment of the read's right end: the same affine
X-drop alignment that extends a hit's span from its last seed to the read's end. (Candidates
whose spans were not extended get this alignment computed for the purpose, except where its
result could not change any quantification score.) From it, two numbers:

- `q`: the position in the read of the last aligned base;
- `d`: the number of transcript bases from that base to the candidate's annotated 3' end.

A candidate **fits** its 3' end if

- `d` is at most `--end-fit-tol`, and
- its alignment reaches as far into the read as the farthest-reaching candidate in the same
  orientation, within `--end-fit-tol`.

The second condition matters for a short isoform whose 3' end lies inside the read: its
alignment ends at its own end, but the read goes on, so it does not fit.

The read end is measured from the alignment, not from the number of read bases left over. A
read with a few net insertions near its end still fits the transcript it came from.

**Measuring from the chain (`--end-fit-chain`, the `drna` default).** The `drna` preset does not
extend spans to the read ends, so the alignments above would be run for the end fit alone. It
measures where the read ends from each candidate's chain instead: the read is taken to continue
along the transcript from the chain's last anchor, base for base, until the read or the
transcript ends. Then

- `q` is the read position where that stops, and
- `d` is what is left of the transcript there (0 if the transcript ends first).

The fit rules are the same. No alignment is run, which makes the end fit's cost small. The
price is that bases the read carries past its transcript (adapter or other unaligned sequence)
and indels between the last anchor and the read's end are not seen, so the measurement is
coarser. Use `--no-end-fit-chain` to align instead, as the other presets do; `--end-fit-chain`
turns the chain measurement on in them.

Candidates whose spans were extended (`--extend-ends`) are always measured from their extension's
alignment, which is already there. So turning span extension on in `drna` also turns the chain
measurement off, and the end fit then aligns every candidate, as in the other presets. Give
`--end-fit-chain` as well to keep the chain for the candidates that were not extended.

**Ends that cannot be told apart get no relative adjustment.** Each candidate's annotated 3'
end is placed on the read (`q + d`, relative to that orientation's read end). A candidate whose
end lies within `--end-fit-tol` of a fitting candidate's end gains as the fitting one does
(and so on, transitively). If every candidate gains, or none fits, nothing changes. So the
adjustment only separates candidates whose ends differ by more than the measurement error,
such as an isoform that ends at an internal poly(A) site of a longer sibling.

## Strand without a tail

When both read orientations align about equally well, spearfish already reports the hits of
both orientations, with the mapping score deciding which orientation is primary. Without a
tail, the end fit is computed over the hits of both orientations together. The orientation
whose transcript's 3' end the read reaches gets the higher quantification score. This
separates antisense transcripts that overlap on opposite strands: their reads fit one
transcript's 3' end, but end inside the other.

The mapping score's choice of orientation does not change, and no hits are added: only the
quantification scores differ.

## Orientation evidence beyond the shared sequence

Two rules weigh evidence about a read's orientation that the mapping score did not use before
0.3.0 (see [Mapping](/spearfish/guide/mapping/)):

- a poly(A) tail followed by up to 120 read bases counts where it starts at the hit transcript's
  3' end (`--orient-tail-after`), as on cDNA reads that kept their primer and sequencing adapter;
- best hits that are one sequence on opposite strands are compared by exact local alignments of
  the whole read (`--orient-same-exact`), so bases one transcript has of its own count.

By default (`--orient-quant`) both act on the quantification score only. The orientation, the
hits, their ranks and their mapping scores are decided as without them. Each hit's
quantification score then moves by what the two rules would have moved its mapping score. So
the read weighs toward the orientation they favour in quantification, while the mapping
output, the best set included, stays as without them. `--no-orient-quant` applies them to the
mapping score instead, so they decide the orientation too.

## Choice of the tolerance

`--end-fit-tol` is the measurement error of where a read ends on a transcript, in parsed
(homopolymer-compressed) bases. Its default, 4, comes from the error profile of nanopore reads,
not from fitting a data set:

- homopolymer-length errors, the most common nanopore error, do not shift positions in
  compressed bases;
- the remaining substitutions and short indels occur at a few percent per compressed base, and
  the alignment's end bonus carries it through an isolated error at the very end, so a read that
  reaches a transcript's end usually lands on it or within 1-2 bases;
- poly(A) trimming can remove the transcript's own terminal A's, one compressed base.

Four bases covers these with one base to spare. It is also small enough that isoforms whose 3'
ends differ by more than a few bases are told apart.

## Settings and presets

| setting (library, Python keyword) | CLI flag | default (`ont`, `drna`, `hifi`) | `-fast` presets |
|---|---|---|---|
| `end_fit` | `--end-fit SCORE` | 5 points (0: off) | 0 (off) |
| `end_fit_tol` | `--end-fit-tol BASES` | 4 | (unused) |
| `end_fit_chain` | `--end-fit-chain` / `--no-end-fit-chain` | off, on (off with `--extend-ends`), off | (unused) |
| `orient_quant` | `--orient-quant` / `--no-orient-quant` | on | off |

In Python the settings are keyword arguments of `spearfish.Mapper` (and `Index.mapper`,
`resolve_settings`), as for every other setting:

```python
spearfish.Mapper(index, "drna", end_fit_chain=False)   # measure by alignment
spearfish.Mapper(index, "ont", end_fit=0, orient_quant=False)  # quant_score == score
```

`Hit.quant_score` is the quantification score, also a column of `hit_columns()` and
`to_pandas()`. `Mapper.tsv_header` is the header of that mapper's TSV, with `qscore` when the end
fit is on.

The `-fast` presets reproduce spearfish 0.1.0's mapping, so the end fit is off there and
their output is unchanged. With `--end-fit 0 --no-orient-quant` (or `--tiebreak false`) the
quantification score equals the mapping score, and the output formats are as before (no
`qscore` column or tag).

Points are on the mapping score's scale: with the default score scale, 2 points weigh a factor
`e` in oarfish, so the default 5 points tilt a read's weight by about `e^2.5`, roughly 12 to 1,
toward the candidates it fits. That is a soft preference: the quantifier still weighs the
transcripts' other reads.

## Output formats

All formats report both scores when the end fit is on. The additions are backward compatible:
readers that do not know them are unaffected.

- **TSV** (`-o`): a 15th column, `qscore`, after `flags`. The header line names it. On the line
  of an unmapped read it is `0`. Columns 1-14 are unchanged. `rank` and the best set follow
  `score`, not `qscore`.
- **PAF** (`--paf`): an optional tag `qs:f:` with the quantification score, after the other
  tags. `AS:i` and `sc:f` stay the mapping score.
- **RAD** (`--rad`): each alignment carries an optional `qscore` (f32) tag next to `score`, and
  the file tag `quant_scores` is 1. A reader that does not know the tag skips it. A reader that
  knows it takes `qscore`, or `score` when the file has none.
- **Library** (`spearfish::Hit`): `quant_score` next to `score`; `Mapper::quant_scores()` says
  whether they can differ.
