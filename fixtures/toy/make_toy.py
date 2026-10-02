#!/usr/bin/env python3
"""Generate the toy transcriptome (deterministic; re-running reproduces the committed files).

Three loci on a synthetic chromosome "toy", plus a paralog pair. See README.md for the
expected splice structure.
"""
import random
from pathlib import Path

rng = random.Random(20261001)
OUT = Path(__file__).parent


def rand_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, rate):
    """Substitutions at the given rate, never to the same base."""
    out = list(seq)
    for i, b in enumerate(out):
        if rng.random() < rate:
            out[i] = rng.choice([c for c in "ACGT" if c != b])
    return "".join(out)


genome = []          # chromosome pieces
pos = 0              # 0-based cursor on the chromosome
gtf = []
fasta = []
tx_counter = [0]


def add_gene(gene_id, gene_name, exons, isoforms):
    """exons: dict name -> sequence; isoforms: dict tx_name -> list of (exon_name, lo, hi)
    slices (lo/hi into the exon sequence, None = full). Exons are laid out in order with
    random introns; each exon's genomic interval covers its full sequence."""
    global pos
    coords = {}
    for name, seq in exons.items():
        intron = rand_seq(rng.randint(300, 1500))
        genome.append(intron)
        pos += len(intron)
        coords[name] = pos  # 0-based start
        genome.append(seq)
        pos += len(seq)
    gene_lo = min(coords.values())
    gene_hi = max(coords[n] + len(s) for n, s in exons.items())
    gtf.append(("toy", "gene", gene_lo + 1, gene_hi, f'gene_id "{gene_id}"; gene_name "{gene_name}";'))
    for tx_name, parts in isoforms.items():
        tx_counter[0] += 1
        tx_id = f"TOYT{tx_counter[0]:07d}"
        seq = ""
        rows = []
        for ex, lo, hi in parts:
            s = exons[ex]
            lo = 0 if lo is None else lo
            hi = len(s) if hi is None else hi
            seq += s[lo:hi]
            g0 = coords[ex] + lo
            rows.append((g0 + 1, g0 + (hi - lo)))
        attrs = f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_name "{gene_name}"; transcript_name "{tx_name}";'
        gtf.append(("toy", "transcript", rows[0][0], rows[-1][1], attrs))
        for i, (a, b) in enumerate(rows, 1):
            gtf.append(("toy", "exon", a, b, attrs + f" exon_number {i};"))
        fasta.append((f"{tx_id}|{gene_id}|-|-|{tx_name}|{gene_name}|{len(seq)}|protein_coding|", seq))


# Gene A: cassette exon A3 (exon skipping).
A = {"A1": rand_seq(400), "A2": rand_seq(250), "A3": rand_seq(120), "A4": rand_seq(500)}
add_gene("TOYG0000001", "GENEA", A, {
    "GENEA-201": [("A1", None, None), ("A2", None, None), ("A3", None, None), ("A4", None, None)],
    "GENEA-202": [("A1", None, None), ("A2", None, None), ("A4", None, None)],
})

# Gene B: alternative 3' splice site on B2 (the short form starts 36 nt later).
B = {"B1": rand_seq(350), "B2": rand_seq(300), "B3": rand_seq(450)}
add_gene("TOYG0000002", "GENEB", B, {
    "GENEB-201": [("B1", None, None), ("B2", None, None), ("B3", None, None)],
    "GENEB-202": [("B1", None, None), ("B2", 36, None), ("B3", None, None)],
})

# Genes C and C': paralogs, C' = C with ~2% substitutions; one isoform each.
C = {"C1": rand_seq(300), "C2": rand_seq(400), "C3": rand_seq(350)}
Cp = {k.replace("C", "D"): mutate(v, 0.02) for k, v in C.items()}
add_gene("TOYG0000003", "GENEC", C, {"GENEC-201": [(k, None, None) for k in C]})
add_gene("TOYG0000004", "GENECP", Cp, {"GENECP-201": [(k, None, None) for k in Cp]})

genome.append(rand_seq(500))


def wrap(s, w=60):
    return "\n".join(s[i:i + w] for i in range(0, len(s), w))


with open(OUT / "transcripts.fa", "w") as f:
    for h, s in fasta:
        f.write(f">{h}\n{wrap(s)}\n")
with open(OUT / "genome.fa", "w") as f:
    f.write(">toy\n" + wrap("".join(genome)) + "\n")
with open(OUT / "annotation.gtf", "w") as f:
    for chrom, feat, a, b, attrs in gtf:
        f.write(f"{chrom}\ttoy\t{feat}\t{a}\t{b}\t.\t+\t.\t{attrs}\n")
