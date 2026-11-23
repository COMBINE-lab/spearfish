"""Python bindings for spearfish, a fast long-read RNA-seq mapper: to a transcriptome, or
spliced to a genome.

>>> import spearfish
>>> index = spearfish.Index.load("gencode_v50.sfi")          # or Index.build("transcripts.fa")
>>> mapper = spearfish.Mapper(index, "ont")                  # preset: "ont", "drna" or "hifi"
>>> read = mapper.map("ACGT...")                             # one read
>>> read.best.transcript, read.best.score
>>> for batch in mapper.map_file("reads.fq.gz"):             # a file, read and mapped in Rust
...     df = batch.to_pandas()                               # or batch.write_tsv(f), rad.write(batch)

Spliced mapping to a genome:

>>> genome = spearfish.GenomeIndex.load("genome.sfg", annotation="genes.gtf")   # annotation optional
>>> smapper = spearfish.SplicedMapper(genome, "ont")
>>> smapper.map("ACGT...").best.exons                        # [(start, end), ...]
>>> smapper.map_file("reads.fq.gz", bam="reads.bam")         # straight to BAM, in Rust

See the README for the full interface.
"""

from ._spearfish import (
    SCORE_DENOM,
    SETTINGS,
    TSV_HEADER,
    FileBatches,
    GenomeIndex,
    Hit,
    Index,
    MappedRead,
    Mapper,
    MappingBatch,
    RadWriter,
    RecordBatches,
    SplicedAlignment,
    SplicedMapper,
    SplicedRead,
    __version__,
    presets,
    resolve_settings,
)

__all__ = [
    "SCORE_DENOM",
    "SETTINGS",
    "TSV_HEADER",
    "FileBatches",
    "GenomeIndex",
    "Hit",
    "Index",
    "MappedRead",
    "Mapper",
    "MappingBatch",
    "RadWriter",
    "RecordBatches",
    "SplicedAlignment",
    "SplicedMapper",
    "SplicedRead",
    "__version__",
    "presets",
    "resolve_settings",
]
