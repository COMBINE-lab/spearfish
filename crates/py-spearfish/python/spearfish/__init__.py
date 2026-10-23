"""Python bindings for spearfish, a fast long-read RNA-seq mapper to the transcriptome.

>>> import spearfish
>>> index = spearfish.Index.load("gencode_v50.sfi")          # or Index.build("transcripts.fa")
>>> mapper = spearfish.Mapper(index, "ont")                  # preset: "ont", "drna" or "hifi"
>>> read = mapper.map("ACGT...")                             # one read
>>> read.best.transcript, read.best.score
>>> for batch in mapper.map_file("reads.fq.gz"):             # a file, read and mapped in Rust
...     df = batch.to_pandas()                               # or batch.write_tsv(f), rad.write(batch)

See the README for the full interface.
"""

from ._spearfish import (
    SCORE_DENOM,
    SETTINGS,
    TSV_HEADER,
    FileBatches,
    Hit,
    Index,
    MappedRead,
    Mapper,
    MappingBatch,
    RadWriter,
    RecordBatches,
    __version__,
    presets,
    resolve_settings,
)

__all__ = [
    "SCORE_DENOM",
    "SETTINGS",
    "TSV_HEADER",
    "FileBatches",
    "Hit",
    "Index",
    "MappedRead",
    "Mapper",
    "MappingBatch",
    "RadWriter",
    "RecordBatches",
    "__version__",
    "presets",
    "resolve_settings",
]
