"""Type stubs for the spearfish Python bindings."""

from __future__ import annotations

import os
from typing import IO, Any, Iterable, Iterator, Literal, Mapping, Tuple, Union

import numpy as np
import numpy.typing as npt

__version__: str
SCORE_DENOM: float
"""Scale of score differences: a hit ``d`` below its read's best is ``exp(-d / SCORE_DENOM)``
times as likely."""
SETTINGS: tuple[str, ...]
"""Every mapping setting a :class:`Mapper` takes as a keyword argument."""
TSV_HEADER: str
"""The header line of ``spearfish map``'s TSV (newline included)."""

PathLike = Union[str, os.PathLike]
Seq = Union[str, bytes, bytearray, memoryview]
Preset = Literal["ont", "drna", "hifi"]

def presets() -> list[str]: ...
def resolve_settings(preset: Preset = "ont", **settings: Any) -> dict[str, Any]:
    """The complete parameters a preset and settings resolve to (``spearfish map
    --print-params``)."""

class Index:
    """A spearfish transcriptome index (immutable; share it across mappers and threads)."""

    @staticmethod
    def load(path: PathLike) -> Index:
        """Open an index file (memory-mapped)."""
    @staticmethod
    def build(
        transcripts: PathLike
        | Mapping[str, Seq]
        | Iterable[Tuple[str, Seq] | Tuple[str, str, Seq]],
        output: PathLike | None = None,
        *,
        t2g: PathLike | Mapping[str, str] | None = None,
        infer_genes: bool | int | Mapping[str, int] | None = None,
        scheme: str = "lcp:K=3",
        hpc: bool = True,
        mask_genes: int = 20,
        collapse: Literal["gene", "global"] = "gene",
        threads: int = 0,
        allow_collisions: bool = False,
    ) -> Index:
        """Build an index from a transcript FASTA or sequences, in memory or to ``output``."""

    path: str | None
    num_transcripts: int
    num_genes: int
    scheme: str
    hpc: bool
    has_raw_sequence: bool
    id: str
    manifest: dict[str, Any] | None
    build_report: dict[str, Any] | None

    def transcript_name(self, i: int) -> str: ...
    def transcript_length(self, i: int) -> int: ...
    def transcript_gene(self, i: int) -> str: ...
    def transcript_names(self) -> list[str]: ...
    def transcript_lengths(self) -> npt.NDArray[np.uint32]: ...
    def transcript_gene_ids(self) -> npt.NDArray[np.uint32]: ...
    def gene_names(self) -> list[str]: ...
    def transcript_id(self, name: str) -> int: ...
    def mapper(self, preset: Preset = "ont", *, threads: int = 0, **settings: Any) -> Mapper: ...
    def __len__(self) -> int: ...

class Hit:
    """One reported transcript for a read (raw coordinates, half-open spans; query spans on the
    read's forward strand)."""

    transcript_id: int
    transcript: str
    gene: str
    transcript_length: int
    rank: int
    reverse: bool
    strand: Literal["+", "-"]
    score: float
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    anchors: int
    anchored: int
    max_discrepancy: float
    extended: bool

class MappedRead:
    """All hits of one read, best first."""

    name: str | None
    length: int
    hits: tuple[Hit, ...]
    best: Hit | None
    is_mapped: bool
    num_hits: int
    low_coverage: bool
    rescued: bool

    def probabilities(self) -> list[float]: ...
    def __iter__(self) -> Iterator[Hit]: ...

class MappingBatch:
    """Results of a batch of reads, in input order, held in Rust until asked for."""

    index: Index
    names: list[str] | None
    num_hits: int
    num_mapped: int

    def hit_columns(self) -> dict[str, np.ndarray]: ...
    def read_columns(self) -> dict[str, np.ndarray]: ...
    def to_pandas(self, *, unmapped: bool = False) -> Any: ...
    def to_tsv(self, *, header: bool = False) -> bytes: ...
    def to_paf(self) -> bytes: ...
    def write_tsv(self, file: IO[Any], *, header: bool = False) -> None: ...
    def write_paf(self, file: IO[Any]) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, i: int) -> MappedRead: ...
    def __iter__(self) -> Iterator[MappedRead]: ...

class FileBatches(Iterator[MappingBatch]):
    def __next__(self) -> MappingBatch: ...

class RecordBatches(Iterator[MappingBatch]):
    def __next__(self) -> MappingBatch: ...

class Mapper:
    """Maps long reads to an index with a preset and settings (immutable, thread-safe)."""

    def __init__(
        self, index: Index, preset: Preset = "ont", *, threads: int = 0, **settings: Any
    ) -> None: ...

    index: Index
    preset: str
    threads: int
    settings: dict[str, Any]
    options: dict[str, Any]
    score_denom: float

    def map(self, seq: Seq, name: str | None = None) -> MappedRead: ...
    def map_batch(
        self, seqs: Iterable[Seq], names: Iterable[str] | None = None
    ) -> MappingBatch: ...
    def map_file(self, path: PathLike, *, chunk_size: int = 20000) -> FileBatches: ...
    def map_records(
        self, records: Iterable[Seq | Tuple[Any, ...]], *, chunk_size: int = 20000
    ) -> RecordBatches: ...

class RadWriter:
    """Writes mapped batches to a RAD file (``bulk_long_read``) for ``oarfish -a``."""

    def __init__(
        self, path: PathLike, mapper: Mapper, *, codec: Literal["lz4", "zstd", "none"] = "lz4"
    ) -> None: ...

    closed: bool
    num_reads: int
    num_mapped: int

    def write(self, batch: MappingBatch) -> None: ...
    def close(self) -> None: ...
    def __enter__(self) -> RadWriter: ...
    def __exit__(self, exc_type: Any = None, exc: Any = None, tb: Any = None) -> bool: ...
