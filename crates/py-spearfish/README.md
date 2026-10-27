# spearfish (Python)

Python bindings for [spearfish](https://github.com/COMBINE-lab/spearfish), a fast long-read
RNA-seq mapper to the transcriptome. Load or build an index, make a mapper from a preset and
settings, and map reads one at a time, in batches, or straight from a FASTA/FASTQ file. The
results are the same as `spearfish map`'s: every reported transcript per read, best first, with
its score and spans.

## Install

```sh
pip install spearfish
```

Wheels are built for Linux (x86-64 and arm64) and macOS (Intel and Apple Silicon). From a checkout
instead (Rust 1.89 or newer and [maturin](https://www.maturin.rs)):

```sh
cd crates/py-spearfish
pip install .                 # or, while developing: maturin develop --release
```

One wheel serves CPython 3.9 and later (the stable ABI). NumPy is required; pandas is
optional, for `to_pandas()`.

## Quick start

```python
import spearfish

index = spearfish.Index.load("transcripts.sfi")      # written by `spearfish index`
mapper = spearfish.Mapper(index, "ont")              # preset: "ont", "drna" or "hifi"

read = mapper.map("ACGT...")                         # one read
if read.is_mapped:
    hit = read.best
    print(hit.transcript, hit.gene, hit.strand, hit.score,
          hit.query_start, hit.query_end, hit.target_start, hit.target_end)

for batch in mapper.map_file("reads.fq.gz"):         # a whole file, read and mapped in Rust
    df = batch.to_pandas()                           # one row per hit, like the CLI's TSV
```

## Indexes

```python
index = spearfish.Index.load("transcripts.sfi")      # memory-mapped

# Build from a FASTA (optionally compressed). GENCODE headers (tx|gene|...) give genes.
index = spearfish.Index.build("gencode.transcripts.fa.gz", threads=16)          # in memory
index = spearfish.Index.build("gencode.transcripts.fa.gz", "gencode.sfi")       # to a file

# ... or from sequences: a dict, or (name, seq) / (name, gene, seq) tuples.
index = spearfish.Index.build([("tx1", "gene1", "ACGT..."), ("tx2", "gene1", "ACGG...")])

# Genes from a table, or inferred from shared sequence when there are none.
index = spearfish.Index.build("transcripts.fa", t2g="t2g.tsv")       # or a dict {tx: gene}
index = spearfish.Index.build("transcripts.fa", infer_genes=True)
```

`Index.build` takes the options of `spearfish index`: `scheme`, `hpc`, `mask_genes`,
`collapse`, `threads`, and `allow_collisions`. An index tells you about itself:
`num_transcripts`, `num_genes`, `transcript_names()`, `transcript_lengths()`,
`transcript_gene(i)`, `gene_names()`, `transcript_id(name)`, `scheme`, `hpc`, `id`, `manifest`
(what built it and from what), and `build_report` (construction statistics, for an index built
in this process).

## Mappers, presets and settings

```python
mapper = spearfish.Mapper(index, "hifi")
mapper = spearfish.Mapper(index, "ont", extend_ends=True, tie_kappa=5.0, threads=8)
mapper = index.mapper("drna")                        # the same thing
```

The preset sets every default. Keyword arguments change single settings, as `spearfish map`'s
flags do: `spearfish.SETTINGS` lists them, and `None` keeps the preset's value. An unknown
setting raises `TypeError`; an invalid value, or settings that contradict each other, raise
`ValueError`. Settings that have no effect, such as tie-break settings with `tiebreak=False`, are
logged as warnings through Python's `logging` module, under the `spearfish` logger.

```python
mapper.settings                                      # what you changed
mapper.options                                       # every parameter the mapper runs with
spearfish.resolve_settings("hifi", extend_ends=True) # the same, without an index
```

Most settings are numbers or booleans. A few take strings:
- `strands`: `"both"`, `"fwd"` or `"rev"`.
- `rescue_window`: `"full"`, `"PAD,DIV"` or `(pad, div)`.
- `chainer`: `"oracle"`, `"refine"` or `("refine", beam)`.

## Mapping

| Call | Input | Returns | Use it for |
| --- | --- | --- | --- |
| `mapper.map(seq, name=None)` | one sequence (`str` or `bytes`) | `MappedRead` | interactive use, a few reads |
| `mapper.map_batch(seqs, names=None)` | an iterable of sequences | `MappingBatch` | reads already in memory |
| `mapper.map_records(records, chunk_size=20000)` | an iterable of sequences or `(name, seq, ...)` tuples, as `mappy.fastx_read` and `pyfastx` give them | iterator of `MappingBatch` | streams from other readers |
| `mapper.map_file(path, chunk_size=20000)` | a FASTA/FASTQ file (plain, gzip, bzip2 or zstd; `"-"` for stdin) | iterator of `MappingBatch` | whole files |

Results come back in input order. A `MappedRead` has the following fields:
- `name` and `length`;
- `hits` (a tuple of `Hit`, best first) and `best`;
- `is_mapped`, `low_coverage` and `rescued`;
- `probabilities()`, the hits' relative likelihoods, `exp(-(best - score) / SCORE_DENOM)` normalized to sum to 1.

A `Hit` has the following fields:
- the transcript: `transcript`, `transcript_id`, `gene` and `transcript_length`;
- the match: `rank`, `strand` (or `reverse`) and `score`;
- the spans: `query_start`, `query_end`, `target_start` and `target_end`, in raw bases, half-open, with query spans on the read's forward strand;
- the chain: `anchors`, `anchored`, `max_discrepancy` and `extended`.

A `MappingBatch` keeps its results in Rust until you ask for them:
- `batch[i]` and iteration give `MappedRead` objects;
- `batch.hit_columns()` and `batch.read_columns()` give whole columns as NumPy arrays, without creating a Python object per hit;
- `batch.to_pandas(unmapped=False)` gives one row per hit, with transcripts as a categorical column.

### Writing spearfish map's formats

A batch writes exactly what `spearfish map` writes, through the same library code:

```python
with open("reads.tsv", "wb") as tsv, open("reads.paf", "wb") as paf, \
        spearfish.RadWriter("reads.rad", mapper) as rad:
    tsv.write(spearfish.TSV_HEADER.encode())
    for batch in mapper.map_file("reads.fq.gz"):
        batch.write_tsv(tsv)          # or batch.to_tsv() -> bytes
        batch.write_paf(paf)          # or batch.to_paf()
        rad.write(batch)              # RAD for `oarfish -a`
```

- **TSV:** one line per hit, plus one `*` line per unmapped read.
- **PAF:** unmapped reads are left out.
- **Read names:** reads without names are named by their position in the batch.
- **Files:** binary and text files both work.
- **RAD:** `RadWriter` records the mapper's preset, settings and index in the file's header. It moves the file into place, with its read counts, when the `with` block ends. On an exception, or if the writer is never closed, no file is left. The codec is `"lz4"` (default), `"zstd"` or `"none"`.

### Why batches

Mapping one long read costs far more than one call from Python into Rust, so calling `map` per
read costs little *on one thread*. The difference is parallelism:
- **Batch calls** (`map_batch`, `map_records`, `map_file`) map on the mapper's own thread pool (`threads=`, default every CPU), with the GIL released, and run at close to `spearfish map`'s speed.
- **The iterators work ahead.** While your code works on one batch from `map_file` or `map_records`, the next chunk is already mapping in the background. `map_file` also reads and decompresses on a thread of its own. Time spent on each batch in Python (`to_pandas`, filtering, writing) is therefore hidden behind mapping, as long as it takes less time than mapping a chunk.
- **`map` from several Python threads** also releases the GIL while it maps, and a mapper can be shared between threads. But building each result's Python objects needs the GIL, so this stops scaling after a few threads.

For more than a few thousand reads, use a batch call. The iterators hold a few chunks at a
time: the one you are using, the one mapping, and for `map_file` one being read. So memory
stays flat as long as you drop each batch when you are done with it; `chunk_size` sets the
scale. Columnar results (`hit_columns`, `to_pandas`) are the cheapest to consume.

## Thread safety

`Index` and `Mapper` are immutable and thread-safe; share one of each across threads. The
iterators returned by `map_file` and `map_records` are meant for one consumer: if a second
thread calls `next()` on a `map_records` iterator while the first is in it, it raises
`RuntimeError` rather than wait.

## Releases

`.github/workflows/publish-python.yml` builds abi3 wheels and an sdist:
- wheels for Linux x86-64 and arm64 (manylinux) and macOS (universal2);
- every wheel is tested on the platforms it serves before anything is published.

To release:
1. Push a tag `py-spearfish-vX.Y.Z`, where `X.Y.Z` is the workspace version.
2. The workflow publishes to PyPI by trusted publishing.

Run by hand, the workflow builds and tests, and publishes to TestPyPI if asked.

## Development

```sh
cd crates/py-spearfish
maturin develop --release
pytest tests
```

The tests use the toy transcriptome in `fixtures/toy`.
