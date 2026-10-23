"""Tests for the spearfish Python bindings, on the toy transcriptome in fixtures/toy."""

import gzip
import os
import subprocess
import threading
from pathlib import Path

import numpy as np
import pytest

import spearfish

TOY = Path(__file__).resolve().parents[3] / "fixtures" / "toy" / "transcripts.fa"
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    recs, name, seq = [], None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(seq)))
            name, seq = line[1:], []
        else:
            seq.append(line)
    recs.append((name, "".join(seq)))
    return recs


@pytest.fixture(scope="module")
def records():
    return read_fasta(TOY)


@pytest.fixture(scope="module")
def index():
    return spearfish.Index.build(TOY, threads=2)


@pytest.fixture(scope="module")
def mapper(index):
    return spearfish.Mapper(index, "ont", threads=2)


@pytest.fixture(scope="module")
def reads(records):
    """Reads cut from each transcript, forward and reverse complemented, with their origin."""
    out = []
    for i, (_, s) in enumerate(records):
        for start in (0, 150):
            r = s[start:start + 600]
            out.append((f"t{i}_{start}_fw", r, i, False))
            out.append((f"t{i}_{start}_rc", revcomp(r), i, True))
    return out


# ─── Index ──────────────────────────────────────────────────────────────────────────────────


def test_index_metadata(index, records):
    assert index.num_transcripts == len(records) == len(index)
    names = [h.split("|")[0] for h, _ in records]
    assert index.transcript_names() == names
    assert list(index.transcript_lengths()) == [len(s) for _, s in records]
    assert index.transcript_lengths().dtype == np.uint32
    genes = [h.split("|")[1] for h, _ in records]
    assert [index.transcript_gene(i) for i in range(len(records))] == genes
    assert index.num_genes == len(set(genes)) == len(index.gene_names())
    gid = index.transcript_gene_ids()
    assert [index.gene_names()[g] for g in gid] == genes
    assert index.transcript_id(names[2]) == 2
    assert index.transcript_name(1) == names[1]
    assert index.scheme == "lcp:K=3" and index.hpc and index.has_raw_sequence
    assert index.path is None
    rep = index.build_report
    assert rep["transcripts"] == len(records) and rep["genes"] == index.num_genes
    assert "seconds" in rep
    assert index.manifest["kind"] == "spearfish-index"
    assert "Index(transcripts=6" in repr(index)


def test_index_errors(index, tmp_path):
    with pytest.raises(IndexError):
        index.transcript_name(index.num_transcripts)
    with pytest.raises(KeyError):
        index.transcript_id("nope")
    with pytest.raises(OSError):
        spearfish.Index.load(tmp_path / "missing.sfi")
    with pytest.raises(ValueError):
        spearfish.Index.build(TOY, collapse="nope")
    with pytest.raises(ValueError):
        spearfish.Index.build([])
    with pytest.raises(ValueError):
        spearfish.Index.build(TOY, t2g={"a": "b"}, infer_genes=True)


def test_build_to_file_and_load(index, tmp_path):
    out = tmp_path / "toy.sfi"
    built = spearfish.Index.build(str(TOY), out, threads=1)
    assert out.exists() and (tmp_path / "toy.sfi.json").exists()
    assert built.path == str(out.resolve())
    loaded = spearfish.Index.load(out)
    assert loaded.id == built.id == index.id
    assert loaded.build_report is None
    assert loaded.manifest["input"]["genes_from"] == "fasta"
    assert loaded.manifest["input"]["path"] == str(TOY)


def test_build_from_sequences(index, records):
    # The same transcripts as tuples or a dict build the same index (the id hashes the input).
    assert spearfish.Index.build(records).id == index.id
    assert spearfish.Index.build(dict(records)).id == index.id
    with_genes = [(h.split("|")[0], h.split("|")[1], s) for h, s in records]
    assert spearfish.Index.build(with_genes).id == index.id
    assert spearfish.Index.build([(h, s.encode()) for h, s in records]).id == index.id


def test_genes(records):
    # A table puts every transcript in one gene.
    names = [h.split("|")[0] for h, _ in records]
    one = spearfish.Index.build(TOY, t2g={n: "G" for n in names})
    assert one.num_genes == 1 and one.manifest["input"]["genes_from"] == "t2g"
    plain = [(n, s) for n, (_, s) in zip(names, records)]
    alone = spearfish.Index.build(plain)
    assert alone.num_genes == len(records)
    inferred = spearfish.Index.build(plain, infer_genes=True)
    assert inferred.num_genes < len(records)
    assert spearfish.Index.build(plain, infer_genes={"span": 12}).num_genes == inferred.num_genes


def test_t2g_file(tmp_path, records):
    t2g = tmp_path / "t2g.tsv"
    t2g.write_text("".join(f"{h.split('|')[0]}\tG\n" for h, _ in records))
    idx = spearfish.Index.build(TOY, t2g=t2g)
    assert idx.num_genes == 1


# ─── Mapping ────────────────────────────────────────────────────────────────────────────────


def test_map_one(mapper, reads, index):
    for name, seq, tx, rev in reads:
        r = mapper.map(seq, name=name)
        assert r.name == name and r.length == len(seq) and r.is_mapped
        best = r.best
        assert best is r.hits[0] and best.rank == 0
        # The read's own transcript is among the best-scoring hits (isoforms share sequence).
        top = [h.transcript_id for h in r.hits if h.score == best.score]
        assert tx in top, (name, r)
        h = next(h for h in r.hits if h.transcript_id == tx)
        assert h.reverse == rev and h.strand == ("-" if rev else "+")
        assert 0 <= h.query_start < h.query_end <= len(seq)
        assert 0 <= h.target_start < h.target_end <= h.transcript_length
        assert h.transcript == index.transcript_name(tx)
        assert h.gene == index.transcript_gene(tx)
        assert [h.rank for h in r] == list(range(r.num_hits))
        p = r.probabilities()
        assert len(p) == r.num_hits and abs(sum(p) - 1) < 1e-9 and p[0] == max(p)


def test_sequence_types(mapper, reads):
    _, seq, _, _ = reads[0]
    a, b, c = mapper.map(seq), mapper.map(seq.encode()), mapper.map(bytearray(seq.encode()))
    key = lambda r: [(h.transcript_id, h.score, h.query_start, h.target_end) for h in r.hits]
    assert key(a) == key(b) == key(c)
    # Direct RNA spells T as U.
    drna = spearfish.Mapper(mapper.index, "drna")
    assert key(drna.map(seq.replace("T", "U"))) == key(drna.map(seq))
    with pytest.raises(TypeError):
        mapper.map(1234)


def test_unmapped(mapper):
    r = mapper.map("ACGT" * 10)
    assert not r.is_mapped and r.best is None and r.hits == () and r.probabilities() == []
    assert mapper.map("").length == 0


def same(a, b):
    return [(h.transcript_id, h.reverse, h.score, h.query_start, h.query_end, h.target_start,
             h.target_end, h.anchors, h.extended) for h in a.hits] == \
           [(h.transcript_id, h.reverse, h.score, h.query_start, h.query_end, h.target_start,
             h.target_end, h.anchors, h.extended) for h in b.hits]


def test_batch_matches_single(mapper, reads):
    seqs = [s for _, s, _, _ in reads] + ["ACGT" * 10]
    names = [n for n, _, _, _ in reads] + ["junk"]
    batch = mapper.map_batch(seqs, names=names)
    assert len(batch) == len(seqs) and batch.names == names
    assert batch.num_mapped == len(seqs) - 1
    for i, s in enumerate(seqs):
        assert batch[i].name == names[i]
        assert same(batch[i], mapper.map(s))
    assert batch[-1].name == "junk" and not batch[-1].is_mapped
    assert [r.name for r in batch] == names
    with pytest.raises(IndexError):
        batch[len(seqs)]
    with pytest.raises(ValueError):
        mapper.map_batch(seqs, names=names[:-1])
    unnamed = mapper.map_batch(iter(seqs))
    assert unnamed.names is None and unnamed[0].name is None
    assert len(mapper.map_batch([])) == 0


def test_columns(mapper, reads):
    seqs = [s for _, s, _, _ in reads] + ["ACGT" * 10]
    batch = mapper.map_batch(seqs)
    h, r = batch.hit_columns(), batch.read_columns()
    assert len(h["read"]) == batch.num_hits == r["hit_offset"][-1]
    assert list(r["num_hits"]) == [x.num_hits for x in batch]
    assert list(r["length"]) == [len(s) for s in seqs]
    for i, read in enumerate(batch):
        lo, hi = r["hit_offset"][i], r["hit_offset"][i + 1]
        assert (h["read"][lo:hi] == i).all()
        assert list(h["transcript_id"][lo:hi]) == [x.transcript_id for x in read.hits]
        assert list(h["rank"][lo:hi]) == list(range(read.num_hits))
        assert np.allclose(h["score"][lo:hi], [x.score for x in read.hits])
        assert list(h["reverse"][lo:hi]) == [x.reverse for x in read.hits]
    assert h["transcript_id"].dtype == np.uint32 and h["reverse"].dtype == np.bool_


def test_to_pandas(mapper, reads, index):
    pd = pytest.importorskip("pandas")
    seqs = [s for _, s, _, _ in reads] + ["ACGT" * 10]
    names = [n for n, _, _, _ in reads] + ["junk"]
    batch = mapper.map_batch(seqs, names=names)
    df = batch.to_pandas()
    assert len(df) == batch.num_hits and "junk" not in set(df["read"])
    first = df.iloc[0]
    hit = batch[0].best
    assert first["read"] == names[0] and first["transcript"] == hit.transcript
    assert first["strand"] == hit.strand and first["query_end"] == hit.query_end
    assert isinstance(df["transcript"].dtype, pd.CategoricalDtype)
    full = batch.to_pandas(unmapped=True)
    assert len(full) == len(df) + 1
    junk = full[full["read"] == "junk"]
    assert len(junk) == 1 and not junk["mapped"].iloc[0] and pd.isna(junk["transcript"].iloc[0])
    assert list(full["read"].drop_duplicates()) == names
    unnamed = mapper.map_batch(seqs).to_pandas()
    assert unnamed["read"].iloc[0] == 0


def write_reads(path, reads, fastq):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as f:
        for name, seq, _, _ in reads:
            if fastq:
                f.write(f"@{name} some description\n{seq}\n+\n{'I' * len(seq)}\n")
            else:
                f.write(f">{name} some description\n{seq}\n")


@pytest.mark.parametrize("fname,fastq", [("r.fq", True), ("r.fa", False), ("r.fq.gz", True)])
def test_map_file(mapper, reads, tmp_path, fname, fastq):
    path = tmp_path / fname
    write_reads(path, reads, fastq)
    batches = list(mapper.map_file(path, chunk_size=3))
    assert [len(b) for b in batches] == [3] * (len(reads) // 3) + ([len(reads) % 3] if len(reads) % 3 else [])
    got = [r for b in batches for r in b]
    assert [r.name for r in got] == [n for n, _, _, _ in reads]
    for r, (_, seq, _, _) in zip(got, reads):
        assert same(r, mapper.map(seq))
    assert sum(len(b) for b in mapper.map_file(str(path))) == len(reads)


def test_map_file_errors(mapper, tmp_path):
    with pytest.raises(OSError):
        mapper.map_file(tmp_path / "missing.fq")
    with pytest.raises(ValueError):
        mapper.map_file(tmp_path / "missing.fq", chunk_size=0)
    bad = tmp_path / "bad.fq"
    bad.write_text("@r1\nACGT\n+\nII\n")
    with pytest.raises(OSError):
        list(mapper.map_file(bad))
    xz = tmp_path / "r.fq.xz"
    xz.write_bytes(b"\xfd7zXZ\x00\x00\x04 not really xz")
    with pytest.raises(OSError, match="xz-compressed input is not supported"):
        mapper.map_file(xz)


def test_map_records(mapper, reads):
    recs = [(n, s, "qual") for n, s, _, _ in reads]
    batches = list(mapper.map_records(iter(recs), chunk_size=5))
    assert sum(len(b) for b in batches) == len(reads)
    assert [r.name for b in batches for r in b] == [n for n, _, _, _ in reads]
    plain = list(mapper.map_records((s for _, s, _, _ in reads), chunk_size=100))
    assert len(plain) == 1 and plain[0].names is None
    for r, (_, s, _, _) in zip(plain[0], reads):
        assert same(r, mapper.map(s))
    with pytest.raises(TypeError):
        list(mapper.map_records([("only-a-name",)]))


def test_threads(index, reads):
    # One mapper, many Python threads: map() releases the GIL and gives the same results.
    m = spearfish.Mapper(index, "ont", threads=1)
    seqs = [s for _, s, _, _ in reads] * 5
    expect = [m.map(s) for s in seqs]
    out = [None] * len(seqs)

    def work(k):
        for i in range(k, len(seqs), 4):
            out[i] = m.map(seqs[i])

    ts = [threading.Thread(target=work, args=(k,)) for k in range(4)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()
    assert all(same(a, b) for a, b in zip(out, expect))
    assert m.threads == 1
    assert spearfish.Mapper(index, threads=3).threads == 3


# ─── Settings ───────────────────────────────────────────────────────────────────────────────


def test_presets_and_settings(index):
    assert spearfish.presets() == ["ont", "drna", "hifi"]
    assert "extend_ends" in spearfish.SETTINGS and "preset" not in spearfish.SETTINGS
    m = spearfish.Mapper(index, "hifi", extend_ends=True, tie_kappa=5, rescue_window="10,50")
    assert m.preset == "hifi"
    assert m.settings == {"extend_ends": True, "tie_kappa": 5.0,
                          "rescue_window": {"window": {"pad": 10, "div": 50}}}
    o = m.options
    assert o["preset"] == "hifi" and o["params"]["ends"] is not None
    assert o["params"]["tiebreak"]["kappa"] == 5.0
    assert spearfish.resolve_settings("hifi")["params"]["rescue"] is None
    assert spearfish.resolve_settings("ont")["params"]["rescue"] is not None
    assert spearfish.resolve_settings("ont", rescue_window=(5, 20)) == \
        spearfish.resolve_settings("ont", rescue_window="5,20")
    full = spearfish.resolve_settings("ont", rescue_window="full")
    assert full["params"]["rescue"]["window"] is None
    assert spearfish.resolve_settings("ont", chainer="refine")["chainer"] == {"refine": {"beam": None}}
    assert spearfish.resolve_settings("ont", chainer=("refine", 8))["chainer"] == {"refine": {"beam": 8}}
    # None keeps the preset's value.
    assert spearfish.resolve_settings("ont", tie_kappa=None) == spearfish.resolve_settings("ont")
    # Every preset maps the toy reads.
    for p in spearfish.presets():
        assert index.mapper(p).map(index.transcript_name(0)) is not None
    assert index.mapper("drna", tiebreak=False).options["params"]["tiebreak"] is None


def test_setting_errors(index):
    with pytest.raises(TypeError, match="unknown field"):
        spearfish.Mapper(index, "ont", no_such_setting=1)
    with pytest.raises(ValueError):
        spearfish.Mapper(index, "nanopore")
    with pytest.raises(ValueError):
        spearfish.Mapper(index, "ont", tie_kappa="high")
    with pytest.raises(ValueError):
        spearfish.Mapper(index, "ont", pad_div=0)
    with pytest.raises(ValueError):
        spearfish.Mapper(index, "ont", rescue_other_always=True, rescue_other_margin=3.0)
    with pytest.raises(ValueError):
        spearfish.Mapper(index, "ont", strands="sideways")
    assert spearfish.Mapper(index, "ont", strands="fwd").options["params"]["strands"] == "fwd"


def test_settings_change_results(index, reads):
    _, seq, tx, _ = reads[1]  # a reverse-complemented read
    assert all(not h.reverse for h in spearfish.Mapper(index, "ont", strands="fwd").map(seq).hits)
    assert any(h.reverse for h in spearfish.Mapper(index, "ont").map(seq).hits)
    ext = spearfish.Mapper(index, "hifi", extend_ends=True).map(seq)
    assert ext.best.extended


# ─── Pipelines (map_file / map_records map one chunk ahead) ────────────────────────────────


def test_records_read_ahead_is_bounded(mapper, reads):
    pulled = []

    def gen():
        for r in reads * 4:
            pulled.append(r)
            yield r[0], r[1]

    it = mapper.map_records(gen(), chunk_size=4)
    first = next(it)
    assert len(first) == 4
    # The returned chunk and the one mapping behind it, no more.
    assert len(pulled) == 8
    rest = list(it)
    assert sum(len(b) for b in rest) + 4 == len(reads) * 4 == len(pulled)
    names = [r.name for b in [first, *rest] for r in b]
    assert names == [r[0] for r in reads * 4]


def test_records_input_error(mapper, reads):
    def gen():
        for i, r in enumerate(reads):
            if i == 5:
                raise RuntimeError("bad input")
            yield r[0], r[1]

    it = mapper.map_records(gen(), chunk_size=2)
    with pytest.raises(RuntimeError, match="bad input"):
        list(it)
    assert list(it) == []


def test_pipelines_can_be_abandoned(mapper, reads, tmp_path):
    path = tmp_path / "r.fq"
    write_reads(path, reads * 20, True)
    for _ in range(5):
        it = mapper.map_file(path, chunk_size=2)
        next(it)
        del it
        it = mapper.map_records(((n, s) for n, s, _, _ in reads * 20), chunk_size=2)
        next(it)
        del it
    # The mapper still works, and full runs still see every read, in order.
    got = [r.name for b in mapper.map_file(path, chunk_size=1) for r in b]
    assert got == [n for n, _, _, _ in reads * 20]


# ─── Output (spearfish map's formats) ──────────────────────────────────────────────────────


def test_tsv_and_paf(mapper, reads, index, tmp_path):
    seqs = [s for _, s, _, _ in reads] + ["ACGT" * 10]
    names = [n for n, _, _, _ in reads] + ["junk"]
    batch = mapper.map_batch(seqs, names=names)
    tsv = batch.to_tsv().decode()
    lines = tsv.splitlines()
    assert len(lines) == batch.num_hits + 1  # the unmapped read's `*` line
    assert lines[-1] == "junk\t40\t*\t0\t*\t0\t0\t0\t0\t0\t0\t0\t0\tunmapped"
    best = batch[0].best
    f = lines[0].split("\t")
    assert f[:5] == [names[0], str(len(seqs[0])), best.strand, "0", best.transcript]
    assert f[5] == f"{best.score:.2f}"
    assert [int(x) for x in f[6:10]] == [best.query_start, best.query_end, best.target_start,
                                         best.target_end]
    assert batch.to_tsv(header=True).decode() == spearfish.TSV_HEADER + tsv
    assert spearfish.TSV_HEADER.startswith("read\tread_len\t") and spearfish.TSV_HEADER.endswith("\n")
    paf = batch.to_paf().decode().splitlines()
    assert len(paf) == batch.num_hits and not any(l.startswith("junk\t") for l in paf)
    p = paf[0].split("\t")
    assert p[0] == names[0] and p[5] == best.transcript and "tp:A:P" in p
    # Files: binary and text, as written by to_tsv / to_paf.
    with open(tmp_path / "o.tsv", "wb") as fb:
        batch.write_tsv(fb, header=True)
    with open(tmp_path / "o.paf", "w") as ft:
        batch.write_paf(ft)
    assert (tmp_path / "o.tsv").read_text() == spearfish.TSV_HEADER + tsv
    assert (tmp_path / "o.paf").read_bytes() == batch.to_paf()
    # Unnamed reads are named by their position in the batch.
    assert mapper.map_batch(seqs[:2]).to_tsv().decode().startswith("0\t")


@pytest.mark.skipif(not os.environ.get("SPEARFISH_BIN"), reason="SPEARFISH_BIN not set")
def test_output_matches_cli(mapper, reads, tmp_path):
    """Byte for byte what `spearfish map` writes (one thread, so reads stay in input order)."""
    fq = tmp_path / "r.fq"
    write_reads(fq, reads, True)
    idx = tmp_path / "toy.sfi"
    spearfish.Index.build(TOY, idx)
    cli_tsv, cli_paf = tmp_path / "cli.tsv", tmp_path / "cli.paf"
    subprocess.run([os.environ["SPEARFISH_BIN"], "map", str(idx), str(fq), "-t", "1", "-q",
                    "-o", str(cli_tsv), "--paf", str(cli_paf)], check=True)
    m = spearfish.Mapper(spearfish.Index.load(idx), "ont")
    batches = list(m.map_file(fq, chunk_size=7))
    tsv = spearfish.TSV_HEADER.encode() + b"".join(b.to_tsv() for b in batches)
    assert tsv == cli_tsv.read_bytes()
    assert b"".join(b.to_paf() for b in batches) == cli_paf.read_bytes()


def test_rad_writer(mapper, reads, tmp_path, index):
    seqs = [s for _, s, _, _ in reads] + ["ACGT" * 10]
    path = tmp_path / "r.rad"
    with spearfish.RadWriter(path, mapper) as rad:
        for b in mapper.map_records(seqs, chunk_size=5):
            rad.write(b)
        assert rad.num_reads == len(seqs) and rad.num_mapped == len(seqs) - 1
        assert not path.exists()  # written under a temporary name until closed
    assert rad.closed and path.exists() and path.read_bytes()[:8] == b"RAD_FILE"
    assert (rad.num_reads, rad.num_mapped) == (len(seqs), len(seqs) - 1)
    with pytest.raises(ValueError, match="closed"):
        rad.write(mapper.map_batch(seqs[:1]))
    rad.close()  # closing twice does nothing
    # An exception inside the block leaves no file.
    bad = tmp_path / "bad.rad"
    with pytest.raises(RuntimeError):
        with spearfish.RadWriter(bad, mapper, codec="zstd") as rad:
            rad.write(mapper.map_batch(seqs))
            raise RuntimeError("stop")
    assert not bad.exists() and not list(tmp_path.glob("bad.rad*"))
    # A batch from another index can't go into this file.
    other = spearfish.Mapper(spearfish.Index.build(TOY), "ont")
    with spearfish.RadWriter(tmp_path / "x.rad", mapper) as rad:
        with pytest.raises(ValueError, match="another index"):
            rad.write(other.map_batch(seqs[:2]))
    with pytest.raises(ValueError, match="codec"):
        spearfish.RadWriter(tmp_path / "y.rad", mapper, codec="gzip")
