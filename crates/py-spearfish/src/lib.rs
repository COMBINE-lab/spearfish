//! Python bindings for spearfish: build or load an index, make a mapper from a preset and
//! settings, and map long reads one at a time, in batches, or straight from a FASTA/FASTQ file.
//!
//! The Python package (`python/spearfish`) re-exports these classes and adds pure-Python
//! helpers (data frames, chunked mapping of any iterable). Everything that runs per read
//! happens here with the GIL released; batches map in parallel on the mapper's own thread pool.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use numpy::IntoPyArray;
use pyo3::exceptions::{PyIOError, PyIndexError, PyKeyError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyString, PyTuple};
use rayon::prelude::*;

use spearfish::provenance::{InputInfo, absolute, file_xxh3};
use spearfish::{BuildReport, Index, IndexOptions, InferGenesOptions, MapSettings, MappedRead, Preset};

/// `anyhow` errors as Python `ValueError`s, with their causes.
fn value_err(e: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn io_err(what: &str, e: impl std::fmt::Display) -> PyErr {
    PyIOError::new_err(format!("{what}: {e:#}"))
}

/// `0` means every available CPU.
fn threads(requested: usize) -> usize {
    if requested == 0 {
        std::thread::available_parallelism().map_or(1, |n| n.get())
    } else {
        requested
    }
}

fn pool(n: usize) -> PyResult<rayon::ThreadPool> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads(n))
        .thread_name(|i| format!("spearfish-{i}"))
        .build()
        .map_err(value_err)
}

/// A path from `str` or `os.PathLike`, or `None` for anything else.
fn as_path(obj: &Bound<'_, PyAny>) -> PyResult<Option<PathBuf>> {
    if obj.is_instance_of::<PyString>() || obj.hasattr("__fspath__")? {
        return Ok(Some(obj.extract::<PathBuf>()?));
    }
    Ok(None)
}

/// A read's bases as `bytes`: `bytes` as is, `str` encoded (ASCII bases), and any other
/// buffer of bytes (`bytearray`, `memoryview`) copied.
fn seq_bytes<'py>(obj: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyBytes>> {
    if let Ok(b) = obj.cast::<PyBytes>() {
        return Ok(b.clone());
    }
    if let Ok(s) = obj.cast::<PyString>() {
        return s.encode_utf8();
    }
    match obj.extract::<Vec<u8>>() {
        Ok(v) => Ok(PyBytes::new(obj.py(), &v)),
        Err(_) => Err(PyTypeError::new_err(format!(
            "a sequence must be str or bytes, not {}",
            obj.get_type().name()?
        ))),
    }
}

/// The read name from a FASTA/FASTQ ID: up to the first whitespace (lossy UTF-8), as
/// `spearfish map` reports it.
fn read_name(id: &[u8]) -> String {
    let id = String::from_utf8_lossy(id);
    id.split_whitespace().next().unwrap_or("").to_string()
}

// ─── Index ──────────────────────────────────────────────────────────────────────────────────

/// A spearfish transcriptome index.
///
/// Open one written by ``spearfish index`` (or :meth:`Index.build` with ``output``) with
/// :meth:`Index.load`, which memory-maps it, or build one from a transcript FASTA or from
/// sequences with :meth:`Index.build`. An index is immutable and can be shared by any number of
/// mappers and threads.
#[pyclass(name = "Index", module = "spearfish", frozen, weakref)]
struct PyIndex {
    inner: Arc<Index>,
    path: Option<String>,
    report: Option<BuildReport>,
}

/// Transcripts as (name, gene, sequence), as `spearfish::Index::build` takes them.
type Transcripts = Vec<(String, String, Vec<u8>)>;

/// Transcripts as (name, gene, sequence) from a FASTA path, a mapping of name to sequence, or
/// an iterable of (name, sequence) or (name, gene, sequence) tuples. Names given without a gene
/// are split as FASTA IDs are (GENCODE `tx|gene|...`).
fn transcripts_from(obj: &Bound<'_, PyAny>) -> PyResult<(Transcripts, Option<PathBuf>)> {
    let py = obj.py();
    if let Some(path) = as_path(obj)? {
        let txs = py
            .detach(|| spearfish::read_transcripts(&path))
            .map_err(|e| io_err(&format!("reading {}", path.display()), e))?;
        return Ok((txs, Some(path)));
    }
    let mut txs = Vec::new();
    let mut push = |name: String, gene: Option<String>, seq: &Bound<'_, PyAny>| -> PyResult<()> {
        let seq = seq_bytes(seq)?.as_bytes().to_vec();
        let (name, gene) = match gene {
            Some(g) => (spearfish::name_and_gene(&name).0, g),
            None => spearfish::name_and_gene(&name),
        };
        txs.push((name, gene, seq));
        Ok(())
    };
    if let Ok(d) = obj.cast::<PyDict>() {
        for (k, v) in d.iter() {
            push(k.extract()?, None, &v)?;
        }
    } else {
        for item in obj.try_iter()? {
            let item = item?;
            let t = item.cast::<PyTuple>().map_err(|_| {
                PyTypeError::new_err("transcripts must be a FASTA path, a dict of name to sequence, or (name, seq) / (name, gene, seq) tuples")
            })?;
            match t.len() {
                2 => push(t.get_item(0)?.extract()?, None, &t.get_item(1)?)?,
                3 => push(
                    t.get_item(0)?.extract()?,
                    Some(t.get_item(1)?.extract()?),
                    &t.get_item(2)?,
                )?,
                n => {
                    return Err(PyValueError::new_err(format!(
                        "a transcript tuple has 2 or 3 fields (name, [gene,] seq), not {n}"
                    )));
                }
            }
        }
    }
    if txs.is_empty() {
        return Err(PyValueError::new_err("no transcripts given"));
    }
    Ok((txs, None))
}

fn infer_genes_from(obj: Option<&Bound<'_, PyAny>>) -> PyResult<Option<InferGenesOptions>> {
    let Some(obj) = obj else { return Ok(None) };
    if obj.is_none() {
        return Ok(None);
    }
    if let Ok(b) = obj.extract::<bool>() {
        return Ok(b.then(InferGenesOptions::default));
    }
    if let Ok(max_tx) = obj.extract::<u32>() {
        return Ok(Some(InferGenesOptions::builder().max_tx(max_tx).build()));
    }
    if obj.is_instance_of::<PyDict>() {
        return pythonize::depythonize::<InferGenesOptions>(obj)
            .map(Some)
            .map_err(|e| PyValueError::new_err(format!("infer_genes: {e}")));
    }
    Err(PyTypeError::new_err(
        "infer_genes must be None, a bool, an int (the max_tx repeat threshold) or a dict of options",
    ))
}

fn report_dict<'py>(py: Python<'py>, r: &BuildReport) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    macro_rules! put {
        ($($f:ident),* $(,)?) => { $( d.set_item(stringify!($f), r.$f)?; )* };
    }
    put!(
        transcripts,
        genes,
        unseeded_transcripts,
        seed_occurrences,
        nodes,
        unitigs,
        cyclic_unitigs,
        unitig_n50_seeds,
        occurrences,
        colors,
        masked_unitigs,
        families,
        largest_family_tx,
        multi_gene_families,
        inconsistent_lengths,
        label_collisions,
        heads_start,
        heads_multi_pred,
        heads_pred_branches,
        largest_group_tx,
    );
    let secs = PyDict::new(py);
    for (k, s) in &r.seconds {
        secs.set_item(*k, *s)?;
    }
    d.set_item("seconds", secs)?;
    Ok(d)
}

impl PyIndex {
    fn check_tx(&self, t: usize) -> PyResult<u32> {
        let n = self.inner.n_transcripts();
        if t >= n {
            return Err(PyIndexError::new_err(format!(
                "transcript {t} out of range ({n} transcripts)"
            )));
        }
        Ok(t as u32)
    }

    fn n_genes(&self) -> usize {
        self.inner.raw().gene_name_off().len().saturating_sub(1)
    }
}

#[pymethods]
impl PyIndex {
    /// Open an index file written by ``spearfish index`` or :meth:`Index.build` (memory-mapped).
    ///
    /// :param path: The index file (``.sfi``).
    /// :raises OSError: If the file can't be read or isn't a compatible index.
    #[staticmethod]
    fn load(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        let inner = py
            .detach(|| Index::open(&path))
            .map_err(|e| io_err(&format!("opening index {}", path.display()), e))?;
        Ok(PyIndex {
            inner: Arc::new(inner),
            path: Some(absolute(&path)),
            report: None,
        })
    }

    /// Build an index, in memory or written to ``output``.
    ///
    /// :param transcripts: A transcript FASTA (path, optionally compressed; GENCODE headers
    ///     ``tx|gene|...`` give gene IDs), a dict of name to sequence, or an iterable of
    ///     ``(name, seq)`` or ``(name, gene, seq)`` tuples.
    /// :param output: Write the index here (with its manifest as ``<output>.json``) and open it
    ///     memory-mapped; without it the index lives in memory only.
    /// :param t2g: Transcript-to-gene table, a path (whitespace-separated ``transcript gene``
    ///     per line) or a dict; overrides genes from the headers.
    /// :param infer_genes: Group transcripts into genes from shared sequence, for
    ///     transcriptomes without gene IDs: ``True`` for the defaults, an int for the repeat
    ///     threshold (``max_tx``), or a dict of ``max_tx``, ``span``, ``attach_span``,
    ///     ``attach_max``.
    /// :param scheme: Base seed scheme, e.g. ``"lcp:K=3"`` or ``"minimizer:k=19,w=23"``.
    /// :param hpc: Homopolymer-compress transcripts (and reads at mapping time).
    /// :param mask_genes: Mask unitigs whose color touches more than this many genes.
    /// :param collapse: Graph nodes per ``"gene"`` (label, gene) or ``"global"`` (label).
    /// :param threads: Build threads (0: every available CPU).
    /// :param allow_collisions: Build even if two seed contents share a label.
    #[staticmethod]
    #[pyo3(signature = (transcripts, output=None, *, t2g=None, infer_genes=None, scheme="lcp:K=3", hpc=true, mask_genes=20, collapse="gene", threads=0, allow_collisions=false))]
    #[allow(clippy::too_many_arguments)]
    fn build(
        py: Python<'_>,
        transcripts: &Bound<'_, PyAny>,
        output: Option<PathBuf>,
        t2g: Option<&Bound<'_, PyAny>>,
        infer_genes: Option<&Bound<'_, PyAny>>,
        scheme: &str,
        hpc: bool,
        mask_genes: u32,
        collapse: &str,
        threads: usize,
        allow_collisions: bool,
    ) -> PyResult<Self> {
        let per_gene = match collapse {
            "gene" => true,
            "global" => false,
            c => {
                return Err(PyValueError::new_err(format!(
                    "collapse must be \"gene\" or \"global\", not {c:?}"
                )));
            }
        };
        let infer = infer_genes_from(infer_genes)?;
        if infer.is_some() && t2g.is_some_and(|t| !t.is_none()) {
            return Err(PyValueError::new_err("give t2g or infer_genes, not both"));
        }
        // The gene table is small: read it first, so a bad path fails before the transcripts.
        let t2g_table: Option<(Option<PathBuf>, HashMap<String, String>)> = match t2g {
            None => None,
            Some(t) if t.is_none() => None,
            Some(t) => match as_path(t)? {
                Some(path) => {
                    let text = std::fs::read_to_string(&path)
                        .map_err(|e| io_err(&format!("reading {}", path.display()), e))?;
                    let map = text
                        .lines()
                        .filter_map(|l| {
                            let mut f = l.split_whitespace();
                            Some((f.next()?.to_string(), f.next()?.to_string()))
                        })
                        .collect();
                    Some((Some(path), map))
                }
                None => Some((None, t.extract::<HashMap<String, String>>()?)),
            },
        };
        let (mut txs, fasta) = transcripts_from(transcripts)?;
        let mut input = InputInfo {
            path: fasta.as_deref().map(absolute).unwrap_or_default(),
            genes_from: if fasta.is_some() { "fasta" } else { "given" }.into(),
            ..Default::default()
        };
        if let Some((path, map)) = &t2g_table {
            let mut missing = 0usize;
            for (name, gene, _) in txs.iter_mut() {
                match map.get(name) {
                    Some(g) => *gene = g.clone(),
                    None => missing += 1,
                }
            }
            if missing > 0 {
                log::warn!("{missing} transcripts not in t2g; they keep the genes they were given");
            }
            input.genes_from = "t2g".into();
            if let Some(p) = path {
                let hash = file_xxh3(p).map_err(|e| io_err(&format!("reading {}", p.display()), e))?;
                input.t2g = Some((absolute(p), hash));
            }
        }
        if infer.is_some() {
            input.genes_from = "inferred".into();
        } else if t2g_table.is_none() && txs.iter().all(|(n, g, _)| n == g) {
            input.genes_from = "none".into();
        }
        if infer.is_none() && per_gene && txs.len() > 1 && txs.iter().all(|(n, g, _)| n == g) {
            log::warn!(
                "every transcript is its own gene: isoforms are not collapsed, so the index is larger and mapping less accurate; give t2g or infer_genes"
            );
        }
        let opts = IndexOptions::builder()
            .scheme(scheme)
            .hpc(hpc)
            .f_max(mask_genes)
            .per_gene(per_gene)
            .allow_collisions(allow_collisions)
            .maybe_infer_genes(infer)
            .build();
        let pool = pool(threads)?;
        let (inner, report, path) = py
            .detach(|| -> anyhow::Result<(Index, BuildReport, Option<String>)> {
                pool.install(|| match &output {
                    Some(out) => {
                        let (_, rep, _) = Index::build_to_file_with(txs, &opts, out, input, Vec::new())?;
                        Ok((Index::open(out)?, rep, Some(absolute(out))))
                    }
                    None => {
                        let (idx, rep) = Index::build_with(txs, &opts, input)?;
                        Ok((idx, rep, None))
                    }
                })
            })
            .map_err(|e| value_err(format!("building the index: {e:#}")))?;
        Ok(PyIndex {
            inner: Arc::new(inner),
            path,
            report: Some(report),
        })
    }

    /// The index file, if it was loaded from or written to one.
    #[getter]
    fn path(&self) -> Option<String> {
        self.path.clone()
    }

    /// Number of transcripts.
    #[getter]
    fn num_transcripts(&self) -> usize {
        self.inner.n_transcripts()
    }

    /// Number of genes (the grouping isoforms are collapsed by).
    #[getter]
    fn num_genes(&self) -> usize {
        self.n_genes()
    }

    /// The base seed scheme the index was built with.
    #[getter]
    fn scheme(&self) -> String {
        self.inner.raw().meta("spec").unwrap_or("").to_string()
    }

    /// Whether transcripts (and reads) are homopolymer-compressed.
    #[getter]
    fn hpc(&self) -> bool {
        self.inner.raw().meta("hpc") == Some("true")
    }

    /// Whether the index holds raw (uncompressed) transcript sequence, for exact transcript
    /// coordinates and raw read ends.
    #[getter]
    fn has_raw_sequence(&self) -> bool {
        self.inner.has_raw_sequence()
    }

    /// A short identifier of the index's contents (stable across rebuilds of the same input).
    #[getter]
    fn id(&self) -> String {
        self.inner.id()
    }

    /// What built the index and from what, as a dict (``None`` for an index without one).
    #[getter]
    fn manifest<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.inner
            .manifest()
            .map(|m| pythonize::pythonize(py, &m).map_err(value_err))
            .transpose()
    }

    /// Construction statistics, for an index built in this process (otherwise ``None``; the
    /// manifest's ``stats`` holds the same for a loaded index).
    #[getter]
    fn build_report<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        self.report.as_ref().map(|r| report_dict(py, r)).transpose()
    }

    /// Name of transcript ``i``.
    fn transcript_name(&self, i: usize) -> PyResult<&str> {
        Ok(self.inner.tx_name(self.check_tx(i)?))
    }

    /// Length of transcript ``i`` in bases.
    fn transcript_length(&self, i: usize) -> PyResult<u32> {
        Ok(self.inner.tx_len(self.check_tx(i)?))
    }

    /// Gene of transcript ``i``.
    fn transcript_gene(&self, i: usize) -> PyResult<&str> {
        let t = self.check_tx(i)?;
        let raw = self.inner.raw();
        Ok(raw.gene_name(raw.tx_gene()[t as usize]))
    }

    /// Every transcript's name, in index order (transcript ids index this list).
    fn transcript_names(&self) -> Vec<&str> {
        (0..self.inner.n_transcripts() as u32)
            .map(|t| self.inner.tx_name(t))
            .collect()
    }

    /// Every transcript's length, as a NumPy ``uint32`` array.
    fn transcript_lengths<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<u32>> {
        self.inner.raw().tx_len().to_vec().into_pyarray(py)
    }

    /// Every transcript's gene id (an index into :meth:`gene_names`), as a NumPy ``uint32``
    /// array.
    fn transcript_gene_ids<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<u32>> {
        self.inner.raw().tx_gene().to_vec().into_pyarray(py)
    }

    /// Every gene's name, in gene id order.
    fn gene_names(&self) -> Vec<&str> {
        let raw = self.inner.raw();
        (0..self.n_genes() as u32).map(|g| raw.gene_name(g)).collect()
    }

    /// The transcript id of ``name``.
    ///
    /// :raises KeyError: If no transcript has that name.
    fn transcript_id(&self, name: &str) -> PyResult<u32> {
        (0..self.inner.n_transcripts() as u32)
            .find(|&t| self.inner.tx_name(t) == name)
            .ok_or_else(|| PyKeyError::new_err(name.to_string()))
    }

    /// A :class:`Mapper` over this index: ``index.mapper("hifi", extend_ends=True)`` is
    /// ``Mapper(index, "hifi", extend_ends=True)``.
    #[pyo3(signature = (preset="ont", *, threads=0, **settings))]
    fn mapper(
        slf: &Bound<'_, Self>,
        preset: &str,
        threads: usize,
        settings: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<PyMapper> {
        PyMapper::new(slf, preset, threads, settings)
    }

    fn __len__(&self) -> usize {
        self.inner.n_transcripts()
    }

    fn __repr__(&self) -> String {
        format!(
            "Index(transcripts={}, genes={}, scheme='{}', hpc={}, id='{}'{})",
            self.inner.n_transcripts(),
            self.n_genes(),
            self.scheme(),
            self.hpc(),
            self.inner.id(),
            self.path
                .as_ref()
                .map_or_else(String::new, |p| format!(", path='{p}'"))
        )
    }
}

// ─── Settings ───────────────────────────────────────────────────────────────────────────────

/// Python-friendly spellings of the settings whose serialized form is nested: `rescue_window`
/// as `"full"`, `"PAD,DIV"` or `(pad, div)`, and `chainer` as `"oracle"`, `"refine"` or
/// `("refine", beam)`.
fn normalize_setting(key: &str, v: serde_json::Value) -> serde_json::Value {
    use serde_json::{Value, json};
    match (key, v) {
        ("rescue_window", Value::String(s)) if s != "full" => match s.split_once(',') {
            Some((p, d)) => match (p.trim().parse::<u64>(), d.trim().parse::<u64>()) {
                (Ok(p), Ok(d)) => json!({"window": {"pad": p, "div": d}}),
                _ => Value::String(s),
            },
            None => Value::String(s),
        },
        ("rescue_window", Value::Array(a)) if a.len() == 2 => {
            json!({"window": {"pad": a[0], "div": a[1]}})
        }
        ("chainer", Value::String(s)) if s == "refine" => json!({"refine": {"beam": null}}),
        ("chainer", Value::Array(a)) if a.len() == 2 && a[0] == "refine" => {
            json!({"refine": {"beam": a[1]}})
        }
        (_, v) => v,
    }
}

/// Settings from a preset name and keyword arguments (each a `MapSettings` field).
fn settings_from(preset: &str, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<MapSettings> {
    let preset: Preset = preset.parse().map_err(value_err)?;
    let mut obj = serde_json::Map::new();
    obj.insert("preset".into(), serde_json::Value::String(preset.name().into()));
    if let Some(kw) = kwargs {
        for (k, v) in kw.iter() {
            let key: String = k.extract()?;
            if key == "preset" {
                return Err(PyTypeError::new_err("preset given twice"));
            }
            // `None` is the preset's value.
            if v.is_none() {
                continue;
            }
            let val: serde_json::Value = pythonize::depythonize(&v)
                .map_err(|e| PyValueError::new_err(format!("setting {key}: {e}")))?;
            let val = normalize_setting(&key, val);
            obj.insert(key, val);
        }
    }
    serde_json::from_value(serde_json::Value::Object(obj)).map_err(|e| {
        let msg = e.to_string();
        if msg.starts_with("unknown field") {
            PyTypeError::new_err(format!(
                "{msg} (spearfish.SETTINGS lists the settings; see spearfish.resolve_settings)"
            ))
        } else {
            PyValueError::new_err(format!("invalid setting: {msg}"))
        }
    })
}

/// The complete parameters a preset and settings resolve to, as a dict (what ``spearfish map
/// --print-params`` prints).
///
/// :param preset: ``"ont"``, ``"drna"`` or ``"hifi"``.
/// :param settings: Settings changed from the preset (see :data:`SETTINGS`).
#[pyfunction]
#[pyo3(signature = (preset="ont", **settings))]
fn resolve_settings<'py>(
    py: Python<'py>,
    preset: &str,
    settings: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyAny>> {
    let opts = settings_from(preset, settings)?.resolve().map_err(value_err)?;
    pythonize::pythonize(py, &opts).map_err(value_err)
}

/// The preset names.
#[pyfunction]
fn presets() -> Vec<&'static str> {
    Preset::ALL.iter().map(|p| p.name()).collect()
}

// ─── Results ────────────────────────────────────────────────────────────────────────────────

/// One reported transcript for a read.
///
/// Query coordinates are raw read bases on the read's forward strand; transcript coordinates
/// are raw transcript bases. Spans are half-open.
#[pyclass(name = "Hit", module = "spearfish", frozen)]
struct PyHit {
    index: Arc<Index>,
    /// Transcript id (index order; see :meth:`Index.transcript_names`).
    #[pyo3(get)]
    transcript_id: u32,
    /// Rank among the read's hits (0: best).
    #[pyo3(get)]
    rank: u32,
    /// The read's reverse complement matched the transcript.
    #[pyo3(get)]
    reverse: bool,
    /// Chain score after the tie-break (minimap2's alignment-score scale; higher is better).
    #[pyo3(get)]
    score: f32,
    /// Start of the query span.
    #[pyo3(get)]
    query_start: u32,
    /// End of the query span (exclusive).
    #[pyo3(get)]
    query_end: u32,
    /// Start of the transcript span.
    #[pyo3(get)]
    target_start: u32,
    /// End of the transcript span (exclusive).
    #[pyo3(get)]
    target_end: u32,
    /// Anchors in the chain.
    #[pyo3(get)]
    anchors: u32,
    /// Read bases the anchors cover (homopolymer-compressed).
    #[pyo3(get)]
    anchored: u32,
    /// Largest gap discrepancy between consecutive anchors.
    #[pyo3(get)]
    max_discrepancy: f32,
    /// The span reaches the read's alignment ends; otherwise it runs between the outermost
    /// anchors.
    #[pyo3(get)]
    extended: bool,
}

#[pymethods]
impl PyHit {
    /// Transcript name.
    #[getter]
    fn transcript(&self) -> &str {
        self.index.tx_name(self.transcript_id)
    }

    /// Gene of the transcript.
    #[getter]
    fn gene(&self) -> &str {
        let raw = self.index.raw();
        raw.gene_name(raw.tx_gene()[self.transcript_id as usize])
    }

    /// Transcript length in bases.
    #[getter]
    fn transcript_length(&self) -> u32 {
        self.index.tx_len(self.transcript_id)
    }

    /// ``"+"`` or ``"-"``.
    #[getter]
    fn strand(&self) -> &'static str {
        if self.reverse { "-" } else { "+" }
    }

    fn __repr__(&self) -> String {
        format!(
            "Hit(transcript='{}', strand='{}', score={:.2}, query={}-{}, target={}-{}, rank={})",
            self.transcript(),
            self.strand(),
            self.score,
            self.query_start,
            self.query_end,
            self.target_start,
            self.target_end,
            self.rank
        )
    }
}

/// All hits of one read, best first.
#[pyclass(name = "MappedRead", module = "spearfish", frozen)]
struct PyMappedRead {
    /// Read name, if known (from the file, or given).
    #[pyo3(get)]
    name: Option<String>,
    /// Read length in bases.
    #[pyo3(get)]
    length: u32,
    /// The reported hits, best first.
    #[pyo3(get)]
    hits: Py<PyTuple>,
    /// The best chain spans less than half the read (spearfish's own EM skips such reads).
    #[pyo3(get)]
    low_coverage: bool,
    /// Hits came from the lower-level rescue.
    #[pyo3(get)]
    rescued: bool,
    scores: Vec<f32>,
}

fn to_py_read(
    py: Python<'_>,
    index: &Arc<Index>,
    name: Option<String>,
    r: &MappedRead,
) -> PyResult<PyMappedRead> {
    let hits: Vec<Bound<'_, PyHit>> = r
        .hits
        .iter()
        .enumerate()
        .map(|(rank, h)| {
            Bound::new(
                py,
                PyHit {
                    index: index.clone(),
                    transcript_id: h.tx,
                    rank: rank as u32,
                    reverse: h.reverse,
                    score: h.score,
                    query_start: h.query_start,
                    query_end: h.query_end,
                    target_start: h.target_start,
                    target_end: h.target_end,
                    anchors: h.anchors,
                    anchored: h.anchored,
                    max_discrepancy: h.max_discrepancy,
                    extended: h.extended,
                },
            )
        })
        .collect::<PyResult<_>>()?;
    Ok(PyMappedRead {
        name,
        length: r.len,
        hits: PyTuple::new(py, hits)?.unbind(),
        low_coverage: r.low_coverage,
        rescued: r.rescued,
        scores: r.hits.iter().map(|h| h.score).collect(),
    })
}

#[pymethods]
impl PyMappedRead {
    /// Whether any transcript was reported.
    #[getter]
    fn is_mapped(&self) -> bool {
        !self.scores.is_empty()
    }

    /// Number of hits.
    #[getter]
    fn num_hits(&self) -> usize {
        self.scores.len()
    }

    /// The best hit, or ``None`` if the read is unmapped.
    #[getter]
    fn best<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        let hits = self.hits.bind(py);
        if hits.is_empty() {
            Ok(None)
        } else {
            hits.get_item(0).map(Some)
        }
    }

    /// The hits' probabilities: ``exp(-(best - score) / SCORE_DENOM)``, normalized to sum to 1
    /// (empty for an unmapped read).
    fn probabilities(&self) -> Vec<f64> {
        let Some(best) = self.scores.first().copied() else {
            return Vec::new();
        };
        let w: Vec<f64> = self
            .scores
            .iter()
            .map(|&s| (-((best - s) as f64) / spearfish::SCORE_DENOM as f64).exp())
            .collect();
        let z: f64 = w.iter().sum();
        w.into_iter().map(|x| x / z).collect()
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        Ok(self.hits.bind(py).as_any().try_iter()?.into_any())
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let best = match self.best(py)? {
            Some(h) => format!(", best={}", h.repr()?),
            None => String::new(),
        };
        Ok(format!(
            "MappedRead(name={}, length={}, num_hits={}{best})",
            self.name.as_deref().map_or("None".into(), |n| format!("'{n}'")),
            self.length,
            self.scores.len()
        ))
    }
}

/// The results of mapping a batch of reads, in input order, held in Rust until asked for.
///
/// Index it or iterate it for :class:`MappedRead` objects, or take whole columns as NumPy
/// arrays with :meth:`hit_columns` and :meth:`read_columns` (no per-hit Python objects); the
/// Python package adds ``to_pandas()``.
#[pyclass(name = "MappingBatch", module = "spearfish", frozen, sequence)]
struct PyBatch {
    index: Arc<Index>,
    index_obj: Py<PyIndex>,
    names: Option<Vec<String>>,
    reads: Vec<MappedRead>,
}

impl PyBatch {
    fn read_at(&self, py: Python<'_>, i: usize) -> PyResult<PyMappedRead> {
        let name = self.names.as_ref().map(|n| n[i].clone());
        to_py_read(py, &self.index, name, &self.reads[i])
    }
}

#[pymethods]
impl PyBatch {
    fn __len__(&self) -> usize {
        self.reads.len()
    }

    fn __getitem__(&self, py: Python<'_>, i: isize) -> PyResult<PyMappedRead> {
        let n = self.reads.len() as isize;
        let j = if i < 0 { i + n } else { i };
        if j < 0 || j >= n {
            return Err(PyIndexError::new_err(format!(
                "read {i} out of range ({n} reads)"
            )));
        }
        self.read_at(py, j as usize)
    }

    fn __iter__(slf: &Bound<'_, Self>) -> PyBatchIter {
        PyBatchIter {
            batch: slf.clone().unbind(),
            next: Mutex::new(0),
        }
    }

    /// The index the reads were mapped against.
    #[getter]
    fn index(&self, py: Python<'_>) -> Py<PyIndex> {
        self.index_obj.clone_ref(py)
    }

    /// Read names (``None`` if the reads were given without names).
    #[getter]
    fn names(&self) -> Option<Vec<String>> {
        self.names.clone()
    }

    /// Total hits over all reads.
    #[getter]
    fn num_hits(&self) -> usize {
        self.reads.iter().map(|r| r.hits.len()).sum()
    }

    /// Reads with at least one hit.
    #[getter]
    fn num_mapped(&self) -> usize {
        self.reads.iter().filter(|r| !r.hits.is_empty()).count()
    }

    /// One row per hit, as a dict of NumPy arrays: ``read`` (the read's position in the
    /// batch), ``rank``, ``transcript_id``, ``reverse``, ``score``, ``query_start``,
    /// ``query_end``, ``target_start``, ``target_end``, ``anchors``, ``anchored``,
    /// ``max_discrepancy`` and ``extended``. A read's hits are contiguous, best first;
    /// unmapped reads have no rows.
    fn hit_columns<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let n = self.num_hits();
        let mut read = Vec::with_capacity(n);
        let mut rank = Vec::with_capacity(n);
        let mut tx = Vec::with_capacity(n);
        let mut reverse = Vec::with_capacity(n);
        let mut score = Vec::with_capacity(n);
        let mut qs = Vec::with_capacity(n);
        let mut qe = Vec::with_capacity(n);
        let mut ts = Vec::with_capacity(n);
        let mut te = Vec::with_capacity(n);
        let mut anchors = Vec::with_capacity(n);
        let mut anchored = Vec::with_capacity(n);
        let mut disc = Vec::with_capacity(n);
        let mut extended = Vec::with_capacity(n);
        for (i, r) in self.reads.iter().enumerate() {
            for (k, h) in r.hits.iter().enumerate() {
                read.push(i as u32);
                rank.push(k as u32);
                tx.push(h.tx);
                reverse.push(h.reverse);
                score.push(h.score);
                qs.push(h.query_start);
                qe.push(h.query_end);
                ts.push(h.target_start);
                te.push(h.target_end);
                anchors.push(h.anchors);
                anchored.push(h.anchored);
                disc.push(h.max_discrepancy);
                extended.push(h.extended);
            }
        }
        let d = PyDict::new(py);
        d.set_item("read", read.into_pyarray(py))?;
        d.set_item("rank", rank.into_pyarray(py))?;
        d.set_item("transcript_id", tx.into_pyarray(py))?;
        d.set_item("reverse", reverse.into_pyarray(py))?;
        d.set_item("score", score.into_pyarray(py))?;
        d.set_item("query_start", qs.into_pyarray(py))?;
        d.set_item("query_end", qe.into_pyarray(py))?;
        d.set_item("target_start", ts.into_pyarray(py))?;
        d.set_item("target_end", te.into_pyarray(py))?;
        d.set_item("anchors", anchors.into_pyarray(py))?;
        d.set_item("anchored", anchored.into_pyarray(py))?;
        d.set_item("max_discrepancy", disc.into_pyarray(py))?;
        d.set_item("extended", extended.into_pyarray(py))?;
        Ok(d)
    }

    /// One row per read, as a dict of NumPy arrays: ``length``, ``num_hits``,
    /// ``low_coverage``, ``rescued``, and ``hit_offset`` (one longer than the batch: read
    /// ``i``'s rows in :meth:`hit_columns` are ``hit_offset[i]:hit_offset[i + 1]``).
    fn read_columns<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let n = self.reads.len();
        let mut len = Vec::with_capacity(n);
        let mut nh = Vec::with_capacity(n);
        let mut low = Vec::with_capacity(n);
        let mut resc = Vec::with_capacity(n);
        let mut off = Vec::with_capacity(n + 1);
        let mut o = 0u64;
        off.push(0u64);
        for r in &self.reads {
            len.push(r.len);
            nh.push(r.hits.len() as u32);
            low.push(r.low_coverage);
            resc.push(r.rescued);
            o += r.hits.len() as u64;
            off.push(o);
        }
        let d = PyDict::new(py);
        d.set_item("length", len.into_pyarray(py))?;
        d.set_item("num_hits", nh.into_pyarray(py))?;
        d.set_item("low_coverage", low.into_pyarray(py))?;
        d.set_item("rescued", resc.into_pyarray(py))?;
        d.set_item("hit_offset", off.into_pyarray(py))?;
        Ok(d)
    }

    /// One row per hit as a :class:`pandas.DataFrame`, like ``spearfish map``'s TSV: read
    /// (name, or position in the batch), read_length, rank, transcript (categorical), strand,
    /// score, spans, anchors, anchored, max_discrepancy, extended, low_coverage, rescued.
    ///
    /// :param unmapped: Also give each unmapped read one row, with no transcript.
    #[pyo3(signature = (*, unmapped=false))]
    fn to_pandas<'py>(slf: &Bound<'py, Self>, unmapped: bool) -> PyResult<Bound<'py, PyAny>> {
        slf.py()
            .import("spearfish._frames")?
            .call_method1("batch_to_pandas", (slf, unmapped))
    }

    fn __repr__(&self) -> String {
        format!(
            "MappingBatch(reads={}, mapped={}, hits={})",
            self.reads.len(),
            self.num_mapped(),
            self.num_hits()
        )
    }
}

#[pyclass(name = "_MappingBatchIter", module = "spearfish", frozen)]
struct PyBatchIter {
    batch: Py<PyBatch>,
    next: Mutex<usize>,
}

#[pymethods]
impl PyBatchIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self, py: Python<'_>) -> PyResult<Option<PyMappedRead>> {
        let b = self.batch.get();
        let mut next = self.next.lock().unwrap();
        if *next >= b.reads.len() {
            return Ok(None);
        }
        let r = b.read_at(py, *next)?;
        *next += 1;
        Ok(Some(r))
    }
}

// ─── Mapper ─────────────────────────────────────────────────────────────────────────────────

/// Maps long reads to an :class:`Index` with a preset and settings.
///
/// ``Mapper(index, "hifi", extend_ends=True, tie_kappa=5.0)``: the preset (``"ont"``,
/// ``"drna"`` or ``"hifi"``) sets every default, and keyword arguments change single settings,
/// as ``spearfish map``'s flags do (:data:`SETTINGS` lists them; ``None`` keeps the preset's
/// value). A mapper is immutable and thread-safe: :meth:`map` releases the GIL, and
/// :meth:`map_batch` and :meth:`map_file` map on the mapper's own pool of ``threads`` threads.
#[pyclass(name = "Mapper", module = "spearfish", frozen)]
struct PyMapper {
    inner: Arc<spearfish::Mapper>,
    index: Py<PyIndex>,
    pool: Arc<rayon::ThreadPool>,
    threads: usize,
    settings: MapSettings,
}

/// Reads `map_batch` maps between checks for Python signals.
const SIGNAL_CHECK_READS: usize = 1 << 16;

/// Map `seqs` in parallel on `pool`, results in input order.
fn map_on(mapper: &spearfish::Mapper, pool: &rayon::ThreadPool, seqs: &[&[u8]]) -> Vec<MappedRead> {
    pool.install(|| seqs.par_iter().map(|s| mapper.map(s)).collect())
}

impl PyMapper {
    fn map_seqs(&self, seqs: &[&[u8]]) -> Vec<MappedRead> {
        map_on(&self.inner, &self.pool, seqs)
    }

    /// A thread that maps chunks in the order they are sent, one at a time on the mapper's
    /// pool, so that a chunk maps while Python works on the one before. Closing the job
    /// channel ends it, as does dropping the results.
    fn spawn_pipeline(&self) -> std::io::Result<(JobSender, ResultReceiver)> {
        let (jobs_tx, jobs_rx) = std::sync::mpsc::sync_channel::<anyhow::Result<Job>>(1);
        let (res_tx, res_rx) = std::sync::mpsc::sync_channel::<Mapped>(1);
        let (mapper, pool) = (self.inner.clone(), self.pool.clone());
        std::thread::Builder::new()
            .name("spearfish-mapper".into())
            .spawn(move || {
                for job in jobs_rx {
                    let out = job.map(|j| {
                        let seqs: Vec<&[u8]> = j.seqs.iter().map(Vec::as_slice).collect();
                        let reads = map_on(&mapper, &pool, &seqs);
                        (j.names, reads)
                    });
                    let failed = out.is_err();
                    if res_tx.send(out).is_err() || failed {
                        break;
                    }
                }
            })?;
        Ok((jobs_tx, res_rx))
    }

    fn batch(&self, py: Python<'_>, names: Option<Vec<String>>, reads: Vec<MappedRead>) -> PyBatch {
        PyBatch {
            index: self.inner.index().clone(),
            index_obj: self.index.clone_ref(py),
            names,
            reads,
        }
    }
}

#[pymethods]
impl PyMapper {
    /// :param index: The :class:`Index` to map against.
    /// :param preset: ``"ont"`` (ONT cDNA, the default), ``"drna"`` (ONT direct RNA) or
    ///     ``"hifi"`` (PacBio HiFi).
    /// :param threads: Threads for :meth:`map_batch` and :meth:`map_file` (0: every
    ///     available CPU).
    /// :param settings: Settings changed from the preset.
    /// :raises TypeError: For an unknown setting.
    /// :raises ValueError: For an invalid value, or settings that contradict each other.
    #[new]
    #[pyo3(signature = (index, preset="ont", *, threads=0, **settings))]
    fn new(
        index: &Bound<'_, PyIndex>,
        preset: &str,
        threads: usize,
        settings: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Self> {
        let settings = settings_from(preset, settings)?;
        let opts = settings.resolve().map_err(value_err)?;
        let inner = spearfish::Mapper::new(index.get().inner.clone(), &opts).map_err(value_err)?;
        let n = self::threads(threads);
        Ok(PyMapper {
            inner: Arc::new(inner),
            index: index.clone().unbind(),
            pool: Arc::new(pool(n)?),
            threads: n,
            settings,
        })
    }

    /// Map one read.
    ///
    /// :param seq: The read's bases (``str`` or ``bytes``; ``U`` reads as ``T``).
    /// :param name: A name to carry on the result.
    /// :returns: A :class:`MappedRead`.
    #[pyo3(signature = (seq, name=None))]
    fn map(&self, py: Python<'_>, seq: &Bound<'_, PyAny>, name: Option<String>) -> PyResult<PyMappedRead> {
        let b = seq_bytes(seq)?;
        let s = b.as_bytes();
        let mapper = &self.inner;
        let r = py.detach(|| mapper.map(s));
        to_py_read(py, mapper.index(), name, &r)
    }

    /// Map a batch of reads in parallel, results in input order.
    ///
    /// :param seqs: The reads' bases (an iterable of ``str`` or ``bytes``).
    /// :param names: Their names, as many as there are reads (optional).
    /// :returns: A :class:`MappingBatch`.
    #[pyo3(signature = (seqs, names=None))]
    fn map_batch(
        &self,
        py: Python<'_>,
        seqs: &Bound<'_, PyAny>,
        names: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PyBatch> {
        let owned: Vec<Bound<'_, PyBytes>> = seqs
            .try_iter()?
            .map(|s| seq_bytes(&s?))
            .collect::<PyResult<_>>()?;
        let names: Option<Vec<String>> = match names {
            Some(n) if !n.is_none() => {
                let v = n
                    .try_iter()?
                    .map(|x| x?.extract::<String>())
                    .collect::<PyResult<Vec<String>>>()?;
                if v.len() != owned.len() {
                    return Err(PyValueError::new_err(format!(
                        "{} names for {} reads",
                        v.len(),
                        owned.len()
                    )));
                }
                Some(v)
            }
            _ => None,
        };
        let slices: Vec<&[u8]> = owned.iter().map(|b| b.as_bytes()).collect();
        // In parts, so Ctrl-C (and other signals) can interrupt a large batch.
        let mut reads = Vec::with_capacity(slices.len());
        for part in slices.chunks(SIGNAL_CHECK_READS) {
            reads.extend(py.detach(|| self.map_seqs(part)));
            py.check_signals()?;
        }
        Ok(self.batch(py, names, reads))
    }

    /// Map a FASTA/FASTQ file (optionally gzip/bzip2/zstd-compressed; ``"-"`` reads stdin)
    /// in chunks, read and mapped in Rust.
    ///
    /// :param path: The reads.
    /// :param chunk_size: Reads per :class:`MappingBatch`.
    /// :returns: An iterator of :class:`MappingBatch`, in file order, with read names.
    #[pyo3(signature = (path, *, chunk_size=20_000))]
    fn map_file(slf: &Bound<'_, Self>, path: PathBuf, chunk_size: usize) -> PyResult<PyFileBatches> {
        if chunk_size == 0 {
            return Err(PyValueError::new_err("chunk_size must be at least 1"));
        }
        // As `spearfish map` opens its reads: compressed input decompresses on a thread of its
        // own (ahead of the reader thread below), and xz input is refused with a clear message.
        let reader = slf
            .py()
            .detach(|| spearfish::open_fastx(Some(&path), None))
            .map_err(|e| io_err(&format!("opening {}", path.display()), e))?;
        let (jobs, results) = slf
            .get()
            .spawn_pipeline()
            .map_err(|e| io_err("starting the mapping thread", e))?;
        spawn_reader(reader, chunk_size, jobs).map_err(|e| io_err("starting the reader thread", e))?;
        Ok(PyFileBatches {
            mapper: slf.clone().unbind(),
            path,
            results: Mutex::new(results),
        })
    }

    /// Map reads from any iterable, in chunks: each item is a sequence (``str`` or ``bytes``)
    /// or a tuple whose first two fields are the name and the sequence, as
    /// ``mappy.fastx_read`` and ``pyfastx`` give them.
    ///
    /// :param records: The reads.
    /// :param chunk_size: Reads per :class:`MappingBatch`.
    /// :returns: An iterator of :class:`MappingBatch`, in input order.
    #[pyo3(signature = (records, *, chunk_size=20_000))]
    fn map_records(
        slf: &Bound<'_, Self>,
        records: &Bound<'_, PyAny>,
        chunk_size: usize,
    ) -> PyResult<PyRecordBatches> {
        if chunk_size == 0 {
            return Err(PyValueError::new_err("chunk_size must be at least 1"));
        }
        let records = records.try_iter()?.unbind();
        let (jobs, results) = slf
            .get()
            .spawn_pipeline()
            .map_err(|e| io_err("starting the mapping thread", e))?;
        Ok(PyRecordBatches {
            mapper: slf.clone().unbind(),
            records,
            chunk_size,
            feed: Mutex::new(Feed {
                jobs: Some(jobs),
                outstanding: 0,
            }),
            results: Mutex::new(results),
        })
    }

    /// The index mapped against.
    #[getter]
    fn index(&self, py: Python<'_>) -> Py<PyIndex> {
        self.index.clone_ref(py)
    }

    /// The preset name.
    #[getter]
    fn preset(&self) -> &'static str {
        self.settings.preset.name()
    }

    /// Threads used by :meth:`map_batch` and :meth:`map_file`.
    #[getter]
    fn threads(&self) -> usize {
        self.threads
    }

    /// The settings changed from the preset, as a dict.
    #[getter]
    fn settings<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let d = pythonize::pythonize(py, &self.settings).map_err(value_err)?;
        d.del_item("preset")?;
        Ok(d)
    }

    /// Every parameter the mapper runs with, as a dict (see :func:`resolve_settings`).
    #[getter]
    fn options<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        pythonize::pythonize(py, self.inner.options()).map_err(value_err)
    }

    /// The scale of score differences: a hit ``d`` below its read's best is
    /// ``exp(-d / score_denom)`` times as likely.
    #[getter]
    fn score_denom(&self) -> f32 {
        spearfish::SCORE_DENOM
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let s = serde_json::to_value(&self.settings).map_err(value_err)?;
        let changed: Vec<String> = s
            .as_object()
            .into_iter()
            .flatten()
            .filter(|(k, _)| *k != "preset")
            .map(|(k, v)| format!(", {k}={v}"))
            .collect();
        Ok(format!(
            "Mapper(index={}, preset='{}', threads={}{})",
            self.index.bind(py).repr()?,
            self.preset(),
            self.threads,
            changed.concat()
        ))
    }
}

/// A chunk of reads to map: their names (if known) and sequences.
struct Job {
    names: Option<Vec<String>>,
    seqs: Vec<Vec<u8>>,
}

/// A mapped chunk: the names it came with and the results, in order.
type Mapped = anyhow::Result<(Option<Vec<String>>, Vec<MappedRead>)>;
type JobSender = std::sync::mpsc::SyncSender<anyhow::Result<Job>>;
type ResultReceiver = std::sync::mpsc::Receiver<Mapped>;

/// Reads `reader` in chunks of `n` records on a thread of its own and sends them to `jobs`,
/// so reading (and decompressing) overlaps mapping. The thread ends at the end of the file,
/// after an error (which it sends), or when the mapping side hangs up.
fn spawn_reader(
    mut reader: paraseq::fastx::Reader<paraseq::BoxedReader>,
    n: usize,
    jobs: JobSender,
) -> std::io::Result<()> {
    std::thread::Builder::new()
        .name("spearfish-reader".into())
        .spawn(move || {
            use paraseq::Record;
            let mut rset = reader.new_record_set();
            let mut pending: std::collections::VecDeque<(String, Vec<u8>)> = Default::default();
            let mut eof = false;
            loop {
                let chunk = (|| -> anyhow::Result<Vec<(String, Vec<u8>)>> {
                    while pending.len() < n && !eof {
                        if !rset.fill(&mut reader)? {
                            eof = true;
                            break;
                        }
                        for rec in rset.iter() {
                            let rec = rec?;
                            pending.push_back((read_name(rec.id()), rec.seq().into_owned()));
                        }
                    }
                    let k = n.min(pending.len());
                    Ok(pending.drain(..k).collect())
                })();
                let job = match chunk {
                    Ok(c) if c.is_empty() => break,
                    Ok(c) => {
                        let (names, seqs) = c.into_iter().unzip();
                        Ok(Job {
                            names: Some(names),
                            seqs,
                        })
                    }
                    Err(e) => Err(e),
                };
                let failed = job.is_err();
                if jobs.send(job).is_err() || failed {
                    break;
                }
            }
        })?;
    Ok(())
}

/// The next mapped chunk, waiting without the GIL but checking for signals (Ctrl-C) every
/// 100 ms; `None` once the pipeline has finished.
fn next_mapped(py: Python<'_>, results: &Mutex<ResultReceiver>) -> PyResult<Option<Mapped>> {
    use std::sync::mpsc::RecvTimeoutError;
    loop {
        let got = py.detach(|| {
            results
                .lock()
                .unwrap()
                .recv_timeout(std::time::Duration::from_millis(100))
        });
        match got {
            Ok(m) => return Ok(Some(m)),
            Err(RecvTimeoutError::Timeout) => py.check_signals()?,
            Err(RecvTimeoutError::Disconnected) => return Ok(None),
        }
    }
}

/// Batches of a mapped read file (see :meth:`Mapper.map_file`).
#[pyclass(name = "FileBatches", module = "spearfish", frozen)]
struct PyFileBatches {
    mapper: Py<PyMapper>,
    path: PathBuf,
    results: Mutex<ResultReceiver>,
}

#[pymethods]
impl PyFileBatches {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self, py: Python<'_>) -> PyResult<Option<PyBatch>> {
        match next_mapped(py, &self.results)? {
            Some(Ok((names, reads))) => Ok(Some(self.mapper.get().batch(py, names, reads))),
            Some(Err(e)) => Err(io_err(&format!("reading {}", self.path.display()), e)),
            None => Ok(None),
        }
    }
}

/// The Python side of a `map_records` pipeline: where to send chunks, and how many are being
/// mapped or wait to be collected.
struct Feed {
    /// `None` once the input is exhausted (or failed).
    jobs: Option<JobSender>,
    outstanding: usize,
}

/// Chunks `map_records` keeps in the pipeline: one returned to Python while the next maps.
const RECORDS_AHEAD: usize = 2;

/// Batches of mapped records from an iterable (see :meth:`Mapper.map_records`).
#[pyclass(name = "RecordBatches", module = "spearfish", frozen)]
struct PyRecordBatches {
    mapper: Py<PyMapper>,
    records: Py<pyo3::types::PyIterator>,
    chunk_size: usize,
    feed: Mutex<Feed>,
    results: Mutex<ResultReceiver>,
}

impl PyRecordBatches {
    /// The next chunk of the input, copied out of Python (it is mapped on another thread);
    /// `None` at the end.
    fn pull(&self, py: Python<'_>) -> PyResult<Option<Job>> {
        let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(self.chunk_size);
        let mut names: Vec<Option<String>> = Vec::with_capacity(self.chunk_size);
        let mut it = self.records.bind(py).clone();
        for item in it.by_ref().take(self.chunk_size) {
            let item = item?;
            match item.cast::<PyTuple>() {
                Ok(t) if t.len() >= 2 => {
                    let name: String = t.get_item(0)?.extract()?;
                    names.push(Some(read_name(name.as_bytes())));
                    seqs.push(seq_bytes(&t.get_item(1)?)?.as_bytes().to_vec());
                }
                Ok(_) => {
                    return Err(PyTypeError::new_err(
                        "a record tuple needs at least (name, sequence)",
                    ));
                }
                Err(_) => {
                    names.push(None);
                    seqs.push(seq_bytes(&item)?.as_bytes().to_vec());
                }
            }
        }
        if seqs.is_empty() {
            return Ok(None);
        }
        // Names only if every read in the chunk has one.
        let names = names.into_iter().collect::<Option<Vec<String>>>();
        Ok(Some(Job { names, seqs }))
    }
}

#[pymethods]
impl PyRecordBatches {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self, py: Python<'_>) -> PyResult<Option<PyBatch>> {
        {
            // Held while Python code runs (the input) and while the GIL is released, so a
            // second consumer must not wait for it: that could deadlock on the GIL.
            let mut feed = self.feed.try_lock().map_err(|_| {
                pyo3::exceptions::PyRuntimeError::new_err(
                    "this RecordBatches is already being iterated (from another thread)",
                )
            })?;
            // Keep the pipeline full: the chunk returned now, and the next one mapping.
            while feed.outstanding < RECORDS_AHEAD && feed.jobs.is_some() {
                match self.pull(py) {
                    Ok(Some(job)) => {
                        let tx = feed.jobs.as_ref().unwrap();
                        // Never blocks for long: at most RECORDS_AHEAD chunks are in flight,
                        // fewer than the pipeline holds.
                        if py.detach(|| tx.send(Ok(job))).is_err() {
                            feed.jobs = None;
                        } else {
                            feed.outstanding += 1;
                        }
                    }
                    Ok(None) => feed.jobs = None,
                    Err(e) => {
                        // A failed input ends the iteration; chunks already sent are dropped.
                        feed.jobs = None;
                        feed.outstanding = 0;
                        return Err(e);
                    }
                }
            }
            if feed.outstanding == 0 {
                return Ok(None);
            }
            feed.outstanding -= 1;
        }
        match next_mapped(py, &self.results)? {
            Some(Ok((names, reads))) => Ok(Some(self.mapper.get().batch(py, names, reads))),
            Some(Err(e)) => Err(value_err(format!("{e:#}"))),
            None => Ok(None),
        }
    }
}

// ─── Module ─────────────────────────────────────────────────────────────────────────────────

#[pymodule]
fn _spearfish(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // spearfish's warnings (e.g. settings that have no effect) go to Python's `logging`.
    let _ = pyo3_log::try_init();
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("SCORE_DENOM", spearfish::SCORE_DENOM)?;
    let fields: Vec<&str> = MapSettings::FIELDS
        .iter()
        .copied()
        .filter(|f| *f != "preset")
        .collect();
    m.add("SETTINGS", PyTuple::new(m.py(), fields)?)?;
    m.add_class::<PyIndex>()?;
    m.add_class::<PyMapper>()?;
    m.add_class::<PyMappedRead>()?;
    m.add_class::<PyHit>()?;
    m.add_class::<PyBatch>()?;
    m.add_class::<PyFileBatches>()?;
    m.add_class::<PyRecordBatches>()?;
    m.add_function(wrap_pyfunction!(resolve_settings, m)?)?;
    m.add_function(wrap_pyfunction!(presets, m)?)?;
    Ok(())
}
