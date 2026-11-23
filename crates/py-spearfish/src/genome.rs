//! Spliced mapping to a genome: `GenomeIndex`, `SplicedMapper`, `SplicedRead`,
//! `SplicedAlignment`.

use std::path::PathBuf;
use std::sync::Arc;

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict};
use rayon::prelude::*;

use spearfish::Preset;
use spearfish::genome::{
    FileOutputs, GenomeIndex, GenomeIndexOptions, SpliceOptions, SplicedMapper, SplicedRead,
    TranscriptStrand, output,
};

use crate::{SIGNAL_CHECK_READS, io_err, pool, seq_bytes, value_err};

/// A spearfish genome index, for spliced mapping.
///
/// Open one written by ``spearfish index --genome`` (or :meth:`GenomeIndex.build` with
/// ``output``) with :meth:`GenomeIndex.load`, or build one from a genome FASTA with
/// :meth:`GenomeIndex.build`. Known transcripts are optional: ``annotation`` names a GTF,
/// GFF3 or BED12 file whose splice sites, junctions and exons are then preferred where a
/// read fits them. The index file stores no annotation. An index is immutable and can be
/// shared by any number of mappers and threads.
#[pyclass(name = "GenomeIndex", module = "spearfish", frozen)]
pub struct PyGenomeIndex {
    inner: GenomeIndex,
    path: Option<PathBuf>,
}

fn with_annotation(mut index: GenomeIndex, annotation: Option<PathBuf>) -> PyResult<GenomeIndex> {
    if let Some(a) = annotation {
        index
            .load_annotation(&a)
            .map_err(|e| PyValueError::new_err(format!("{e:#}")))?;
    }
    Ok(index)
}

#[pymethods]
impl PyGenomeIndex {
    /// Open a genome index file (memory-mapped).
    ///
    /// :param path: The index.
    /// :param annotation: Known transcripts for this run (GTF, GFF3 or BED12, optionally
    ///     compressed); sequence names must match the genome's.
    /// :raises OSError: If the file can't be read or is not a genome index.
    /// :raises ValueError: If no transcript of the annotation is on a sequence of the index.
    #[staticmethod]
    #[pyo3(signature = (path, *, annotation=None))]
    fn load(py: Python<'_>, path: PathBuf, annotation: Option<PathBuf>) -> PyResult<Self> {
        let inner = py
            .detach(|| GenomeIndex::open(&path))
            .map_err(|e| io_err(&format!("index {}", path.display()), e))?;
        let inner = py.detach(|| with_annotation(inner, annotation))?;
        Ok(PyGenomeIndex {
            inner,
            path: Some(path),
        })
    }

    /// Build a genome index from a FASTA file (optionally gzip, bzip2 or zstd compressed).
    ///
    /// :param fasta: The genome.
    /// :param output: Also write the index (and its manifest, ``<output>.json``) there.
    /// :param annotation: As for :meth:`load`.
    /// :param max_occ: Seeds that occur more often than this (repeats) are left out.
    /// :param threads: Threads for the build (0: every available CPU).
    #[staticmethod]
    #[pyo3(signature = (fasta, *, output=None, annotation=None, max_occ=1000, threads=0))]
    fn build(
        py: Python<'_>,
        fasta: PathBuf,
        output: Option<PathBuf>,
        annotation: Option<PathBuf>,
        max_occ: u32,
        threads: usize,
    ) -> PyResult<Self> {
        let opts = GenomeIndexOptions::builder().max_occ(max_occ).build();
        let pool = pool(threads)?;
        let inner = py
            .detach(|| {
                pool.install(|| match &output {
                    Some(out) => GenomeIndex::build_to_file(&fasta, &opts, out, Vec::new()).map(|r| r.0),
                    None => GenomeIndex::build_from_fasta(&fasta, &opts),
                })
            })
            .map_err(|e| io_err(&format!("building from {}", fasta.display()), e))?;
        let inner = py.detach(|| with_annotation(inner, annotation))?;
        Ok(PyGenomeIndex { inner, path: output })
    }

    /// The file the index was loaded from or written to, if any.
    #[getter]
    fn path(&self) -> Option<String> {
        self.path.as_ref().map(|p| p.display().to_string())
    }

    /// Number of reference sequences.
    #[getter]
    fn num_sequences(&self) -> usize {
        self.inner.n_sequences()
    }

    /// The sequences' names, in index order.
    #[getter]
    fn sequence_names(&self) -> Vec<&str> {
        (0..self.inner.n_sequences() as u32)
            .map(|i| self.inner.seq_name(i))
            .collect()
    }

    /// The sequences' lengths in bases, in index order.
    #[getter]
    fn sequence_lengths(&self) -> Vec<u32> {
        (0..self.inner.n_sequences() as u32)
            .map(|i| self.inner.seq_len(i))
            .collect()
    }

    /// Bases in all sequences.
    #[getter]
    fn total_length(&self) -> u64 {
        self.inner.total_len()
    }

    /// Identifier of the index's contents (``"unknown"`` for an index built in memory).
    #[getter]
    fn id(&self) -> String {
        self.inner.id()
    }

    /// What built the index and from what (a ``dict``), if the index records it.
    #[getter]
    fn manifest<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.inner
            .manifest()
            .map(|m| pythonize::pythonize(py, &m).map_err(value_err))
            .transpose()
    }

    /// What the loaded annotation held (``format``, ``transcripts``, ``introns``, ``exons``,
    /// ``skipped_records``), or ``None`` without one.
    #[getter]
    fn annotation<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        let Some(a) = self.inner.annotation() else {
            return Ok(None);
        };
        let d = PyDict::new(py);
        d.set_item("format", a.format)?;
        d.set_item("transcripts", a.transcripts)?;
        d.set_item("introns", a.introns)?;
        d.set_item("exons", a.exons)?;
        d.set_item("skipped_records", a.skipped_records)?;
        Ok(Some(d))
    }

    fn __len__(&self) -> usize {
        self.inner.n_sequences()
    }

    fn __repr__(&self) -> String {
        format!(
            "GenomeIndex(sequences={}, bases={}, annotation={})",
            self.inner.n_sequences(),
            self.inner.total_len(),
            if self.inner.annotation().is_some() {
                "yes"
            } else {
                "no"
            }
        )
    }
}

/// One alignment of a read to the genome: a locus, a strand and exon blocks. Coordinates are
/// 0-based and half-open, on the named sequence.
#[pyclass(name = "SplicedAlignment", module = "spearfish", frozen)]
pub struct PySplicedAlignment {
    /// Name of the reference sequence.
    #[pyo3(get)]
    sequence: String,
    /// Index of the reference sequence.
    #[pyo3(get)]
    sequence_id: u32,
    /// Start of the first exon block on the sequence.
    #[pyo3(get)]
    start: u32,
    /// End of the last exon block (exclusive).
    #[pyo3(get)]
    end: u32,
    /// Start of the aligned span on the read as sequenced.
    #[pyo3(get)]
    query_start: u32,
    /// End of the aligned span on the read (exclusive).
    #[pyo3(get)]
    query_end: u32,
    /// ``"+"`` if the read maps as sequenced, ``"-"`` as its reverse complement.
    #[pyo3(get)]
    strand: &'static str,
    /// The transcript's strand from the splice sites: ``"+"``, ``"-"``, or ``"."`` when
    /// there is no evidence.
    #[pyo3(get)]
    transcript_strand: &'static str,
    /// Mapping quality, 0 to 60.
    #[pyo3(get)]
    mapq: u8,
    /// Score; higher is better, comparable between the alignments of one read.
    #[pyo3(get)]
    score: i32,
    /// Whether this is the read's best alignment.
    #[pyo3(get)]
    is_primary: bool,
    /// Seed matches in the chain.
    #[pyo3(get)]
    anchors: u32,
    /// Edits in the base-level alignment, or ``None`` when none was made.
    #[pyo3(get)]
    edits: Option<u32>,
    /// The exon blocks as ``(start, end)`` on the sequence, in genome order.
    #[pyo3(get)]
    exons: Vec<(u32, u32)>,
    /// The introns as ``(start, end)`` on the sequence, in genome order.
    #[pyo3(get)]
    introns: Vec<(u32, u32)>,
    /// The blocks as ``(query_start, query_end, start, end)``; read coordinates are in the
    /// orientation the read maps in.
    #[pyo3(get)]
    blocks: Vec<(u32, u32, u32, u32)>,
    /// The CIGAR string (introns as ``N``): base by base when a base-level alignment was
    /// made, otherwise one ``M`` run per block.
    #[pyo3(get)]
    cigar: String,
}

#[pymethods]
impl PySplicedAlignment {
    /// Whether the alignment has an intron.
    #[getter]
    fn is_spliced(&self) -> bool {
        !self.introns.is_empty()
    }

    fn __repr__(&self) -> String {
        format!(
            "SplicedAlignment({}:{}-{} {}, exons={}, mapq={}, score={})",
            self.sequence,
            self.start,
            self.end,
            self.strand,
            self.exons.len(),
            self.mapq,
            self.score
        )
    }
}

/// The alignments of one read, best first; none when the read does not map.
#[pyclass(name = "SplicedRead", module = "spearfish", frozen)]
pub struct PySplicedRead {
    /// The read's name, if one was given.
    #[pyo3(get)]
    name: Option<String>,
    /// Read length in bases.
    #[pyo3(get)]
    length: u32,
    /// The :class:`SplicedAlignment` objects: the primary, then the secondary ones.
    #[pyo3(get)]
    alignments: Vec<Py<PySplicedAlignment>>,
    paf: Vec<u8>,
}

#[pymethods]
impl PySplicedRead {
    /// Whether the read has an alignment.
    #[getter]
    fn is_mapped(&self) -> bool {
        !self.alignments.is_empty()
    }

    /// The primary alignment, or ``None``.
    #[getter]
    fn best(&self, py: Python<'_>) -> Option<Py<PySplicedAlignment>> {
        self.alignments.first().map(|a| a.clone_ref(py))
    }

    /// The read's PAF lines, as ``spearfish map-spliced`` writes them (empty if unmapped).
    fn to_paf<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.paf)
    }

    fn __len__(&self) -> usize {
        self.alignments.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "SplicedRead(name={:?}, length={}, alignments={})",
            self.name.as_deref().unwrap_or(""),
            self.length,
            self.alignments.len()
        )
    }
}

fn to_py(
    py: Python<'_>,
    index: &GenomeIndex,
    name: Option<String>,
    r: &SplicedRead,
) -> PyResult<PySplicedRead> {
    let mut paf = Vec::new();
    output::paf_lines(&mut paf, name.as_deref().unwrap_or("*"), r, index);
    let alignments = r
        .alignments
        .iter()
        .map(|a| {
            let cigar: String = a
                .cigar(r.len)
                .into_iter()
                .map(|op| {
                    let (n, c) = op.len_and_letter();
                    format!("{n}{c}")
                })
                .collect();
            Py::new(
                py,
                PySplicedAlignment {
                    sequence: index.seq_name(a.seq).to_string(),
                    sequence_id: a.seq,
                    start: a.target_start,
                    end: a.target_end,
                    query_start: a.query_start,
                    query_end: a.query_end,
                    strand: if a.reverse { "-" } else { "+" },
                    transcript_strand: match a.transcript_strand {
                        TranscriptStrand::Forward => "+",
                        TranscriptStrand::Reverse => "-",
                        TranscriptStrand::Unknown => ".",
                    },
                    mapq: a.mapq,
                    score: a.score,
                    is_primary: a.primary,
                    anchors: a.anchors,
                    edits: a.edits,
                    exons: a.blocks().map(|b| (b.target_start, b.target_end)).collect(),
                    introns: a.introns().collect(),
                    blocks: a
                        .blocks()
                        .map(|b| (b.query_start, b.query_end, b.target_start, b.target_end))
                        .collect(),
                    cigar,
                },
            )
        })
        .collect::<PyResult<_>>()?;
    Ok(PySplicedRead {
        name,
        length: r.len,
        alignments,
        paf,
    })
}

/// Maps long RNA-seq reads, spliced, to a :class:`GenomeIndex`.
///
/// ``SplicedMapper(index, "drna", secondary=2)``. A mapper is immutable and thread-safe:
/// :meth:`map` releases the GIL, and :meth:`map_batch` and :meth:`map_file` map on the
/// mapper's own pool of ``threads`` threads.
#[pyclass(name = "SplicedMapper", module = "spearfish", frozen)]
pub struct PySplicedMapper {
    inner: Arc<SplicedMapper>,
    index: Py<PyGenomeIndex>,
    pool: Arc<rayon::ThreadPool>,
    threads: usize,
    preset: Preset,
}

#[pymethods]
impl PySplicedMapper {
    /// :param index: The :class:`GenomeIndex` to map against (with its annotation, if any).
    /// :param preset: ``"ont"`` (ONT cDNA, the default), ``"drna"`` (ONT direct RNA: reads
    ///     are the transcripts' sense strand) or ``"hifi"`` (PacBio).
    /// :param threads: Threads for :meth:`map_batch` and :meth:`map_file` (0: every
    ///     available CPU).
    /// :param settings: Settings changed from the preset, by name (:data:`SPLICE_SETTINGS`
    ///     lists them; :attr:`options` shows their values). The common ones: ``secondary``
    ///     (secondary alignments reported per read, at most; the same as ``max_hits`` minus
    ///     one), ``base_level`` (``"blocks"``, the default: exon blocks only; ``"edits"``:
    ///     each block aligned with the fewest edits; ``"affine"``: with affine gap costs),
    ///     ``max_intron`` and ``min_score``.
    /// :raises TypeError: For an unknown setting.
    /// :raises ValueError: For an unknown preset or an invalid value.
    #[new]
    #[pyo3(signature = (index, preset="ont", *, threads=0, **settings))]
    fn new(
        index: &Bound<'_, PyGenomeIndex>,
        preset: &str,
        threads: usize,
        settings: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Self> {
        let preset: Preset = preset.parse().map_err(value_err)?;
        // The preset's options as a map, the caller's settings over them, and back: every
        // field of the options is a setting, under its own name.
        let mut opts = serde_json::to_value(SpliceOptions::preset(preset)).map_err(value_err)?;
        let fields = opts.as_object_mut().expect("options are a map");
        for (k, v) in settings.into_iter().flat_map(|d| d.iter()) {
            let key: String = k.extract()?;
            if v.is_none() {
                continue;
            }
            let mut value: serde_json::Value = pythonize::depythonize(&v).map_err(value_err)?;
            let name = match key.as_str() {
                "secondary" => {
                    let n = value
                        .as_u64()
                        .ok_or_else(|| PyValueError::new_err("secondary is a count of alignments"))?;
                    value = (n + 1).into();
                    "max_hits"
                }
                k if fields.contains_key(k) => k,
                _ => {
                    return Err(pyo3::exceptions::PyTypeError::new_err(format!(
                        "unknown setting {key:?}; the settings are secondary, {}",
                        fields.keys().cloned().collect::<Vec<_>>().join(", ")
                    )));
                }
            };
            fields.insert(name.to_string(), value);
        }
        let opts: SpliceOptions = serde_json::from_value(opts)
            .map_err(|e| PyValueError::new_err(format!("invalid setting: {e}")))?;
        let n = crate::threads(threads);
        Ok(PySplicedMapper {
            inner: Arc::new(SplicedMapper::new(index.get().inner.clone(), &opts)),
            index: index.clone().unbind(),
            pool: Arc::new(pool(n)?),
            threads: n,
            preset,
        })
    }

    /// Map one read.
    ///
    /// :param seq: The read's bases (``str`` or ``bytes``; ``U`` reads as ``T``).
    /// :param name: A name to carry on the result (and in its PAF lines).
    /// :returns: A :class:`SplicedRead`.
    #[pyo3(signature = (seq, name=None))]
    fn map(&self, py: Python<'_>, seq: &Bound<'_, PyAny>, name: Option<String>) -> PyResult<PySplicedRead> {
        let b = seq_bytes(seq)?;
        let s = b.as_bytes();
        let mapper = &self.inner;
        let r = py.detach(|| mapper.map(s));
        to_py(py, mapper.index(), name, &r)
    }

    /// Map a batch of reads in parallel, results in input order.
    ///
    /// :param seqs: The reads' bases (an iterable of ``str`` or ``bytes``).
    /// :param names: Their names, as many as there are reads (optional).
    /// :returns: A list of :class:`SplicedRead`.
    #[pyo3(signature = (seqs, names=None))]
    fn map_batch(
        &self,
        py: Python<'_>,
        seqs: &Bound<'_, PyAny>,
        names: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Vec<PySplicedRead>> {
        let owned: Vec<Bound<'_, PyBytes>> = seqs
            .try_iter()?
            .map(|s| seq_bytes(&s?))
            .collect::<PyResult<_>>()?;
        let mut names: Vec<Option<String>> = match names {
            Some(n) if !n.is_none() => n
                .try_iter()?
                .map(|x| x?.extract::<String>().map(Some))
                .collect::<PyResult<_>>()?,
            _ => vec![None; owned.len()],
        };
        if names.len() != owned.len() {
            return Err(PyValueError::new_err(format!(
                "{} names for {} reads",
                names.len(),
                owned.len()
            )));
        }
        let slices: Vec<&[u8]> = owned.iter().map(|b| b.as_bytes()).collect();
        let (mapper, pool) = (&self.inner, &self.pool);
        let mut out = Vec::with_capacity(slices.len());
        // In parts, so Ctrl-C (and other signals) can interrupt a large batch.
        for part in slices.chunks(SIGNAL_CHECK_READS) {
            let reads: Vec<SplicedRead> =
                py.detach(|| pool.install(|| part.par_iter().map(|s| mapper.map(s)).collect()));
            py.check_signals()?;
            for r in &reads {
                let name = names[out.len()].take();
                out.push(to_py(py, mapper.index(), name, r)?);
            }
        }
        Ok(out)
    }

    /// Map a FASTA/FASTQ file (optionally gzip/bzip2/zstd-compressed) straight to alignment
    /// files, read, mapped and written in Rust, as ``spearfish map-spliced`` does. Reads
    /// are written as they finish, not in input order.
    ///
    /// :param reads: The reads.
    /// :param paf: PAF file to write (mapped reads only).
    /// :param sam: SAM file to write.
    /// :param bam: BAM file to write (unsorted).
    /// :param bam_level: BAM compression level, 0 to 9.
    /// :param unmapped: Write a record for each read that does not map (SAM and BAM).
    /// :param qualities: Write base qualities when the reads have them (SAM and BAM).
    /// :returns: A ``dict`` of counts: ``reads``, ``mapped``, ``unmapped``, ``spliced``,
    ///     ``junctions``, ``multi_mapped``, ``forward``, ``alignments``, ``bases``.
    /// :raises ValueError: If no output is named.
    #[pyo3(signature = (reads, *, paf=None, sam=None, bam=None, bam_level=1, unmapped=true, qualities=true))]
    #[allow(clippy::too_many_arguments)]
    fn map_file<'py>(
        &self,
        py: Python<'py>,
        reads: PathBuf,
        paf: Option<PathBuf>,
        sam: Option<PathBuf>,
        bam: Option<PathBuf>,
        bam_level: u8,
        unmapped: bool,
        qualities: bool,
    ) -> PyResult<Bound<'py, PyAny>> {
        if paf.is_none() && sam.is_none() && bam.is_none() {
            return Err(PyValueError::new_err("name at least one of paf, sam and bam"));
        }
        if bam_level > 9 {
            return Err(PyValueError::new_err("bam_level is 0 to 9"));
        }
        let create = |p: &Option<PathBuf>| -> PyResult<Option<Box<dyn std::io::Write + Send>>> {
            match p {
                Some(p) => {
                    let f = std::fs::File::create(p)
                        .map_err(|e| PyIOError::new_err(format!("creating {}: {e}", p.display())))?;
                    Ok(Some(Box::new(f)))
                }
                None => Ok(None),
            }
        };
        let outputs = FileOutputs {
            paf: create(&paf)?,
            sam: create(&sam)?,
            bam: create(&bam)?,
            bam_level,
            unmapped,
            qualities,
            command_line: String::new(),
        };
        let (mapper, threads) = (&self.inner, self.threads);
        let (counts, _) = py
            .detach(|| spearfish::genome::map_file(mapper, Some(&reads), outputs, threads, 32, None))
            .map_err(|e| io_err(&format!("mapping {}", reads.display()), e))?;
        pythonize::pythonize(py, &counts).map_err(value_err)
    }

    /// The :class:`GenomeIndex` mapped against.
    #[getter]
    fn index(&self, py: Python<'_>) -> Py<PyGenomeIndex> {
        self.index.clone_ref(py)
    }

    /// The preset's name.
    #[getter]
    fn preset(&self) -> &'static str {
        self.preset.name()
    }

    /// Threads :meth:`map_batch` and :meth:`map_file` use.
    #[getter]
    fn threads(&self) -> usize {
        self.threads
    }

    /// Every parameter the mapper runs with (a ``dict``).
    #[getter]
    fn options<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        pythonize::pythonize(py, self.inner.options()).map_err(value_err)
    }

    fn __repr__(&self) -> String {
        format!(
            "SplicedMapper(preset='{}', threads={})",
            self.preset.name(),
            self.threads
        )
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let opts = serde_json::to_value(SpliceOptions::default()).map_err(value_err)?;
    let mut names: Vec<String> = vec!["secondary".into()];
    names.extend(opts.as_object().expect("options are a map").keys().cloned());
    m.add("SPLICE_SETTINGS", pyo3::types::PyTuple::new(m.py(), names)?)?;
    m.add_class::<PyGenomeIndex>()?;
    m.add_class::<PySplicedMapper>()?;
    m.add_class::<PySplicedRead>()?;
    m.add_class::<PySplicedAlignment>()?;
    Ok(())
}
