//! RAD records for bulk long-read mappings, shared by spearfish (writer) and oarfish (reader).
//!
//! A file holds one record per mapped read with all of the read's alignments, in chunks that
//! may each be compressed with libradicl's per-chunk codec (LZ4 or zstd). The schema does not
//! depend on the mapper: any long-read mapper can write it, and file-level tags carry what a
//! quantifier needs to interpret the scores and positions ([`Meta`]).
//!
//! Schema (`known_rad_type = "bulk_long_read"`, `schema_version = 1`):
//!
//! | level     | tag                          | type | meaning                                      |
//! |-----------|------------------------------|------|----------------------------------------------|
//! | read      | `read_len`                   | u32  | read length in bases                         |
//! | read      | `read_flags`                 | u8   | [`read_flags`] bits                          |
//! | alignment | `ref_ori`                    | u32  | target id in the low 31 bits; bit 31 set if the read aligns forward |
//! | alignment | `ref_start`, `ref_end`       | u32  | aligned span on the target, 0-based, end exclusive |
//! | alignment | `query_start`, `query_end`   | u32  | aligned span on the read's forward strand    |
//! | alignment | `score`                      | f32  | alignment score, higher is better            |
//! | alignment | `aln_flags`                  | u8   | [`aln_flags`] bits                           |
//!
//! A reader finds tags by name and skips the ones it does not know, so a producer can add
//! tags. `ref_ori`, `ref_start`, `ref_end` and `score` are required; the others read as 0 when
//! absent. Reads without alignments are not written; [`Counts`] records how many there were.

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context as _, Result, bail};
use bio_types::strand::Strand;
use libradicl::header::{RadHeader, RadPrelude, SpecVersion};
use libradicl::rad_types::{
    RadFloatId, RadIntId, RadType, TagDesc, TagMap, TagRole, TagSection, TagSectionLabel, TagValue,
};
use libradicl::readers::ParallelRadReader;
use libradicl::record::{MappedRecord, RecordContext};
use libradicl::{CHUNK_CODEC_TAG, ChunkBuf, ConcurrentChunkWriter, RadFileWriter};

pub use libradicl;
pub use libradicl::ChunkCodec;

pub const KNOWN_RAD_TYPE: &str = "bulk_long_read";
pub const SCHEMA_VERSION: u32 = 1;

/// Bits of the `read_flags` tag.
pub mod read_flags {
    /// The read's best alignment covers little of it.
    pub const LOW_COVERAGE: u8 = 1;
    /// The read was mapped by a rescue pass.
    pub const RESCUED: u8 = 2;
}

/// Bits of the `aln_flags` tag.
pub mod aln_flags {
    /// The span runs between the outermost anchors (seeds) only and was not extended to the
    /// read's ends, so its endpoints are not reliable read-end positions.
    pub const ANCHOR_SPAN: u8 = 1;
}

/// One alignment of a read.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Aln {
    pub tid: u32,
    pub fw: bool,
    pub ref_start: u32,
    pub ref_end: u32,
    pub query_start: u32,
    pub query_end: u32,
    pub score: f32,
    pub flags: u8,
}

/// A read and its alignments, stored by field.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadRecord {
    pub read_len: u32,
    pub flags: u8,
    pub refs: Vec<u32>,
    pub fw: Vec<bool>,
    pub ref_start: Vec<u32>,
    pub ref_end: Vec<u32>,
    pub query_start: Vec<u32>,
    pub query_end: Vec<u32>,
    pub score: Vec<f32>,
    pub aln_flags: Vec<u8>,
}

impl ReadRecord {
    /// Empty the record, keeping its allocations.
    pub fn clear(&mut self) {
        self.read_len = 0;
        self.flags = 0;
        self.refs.clear();
        self.fw.clear();
        self.ref_start.clear();
        self.ref_end.clear();
        self.query_start.clear();
        self.query_end.clear();
        self.score.clear();
        self.aln_flags.clear();
    }

    pub fn push(&mut self, a: Aln) {
        self.refs.push(a.tid);
        self.fw.push(a.fw);
        self.ref_start.push(a.ref_start);
        self.ref_end.push(a.ref_end);
        self.query_start.push(a.query_start);
        self.query_end.push(a.query_end);
        self.score.push(a.score);
        self.aln_flags.push(a.flags);
    }

    pub fn len(&self) -> usize {
        self.refs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.refs.is_empty()
    }

    pub fn aln(&self, i: usize) -> Aln {
        Aln {
            tid: self.refs[i],
            fw: self.fw[i],
            ref_start: self.ref_start[i],
            ref_end: self.ref_end[i],
            query_start: self.query_start[i],
            query_end: self.query_end[i],
            score: self.score[i],
            flags: self.aln_flags[i],
        }
    }

    pub fn alns(&self) -> impl ExactSizeIterator<Item = Aln> + '_ {
        (0..self.len()).map(|i| self.aln(i))
    }
}

/// Which record field a tag holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Field {
    ReadLen,
    ReadFlags,
    RefOri,
    RefStart,
    RefEnd,
    QueryStart,
    QueryEnd,
    Score,
    AlnFlags,
    /// A tag this crate does not know: skipped when reading, written as zero or empty.
    Other,
}

const ORI_FW: u32 = 1 << 31;

/// The read and alignment tag layouts of a file.
#[derive(Clone, Debug)]
pub struct Context {
    read: Vec<(Field, TagDesc)>,
    aln: Vec<(Field, TagDesc)>,
}

fn u32_tag(name: &str) -> TagDesc {
    TagDesc::new(name, RadType::Int(RadIntId::U32))
}

fn u8_tag(name: &str) -> TagDesc {
    TagDesc::new(name, RadType::Int(RadIntId::U8))
}

impl Context {
    /// The schema this crate writes.
    pub fn canonical() -> Self {
        Context {
            read: vec![
                (Field::ReadLen, u32_tag("read_len")),
                (Field::ReadFlags, u8_tag("read_flags")),
            ],
            aln: vec![
                (Field::RefOri, u32_tag("ref_ori").with_role(TagRole::Orientation)),
                (
                    Field::RefStart,
                    u32_tag("ref_start").with_role(TagRole::MappingPosition),
                ),
                (Field::RefEnd, u32_tag("ref_end")),
                (Field::QueryStart, u32_tag("query_start")),
                (Field::QueryEnd, u32_tag("query_end")),
                (
                    Field::Score,
                    TagDesc::new("score", RadType::Float(RadFloatId::F32)),
                ),
                (Field::AlnFlags, u8_tag("aln_flags")),
            ],
        }
    }

    fn sections(&self) -> (TagSection, TagSection) {
        let mut rt = TagSection::new_with_label(TagSectionLabel::ReadTags);
        for (_, d) in &self.read {
            rt.add_tag_desc(d.clone());
        }
        let mut at = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
        for (_, d) in &self.aln {
            at.add_tag_desc(d.clone());
        }
        (rt, at)
    }
}

/// Integer tags may use any integer width; the score may be an integer or a float.
fn check_type(field: Field, d: &TagDesc) -> Result<()> {
    let ok = match field {
        Field::Other => true,
        Field::Score => matches!(d.typeid, RadType::Int(_) | RadType::Float(_)),
        _ => matches!(d.typeid, RadType::Int(_)),
    };
    if !ok {
        bail!("RAD tag {:?} has unsupported type {:?}", d.name, d.typeid);
    }
    Ok(())
}

impl RecordContext for Context {
    fn get_context_from_tag_section(_ft: &TagSection, rt: &TagSection, at: &TagSection) -> Result<Self> {
        let read = rt
            .iter_desc()
            .map(|d| {
                let f = match d.name.as_str() {
                    "read_len" => Field::ReadLen,
                    "read_flags" => Field::ReadFlags,
                    _ => Field::Other,
                };
                check_type(f, d).map(|_| (f, d.clone()))
            })
            .collect::<Result<Vec<_>>>()?;
        let aln = at
            .iter_desc()
            .map(|d| {
                let f = match d.name.as_str() {
                    "ref_ori" => Field::RefOri,
                    "ref_start" => Field::RefStart,
                    "ref_end" => Field::RefEnd,
                    "query_start" => Field::QueryStart,
                    "query_end" => Field::QueryEnd,
                    "score" => Field::Score,
                    "aln_flags" => Field::AlnFlags,
                    _ => Field::Other,
                };
                check_type(f, d).map(|_| (f, d.clone()))
            })
            .collect::<Result<Vec<_>>>()?;
        for need in [Field::RefOri, Field::RefStart, Field::RefEnd, Field::Score] {
            if !aln.iter().any(|(f, _)| *f == need) {
                bail!("RAD file lacks the alignment tag for {need:?} (not a {KNOWN_RAD_TYPE} file?)");
            }
        }
        Ok(Context { read, aln })
    }
}

fn read_exact<const N: usize, R: Read>(r: &mut R) -> [u8; N] {
    let mut b = [0u8; N];
    r.read_exact(&mut b).expect("truncated RAD record");
    b
}

/// An integer tag's value. Signed values are read as their two's-complement bits, which is
/// what the unsigned fields this crate reads should hold anyway.
fn read_int<R: Read>(id: RadIntId, r: &mut R) -> u64 {
    match id {
        RadIntId::U8 | RadIntId::I8 => read_exact::<1, _>(r)[0] as u64,
        RadIntId::U16 | RadIntId::I16 => u16::from_le_bytes(read_exact(r)) as u64,
        RadIntId::U32 | RadIntId::I32 => u32::from_le_bytes(read_exact(r)) as u64,
        RadIntId::U64 | RadIntId::I64 => u64::from_le_bytes(read_exact(r)),
        RadIntId::U128 | RadIntId::I128 => u128::from_le_bytes(read_exact(r)) as u64,
    }
}

fn read_signed<R: Read>(id: RadIntId, r: &mut R) -> i64 {
    match id {
        RadIntId::I8 => read_exact::<1, _>(r)[0] as i8 as i64,
        RadIntId::I16 => i16::from_le_bytes(read_exact(r)) as i64,
        RadIntId::I32 => i32::from_le_bytes(read_exact(r)) as i64,
        RadIntId::I64 => i64::from_le_bytes(read_exact(r)),
        RadIntId::I128 => i128::from_le_bytes(read_exact(r)) as i64,
        u => read_int(u, r) as i64,
    }
}

fn read_score<R: Read>(t: RadType, r: &mut R) -> f32 {
    match t {
        RadType::Float(RadFloatId::F32) => f32::from_le_bytes(read_exact(r)),
        RadType::Float(RadFloatId::F64) => f64::from_le_bytes(read_exact(r)) as f32,
        RadType::Int(id) => read_signed(id, r) as f32,
        _ => unreachable!("checked in get_context_from_tag_section"),
    }
}

fn write_int<W: Write>(id: RadIntId, v: u64, w: &mut W) -> std::io::Result<()> {
    match id {
        RadIntId::U8 | RadIntId::I8 => w.write_all(&[v as u8]),
        RadIntId::U16 | RadIntId::I16 => w.write_all(&(v as u16).to_le_bytes()),
        RadIntId::U32 | RadIntId::I32 => w.write_all(&(v as u32).to_le_bytes()),
        RadIntId::U64 | RadIntId::I64 => w.write_all(&v.to_le_bytes()),
        RadIntId::U128 | RadIntId::I128 => w.write_all(&(v as u128).to_le_bytes()),
    }
}

/// A zero (or empty) value of a tag's type, for tags the record does not hold.
fn write_zero<W: Write>(d: &TagDesc, w: &mut W) -> Result<()> {
    let zero = match d.typeid {
        RadType::Bool => TagValue::U8(0),
        RadType::Int(id) => return Ok(write_int(id, 0, w)?),
        RadType::Float(RadFloatId::F32) => TagValue::F32(0.0),
        RadType::Float(RadFloatId::F64) => TagValue::F64(0.0),
        RadType::String => TagValue::String(String::new()),
        RadType::Array(len, _) => return Ok(write_int(len, 0, w)?),
    };
    zero.write_with_type(&d.typeid, w)?;
    Ok(())
}

impl MappedRecord for ReadRecord {
    type ParsingContext = Context;
    /// The number of alignments.
    type PeekResult = u32;

    fn peek_record(buf: &[u8], _ctx: &Context) -> u32 {
        u32::from_le_bytes(buf[..4].try_into().expect("RAD record shorter than 4 bytes"))
    }

    fn from_bytes_with_context<T: Read>(r: &mut T, ctx: &Context) -> Self {
        let n = u32::from_le_bytes(read_exact(r)) as usize;
        let mut rec = ReadRecord {
            refs: Vec::with_capacity(n),
            fw: Vec::with_capacity(n),
            ref_start: Vec::with_capacity(n),
            ref_end: Vec::with_capacity(n),
            query_start: Vec::with_capacity(n),
            query_end: Vec::with_capacity(n),
            score: Vec::with_capacity(n),
            aln_flags: Vec::with_capacity(n),
            ..Default::default()
        };
        for (f, d) in &ctx.read {
            match (f, d.typeid) {
                (Field::ReadLen, RadType::Int(id)) => rec.read_len = read_int(id, r) as u32,
                (Field::ReadFlags, RadType::Int(id)) => rec.flags = read_int(id, r) as u8,
                _ => {
                    d.value_from_bytes(r);
                }
            }
        }
        for _ in 0..n {
            let mut a = Aln::default();
            for (f, d) in &ctx.aln {
                let int = |r: &mut T| match d.typeid {
                    RadType::Int(id) => read_int(id, r),
                    _ => unreachable!("checked in get_context_from_tag_section"),
                };
                match f {
                    Field::RefOri => {
                        let v = int(r) as u32;
                        a.tid = v & !ORI_FW;
                        a.fw = v & ORI_FW != 0;
                    }
                    Field::RefStart => a.ref_start = int(r) as u32,
                    Field::RefEnd => a.ref_end = int(r) as u32,
                    Field::QueryStart => a.query_start = int(r) as u32,
                    Field::QueryEnd => a.query_end = int(r) as u32,
                    Field::Score => a.score = read_score(d.typeid, r),
                    Field::AlnFlags => a.flags = int(r) as u8,
                    _ => {
                        d.value_from_bytes(r);
                    }
                }
            }
            rec.push(a);
        }
        rec
    }

    fn write<W: Write>(&self, w: &mut W, ctx: &Context) -> Result<()> {
        w.write_all(&(self.len() as u32).to_le_bytes())?;
        for (f, d) in &ctx.read {
            match (f, d.typeid) {
                (Field::ReadLen, RadType::Int(id)) => write_int(id, self.read_len as u64, w)?,
                (Field::ReadFlags, RadType::Int(id)) => write_int(id, self.flags as u64, w)?,
                _ => write_zero(d, w)?,
            }
        }
        for a in self.alns() {
            for (f, d) in &ctx.aln {
                let v = match f {
                    Field::RefOri => (a.tid & !ORI_FW | if a.fw { ORI_FW } else { 0 }) as u64,
                    Field::RefStart => a.ref_start as u64,
                    Field::RefEnd => a.ref_end as u64,
                    Field::QueryStart => a.query_start as u64,
                    Field::QueryEnd => a.query_end as u64,
                    Field::AlnFlags => a.flags as u64,
                    Field::Score => {
                        match d.typeid {
                            RadType::Float(RadFloatId::F32) => w.write_all(&a.score.to_le_bytes())?,
                            RadType::Float(RadFloatId::F64) => {
                                w.write_all(&(a.score as f64).to_le_bytes())?
                            }
                            RadType::Int(id) => write_int(id, a.score.round() as i64 as u64, w)?,
                            _ => unreachable!("checked in get_context_from_tag_section"),
                        }
                        continue;
                    }
                    _ => {
                        write_zero(d, w)?;
                        continue;
                    }
                };
                match d.typeid {
                    RadType::Int(id) => write_int(id, v, w)?,
                    _ => unreachable!("checked in get_context_from_tag_section"),
                }
            }
        }
        Ok(())
    }

    fn is_empty(&self) -> bool {
        self.refs.is_empty()
    }

    fn num_aln(&self) -> usize {
        self.refs.len()
    }

    fn refs(&self) -> &[u32] {
        &self.refs
    }

    fn has_alignment_on_strand(&self, s: Strand) -> bool {
        match s {
            Strand::Forward => self.fw.iter().any(|&f| f),
            Strand::Reverse => self.fw.iter().any(|&f| !f),
            Strand::Unknown => true,
        }
    }
}

/// What a quantifier needs to know about how the reads were mapped, from file-level tags.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Meta {
    /// The mapper's name, version, and command line.
    pub mapper: String,
    pub mapper_version: String,
    pub cmdline: String,
    /// The mapper's settings, free form (spearfish writes JSON).
    pub params: String,
    /// The mapper's preset, if any.
    pub preset: String,
    /// The sequencing technology as oarfish's `--seq-tech` names it (`ont-cdna`, `ont-drna`,
    /// `pac-bio`, `pac-bio-hifi`), or empty if unknown.
    pub seq_tech: String,
    /// The mapper's suggested scale for score differences: an alignment scoring `d` below
    /// the read's best is `exp(-d / score_denom)` times as likely. 0 if none is suggested.
    pub score_denom: f32,
    /// Whether the spans of the best alignments were extended to the read's ends, so that
    /// their endpoints are read-end positions (as a truncation model needs). Alignments whose
    /// spans were not extended carry [`aln_flags::ANCHOR_SPAN`].
    pub ends_extended: bool,
    /// Which strands of the targets reads were mapped to: `both`, `fwd` or `rev`.
    pub strands: String,
}

/// Read counts, patched into the file when it is finished.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    /// Reads processed, mapped or not.
    pub reads: u64,
    /// Reads with at least one alignment (the records in the file).
    pub mapped: u64,
    /// Bases in all processed reads.
    pub bases: u64,
}

fn file_tags(meta: &Meta, ref_lengths: &[u32], codec: ChunkCodec) -> (TagSection, TagMap) {
    let s = |x: &str| TagValue::String(x.to_string());
    TagSection::from_tag_values(
        TagSectionLabel::FileTags,
        &[
            ("known_rad_type", s(KNOWN_RAD_TYPE)),
            ("schema_version", TagValue::U32(SCHEMA_VERSION)),
            (CHUNK_CODEC_TAG, TagValue::U8(codec.as_u8())),
            ("ref_lengths", TagValue::ArrayU32(ref_lengths.to_vec())),
            ("mapper", s(&meta.mapper)),
            ("mapper_version", s(&meta.mapper_version)),
            ("cmdline", s(&meta.cmdline)),
            ("params", s(&meta.params)),
            ("preset", s(&meta.preset)),
            ("seq_tech", s(&meta.seq_tech)),
            ("score_denom", TagValue::F32(meta.score_denom)),
            ("ends_extended", TagValue::U8(meta.ends_extended as u8)),
            ("strands", s(&meta.strands)),
            ("num_reads", TagValue::U64(0)),
            ("num_mapped", TagValue::U64(0)),
            ("num_bases", TagValue::U64(0)),
        ],
    )
}

/// Records are gathered into chunks of about this many bytes (before compression).
const CHUNK_BYTES: usize = 128 * 1024;

type Out = ConcurrentChunkWriter<BufWriter<File>>;

/// Writes a RAD file. The file is written under a temporary name and moved into place by
/// [`RadWriter::finish`], so an interrupted run leaves no file that looks complete.
pub struct RadWriter {
    out: Arc<Out>,
    ctx: Arc<Context>,
    codec: ChunkCodec,
    path: PathBuf,
    tmp: PathBuf,
}

impl RadWriter {
    /// Start a file over the targets `refs` (names and lengths, in target-id order).
    pub fn create<S: AsRef<str>>(
        path: &Path,
        refs: &[(S, u32)],
        meta: &Meta,
        codec: ChunkCodec,
    ) -> Result<Self> {
        let mut hdr = RadHeader::from_ref_names(refs.iter().map(|(n, _)| n.as_ref()));
        hdr.version = SpecVersion::current();
        let lens: Vec<u32> = refs.iter().map(|(_, l)| *l).collect();
        let (ft, values) = file_tags(meta, &lens, codec);
        let ctx = Context::canonical();
        let (rt, at) = ctx.sections();
        let prelude = RadPrelude::from_header_and_tag_sections(hdr, ft, rt, at);
        let mut tmp = path.as_os_str().to_owned();
        tmp.push(format!(".partial-{}", std::process::id()));
        let tmp = PathBuf::from(tmp);
        let f = File::create(&tmp).with_context(|| format!("creating {}", tmp.display()))?;
        let w = RadFileWriter::new(BufWriter::new(f), &prelude, &values)?;
        Ok(RadWriter {
            out: Arc::new(ConcurrentChunkWriter::new(w)),
            ctx: Arc::new(ctx),
            codec,
            path: path.to_path_buf(),
            tmp,
        })
    }

    /// A per-thread record buffer that appends whole chunks to the file.
    pub fn sink(&self) -> Sink {
        Sink {
            out: self.out.clone(),
            ctx: self.ctx.clone(),
            codec: self.codec,
            buf: ChunkBuf::with_capacity(CHUNK_BYTES + CHUNK_BYTES / 4),
        }
    }

    /// Patch in the counts and move the file into place. Every [`Sink`] must have been
    /// flushed and dropped.
    pub fn finish(self, counts: Counts) -> Result<()> {
        let out = Arc::try_unwrap(self.out)
            .map_err(|_| anyhow::anyhow!("RAD writer finished while a sink is still alive"))?;
        let mut w = out.into_writer()?;
        w.backpatch_file_tag_value("num_reads", &TagValue::U64(counts.reads))?;
        w.backpatch_file_tag_value("num_mapped", &TagValue::U64(counts.mapped))?;
        w.backpatch_file_tag_value("num_bases", &TagValue::U64(counts.bases))?;
        let f = w.finalize()?.into_inner().map_err(|e| e.into_error())?;
        f.sync_all()?;
        std::fs::rename(&self.tmp, &self.path)
            .with_context(|| format!("moving {} into place", self.path.display()))?;
        Ok(())
    }
}

/// Gathers one thread's records into chunks. Call [`Sink::flush`] when the thread is done;
/// records still buffered when a sink is dropped are lost.
pub struct Sink {
    out: Arc<Out>,
    ctx: Arc<Context>,
    codec: ChunkCodec,
    buf: ChunkBuf,
}

impl Clone for Sink {
    /// A sink to the same file with an empty buffer.
    fn clone(&self) -> Self {
        Sink {
            out: self.out.clone(),
            ctx: self.ctx.clone(),
            codec: self.codec,
            buf: ChunkBuf::with_capacity(CHUNK_BYTES + CHUNK_BYTES / 4),
        }
    }
}

impl Sink {
    pub fn push(&mut self, rec: &ReadRecord) -> Result<()> {
        self.buf.write_record(rec, &self.ctx)?;
        if self.buf.byte_len() >= CHUNK_BYTES {
            self.flush()?;
        }
        Ok(())
    }

    /// Write out the buffered records as a chunk.
    pub fn flush(&mut self) -> Result<()> {
        if self.buf.nrec() == 0 {
            return Ok(());
        }
        let full = std::mem::replace(
            &mut self.buf,
            ChunkBuf::with_capacity(CHUNK_BYTES + CHUNK_BYTES / 4),
        );
        let bytes = full.into_bytes_with_codec(self.codec)?;
        self.out.append_chunk_bytes(&bytes)
    }
}

/// The start of a versioned (spec 2) RAD file.
const MAGIC: &[u8; 8] = b"RAD_FILE";

/// Whether `path` starts like a versioned RAD file.
pub fn is_rad(path: &Path) -> Result<bool> {
    let mut b = [0u8; 8];
    let mut f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    Ok(f.read_exact(&mut b).is_ok() && &b == MAGIC)
}

/// A file's header: target names and lengths, mapping metadata, and counts.
#[derive(Clone, Debug)]
pub struct Header {
    pub ref_names: Vec<String>,
    pub ref_lengths: Vec<u32>,
    pub meta: Meta,
    pub counts: Counts,
    pub codec: ChunkCodec,
}

impl Header {
    fn from_parts(prelude: &RadPrelude, tags: &TagMap) -> Result<Self> {
        let string = |k: &str| match tags.get(k) {
            Some(TagValue::String(s)) => s.clone(),
            _ => String::new(),
        };
        let int = |k: &str| match tags.get(k) {
            Some(TagValue::U8(v)) => *v as u64,
            Some(TagValue::U32(v)) => *v as u64,
            Some(TagValue::U64(v)) => *v,
            _ => 0,
        };
        let kind = string("known_rad_type");
        if kind != KNOWN_RAD_TYPE {
            bail!("RAD file is of type {kind:?}, not {KNOWN_RAD_TYPE:?}");
        }
        let ref_lengths = match tags.get("ref_lengths") {
            Some(TagValue::ArrayU32(v)) => v.clone(),
            Some(TagValue::ArrayU64(v)) => v.iter().map(|&x| x as u32).collect(),
            _ => bail!("RAD file has no ref_lengths file tag"),
        };
        let ref_names = prelude.hdr.ref_names.clone();
        if ref_names.len() != ref_lengths.len() {
            bail!(
                "RAD file names {} targets but gives {} lengths",
                ref_names.len(),
                ref_lengths.len()
            );
        }
        Ok(Header {
            ref_names,
            ref_lengths,
            meta: Meta {
                mapper: string("mapper"),
                mapper_version: string("mapper_version"),
                cmdline: string("cmdline"),
                params: string("params"),
                preset: string("preset"),
                seq_tech: string("seq_tech"),
                score_denom: match tags.get("score_denom") {
                    Some(TagValue::F32(v)) => *v,
                    Some(TagValue::F64(v)) => *v as f32,
                    _ => 0.0,
                },
                ends_extended: int("ends_extended") != 0,
                strands: string("strands"),
            },
            counts: Counts {
                reads: int("num_reads"),
                mapped: int("num_mapped"),
                bases: int("num_bases"),
            },
            codec: libradicl::codec::chunk_codec_from_tag_map(tags)?,
        })
    }
}

/// Read a file's header without its records.
pub fn read_header(path: &Path) -> Result<Header> {
    let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut r = BufReader::new(f);
    let prelude = RadPrelude::from_bytes(&mut r).context("reading the RAD prelude")?;
    let tags = prelude
        .file_tags
        .try_parse_tags_from_bytes(&mut r)
        .context("reading the RAD file tags")?;
    Header::from_parts(&prelude, &tags)
}

/// A parallel reader over a file's records, and its header. Records arrive in chunks, in
/// file order, decompressed; see [`ParallelRadReader::process_parallel`].
pub fn open(
    path: &Path,
    nworkers: std::num::NonZeroUsize,
) -> Result<(Header, ParallelRadReader<ReadRecord, BufReader<File>>)> {
    let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let r = ParallelRadReader::<ReadRecord, _>::try_new(BufReader::new(f), nworkers)?;
    let hdr = Header::from_parts(&r.prelude, &r.file_tag_map)?;
    Ok((hdr, r))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    fn rec(i: u32) -> ReadRecord {
        let mut r = ReadRecord {
            read_len: 1000 + i,
            flags: (i % 4) as u8,
            ..Default::default()
        };
        for k in 0..(i % 5 + 1) {
            r.push(Aln {
                tid: (i + k) % 3,
                fw: (i + k).is_multiple_of(2),
                ref_start: i,
                ref_end: i + 900,
                query_start: k,
                query_end: 950 + k,
                score: 1234.5 - k as f32 * 3.25,
                flags: (k % 2) as u8,
            });
        }
        r
    }

    fn roundtrip(codec: ChunkCodec) {
        let dir = std::env::temp_dir().join(format!("lrrad-{}-{}", std::process::id(), codec.as_u8()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.rad");
        let refs = [("tA", 2000u32), ("tB", 1500), ("tC", 3000)];
        let meta = Meta {
            mapper: "test".into(),
            mapper_version: "1".into(),
            preset: "ont".into(),
            seq_tech: "ont-cdna".into(),
            score_denom: 2.0,
            ends_extended: true,
            strands: "both".into(),
            ..Default::default()
        };
        let w = RadWriter::create(&path, &refs, &meta, codec).unwrap();
        let n = 20_000u32;
        std::thread::scope(|s| {
            for t in 0..4 {
                let mut sink = w.sink();
                s.spawn(move || {
                    for i in (t..n).step_by(4) {
                        sink.push(&rec(i)).unwrap();
                    }
                    sink.flush().unwrap();
                });
            }
        });
        let counts = Counts {
            reads: n as u64 + 7,
            mapped: n as u64,
            bases: 42,
        };
        w.finish(counts).unwrap();
        assert!(is_rad(&path).unwrap());
        let h = read_header(&path).unwrap();
        assert_eq!(h.meta, meta);
        assert_eq!(h.counts, counts);
        assert_eq!(h.codec, codec);
        assert_eq!(h.ref_names, vec!["tA", "tB", "tC"]);
        assert_eq!(h.ref_lengths, vec![2000, 1500, 3000]);
        let (_, mut r) = open(&path, std::num::NonZeroUsize::new(3).unwrap()).unwrap();
        let got = Mutex::new(Vec::new());
        r.process_parallel(std::num::NonZeroUsize::new(3).unwrap(), |mc| {
            let mut v = Vec::new();
            for c in mc.iter() {
                v.extend(c.reads);
            }
            got.lock().unwrap().extend(v);
        })
        .unwrap();
        let mut got = got.into_inner().unwrap();
        got.sort_by_key(|r| r.read_len);
        assert_eq!(got.len(), n as usize);
        for (i, r) in got.iter().enumerate() {
            assert_eq!(*r, rec(i as u32));
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn roundtrip_uncompressed() {
        roundtrip(ChunkCodec::None);
    }

    #[test]
    fn roundtrip_lz4() {
        roundtrip(ChunkCodec::Lz4);
    }

    #[test]
    fn roundtrip_zstd() {
        roundtrip(ChunkCodec::Zstd);
    }

    /// A file with an extra tag at each level, and the score as an integer, still reads.
    #[test]
    fn unknown_tags_are_skipped() {
        let mut ctx = Context::canonical();
        ctx.read
            .insert(0, (Field::Other, TagDesc::new("name", RadType::String)));
        ctx.aln.push((Field::Other, u32_tag("nm")));
        ctx.aln.iter_mut().for_each(|(f, d)| {
            if *f == Field::Score {
                *d = TagDesc::new("score", RadType::Int(RadIntId::I32));
            }
        });
        let (rt, at) = ctx.sections();
        let ft = TagSection::new_with_label(TagSectionLabel::FileTags);
        let read = Context::get_context_from_tag_section(&ft, &rt, &at).unwrap();
        let mut r = rec(3);
        r.score.iter_mut().for_each(|s| *s = s.round());
        let mut bytes = Vec::new();
        r.write(&mut bytes, &ctx).unwrap();
        let back = ReadRecord::from_bytes_with_context(&mut bytes.as_slice(), &read);
        assert_eq!(back, r);
    }
}
