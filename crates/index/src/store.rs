//! A flat container of named, typed arrays: written once, memory-mapped and read zero-copy.
//!
//! Layout (little-endian): magic `SPFISH01`, `u32` format version, `u32` section count, a
//! table of `(name: [u8; 24], offset: u64, len_bytes: u64)` entries, then each section's
//! bytes at a 64-byte-aligned offset.

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail, ensure};
use bytemuck::Pod;
use memmap2::Mmap;

const MAGIC: &[u8; 8] = b"SPFISH01";
pub const FORMAT_VERSION: u32 = 3;
const NAME_LEN: usize = 24;
const ALIGN: usize = 64;

/// Sections to write, in order. Sections borrow their data where they can, so writing an
/// index to disk streams it without another copy.
#[derive(Default)]
pub struct Writer<'a> {
    sections: Vec<(String, Section<'a>)>,
}

/// A section's bytes: borrowed, or a vector kept as it was given (no copy to bytes).
enum Section<'a> {
    Borrowed(&'a [u8]),
    Owned(Box<dyn Bytes>),
}

trait Bytes {
    fn bytes(&self) -> &[u8];
}

impl<T: Pod> Bytes for Vec<T> {
    fn bytes(&self) -> &[u8] {
        bytemuck::cast_slice(self)
    }
}

impl std::ops::Deref for Section<'_> {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        match self {
            Section::Borrowed(b) => b,
            Section::Owned(v) => v.bytes(),
        }
    }
}

impl<'a> Writer<'a> {
    /// Add a section borrowing `data`.
    pub fn add<T: Pod>(&mut self, name: &str, data: &'a [T]) {
        assert!(name.len() <= NAME_LEN, "section name too long: {name}");
        self.sections
            .push((name.to_string(), Section::Borrowed(bytemuck::cast_slice(data))));
    }

    /// Add a section owning its data.
    pub fn add_owned<T: Pod>(&mut self, name: &str, data: Vec<T>) {
        assert!(name.len() <= NAME_LEN, "section name too long: {name}");
        self.sections
            .push((name.to_string(), Section::Owned(Box::new(data))));
    }

    /// The header (magic, version, section table) and each section's offset; the file is
    /// `total` bytes.
    fn layout(&self) -> (Vec<u8>, Vec<usize>, usize) {
        let header = 16 + self.sections.len() * (NAME_LEN + 16);
        let mut offsets = Vec::with_capacity(self.sections.len());
        let mut pos = header.next_multiple_of(ALIGN);
        for (_, d) in &self.sections {
            offsets.push(pos);
            pos = (pos + d.len()).next_multiple_of(ALIGN);
        }
        let mut h = vec![0u8; header];
        h[..8].copy_from_slice(MAGIC);
        h[8..12].copy_from_slice(&FORMAT_VERSION.to_le_bytes());
        h[12..16].copy_from_slice(&(self.sections.len() as u32).to_le_bytes());
        for (i, ((name, d), &off)) in self.sections.iter().zip(&offsets).enumerate() {
            let e = 16 + i * (NAME_LEN + 16);
            h[e..e + name.len()].copy_from_slice(name.as_bytes());
            h[e + NAME_LEN..e + NAME_LEN + 8].copy_from_slice(&(off as u64).to_le_bytes());
            h[e + NAME_LEN + 8..e + NAME_LEN + 16].copy_from_slice(&(d.len() as u64).to_le_bytes());
        }
        (h, offsets, pos)
    }

    /// Serialize to bytes (as written to disk).
    pub fn to_bytes(&self) -> Vec<u8> {
        let (h, offsets, total) = self.layout();
        let mut out = vec![0u8; total];
        self.fill(&mut out, &h, &offsets);
        out
    }

    /// [`Writer::to_bytes`] into 8-byte-aligned words, with the length in bytes: what
    /// [`Store::from_words`] takes without another copy.
    pub(crate) fn to_words(&self) -> (Vec<u64>, usize) {
        let (h, offsets, total) = self.layout();
        let mut words = vec![0u64; total.div_ceil(8)];
        self.fill(&mut bytemuck::cast_slice_mut(&mut words)[..total], &h, &offsets);
        (words, total)
    }

    /// The container in memory, as [`Store::open`] would load it from the written file. Each
    /// owned section is freed as soon as it is copied, so the peak is one copy of the data
    /// plus the largest section.
    pub fn into_store(self) -> Result<Store> {
        let (h, offsets, total) = self.layout();
        let mut words = vec![0u64; total.div_ceil(8)];
        let out = &mut bytemuck::cast_slice_mut::<u64, u8>(&mut words)[..total];
        out[..h.len()].copy_from_slice(&h);
        for ((_, d), off) in self.sections.into_iter().zip(offsets) {
            out[off..off + d.len()].copy_from_slice(&d);
        }
        Store::from_words(words, total)
    }

    fn fill(&self, out: &mut [u8], h: &[u8], offsets: &[usize]) {
        out[..h.len()].copy_from_slice(h);
        for ((_, d), &off) in self.sections.iter().zip(offsets) {
            out[off..off + d.len()].copy_from_slice(d);
        }
    }

    /// Stream to `w`, as [`Writer::to_bytes`] lays the bytes out; returns bytes written.
    pub fn write_to(&self, w: &mut impl Write) -> Result<u64> {
        let (h, offsets, total) = self.layout();
        w.write_all(&h)?;
        let mut pos = h.len();
        let zeros = [0u8; ALIGN];
        for ((_, d), &off) in self.sections.iter().zip(&offsets) {
            w.write_all(&zeros[..off - pos])?;
            w.write_all(d)?;
            pos = off + d.len();
        }
        w.write_all(&zeros[..total - pos])?;
        Ok(total as u64)
    }

    /// Write to `path` (through `<path>.partial-<pid>`, renamed into place, or removed if the
    /// write fails); returns bytes written.
    pub fn write(&self, path: &Path) -> Result<u64> {
        let mut tmp = path.as_os_str().to_owned();
        tmp.push(format!(".partial-{}", std::process::id()));
        let tmp = PathBuf::from(tmp);
        let res = self.write_via(&tmp, path);
        if res.is_err() {
            let _ = std::fs::remove_file(&tmp);
        }
        res
    }

    fn write_via(&self, tmp: &Path, path: &Path) -> Result<u64> {
        let f = std::fs::File::create(tmp).with_context(|| format!("creating {}", tmp.display()))?;
        let mut w = std::io::BufWriter::with_capacity(1 << 20, f);
        let n = self.write_to(&mut w)?;
        let f = w.into_inner().map_err(|e| e.into_error())?;
        f.sync_all()?;
        std::fs::rename(tmp, path).with_context(|| format!("moving {} into place", path.display()))?;
        Ok(n)
    }
}

/// Somewhere finished sections go: a [`Writer`], which keeps them, or a [`FileWriter`], which
/// writes each one out and frees it.
pub trait SectionSink {
    fn put<T: Pod>(&mut self, name: &str, data: Vec<T>) -> Result<()>;
}

impl SectionSink for Writer<'_> {
    fn put<T: Pod>(&mut self, name: &str, data: Vec<T>) -> Result<()> {
        self.add_owned(name, data);
        Ok(())
    }
}

/// A container written section by section, for an index too large to hold twice: the space
/// for the section table is left at the start and filled in by [`FileWriter::finish`]. The
/// file is `<path>.partial-<pid>` until then, and is removed if the writer is dropped first.
pub struct FileWriter {
    tmp: PathBuf,
    path: PathBuf,
    file: Option<std::fs::File>,
    pos: usize,
    max_sections: usize,
    table: Vec<(String, usize, usize)>,
}

impl FileWriter {
    /// A container at `path` of at most `max_sections` sections.
    pub fn create(path: &Path, max_sections: usize) -> Result<Self> {
        let mut tmp = path.as_os_str().to_owned();
        tmp.push(format!(".partial-{}", std::process::id()));
        let tmp = PathBuf::from(tmp);
        let mut file = std::fs::File::create(&tmp).with_context(|| format!("creating {}", tmp.display()))?;
        let pos = (16 + max_sections * (NAME_LEN + 16)).next_multiple_of(ALIGN);
        file.write_all(&vec![0u8; pos])?;
        Ok(FileWriter {
            tmp,
            path: path.to_path_buf(),
            file: Some(file),
            pos,
            max_sections,
            table: Vec::new(),
        })
    }

    /// Write the section table and move the file into place; returns bytes written.
    pub fn finish(mut self) -> Result<u64> {
        use std::io::{Seek, SeekFrom};
        let mut file = self.file.take().expect("open until finished");
        let mut h = vec![0u8; 16 + self.table.len() * (NAME_LEN + 16)];
        h[..8].copy_from_slice(MAGIC);
        h[8..12].copy_from_slice(&FORMAT_VERSION.to_le_bytes());
        h[12..16].copy_from_slice(&(self.table.len() as u32).to_le_bytes());
        for (i, (name, off, len)) in self.table.iter().enumerate() {
            let e = 16 + i * (NAME_LEN + 16);
            h[e..e + name.len()].copy_from_slice(name.as_bytes());
            h[e + NAME_LEN..e + NAME_LEN + 8].copy_from_slice(&(*off as u64).to_le_bytes());
            h[e + NAME_LEN + 8..e + NAME_LEN + 16].copy_from_slice(&(*len as u64).to_le_bytes());
        }
        let res = (|| -> Result<()> {
            file.seek(SeekFrom::Start(0))?;
            file.write_all(&h)?;
            file.sync_all()?;
            drop(file);
            std::fs::rename(&self.tmp, &self.path)
                .with_context(|| format!("moving {} into place", self.path.display()))
        })();
        if res.is_err() {
            let _ = std::fs::remove_file(&self.tmp);
        }
        res.map(|()| self.pos as u64)
    }
}

impl SectionSink for FileWriter {
    fn put<T: Pod>(&mut self, name: &str, data: Vec<T>) -> Result<()> {
        assert!(name.len() <= NAME_LEN, "section name too long: {name}");
        assert!(
            self.table.len() < self.max_sections,
            "more sections than declared"
        );
        let bytes: &[u8] = bytemuck::cast_slice(&data);
        let file = self.file.as_mut().expect("open until finished");
        file.write_all(bytes)?;
        self.table.push((name.to_string(), self.pos, bytes.len()));
        let end = self.pos + bytes.len();
        self.pos = end.next_multiple_of(ALIGN);
        file.write_all(&[0u8; ALIGN][..self.pos - end])?;
        Ok(())
    }
}

impl Drop for FileWriter {
    fn drop(&mut self) {
        if self.file.take().is_some() {
            let _ = std::fs::remove_file(&self.tmp);
        }
    }
}

enum Backing {
    Mmap(Mmap),
    /// In-memory copy, held as `u64`s so it is 8-byte aligned.
    Owned(Vec<u64>, usize),
}

/// A loaded container.
pub struct Store {
    backing: Backing,
    sections: HashMap<String, (usize, usize)>,
}

impl Store {
    pub fn open(path: &Path) -> Result<Self> {
        let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        // SAFETY: the file is opened read-only and not modified while mapped.
        let m = unsafe { Mmap::map(&f)? };
        Self::new(Backing::Mmap(m))
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let mut words = vec![0u64; bytes.len().div_ceil(8)];
        bytemuck::cast_slice_mut::<u64, u8>(&mut words)[..bytes.len()].copy_from_slice(bytes);
        Self::from_words(words, bytes.len())
    }

    /// A container over the first `len` bytes of `words` (no copy).
    pub(crate) fn from_words(words: Vec<u64>, len: usize) -> Result<Self> {
        Self::new(Backing::Owned(words, len))
    }

    fn bytes_of(b: &Backing) -> &[u8] {
        match b {
            Backing::Mmap(m) => m,
            Backing::Owned(w, n) => &bytemuck::cast_slice(w)[..*n],
        }
    }

    fn new(backing: Backing) -> Result<Self> {
        let b = Self::bytes_of(&backing);
        ensure!(
            b.len() >= 16 && &b[..8] == MAGIC,
            "not a spearfish index (bad magic)"
        );
        let version = u32::from_le_bytes(b[8..12].try_into().unwrap());
        if version != FORMAT_VERSION {
            bail!("index format version {version}, this build reads {FORMAT_VERSION}; rebuild the index");
        }
        let n = u32::from_le_bytes(b[12..16].try_into().unwrap()) as usize;
        let mut sections = HashMap::with_capacity(n);
        for i in 0..n {
            let e = 16 + i * (NAME_LEN + 16);
            ensure!(e + NAME_LEN + 16 <= b.len(), "truncated section table");
            let name = std::str::from_utf8(&b[e..e + NAME_LEN])?
                .trim_end_matches('\0')
                .to_string();
            let off = u64::from_le_bytes(b[e + NAME_LEN..e + NAME_LEN + 8].try_into().unwrap()) as usize;
            let len = u64::from_le_bytes(b[e + NAME_LEN + 8..e + NAME_LEN + 16].try_into().unwrap()) as usize;
            ensure!(
                off.is_multiple_of(ALIGN) && off.checked_add(len).is_some_and(|e| e <= b.len()),
                "section {name} out of bounds"
            );
            sections.insert(name, (off, len));
        }
        Ok(Store { backing, sections })
    }

    pub fn get<T: Pod>(&self, name: &str) -> Result<&[T]> {
        let &(off, len) = self
            .sections
            .get(name)
            .with_context(|| format!("index has no section {name}"))?;
        let bytes = &Self::bytes_of(&self.backing)[off..off + len];
        bytemuck::try_cast_slice(bytes).map_err(|e| anyhow::anyhow!("section {name}: {e:?}"))
    }

    /// Whether the container has section `name`.
    pub fn has(&self, name: &str) -> bool {
        self.sections.contains_key(name)
    }

    /// Section sizes in bytes, by name.
    pub fn sizes(&self) -> Vec<(String, usize)> {
        let mut v: Vec<_> = self.sections.iter().map(|(k, &(_, l))| (k.clone(), l)).collect();
        v.sort();
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        let mut w = Writer::default();
        w.add("a", &[1u64, 2, 3]);
        w.add("b", &[7u32, 8]);
        w.add("text", b"hello");
        w.add_owned("owned", vec![9u32, 10, 11]);
        let bytes = w.to_bytes();
        let mut streamed = Vec::new();
        assert_eq!(w.write_to(&mut streamed).unwrap(), bytes.len() as u64);
        assert_eq!(streamed, bytes);
        let (words, len) = w.to_words();
        assert_eq!(&bytemuck::cast_slice::<u64, u8>(&words)[..len], &bytes[..]);
        let s = Store::from_bytes(&bytes).unwrap();
        assert_eq!(s.get::<u32>("owned").unwrap(), &[9, 10, 11]);
        assert!(s.has("a") && !s.has("missing"));
        assert_eq!(s.get::<u64>("a").unwrap(), &[1, 2, 3]);
        assert_eq!(s.get::<u32>("b").unwrap(), &[7, 8]);
        assert_eq!(s.get::<u8>("text").unwrap(), b"hello");
        assert!(s.get::<u8>("missing").is_err());
        let (a, b) = (vec![1u64, 2, 3], vec![7u32, 8]);
        let mut w = Writer::default();
        w.add("a", &a);
        w.add("b", &b);
        w.add("text", b"hello");
        w.add_owned("owned", vec![9u32, 10, 11]);
        let s = w.into_store().unwrap();
        assert_eq!(Store::bytes_of(&s.backing), &bytes[..]);

        // Section by section to a file: the same sections, and no file if not finished.
        let dir = std::env::temp_dir().join(format!("spearfish-store-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("c.idx");
        let mut f = FileWriter::create(&path, 5).unwrap();
        f.put("a", a.clone()).unwrap();
        f.put("text", b"hello".to_vec()).unwrap();
        f.put("empty", Vec::<u32>::new()).unwrap();
        f.put("owned", vec![9u32, 10, 11]).unwrap();
        let n = f.finish().unwrap();
        assert_eq!(n, std::fs::metadata(&path).unwrap().len());
        let s = Store::open(&path).unwrap();
        assert_eq!(s.get::<u64>("a").unwrap(), &[1, 2, 3]);
        assert_eq!(s.get::<u8>("text").unwrap(), b"hello");
        assert!(s.get::<u32>("empty").unwrap().is_empty());
        assert_eq!(s.get::<u32>("owned").unwrap(), &[9, 10, 11]);
        drop(s);
        let other = dir.join("d.idx");
        let mut f = FileWriter::create(&other, 1).unwrap();
        f.put("a", a.clone()).unwrap();
        drop(f);
        assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 1);
        std::fs::remove_dir_all(&dir).unwrap();
        assert!(Store::from_bytes(b"not an index at all").is_err());
        // A section length that wraps around past the end is an error, not a slice panic.
        let mut bad = bytes.clone();
        let e = 16 + NAME_LEN + 8;
        bad[e..e + 8].copy_from_slice(&(u64::MAX - 63).to_le_bytes());
        assert!(Store::from_bytes(&bad).is_err());
    }
}
