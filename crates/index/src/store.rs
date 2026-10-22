//! A flat container of named, typed arrays: written once, memory-mapped and read zero-copy.
//!
//! Layout (little-endian): magic `SPFISH01`, `u32` format version, `u32` section count, a
//! table of `(name: [u8; 24], offset: u64, len_bytes: u64)` entries, then each section's
//! bytes at a 64-byte-aligned offset.

use std::borrow::Cow;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

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
    sections: Vec<(String, Cow<'a, [u8]>)>,
}

impl<'a> Writer<'a> {
    /// Add a section borrowing `data`.
    pub fn add<T: Pod>(&mut self, name: &str, data: &'a [T]) {
        assert!(name.len() <= NAME_LEN, "section name too long: {name}");
        self.sections
            .push((name.to_string(), Cow::Borrowed(bytemuck::cast_slice(data))));
    }

    /// Add a section owning its data.
    pub fn add_owned<T: Pod>(&mut self, name: &str, data: Vec<T>) {
        assert!(name.len() <= NAME_LEN, "section name too long: {name}");
        let bytes = bytemuck::cast_slice(&data).to_vec();
        self.sections.push((name.to_string(), Cow::Owned(bytes)));
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

    /// Write to `path` (through a temporary file renamed into place); returns bytes written.
    pub fn write(&self, path: &Path) -> Result<u64> {
        let tmp = path.with_extension("tmp");
        let f = std::fs::File::create(&tmp).with_context(|| format!("creating {}", tmp.display()))?;
        let mut w = std::io::BufWriter::with_capacity(1 << 20, f);
        let n = self.write_to(&mut w)?;
        let f = w.into_inner().map_err(|e| e.into_error())?;
        f.sync_all()?;
        std::fs::rename(&tmp, path)?;
        Ok(n)
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
                off.is_multiple_of(ALIGN) && off + len <= b.len(),
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
        assert!(Store::from_bytes(b"not an index at all").is_err());
    }
}
