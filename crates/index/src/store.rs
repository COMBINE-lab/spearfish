//! A flat container of named, typed arrays: written once, memory-mapped and read zero-copy.
//!
//! Layout (little-endian): magic `SPFISH01`, `u32` format version, `u32` section count, a
//! table of `(name: [u8; 24], offset: u64, len_bytes: u64)` entries, then each section's
//! bytes at a 64-byte-aligned offset.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result, bail, ensure};
use bytemuck::Pod;
use memmap2::Mmap;

const MAGIC: &[u8; 8] = b"SPFISH01";
pub const FORMAT_VERSION: u32 = 2;
const NAME_LEN: usize = 24;
const ALIGN: usize = 64;

/// Sections to write, in order.
#[derive(Default)]
pub struct Writer {
    sections: Vec<(String, Vec<u8>)>,
}

impl Writer {
    pub fn add<T: Pod>(&mut self, name: &str, data: &[T]) {
        assert!(name.len() <= NAME_LEN, "section name too long: {name}");
        self.sections
            .push((name.to_string(), bytemuck::cast_slice(data).to_vec()));
    }

    /// Serialize to bytes (as written to disk).
    pub fn to_bytes(&self) -> Vec<u8> {
        let header = 16 + self.sections.len() * (NAME_LEN + 16);
        let mut offsets = Vec::with_capacity(self.sections.len());
        let mut pos = header.next_multiple_of(ALIGN);
        for (_, d) in &self.sections {
            offsets.push(pos);
            pos = (pos + d.len()).next_multiple_of(ALIGN);
        }
        let mut out = vec![0u8; pos];
        out[..8].copy_from_slice(MAGIC);
        out[8..12].copy_from_slice(&FORMAT_VERSION.to_le_bytes());
        out[12..16].copy_from_slice(&(self.sections.len() as u32).to_le_bytes());
        for (i, ((name, d), &off)) in self.sections.iter().zip(&offsets).enumerate() {
            let e = 16 + i * (NAME_LEN + 16);
            out[e..e + name.len()].copy_from_slice(name.as_bytes());
            out[e + NAME_LEN..e + NAME_LEN + 8].copy_from_slice(&(off as u64).to_le_bytes());
            out[e + NAME_LEN + 8..e + NAME_LEN + 16].copy_from_slice(&(d.len() as u64).to_le_bytes());
            out[off..off + d.len()].copy_from_slice(d);
        }
        out
    }

    pub fn write(&self, path: &Path) -> Result<u64> {
        let bytes = self.to_bytes();
        let tmp = path.with_extension("tmp");
        let mut f = std::fs::File::create(&tmp).with_context(|| format!("creating {}", tmp.display()))?;
        f.write_all(&bytes)?;
        f.sync_all()?;
        std::fs::rename(&tmp, path)?;
        Ok(bytes.len() as u64)
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
        Self::new(Backing::Owned(words, bytes.len()))
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
        let s = Store::from_bytes(&w.to_bytes()).unwrap();
        assert_eq!(s.get::<u64>("a").unwrap(), &[1, 2, 3]);
        assert_eq!(s.get::<u32>("b").unwrap(), &[7, 8]);
        assert_eq!(s.get::<u8>("text").unwrap(), b"hello");
        assert!(s.get::<u8>("missing").is_err());
        assert!(Store::from_bytes(b"not an index at all").is_err());
    }
}
