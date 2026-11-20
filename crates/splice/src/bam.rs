//! BAM output: records encoded here, BGZF blocks compressed by whichever thread made them.
//!
//! A BGZF file is a series of independent gzip members, so worker threads can each turn their
//! records into complete blocks and the writer only has to append bytes in order. The DEFLATE
//! engine is chosen at build time: `zflate` (the default; a Rust port of libdeflate),
//! `libdeflate` (the C library) or `zlib-rs` (through flate2). The record encoder does not
//! depend on the engine, so their outputs decompress to the same bytes.

use crate::index::GenomeIndex;
use crate::map::{Block, SplicedHit};

/// Uncompressed bytes per BGZF block (the format's limit is 64 KiB per member).
const BLOCK: usize = 0xff00;

/// Compression level: 0 stores, 1 is fast, higher is smaller.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Level(pub u8);

/// The name of the DEFLATE engine this build uses.
pub fn engine() -> &'static str {
    if cfg!(feature = "libdeflate") {
        "libdeflate"
    } else if cfg!(feature = "zlib-rs") {
        "zlib-rs"
    } else {
        "zflate"
    }
}

#[cfg(not(any(feature = "libdeflate", feature = "zlib-rs")))]
fn deflate(src: &[u8], level: Level, dst: &mut [u8]) -> (usize, u32) {
    thread_local! {
        static C: std::cell::RefCell<Option<(u8, zflate::Compressor)>> = const { std::cell::RefCell::new(None) };
    }
    C.with_borrow_mut(|c| {
        if c.as_ref().is_none_or(|x| x.0 != level.0) {
            let l = zflate::CompressionLevel::new(level.0).expect("level 0 to 12");
            *c = Some((
                level.0,
                zflate::Compressor::try_new(l).expect("allocating the compressor"),
            ));
        }
        let n = c
            .as_mut()
            .unwrap()
            .1
            .deflate_compress(src, dst)
            .expect("output buffer large enough");
        (n, zflate::crc32(0, src))
    })
}

#[cfg(feature = "libdeflate")]
fn deflate(src: &[u8], level: Level, dst: &mut [u8]) -> (usize, u32) {
    let l = libdeflater::CompressionLvl::new(level.0 as i32).expect("level 0 to 12");
    let n = libdeflater::Compressor::new(l)
        .deflate_compress(src, dst)
        .expect("output buffer large enough");
    (n, libdeflater::crc32(src))
}

#[cfg(all(feature = "zlib-rs", not(feature = "libdeflate")))]
fn deflate(src: &[u8], level: Level, dst: &mut [u8]) -> (usize, u32) {
    let mut c = flate2::Compress::new(flate2::Compression::new(level.0.min(9) as u32), false);
    c.compress(src, dst, flate2::FlushCompress::Finish)
        .expect("output buffer large enough");
    let mut h = flate2::Crc::new();
    h.update(src);
    (c.total_out() as usize, h.sum())
}

/// Append `data` to `out` as complete BGZF blocks.
pub fn bgzf(data: &[u8], level: Level, out: &mut Vec<u8>) {
    let mut buf = vec![0u8; BLOCK + BLOCK / 8 + 1024];
    for blk in data.chunks(BLOCK) {
        let (n, crc) = deflate(blk, level, &mut buf);
        out.extend_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0]);
        out.extend_from_slice(&((n + 25) as u16).to_le_bytes());
        out.extend_from_slice(&buf[..n]);
        out.extend_from_slice(&crc.to_le_bytes());
        out.extend_from_slice(&(blk.len() as u32).to_le_bytes());
    }
}

/// The empty block that ends a BGZF file.
pub const EOF: [u8; 28] = [
    0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0x1b, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

/// The uncompressed BAM header (magic, SAM text, reference list) for the index's sequences.
pub fn header(ix: &GenomeIndex, out: &mut Vec<u8>) {
    let mut text = String::new();
    crate::map::sam_header(ix, &mut text);
    out.extend_from_slice(b"BAM\x01");
    out.extend_from_slice(&(text.len() as u32).to_le_bytes());
    out.extend_from_slice(text.as_bytes());
    out.extend_from_slice(&(ix.names().len() as u32).to_le_bytes());
    for (i, n) in ix.names().iter().enumerate() {
        out.extend_from_slice(&(n.len() as u32 + 1).to_le_bytes());
        out.extend_from_slice(n.as_bytes());
        out.push(0);
        out.extend_from_slice(&(ix.offsets()[i + 1] - ix.offsets()[i]).to_le_bytes());
    }
}

/// The UCSC bin of `[beg, end)`.
fn reg2bin(beg: u32, end: u32) -> u16 {
    let end = end.saturating_sub(1);
    for (shift, off) in [(14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)] {
        if beg >> shift == end >> shift {
            return (off + (beg >> shift)) as u16;
        }
    }
    0
}

/// CIGAR operations of a hit, as `len << 4 | op`: `M` per block, `I`/`D` for a length
/// difference inside it, `N` for introns, `S` for unaligned read ends.
fn cigar(blocks: &[Block], qlen: u32, out: &mut Vec<u32>) {
    let (m, ins, del, skip, clip) = (0u32, 1u32, 2u32, 3u32, 4u32);
    out.clear();
    let (first, last) = (blocks[0], blocks[blocks.len() - 1]);
    if first.qs > 0 {
        out.push(first.qs << 4 | clip);
    }
    for (i, b) in blocks.iter().enumerate() {
        if i > 0 {
            out.push((b.ts - blocks[i - 1].te) << 4 | skip);
        }
        let (lq, lt) = (b.qe - b.qs, b.te - b.ts);
        let mm = lq.min(lt);
        if mm > 0 {
            out.push(mm << 4 | m);
        }
        if lq > lt {
            out.push((lq - lt) << 4 | ins);
        } else if lt > lq {
            out.push((lt - lq) << 4 | del);
        }
    }
    if last.qe < qlen {
        out.push((qlen - last.qe) << 4 | clip);
    }
}

/// Append the BAM record of a hit. `read` is the read as sequenced; secondary records carry
/// no sequence. Tags: `AS`, `cm`, and `XS` (transcript strand) when known.
pub fn record(
    ix: &GenomeIndex,
    name: &str,
    read: &[u8],
    h: &SplicedHit,
    mapq: u8,
    ops: &mut Vec<u32>,
    out: &mut Vec<u8>,
) {
    let qlen = read.len() as u32;
    let (first, last) = (h.blocks[0], h.blocks[h.blocks.len() - 1]);
    let (rid, pos) = ix.locate(first.ts);
    let ref_end = pos + (last.te - first.ts);
    cigar(&h.blocks, qlen, ops);
    let seq: &[u8] = if h.primary { read } else { b"" };
    let flag: u16 = if h.reverse { 16 } else { 0 } | if h.primary { 0 } else { 256 };
    let start = out.len();
    out.extend_from_slice(&[0; 4]);
    out.extend_from_slice(&(rid as i32).to_le_bytes());
    out.extend_from_slice(&(pos as i32).to_le_bytes());
    out.push(name.len().min(254) as u8 + 1);
    out.push(mapq);
    out.extend_from_slice(&reg2bin(pos, ref_end).to_le_bytes());
    out.extend_from_slice(&(ops.len().min(65535) as u16).to_le_bytes());
    out.extend_from_slice(&flag.to_le_bytes());
    out.extend_from_slice(&(seq.len() as u32).to_le_bytes());
    out.extend_from_slice(&(-1i32).to_le_bytes());
    out.extend_from_slice(&(-1i32).to_le_bytes());
    out.extend_from_slice(&0i32.to_le_bytes());
    out.extend_from_slice(&name.as_bytes()[..name.len().min(254)]);
    out.push(0);
    for &c in ops.iter().take(65535) {
        out.extend_from_slice(&c.to_le_bytes());
    }
    // 4-bit bases, in the orientation the read maps in
    let code = |b: u8| match b.to_ascii_uppercase() {
        b'A' => 1u8,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        _ => 15,
    };
    let comp = |c: u8| match c {
        1 => 8,
        2 => 4,
        4 => 2,
        8 => 1,
        x => x,
    };
    let at = |i: usize| {
        if h.reverse {
            comp(code(seq[seq.len() - 1 - i]))
        } else {
            code(seq[i])
        }
    };
    for i in (0..seq.len()).step_by(2) {
        out.push(at(i) << 4 | if i + 1 < seq.len() { at(i + 1) } else { 0 });
    }
    out.extend(std::iter::repeat_n(0xff, seq.len()));
    out.extend_from_slice(b"ASi");
    out.extend_from_slice(&h.score.to_le_bytes());
    out.extend_from_slice(b"cmi");
    out.extend_from_slice(&(h.anchors as i32).to_le_bytes());
    if h.tx_strand != b'.' {
        out.extend_from_slice(b"XSA");
        out.push(h.tx_strand);
    }
    let n = (out.len() - start - 4) as u32;
    out[start..start + 4].copy_from_slice(&n.to_le_bytes());
}
