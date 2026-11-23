//! BAM output: records encoded here, BGZF blocks compressed by whichever thread made them.
//!
//! A BGZF file is a series of independent gzip members, so worker threads can each turn their
//! records into complete blocks and the writer only has to append bytes in order. The DEFLATE
//! engine is chosen at build time: `zlib-rs` through flate2 (the default), or with a build
//! feature `zflate` (a Rust port of libdeflate; fastest, needs Rust 1.91) or `libdeflate` (the C
//! library). The record encoder does not
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
    if cfg!(feature = "zflate") {
        "zflate"
    } else if cfg!(feature = "libdeflate") {
        "libdeflate"
    } else {
        "zlib-rs"
    }
}

#[cfg(feature = "zflate")]
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

#[cfg(all(feature = "libdeflate", not(feature = "zflate")))]
fn deflate(src: &[u8], level: Level, dst: &mut [u8]) -> (usize, u32) {
    let l = libdeflater::CompressionLvl::new(level.0 as i32).expect("level 0 to 12");
    let n = libdeflater::Compressor::new(l)
        .deflate_compress(src, dst)
        .expect("output buffer large enough");
    (n, libdeflater::crc32(src))
}

#[cfg(not(any(feature = "zflate", feature = "libdeflate")))]
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
/// `program` is as for [`crate::map::sam_header`].
pub fn header(ix: &GenomeIndex, program: Option<(&str, &str)>, out: &mut Vec<u8>) {
    let mut text = Vec::new();
    crate::map::sam_header(ix, program, &mut text);
    out.extend_from_slice(b"BAM\x01");
    out.extend_from_slice(&(text.len() as u32).to_le_bytes());
    out.extend_from_slice(&text);
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

/// 4-bit code of a base.
fn code4(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        _ => 15,
    }
}

/// The fixed fields, name, CIGAR, bases and qualities of a record; returns where the record
/// starts, for [`finish`]. `seq` and `qual` are as sequenced, and reversed here for a
/// reverse-strand record.
#[allow(clippy::too_many_arguments)]
fn start(
    name: &str,
    rid: i32,
    pos: i32,
    mapq: u8,
    bin: u16,
    flag: u16,
    ops: &[u32],
    seq: &[u8],
    qual: Option<&[u8]>,
    reverse: bool,
    out: &mut Vec<u8>,
) -> usize {
    let start = out.len();
    out.extend_from_slice(&[0; 4]);
    out.extend_from_slice(&rid.to_le_bytes());
    out.extend_from_slice(&pos.to_le_bytes());
    out.push(name.len().min(254) as u8 + 1);
    out.push(mapq);
    out.extend_from_slice(&bin.to_le_bytes());
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
    let comp = |c: u8| match c {
        1 => 8,
        2 => 4,
        4 => 2,
        8 => 1,
        x => x,
    };
    let at = |i: usize| {
        if reverse {
            comp(code4(seq[seq.len() - 1 - i]))
        } else {
            code4(seq[i])
        }
    };
    for i in (0..seq.len()).step_by(2) {
        out.push(at(i) << 4 | if i + 1 < seq.len() { at(i + 1) } else { 0 });
    }
    // Phred qualities (0xff: none)
    match qual.filter(|q| q.len() == seq.len()) {
        Some(q) if reverse => out.extend(q.iter().rev().map(|&c| c.saturating_sub(33))),
        Some(q) => out.extend(q.iter().map(|&c| c.saturating_sub(33))),
        None => out.extend(std::iter::repeat_n(0xff, seq.len())),
    }
    start
}

/// Write the length of the record that began at `start`.
fn finish(start: usize, out: &mut [u8]) {
    let n = (out.len() - start - 4) as u32;
    out[start..start + 4].copy_from_slice(&n.to_le_bytes());
}

/// Append the BAM record of a hit. `read` is the read as sequenced and `qual` its base
/// qualities (Phred+33), if it has them; secondary records carry neither. Tags: `AS`, `cm`,
/// and `XS` (transcript strand) when known. With `edits`, `ops` already holds the hit's
/// base-level CIGAR (from [`crate::base::Aligner::cigar`]) and an `NM` tag is written;
/// otherwise the block-level CIGAR is made here.
#[allow(clippy::too_many_arguments)]
pub fn record(
    ix: &GenomeIndex,
    name: &str,
    read: &[u8],
    qual: Option<&[u8]>,
    h: &SplicedHit,
    mapq: u8,
    ops: &mut Vec<u32>,
    edits: Option<u32>,
    out: &mut Vec<u8>,
) {
    let qlen = read.len() as u32;
    let (first, last) = (h.blocks[0], h.blocks[h.blocks.len() - 1]);
    let (rid, pos) = ix.locate(first.ts);
    let ref_end = pos + (last.te - first.ts);
    // (a record holds 65535 operations at most; a longer alignment falls back to blocks)
    let edits = edits.filter(|_| ops.len() <= 65535);
    if edits.is_none() {
        cigar(&h.blocks, qlen, ops);
    }
    let (seq, qual): (&[u8], Option<&[u8]>) = if h.primary { (read, qual) } else { (b"", None) };
    let flag: u16 = if h.reverse { 16 } else { 0 } | if h.primary { 0 } else { 256 };
    let start = start(
        name,
        rid as i32,
        pos as i32,
        mapq,
        reg2bin(pos, ref_end),
        flag,
        ops,
        seq,
        qual,
        h.reverse,
        out,
    );
    out.extend_from_slice(b"ASi");
    out.extend_from_slice(&h.score.to_le_bytes());
    out.extend_from_slice(b"cmi");
    out.extend_from_slice(&(h.anchors as i32).to_le_bytes());
    if h.tx_strand != b'.' {
        out.extend_from_slice(b"XSA");
        out.push(h.tx_strand);
    }
    if let Some(e) = edits {
        out.extend_from_slice(b"NMi");
        out.extend_from_slice(&(e as i32).to_le_bytes());
    }
    finish(start, out);
}

/// Append the BAM record of a read that did not map (flag 4), with its sequence and
/// qualities.
pub fn unmapped(name: &str, read: &[u8], qual: Option<&[u8]>, out: &mut Vec<u8>) {
    // (bin 4680 is that of a record without a position)
    let start = start(name, -1, -1, 0, 4680, 4, &[], read, qual, false, out);
    finish(start, out);
}
