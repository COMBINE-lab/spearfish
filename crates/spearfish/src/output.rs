//! Output formats of mapped reads, as `spearfish map` writes them.
//!
//! - **TSV** ([`tsv_lines`], after [`TSV_HEADER`]): one line per reported (read, transcript),
//!   best first; a read's lines are contiguous. The strand is per line: a read whose
//!   orientations nearly tie (antisense transcripts) is reported on both. Columns: read,
//!   read_len, strand, rank, transcript, score, qstart, qend, tstart, tend, anchors, anchored,
//!   max_disc, flags. qstart/qend are raw read coordinates on the read's forward strand and
//!   tstart/tend raw transcript coordinates of the alignment, extended past the outermost
//!   anchors for hits within `extend_within` of the best score (with `extend_ends`); other
//!   hits' spans run between their outermost anchors and carry the flag `anchor_span`.
//!   anchored counts compressed bases. An unmapped read gets one line with transcript `*`.
//! - **PAF** ([`paf_lines`]): the same hits as PAF lines (no CIGAR; unmapped reads omitted).
//! - **RAD** (feature `rad`; `rad_writer`, `rad_record`): `bulk_long_read` records (see
//!   the `long-read-rad` crate), one per mapped read, which oarfish quantifies directly. The
//!   file's tags carry the preset, the score scale, the index identifier, and whether spans
//!   were extended to the read ends.
//!
//! The text writers append to a byte buffer, so a caller can format a batch of reads on each
//! thread and write the buffers whole. Numbers are written digit for digit as Rust's `{}` and
//! `{:.N}` would, without going through `fmt`.

use std::io::Write;

use crate::{Index, MappedRead};

/// The TSV header line, newline included.
pub const TSV_HEADER: &str = "read\tread_len\tstrand\trank\ttranscript\tscore\tqstart\tqend\ttstart\ttend\tanchors\tanchored\tmax_disc\tflags\n";

/// Append `x` in decimal (as `{}` would).
fn push_uint(buf: &mut Vec<u8>, mut x: u64) {
    let mut tmp = [0u8; 20];
    let mut i = tmp.len();
    loop {
        i -= 1;
        tmp[i] = b'0' + (x % 10) as u8;
        x /= 10;
        if x == 0 {
            break;
        }
    }
    buf.extend_from_slice(&tmp[i..]);
}

/// Append `v` with `dec` (at most 4) decimals, digit for digit as `{:.dec$}` would: `v * 10^dec`
/// is exact in an `f64`, and std rounds the exact value half to even. Values that are not finite
/// or not below 1e9 go through `write!`.
fn push_fixed(buf: &mut Vec<u8>, v: f32, dec: u32) {
    debug_assert!(dec <= 4);
    if !v.is_finite() || v.abs() >= 1e9 {
        let _ = write!(buf, "{:.*}", dec as usize, v);
        return;
    }
    if v.is_sign_negative() {
        buf.push(b'-');
    }
    let p = 10u64.pow(dec);
    let a = ((v as f64).abs() * p as f64).round_ties_even() as u64;
    push_uint(buf, a / p);
    if dec > 0 {
        buf.push(b'.');
        let f = a % p;
        let mut d = p / 10;
        while d > 0 {
            buf.push(b'0' + ((f / d) % 10) as u8);
            d /= 10;
        }
    }
}

/// Append the TSV lines of one read named `name` (see the [module docs](self)): one per hit,
/// best first, or one `*` line if it is unmapped.
pub fn tsv_lines(buf: &mut Vec<u8>, name: &str, r: &MappedRead, idx: &Index) {
    if r.hits.is_empty() {
        let _ = writeln!(
            buf,
            "{name}\t{}\t*\t0\t*\t0\t0\t0\t0\t0\t0\t0\t0\tunmapped",
            r.len
        );
        return;
    }
    let read_flags = read_flags(r);
    for (rank, h) in r.hits.iter().enumerate() {
        buf.extend_from_slice(name.as_bytes());
        buf.push(b'\t');
        push_uint(buf, r.len as u64);
        buf.extend_from_slice(if h.reverse { b"\t-\t" } else { b"\t+\t" });
        push_uint(buf, rank as u64);
        buf.push(b'\t');
        buf.extend_from_slice(idx.tx_name(h.tx).as_bytes());
        buf.push(b'\t');
        push_fixed(buf, h.score, 2);
        for x in [
            h.query_start,
            h.query_end,
            h.target_start,
            h.target_end,
            h.anchors,
            h.anchored,
        ] {
            buf.push(b'\t');
            push_uint(buf, x as u64);
        }
        buf.push(b'\t');
        push_fixed(buf, h.max_discrepancy, 0);
        buf.push(b'\t');
        if !push_hit_flags(buf, read_flags, h.extended) {
            buf.push(b'-');
        }
        buf.push(b'\n');
    }
}

/// The flags every TSV line of a read carries: `rescued`, `low_coverage`, both, or `-`.
pub fn read_flags(r: &MappedRead) -> &'static str {
    match (r.low_coverage, r.rescued) {
        (true, true) => "rescued,low_coverage",
        (true, false) => "low_coverage",
        (false, true) => "rescued",
        (false, false) => "-",
    }
}

/// Append a hit's flags (the read's, then `anchor_span` if the span was not extended); false
/// (nothing written) if it has none.
fn push_hit_flags(buf: &mut Vec<u8>, read_flags: &str, extended: bool) -> bool {
    let has_read = read_flags != "-";
    if has_read {
        buf.extend_from_slice(read_flags.as_bytes());
    }
    if !extended {
        buf.extend_from_slice(if has_read { b",anchor_span" } else { b"anchor_span" });
    }
    has_read || !extended
}

/// PAF lines of one read (none if unmapped). Columns as minimap2 writes them without `-c`:
/// residue matches are the read bases the hit's seeds cover (homopolymer-compressed), the
/// block length the longer of the two spans, and the mapping quality 255 (not computed: tied
/// isoforms are reported as hits instead). Tags: `tp:A` P for the best hit, S for the others;
/// `AS:i` the score rounded; `sc:f` the score; `cm:i` anchors in the chain; `rk:i` the rank;
/// `fl:Z` the TSV flags, when any.
pub fn paf_lines(buf: &mut Vec<u8>, name: &str, r: &MappedRead, idx: &Index) {
    let read_flags = read_flags(r);
    for (rank, h) in r.hits.iter().enumerate() {
        let block = (h.query_end - h.query_start).max(h.target_end - h.target_start);
        buf.extend_from_slice(name.as_bytes());
        for x in [r.len, h.query_start, h.query_end] {
            buf.push(b'\t');
            push_uint(buf, x as u64);
        }
        buf.extend_from_slice(if h.reverse { b"\t-\t" } else { b"\t+\t" });
        buf.extend_from_slice(idx.tx_name(h.tx).as_bytes());
        for x in [
            idx.tx_len(h.tx),
            h.target_start,
            h.target_end,
            h.anchored.min(block),
            block,
        ] {
            buf.push(b'\t');
            push_uint(buf, x as u64);
        }
        buf.extend_from_slice(if rank == 0 {
            b"\t255\ttp:A:P\tAS:i:"
        } else {
            b"\t255\ttp:A:S\tAS:i:"
        });
        push_fixed(buf, h.score, 0);
        buf.extend_from_slice(b"\tsc:f:");
        push_fixed(buf, h.score, 2);
        buf.extend_from_slice(b"\tcm:i:");
        push_uint(buf, h.anchors as u64);
        buf.extend_from_slice(b"\trk:i:");
        push_uint(buf, rank as u64);
        let mark = buf.len();
        buf.extend_from_slice(b"\tfl:Z:");
        if !push_hit_flags(buf, read_flags, h.extended) {
            buf.truncate(mark);
        }
        buf.push(b'\n');
    }
}

/// The header of a RAD file for reads mapped by `mapper`: the preset, sequencing technology,
/// score scale, mapped strands, whether spans were extended to the read ends, and the index
/// identifier, with the settings changed from the preset (`None` if unknown) and every
/// parameter as JSON. `command_line` is recorded as given.
#[cfg(feature = "rad")]
pub fn rad_meta(
    mapper: &crate::Mapper,
    settings: Option<&crate::MapSettings>,
    command_line: &str,
) -> anyhow::Result<long_read_rad::Meta> {
    let opts = mapper.options();
    Ok(long_read_rad::Meta {
        mapper: "spearfish".into(),
        mapper_version: crate::provenance::BuildInfo::current().short(),
        cmdline: command_line.into(),
        params: serde_json::to_string(&serde_json::json!({
            "settings": settings,
            "parameters": opts,
        }))?,
        preset: opts.preset.name().into(),
        seq_tech: opts.preset.seq_tech().into(),
        score_denom: mapper.score_denom(),
        ends_extended: opts.params.ends.is_some(),
        strands: opts.params.strands.name().into(),
        index_id: mapper.index().id(),
    })
}

/// Start a RAD file at `path` for reads mapped by `mapper` (see [`rad_meta`]), over every
/// transcript of its index. Push records through [`long_read_rad::RadWriter::sink`]s and
/// [`finish`](long_read_rad::RadWriter::finish) it with the read counts.
#[cfg(feature = "rad")]
pub fn rad_writer(
    path: &std::path::Path,
    mapper: &crate::Mapper,
    settings: Option<&crate::MapSettings>,
    command_line: &str,
    codec: long_read_rad::ChunkCodec,
) -> anyhow::Result<long_read_rad::RadWriter> {
    let idx = mapper.index();
    let refs: Vec<(&str, u32)> = (0..idx.n_transcripts() as u32)
        .map(|t| (idx.tx_name(t), idx.tx_len(t)))
        .collect();
    long_read_rad::RadWriter::create(path, &refs, &rad_meta(mapper, settings, command_line)?, codec)
}

/// Fill `rec` with the read `r`; false (and `rec` left empty) if it is unmapped, since RAD
/// files hold mapped reads only.
#[cfg(feature = "rad")]
pub fn rad_record(rec: &mut long_read_rad::ReadRecord, r: &MappedRead) -> bool {
    use long_read_rad::{Aln, aln_flags, read_flags};
    rec.clear();
    if r.hits.is_empty() {
        return false;
    }
    rec.read_len = r.len;
    rec.flags = if r.low_coverage {
        read_flags::LOW_COVERAGE
    } else {
        0
    } | if r.rescued { read_flags::RESCUED } else { 0 };
    for h in &r.hits {
        rec.push(Aln {
            tid: h.tx,
            fw: !h.reverse,
            ref_start: h.target_start,
            ref_end: h.target_end,
            query_start: h.query_start,
            query_end: h.query_end,
            score: h.score,
            flags: if h.extended { 0 } else { aln_flags::ANCHOR_SPAN },
        });
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    /// SplitMix64, seeded per record (as the CLI's simulations use it).
    struct SplitMix64(u64);

    impl SplitMix64 {
        fn for_record(seed: u64, index: u64) -> Self {
            let mut r = SplitMix64(seed ^ index.wrapping_mul(0xd1b5_4a32_d192_ed03));
            r.next_u64();
            r
        }

        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^ (z >> 31)
        }

        fn below(&mut self, n: u64) -> u64 {
            ((self.next_u64() as u128 * n as u128) >> 64) as u64
        }

        fn unit(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// `push_fixed` writes what `{:.0}` and `{:.2}` write (checked offline on every f32), and
    /// `push_uint` what `{}` writes.
    #[test]
    fn digits_match_std_formatting() {
        let mut buf = Vec::new();
        let mut check = |v: f32| {
            for dec in [0, 2] {
                buf.clear();
                push_fixed(&mut buf, v, dec);
                let want = format!("{:.*}", dec as usize, v);
                assert_eq!(
                    std::str::from_utf8(&buf).unwrap(),
                    want,
                    "{v:e} ({:#x}) to {dec} decimals",
                    v.to_bits()
                );
            }
        };
        let mut rng = SplitMix64::for_record(7, 0);
        // Any bit pattern, scores and discrepancies in their usual ranges, every 4093rd bit
        // pattern, and exact ties at both precisions.
        for _ in 0..200_000 {
            check(f32::from_bits(rng.next_u64() as u32));
            check((rng.unit() * 2e5 - 1e5) as f32);
            check((rng.unit() * 2.0 - 1.0) as f32);
        }
        for b in (0..=u32::MAX).step_by(4093) {
            check(f32::from_bits(b));
        }
        for k in 0..20_000 {
            for f in [0.125, 0.375, 0.5, 0.625, 0.875] {
                check(k as f32 + f);
                check(-(k as f32 + f));
            }
        }
        for v in [
            0.0,
            -0.0,
            -0.001,
            -0.4,
            0.005,
            9.995,
            f32::MIN_POSITIVE,
            1e-40,
            999_999_940.0,
            1e9,
            f32::MAX,
            f32::NAN,
            f32::INFINITY,
            f32::NEG_INFINITY,
        ] {
            check(v);
        }
        let mut buf = Vec::new();
        for x in [0, 1, 9, 10, 99, 100, u32::MAX as u64, u64::MAX]
            .into_iter()
            .chain((0..10_000).map(|_| rng.next_u64() >> rng.below(64)))
        {
            buf.clear();
            push_uint(&mut buf, x);
            assert_eq!(buf, x.to_string().as_bytes());
        }
    }
}
