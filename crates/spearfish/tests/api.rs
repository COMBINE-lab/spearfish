//! Library API: build an index in memory, map reads from its transcripts on both strands.

use std::path::PathBuf;
use std::sync::Arc;

use lcp_core::encode::revcomp;
use spearfish::{Index, IndexOptions, MapOptions, Mapper, Preset, read_transcripts};

fn assert_send_sync<T: Send + Sync>() {}

#[test]
fn map_toy_reads() {
    assert_send_sync::<Mapper>();
    let fa = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/toy/transcripts.fa");
    let txs = read_transcripts(&fa).unwrap();
    let (index, rep) = Index::build(txs.clone(), &IndexOptions::default()).unwrap();
    assert_eq!(rep.transcripts, txs.len());
    let index = Arc::new(index);
    let mapper = Mapper::new(index.clone(), &MapOptions::preset(Preset::Ont)).unwrap();
    for (t, (name, _, seq)) in txs.iter().enumerate() {
        assert_eq!(index.tx_name(t as u32), name);
        assert_eq!(index.tx_len(t as u32) as usize, seq.len());
        for rev in [false, true] {
            let read = if rev { revcomp(seq) } else { seq.clone() };
            let r = mapper.map(&read);
            let top = &r.hits[0];
            assert_eq!(index.tx_name(top.tx), name, "{name} rev {rev}");
            assert_eq!(top.reverse, rev);
            // Raw coordinates: the chain spans most of the read and of the transcript.
            assert!(top.query_end <= read.len() as u32 && top.query_start < top.query_end);
            assert!((top.query_end - top.query_start) as usize * 10 >= read.len() * 8);
            assert!(top.target_end <= seq.len() as u32);
            assert!((top.target_end - top.target_start) as usize * 10 >= seq.len() * 8);
        }
    }
    // Direct RNA preset: sense orientation only.
    let drna = Mapper::new(index.clone(), &MapOptions::preset(Preset::Drna)).unwrap();
    let (_, _, seq) = &txs[0];
    assert!(drna.map(&revcomp(seq)).hits.iter().all(|h| !h.reverse));
    // Chain limits set directly, without settings resolution, are checked too.
    let mut opts = MapOptions::preset(Preset::Ont);
    opts.params.chain.g_max = 3_000_000_000;
    assert!(Mapper::new(index.clone(), &opts).is_err());
}

/// A read that spells T partly as U (either case) maps as its DNA spelling does, query span
/// included: U is read as T, so a T next to a U is one homopolymer run.
#[test]
fn mixed_t_and_u_map_as_t() {
    let fa = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/toy/transcripts.fa");
    let txs = read_transcripts(&fa).unwrap();
    let index = Arc::new(Index::build(txs.clone(), &IndexOptions::default()).unwrap().0);
    for p in [Preset::Ont, Preset::Hifi, Preset::Drna] {
        let mapper = Mapper::new(index.clone(), &MapOptions::preset(p)).unwrap();
        for (_, _, seq) in &txs {
            for (rev, read) in [(false, seq.clone()), (true, revcomp(seq))] {
                // Every second base of each T run as U, alternating U and u.
                let mut mixed = read.clone();
                let mut k = 0;
                for i in 1..mixed.len() {
                    if read[i] == b'T' && read[i - 1] == b'T' && mixed[i - 1] == b'T' {
                        mixed[i] = if k % 2 == 0 { b'U' } else { b'u' };
                        k += 1;
                    }
                }
                assert!(k > 0);
                let want = mapper.map(&read);
                // Direct RNA maps the sense orientation only.
                assert!(rev && p == Preset::Drna || !want.hits.is_empty());
                assert_eq!(mapper.map(&mixed), want, "{p:?} rev {rev}");
            }
        }
    }
}

/// The settings builder gives the same mapper as the preset's options, and records what it
/// resolved to.
#[test]
fn builder_matches_presets() {
    let fa = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/toy/transcripts.fa");
    let txs = read_transcripts(&fa).unwrap();
    let index = Arc::new(
        Index::build(txs.clone(), &IndexOptions::builder().scheme("lcp:K=3").build())
            .unwrap()
            .0,
    );
    assert!(index.manifest().is_some());
    for p in Preset::ALL {
        let a = spearfish::MapSettings::builder(p).mapper(index.clone()).unwrap();
        let b = Mapper::new(index.clone(), &MapOptions::preset(p)).unwrap();
        assert_eq!(a.options(), b.options());
        for (_, _, seq) in &txs {
            assert_eq!(a.map(seq), b.map(seq));
        }
    }
    let m = spearfish::MapSettings::builder(Preset::Hifi)
        .extend_ends(true)
        .mapper(index.clone())
        .unwrap();
    assert!(m.options().params.ends.is_some());
    assert_eq!(m.score_denom(), spearfish::SCORE_DENOM);
}
