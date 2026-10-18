---
title: Using the library
description: Build an index, set up a mapper and map reads from Rust.
---

The `spearfish` crate is the library behind the command line.

```rust
use std::sync::Arc;
use spearfish::{Index, IndexOptions, MapSettings, Preset};

// Build an index in memory (or open one written by `spearfish index`).
let txs = spearfish::read_transcripts("transcripts.fa.gz".as_ref())?;
let opts = IndexOptions::builder().scheme("lcp:K=3").build();
let (index, report) = Index::build(txs, &opts)?;
let index = Arc::new(index);

// A preset, and the settings to change from it.
let mapper = MapSettings::builder(Preset::Hifi)
    .extend_ends(true)
    .tie_kappa(5.0)
    .mapper(index.clone())?;

let read = mapper.map(b"ACGT...");
for hit in &read.hits {
    println!("{} {} {}", index.tx_name(hit.tx), hit.score, hit.reverse);
}
```

- `MapSettings` holds a preset and one optional field per setting; unset fields keep the
  preset's value. `resolve()` gives the complete `MapOptions` (which serialize to JSON), and
  rejects settings that contradict each other.
- `Mapper` is `Send + Sync`: share one across threads.
- `Hit` coordinates are raw bases: query spans on the read's forward strand, transcript spans
  on the transcript. Score differences convert to probabilities with
  `exp(-d / spearfish::SCORE_DENOM)`.
- `Index::manifest()` and `Index::id()` give an index's provenance and identifier;
  `spearfish::provenance::BuildInfo::current()` the library's version and git revision.

The API documentation is `cargo doc -p spearfish --open`.
