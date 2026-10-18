//! Provenance: what built this library, what built an index, and identifiers that tie
//! outputs to the index they came from.

use serde::{Deserialize, Serialize};

/// `version (git revision)` of this build.
pub const VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), " (", env!("SPEARFISH_GIT"), ")");

/// The library build: version, git revision, target and SIMD level.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct BuildInfo {
    /// Crate version.
    pub version: String,
    /// `git describe` of the source tree (`unknown` outside a checkout).
    pub git: String,
    /// Target triple.
    pub target: String,
    /// Cargo profile (`release`, `debug`).
    pub profile: String,
}

impl BuildInfo {
    /// This build.
    pub fn current() -> Self {
        BuildInfo {
            version: env!("CARGO_PKG_VERSION").into(),
            git: env!("SPEARFISH_GIT").into(),
            target: env!("SPEARFISH_TARGET").into(),
            profile: env!("SPEARFISH_PROFILE").into(),
        }
    }

    /// `version (git)`, e.g. for `--version` and output headers.
    pub fn short(&self) -> String {
        format!("{} ({})", self.version, self.git)
    }
}

/// The input transcripts of an index.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub struct InputInfo {
    /// Transcript FASTA as given (empty when built in memory).
    pub path: String,
    /// Transcripts and genes.
    pub transcripts: usize,
    /// Genes (given, from a table, or inferred).
    pub genes: usize,
    /// xxh3-64 over every transcript's name and sequence, in order.
    pub xxh3: String,
    /// Where the genes came from: `fasta` (headers), `t2g` (a table), `inferred`, or `none`
    /// (every transcript its own gene).
    pub genes_from: String,
    /// Transcript-to-gene table, if one was given, and its xxh3-64.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub t2g: Option<(String, String)>,
}

/// What built an index and from what, stored in the index (section `manifest`) and written
/// next to it as `<index>.json`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct IndexManifest {
    /// Always `spearfish-index`.
    pub kind: String,
    /// Container format version.
    pub format_version: u32,
    /// Identifier of the contents (see [`index_id`]).
    pub index_id: String,
    /// The build that wrote the index.
    pub built_by: BuildInfo,
    /// Seconds since the Unix epoch when the build finished.
    pub created_unix: u64,
    /// Command line, if built by the `spearfish` binary.
    pub command_line: Vec<String>,
    /// Build options.
    pub options: crate::IndexOptions,
    /// LCP parser version and hash seed the index was parsed with.
    pub lcp_parser_version: u32,
    /// LCP hash seed.
    pub lcp_seed: String,
    /// The input.
    pub input: InputInfo,
    /// Construction statistics (nodes, unitigs, families, ...).
    pub stats: serde_json::Map<String, serde_json::Value>,
}

impl Default for IndexManifest {
    fn default() -> Self {
        IndexManifest {
            kind: "spearfish-index".into(),
            format_version: index::store::FORMAT_VERSION,
            index_id: String::new(),
            built_by: BuildInfo::current(),
            created_unix: 0,
            command_line: Vec::new(),
            options: Default::default(),
            lcp_parser_version: lcp_core::LCP_PARSER_VERSION,
            lcp_seed: format!("{:#x}", lcp_core::LcpConfig::default().seed),
            input: Default::default(),
            stats: Default::default(),
        }
    }
}

impl IndexManifest {
    /// The manifest of a loaded index, if it has one.
    pub fn of(idx: &index::Index) -> Option<Self> {
        serde_json::from_str(idx.manifest()?).ok()
    }

    /// The manifest of a freshly built index.
    pub fn new(
        data: &index::IndexData,
        opts: &crate::IndexOptions,
        rep: &index::BuildReport,
        mut input: InputInfo,
    ) -> Self {
        let created_unix = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_secs());
        let meta = |k: &str| {
            data.meta
                .iter()
                .find(|(key, _)| key == k)
                .map(|(_, v)| v.clone())
                .unwrap_or_default()
        };
        input.transcripts = rep.transcripts;
        input.genes = rep.genes;
        input.xxh3 = meta("input_xxh3");
        let mut stats = serde_json::Map::new();
        for (k, v) in [
            ("unseeded_transcripts", rep.unseeded_transcripts),
            ("seed_occurrences", rep.seed_occurrences as usize),
            ("nodes", rep.nodes),
            ("unitigs", rep.unitigs),
            ("cyclic_unitigs", rep.cyclic_unitigs),
            ("occurrences", rep.occurrences),
            ("colors", rep.colors),
            ("masked_unitigs", rep.masked_unitigs),
            ("families", rep.families),
            ("largest_family_tx", rep.largest_family_tx),
            ("multi_gene_families", rep.multi_gene_families),
            ("label_collisions", rep.label_collisions),
        ] {
            stats.insert(k.into(), (v as u64).into());
        }
        let secs: serde_json::Map<String, serde_json::Value> = rep
            .seconds
            .iter()
            .map(|(k, s)| (k.to_string(), (*s).into()))
            .collect();
        stats.insert("seconds".into(), secs.into());
        IndexManifest {
            index_id: id_of_meta(&data.meta_text()),
            created_unix,
            options: opts.clone(),
            input,
            stats,
            ..Default::default()
        }
    }
}

/// A short identifier of an index's contents: xxh3-64 of its `meta` section, which holds the
/// input transcripts' hash and every build parameter. Rebuilding from the same inputs with the
/// same parameters gives the same identifier; outputs record it to name their index.
pub fn index_id(idx: &index::Index) -> String {
    id_of_meta(idx.meta_text())
}

fn id_of_meta(meta: &str) -> String {
    format!("{:016x}", xxhash_rust::xxh3::xxh3_64(meta.as_bytes()))
}

/// Where an index's manifest is written: `<index>.json`.
pub fn sidecar(index: &std::path::Path) -> std::path::PathBuf {
    let mut p = index.as_os_str().to_owned();
    p.push(".json");
    p.into()
}

/// xxh3-64 of a file's bytes, as 16 hex digits.
pub fn file_xxh3(path: &std::path::Path) -> std::io::Result<String> {
    use std::io::Read;
    let mut f = std::fs::File::open(path)?;
    let mut h = xxhash_rust::xxh3::Xxh3::new();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    Ok(format!("{:016x}", h.digest()))
}

/// The SIMD kernels this process runs, as a short name.
pub fn simd_level() -> &'static str {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx512bw") && std::arch::is_x86_feature_detected!("avx512vl")
        {
            return "avx512";
        }
        if std::arch::is_x86_feature_detected!("avx2") {
            return "avx2";
        }
    }
    "scalar"
}
