//! Run provenance: what ran, on what, with which settings, and how it went.

use serde::Serialize;
use spearfish::provenance::BuildInfo;

/// Process resource usage so far.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct Usage {
    /// CPU seconds in user and system mode.
    pub cpu_user_seconds: f64,
    pub cpu_system_seconds: f64,
    /// Peak resident set size.
    pub peak_rss_bytes: u64,
}

impl Usage {
    pub fn now() -> Self {
        // SAFETY: getrusage writes the struct it is given and nothing else.
        let mut ru: libc::rusage = unsafe { std::mem::zeroed() };
        if unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut ru) } != 0 {
            return Usage::default();
        }
        let secs = |t: libc::timeval| t.tv_sec as f64 + t.tv_usec as f64 * 1e-6;
        // Linux reports kilobytes, macOS bytes.
        let rss_unit = if cfg!(target_os = "macos") { 1 } else { 1024 };
        Usage {
            cpu_user_seconds: secs(ru.ru_utime),
            cpu_system_seconds: secs(ru.ru_stime),
            peak_rss_bytes: ru.ru_maxrss as u64 * rss_unit,
        }
    }
}

/// What ran: the build, the command line, the host and the start time.
#[derive(Clone, Debug, Serialize)]
pub struct Invocation {
    pub spearfish: BuildInfo,
    /// SIMD kernels used.
    pub simd: &'static str,
    pub command_line: Vec<String>,
    pub working_directory: String,
    pub started_unix: u64,
    /// Worker threads used, and the CPUs available to the process.
    pub threads: usize,
    pub cpus_available: usize,
}

impl Invocation {
    pub fn new(threads: usize) -> Self {
        Invocation {
            spearfish: BuildInfo::current(),
            simd: spearfish::provenance::simd_level(),
            command_line: std::env::args().collect(),
            working_directory: std::env::current_dir()
                .map(|d| d.display().to_string())
                .unwrap_or_default(),
            started_unix: unix_now(),
            threads,
            cpus_available: available_cpus(),
        }
    }
}

pub fn unix_now() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs())
}

/// CPUs this process may use (respects affinity masks and cgroup quotas).
pub fn available_cpus() -> usize {
    std::thread::available_parallelism().map_or(1, |n| n.get())
}

/// `0` means every available CPU.
pub fn threads(requested: usize) -> usize {
    if requested == 0 {
        available_cpus()
    } else {
        requested
    }
}

/// Write `value` as pretty JSON to `path`.
pub fn write_json(path: &std::path::Path, value: &impl Serialize) -> anyhow::Result<()> {
    let mut s = serde_json::to_string_pretty(value)?;
    s.push('\n');
    std::fs::write(path, s).map_err(|e| anyhow::anyhow!("writing {}: {e}", path.display()))
}

/// The index part of a run record: where it is, what it is and what built it.
#[derive(Clone, Debug, Serialize)]
pub struct IndexInfo {
    pub path: String,
    pub index_id: String,
    pub spec: String,
    pub transcripts: usize,
    /// The index's manifest (absent for indexes written before manifests were recorded).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub manifest: Option<spearfish::provenance::IndexManifest>,
}

impl IndexInfo {
    pub fn of(path: &std::path::Path, idx: &spearfish::Index) -> Self {
        IndexInfo {
            path: path.display().to_string(),
            index_id: idx.id(),
            spec: idx.raw().meta("spec").unwrap_or("").to_string(),
            transcripts: idx.n_transcripts(),
            manifest: idx.manifest(),
        }
    }
}
