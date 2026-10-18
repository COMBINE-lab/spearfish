//! `spearfish quant`: reads to transcript abundances, by mapping to a RAD file and
//! quantifying it with oarfish (`oarfish -a`), which takes the score scale, preset and span
//! settings from the file. Writes oarfish's outputs under the prefix, the mapping's run record
//! as `<prefix>.map.meta.json`, and oarfish's log as `<prefix>.oarfish.log` (on stderr with
//! `-v`).

use std::io::IsTerminal;
use std::path::PathBuf;
use std::process::Command;

use anyhow::{Context, Result, bail};

use crate::map_args::MapArgs;
use crate::map_cmd::{Codec, Job, execute};

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Index built by `spearfish index`.
    index: PathBuf,
    /// Reads (FASTA/FASTQ, optionally compressed); `-` reads stdin.
    reads: PathBuf,
    /// Output prefix (oarfish writes <PREFIX>.quant, <PREFIX>.meta_info.json, ...).
    #[arg(short, long, value_name = "PREFIX", help_heading = "Input and output")]
    output: PathBuf,
    /// Worker threads for mapping and quantification (0: every available CPU).
    #[arg(short, long, default_value_t = 0, help_heading = "Input and output")]
    threads: usize,
    /// The oarfish binary [default: $OARFISH, else `oarfish` on the PATH].
    #[arg(long, value_name = "PATH", help_heading = "Input and output")]
    oarfish: Option<PathBuf>,
    /// Keep the intermediate RAD file (<PREFIX>.rad).
    #[arg(long, help_heading = "Input and output")]
    keep_rad: bool,
    /// Further oarfish arguments, after `--` (e.g. `-- --model-truncation`).
    #[arg(last = true, value_name = "OARFISH ARGS")]
    oarfish_args: Vec<String>,
    #[command(flatten)]
    settings: MapArgs,
}

/// The oarfish to run, checked to read spearfish's RAD files.
fn find_oarfish(given: Option<PathBuf>) -> Result<PathBuf> {
    let bin = given
        .or_else(|| std::env::var_os("OARFISH").map(PathBuf::from))
        .unwrap_or_else(|| "oarfish".into());
    let out = Command::new(&bin).arg("--help").output().with_context(|| {
        format!(
            "running {} (install oarfish, or give its path with --oarfish or $OARFISH)",
            bin.display()
        )
    })?;
    let help = String::from_utf8_lossy(&out.stdout);
    if !help.contains("spearfish") {
        bail!(
            "{} cannot read spearfish's RAD files: an oarfish with the spearfish backend is needed",
            bin.display()
        );
    }
    let v = Command::new(&bin).arg("--version").output()?;
    log::debug!(
        "using {} ({})",
        bin.display(),
        String::from_utf8_lossy(&v.stdout).trim()
    );
    Ok(bin)
}

fn with_suffix(prefix: &std::path::Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(suffix);
    s.into()
}

pub fn run(args: Args) -> Result<()> {
    let oarfish = find_oarfish(args.oarfish.clone())?;
    if let Some(dir) = args.output.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir)?;
    }
    let threads = crate::runinfo::threads(args.threads);
    let rad = with_suffix(&args.output, ".rad");
    execute(Job {
        index: args.index,
        reads: args.reads,
        threads,
        output: None,
        paf: None,
        rad: Some(rad.clone()),
        rad_codec: Codec::Lz4,
        meta: Some(with_suffix(&args.output, ".map.meta.json")),
        batch_size: 32,
        settings: args.settings.settings(),
    })?;
    let mut cmd = Command::new(&oarfish);
    cmd.arg("-a")
        .arg(&rad)
        .arg("-o")
        .arg(&args.output)
        .arg("-j")
        .arg(threads.max(2).to_string());
    if log::max_level() < log::LevelFilter::Info {
        cmd.arg("--quiet");
    }
    cmd.args(&args.oarfish_args);
    // oarfish's log goes to a file unless -v asks for it on stderr; files get no colour codes.
    let log_path = with_suffix(&args.output, ".oarfish.log");
    let to_file = log::max_level() < log::LevelFilter::Debug;
    if to_file {
        cmd.stderr(std::fs::File::create(&log_path)?);
    }
    if to_file || !std::io::stderr().is_terminal() {
        cmd.env("NO_COLOR", "1");
    }
    log::debug!("running {cmd:?}");
    let status = cmd
        .status()
        .with_context(|| format!("running {}", oarfish.display()))?;
    if !args.keep_rad {
        let _ = std::fs::remove_file(&rad);
    }
    if !status.success() {
        if to_file {
            let text = std::fs::read_to_string(&log_path).unwrap_or_default();
            let tail: Vec<&str> = text.lines().rev().take(20).collect();
            for l in tail.into_iter().rev() {
                eprintln!("  oarfish: {l}");
            }
        }
        bail!("oarfish failed ({status}); its log is {}", log_path.display());
    }
    log::info!(
        "quant: estimates in {}",
        with_suffix(&args.output, ".quant").display()
    );
    Ok(())
}
