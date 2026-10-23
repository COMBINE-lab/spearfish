//! `spearfish`: fast long-read RNA-seq mapping with locally consistent parsing.

mod annot;
mod doctor_cmd;
mod em;
mod em_cmd;
mod halo_cmd;
mod hist;
mod index_cmd;
mod lcp_stats;
mod locality_cmd;
mod logging;
mod map_args;
mod map_cmd;
mod rng;
mod runinfo;
mod seeding;
mod sim;
mod simulate_cmd;
mod survival_cmd;

use std::path::PathBuf;

use clap::{CommandFactory, Parser, Subcommand};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// `--version`: the version and git revision; `--version` (long) also the target and SIMD.
fn long_version() -> &'static str {
    let b = spearfish::provenance::BuildInfo::current();
    Box::leak(
        format!(
            "{}\ngit {}\ntarget {} ({})\nsimd {}",
            b.version,
            b.git,
            b.target,
            b.profile,
            spearfish::provenance::simd_level()
        )
        .into_boxed_str(),
    )
}

#[derive(Parser)]
#[command(
    name = "spearfish",
    version = spearfish::provenance::VERSION,
    long_version = long_version(),
    about,
    propagate_version = true
)]
struct Cli {
    /// Only report errors.
    #[arg(short, long, global = true, conflicts_with = "verbose")]
    quiet: bool,
    /// More detail on stderr (repeat for more).
    #[arg(short, long, global = true, action = clap::ArgAction::Count)]
    verbose: u8,
    #[command(subcommand)]
    command: Command,
}

// Parsed once; the size difference between subcommands does not matter.
#[allow(clippy::large_enum_variant)]
#[derive(Subcommand)]
enum Command {
    /// Build the index from a transcript FASTA.
    Index(index_cmd::Args),
    /// Map long reads to the indexed transcriptome (TSV, PAF and/or RAD output).
    Map(map_cmd::Args),
    /// Show what an index holds and what built it.
    #[command(alias = "index-stats")]
    Inspect(index_cmd::InspectArgs),
    /// Check this build and machine (SIMD self-test) and, optionally, an index.
    Doctor(doctor_cmd::Args),
    /// Print a shell completion script.
    Completions {
        /// Shell.
        shell: clap_complete::Shell,
    },
    /// Write man pages for spearfish and its commands.
    Man {
        /// Directory to write them to.
        #[arg(short, long, default_value = ".")]
        output: PathBuf,
    },
    /// spearfish's own EM over a mapping TSV (a baseline for method development).
    #[command(hide = true)]
    Em(em_cmd::Args),
    /// Simulate long reads with known origin.
    #[command(hide = true)]
    Simulate(simulate_cmd::Args),
    /// Per-level LCP core statistics for a FASTA/FASTQ file.
    #[command(hide = true)]
    LcpStats(lcp_stats::Args),
    /// Measure LCP context radii by random single edits.
    #[command(hide = true)]
    Locality(locality_cmd::Args),
    /// Seed survival under simulated sequencing errors, per seed scheme.
    #[command(hide = true)]
    Survival(survival_cmd::Args),
    /// Junction halos and isoform resolution per seed scheme.
    #[command(hide = true)]
    Halo(halo_cmd::Args),
}

/// Man pages for the command and each visible subcommand.
fn write_man(dir: &std::path::Path) -> anyhow::Result<()> {
    std::fs::create_dir_all(dir)?;
    let cmd = Cli::command();
    let write = |name: String, c: clap::Command| -> anyhow::Result<()> {
        let path = dir.join(format!("{name}.1"));
        let mut f = std::fs::File::create(&path)?;
        clap_mangen::Man::new(c).render(&mut f)?;
        log::info!("wrote {}", path.display());
        Ok(())
    };
    write("spearfish".into(), cmd.clone())?;
    for sub in cmd.get_subcommands().filter(|s| !s.is_hide_set()) {
        let name = format!("spearfish-{}", sub.get_name());
        let leaked: &'static str = Box::leak(name.clone().into_boxed_str());
        write(name, sub.clone().name(leaked))?;
    }
    Ok(())
}

fn main() -> std::process::ExitCode {
    let cli = Cli::parse();
    logging::init(cli.quiet, cli.verbose);
    match run(cli.command) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        // Output piped into a reader that stopped early (`| head`): not an error.
        Err(e)
            if e.chain().any(|c| {
                c.downcast_ref::<std::io::Error>()
                    .is_some_and(|io| io.kind() == std::io::ErrorKind::BrokenPipe)
            }) =>
        {
            std::process::ExitCode::SUCCESS
        }
        Err(e) => {
            log::error!("{e:#}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(command: Command) -> anyhow::Result<()> {
    match command {
        Command::Index(args) => index_cmd::run(args),
        Command::Map(args) => map_cmd::run(args),
        Command::Inspect(args) => index_cmd::inspect(args),
        Command::Doctor(args) => doctor_cmd::run(args),
        Command::Completions { shell } => {
            clap_complete::generate(shell, &mut Cli::command(), "spearfish", &mut std::io::stdout());
            Ok(())
        }
        Command::Man { output } => write_man(&output),
        Command::Em(args) => em_cmd::run(args),
        Command::Simulate(args) => simulate_cmd::run(args),
        Command::LcpStats(args) => lcp_stats::run(args),
        Command::Locality(args) => locality_cmd::run(args),
        Command::Survival(args) => survival_cmd::run(args),
        Command::Halo(args) => halo_cmd::run(args),
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn cli_is_consistent() {
        use clap::CommandFactory;
        super::Cli::command().debug_assert();
    }
}
