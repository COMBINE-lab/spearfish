//! `spearfish`: fast long-read RNA-seq mapping with locally consistent parsing.
//! Subcommands arrive as their milestones complete (see PLAN.md).

mod annot;
mod halo_cmd;
mod hist;
mod index_cmd;
mod lcp_stats;
mod locality_cmd;
mod map_cmd;
mod quant_cmd;
mod rng;
mod seeding;
mod sim;
mod simulate_cmd;
mod survival_cmd;

use clap::{Parser, Subcommand};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[derive(Parser)]
#[command(name = "spearfish", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

// Parsed once; the size difference between subcommands does not matter.
#[allow(clippy::large_enum_variant)]
#[derive(Subcommand)]
enum Command {
    /// Build the index from a transcript FASTA.
    Index(index_cmd::Args),
    /// Report on a built index.
    IndexStats(index_cmd::StatsArgs),
    /// Map long reads to the indexed transcriptome.
    Map(map_cmd::Args),
    /// Quantify transcripts by EM from a mapping TSV.
    Quant(quant_cmd::Args),
    /// Simulate long reads with known origin.
    Simulate(simulate_cmd::Args),
    /// Per-level LCP core statistics for a FASTA/FASTQ file.
    LcpStats(lcp_stats::Args),
    /// Measure LCP context radii by random single edits.
    Locality(locality_cmd::Args),
    /// Seed survival under simulated sequencing errors, per seed scheme (M2).
    Survival(survival_cmd::Args),
    /// Junction halos and isoform resolution per seed scheme (M4 analysis).
    Halo(halo_cmd::Args),
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        Command::Index(args) => index_cmd::run(args),
        Command::IndexStats(args) => index_cmd::stats(args),
        Command::Map(args) => map_cmd::run(args),
        Command::Quant(args) => quant_cmd::run(args),
        Command::Simulate(args) => simulate_cmd::run(args),
        Command::LcpStats(args) => lcp_stats::run(args),
        Command::Locality(args) => locality_cmd::run(args),
        Command::Survival(args) => survival_cmd::run(args),
        Command::Halo(args) => halo_cmd::run(args),
    }
}
