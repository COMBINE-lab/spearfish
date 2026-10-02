//! `spearfish`: fast long-read RNA-seq mapping with locally consistent parsing.
//! Subcommands arrive as their milestones complete (see PLAN.md).

mod annot;
mod halo_cmd;
mod hist;
mod lcp_stats;
mod locality_cmd;
mod rng;
mod seeding;
mod sim;
mod survival_cmd;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "spearfish", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
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
        Command::LcpStats(args) => lcp_stats::run(args),
        Command::Locality(args) => locality_cmd::run(args),
        Command::Survival(args) => survival_cmd::run(args),
        Command::Halo(args) => halo_cmd::run(args),
    }
}
