//! `spearfish`: fast long-read RNA-seq mapping with locally consistent parsing.
//! Subcommands arrive as their milestones complete (see PLAN.md).

mod hist;
mod lcp_stats;

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
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        Command::LcpStats(args) => lcp_stats::run(args),
    }
}
