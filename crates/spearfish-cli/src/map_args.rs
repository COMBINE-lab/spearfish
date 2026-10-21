//! Mapping settings on the command line: one flag per [`MapSettings`] field. Every flag is
//! optional; an unset flag keeps the preset's value (`spearfish map --print-params` shows
//! them all). On/off settings take `--x` / `--no-x`. `-h` lists the common flags, `--help`
//! all of them; experimental flags are hidden (see the documentation's developer notes).
//! Old flag names are kept as hidden aliases.

use clap::{Args, ValueEnum};
use spearfish::{ChainerKind, MapSettings, Preset, RescueWindow, Strands};

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
pub enum PresetArg {
    /// ONT cDNA (both orientations, minimap2 map-ont scores).
    #[default]
    Ont,
    /// ONT direct RNA (sense orientation only).
    Drna,
    /// PacBio HiFi (minimap2 map-hifi scores).
    Hifi,
}

impl From<PresetArg> for Preset {
    fn from(p: PresetArg) -> Self {
        match p {
            PresetArg::Ont => Preset::Ont,
            PresetArg::Drna => Preset::Drna,
            PresetArg::Hifi => Preset::Hifi,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum StrandsArg {
    /// Both orientations.
    Both,
    /// The read as given.
    Fwd,
    /// The read's reverse complement only.
    Rev,
}

impl From<StrandsArg> for Strands {
    fn from(s: StrandsArg) -> Self {
        match s {
            StrandsArg::Both => Strands::Both,
            StrandsArg::Fwd => Strands::Forward,
            StrandsArg::Rev => Strands::Reverse,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum ChainerArg {
    Oracle,
    Refine,
}

/// `Some(true)` for `--x`, `Some(false)` for `--no-x`, `None` for neither.
fn tri(on: bool, off: bool) -> Option<bool> {
    match (on, off) {
        (true, _) => Some(true),
        (_, true) => Some(false),
        _ => None,
    }
}

/// An on/off pair of flags.
macro_rules! onoff {
    ($self:ident . $on:ident / $off:ident) => {
        tri($self.$on, $self.$off)
    };
}

#[derive(Args, Debug, Default)]
pub struct MapArgs {
    /// Read type; sets the defaults of every other setting.
    #[arg(long, value_enum, default_value = "ont", help_heading = "Settings")]
    pub preset: PresetArg,
    /// Read orientations to map [default: drna fwd, others both].
    #[arg(long, value_enum, help_heading = "Settings")]
    pub strands: Option<StrandsArg>,
    /// Report transcripts scoring within this fraction of the read's best [default: 0.3].
    #[arg(long, alias = "delta-rep", value_name = "FRAC", help_heading = "Settings")]
    pub report_within: Option<f32>,
    /// Extend reported spans to where the read's alignment ends [default: on for ont].
    #[arg(long, overrides_with = "no_extend_ends", help_heading = "Settings")]
    pub extend_ends: bool,
    /// Report spans between the outermost seeds.
    #[arg(long, overrides_with = "extend_ends", help_heading = "Settings")]
    pub no_extend_ends: bool,
    /// Base-level tie-break among near-best transcripts [default: on].
    #[arg(long, overrides_with = "no_tiebreak", help_heading = "Settings")]
    pub tiebreak: bool,
    /// Turn the tie-break off (chain scores only).
    #[arg(long, overrides_with = "tiebreak", help_heading = "Settings")]
    pub no_tiebreak: bool,
    /// Re-seed poorly mapped reads at a lower LCP level [default: on except hifi].
    #[arg(long, overrides_with = "no_rescue", help_heading = "Settings")]
    pub rescue: bool,
    /// Turn the rescue off.
    #[arg(long, overrides_with = "rescue", help_heading = "Settings")]
    pub no_rescue: bool,
    /// Trim a poly(A) tail from a read's 3' end before aligning it [default: on].
    #[arg(long, overrides_with = "no_tail_trim", help_heading = "Settings")]
    pub tail_trim: bool,
    /// Align poly(A) tails like other read bases.
    #[arg(long, overrides_with = "tail_trim", help_heading = "Settings")]
    pub no_tail_trim: bool,

    // Reported hits.
    /// Extend only hits scoring within this fraction of the read's best [default: 0.05].
    #[arg(
        long,
        value_name = "FRAC",
        hide_short_help = true,
        help_heading = "Reported hits"
    )]
    pub extend_within: Option<f32>,
    /// Flag reads whose best hit spans less than this fraction of them [default: 0.5].
    #[arg(
        long,
        value_name = "FRAC",
        hide_short_help = true,
        help_heading = "Reported hits"
    )]
    pub min_coverage: Option<f32>,

    // Candidate genes.
    /// Candidate genes kept per read, best first [default: 8].
    #[arg(
        long,
        alias = "f-keep",
        hide_short_help = true,
        help_heading = "Candidate genes"
    )]
    pub keep_genes: Option<usize>,
    /// Also keep genes within this relative score of the last kept one [default: 0.01]...
    #[arg(
        long,
        alias = "f-tie",
        value_name = "FRAC",
        hide_short_help = true,
        help_heading = "Candidate genes"
    )]
    pub keep_genes_tie: Option<f32>,
    /// ...up to this many in total [default: 64].
    #[arg(
        long,
        alias = "f-max",
        hide_short_help = true,
        help_heading = "Candidate genes"
    )]
    pub max_genes: Option<usize>,
    /// Ignore seeds found in more than this many genes [default: 1000].
    #[arg(
        long,
        alias = "h-max",
        hide_short_help = true,
        help_heading = "Candidate genes"
    )]
    pub max_seed_genes: Option<usize>,
    /// Keep genes scoring at least this fraction of the best gene [default: 0.5].
    #[arg(
        long,
        value_name = "FRAC",
        hide_short_help = true,
        help_heading = "Candidate genes"
    )]
    pub min_gene_frac: Option<f32>,
    /// Cap on a seed's span in the gene score [default: 100].
    #[arg(long, hide_short_help = true, help_heading = "Candidate genes")]
    pub seed_span_cap: Option<u32>,
    #[arg(long, alias = "g-max", hide = true)]
    pub seed_gene_weight_max: Option<f32>,

    // Chaining.
    /// Chaining cost per base of gap discrepancy [default: 1.5].
    #[arg(long, alias = "alpha", hide_short_help = true, help_heading = "Chaining")]
    pub chain_alpha: Option<f32>,
    /// Predecessors considered per anchor [default: 50].
    #[arg(long, alias = "window", hide_short_help = true, help_heading = "Chaining")]
    pub chain_window: Option<usize>,
    /// Longest gap between chained anchors, in bases [default: 5000].
    #[arg(long, hide_short_help = true, help_heading = "Chaining")]
    pub chain_max_gap: Option<i64>,
    /// Largest overlap between chained anchors, in bases [default: 100].
    #[arg(long, hide_short_help = true, help_heading = "Chaining")]
    pub chain_max_overlap: Option<i64>,
    #[arg(long, alias = "beta", hide = true)]
    pub chain_beta: Option<f32>,
    #[arg(long = "no-merge-anchors", alias = "no-merge", hide = true)]
    pub no_merge_anchors: bool,
    #[arg(long, value_enum, hide = true)]
    pub chainer: Option<ChainerArg>,
    #[arg(long, hide = true, requires = "chainer")]
    pub beam: Option<usize>,

    // Tie-break.
    /// Absolute chain-score margin of the tie set [default: 30].
    #[arg(long, hide_short_help = true, help_heading = "Tie-break")]
    pub tie_margin: Option<f32>,
    /// Relative margin of the tie set; the larger margin applies [default: 0.05].
    #[arg(long, value_name = "FRAC", hide_short_help = true, help_heading = "Tie-break")]
    pub tie_margin_rel: Option<f32>,
    /// Align at most this many tied transcripts [default: 256].
    #[arg(long, hide_short_help = true, help_heading = "Tie-break")]
    pub tie_cands: Option<usize>,
    /// Score charged per edit [default: 6].
    #[arg(long, alias = "kappa", hide_short_help = true, help_heading = "Tie-break")]
    pub tie_kappa: Option<f32>,
    /// Score charged per affine-score point below the best [default: 1].
    #[arg(long, alias = "lambda", hide_short_help = true, help_heading = "Tie-break")]
    pub tie_lambda: Option<f32>,
    /// Affine-score only texts within this many edits of the best [default: drna 4, others
    /// all].
    #[arg(
        long,
        alias = "affine-within",
        hide_short_help = true,
        help_heading = "Tie-break"
    )]
    pub tie_affine_within: Option<u32>,
    /// Group tied hits without a shared anchor by the anchors they share [default: hifi].
    #[arg(
        long,
        overrides_with = "no_tie_group",
        hide_short_help = true,
        help_heading = "Tie-break"
    )]
    pub tie_group: bool,
    /// Align each such hit's whole read.
    #[arg(
        long,
        overrides_with = "tie_group",
        hide_short_help = true,
        help_heading = "Tie-break"
    )]
    pub no_tie_group: bool,
    /// Distances beyond the best candidate's plus this are not computed exactly [default: 16].
    #[arg(long, hide_short_help = true, help_heading = "Tie-break")]
    pub tie_cap: Option<u32>,
    #[arg(long, alias = "unit-tiebreak", hide = true)]
    pub tie_unit_cost: bool,
    #[arg(long, alias = "whole-read-tiebreak", hide = true)]
    pub tie_whole_read: bool,
    #[arg(long, hide = true)]
    pub no_tie_whole_trim: bool,
    #[arg(long, alias = "orient-edit", hide = true)]
    pub tie_orient_edit: bool,
    #[arg(long, alias = "flank-min", hide = true)]
    pub tie_flank_min: Option<i64>,

    // Read ends.
    /// Clip read ends that don't align (adapters, primers) [default: on].
    #[arg(
        long,
        overrides_with = "no_clip_ends",
        hide_short_help = true,
        help_heading = "Read ends"
    )]
    pub clip_ends: bool,
    /// Force whole read ends into the alignment.
    #[arg(
        long,
        alias = "no-tie-clip",
        overrides_with = "clip_ends",
        hide_short_help = true,
        help_heading = "Read ends"
    )]
    pub no_clip_ends: bool,
    /// Bonus for aligning through to a read's end [default: 5].
    #[arg(long, hide_short_help = true, help_heading = "Read ends")]
    pub end_bonus: Option<i32>,
    /// X-drop of read-end alignments [default: 50].
    #[arg(long, hide_short_help = true, help_heading = "Read ends")]
    pub xdrop: Option<i32>,
    /// Band of the affine alignments, in cells off the diagonals [default: 16].
    #[arg(long, alias = "tie-band", hide_short_help = true, help_heading = "Read ends")]
    pub band: Option<usize>,
    /// Transcript bases searched beyond a read end's projection [default: 20]...
    #[arg(long, alias = "tie-pad", hide_short_help = true, help_heading = "Read ends")]
    pub pad: Option<u32>,
    /// ...plus 1/DIV of the unanchored read end [default: 10].
    #[arg(
        long,
        alias = "tie-pad-div",
        value_name = "DIV",
        hide_short_help = true,
        help_heading = "Read ends"
    )]
    pub pad_div: Option<u32>,
    /// Keep read-end bands on the last shared anchor's diagonal [default: hifi].
    #[arg(
        long,
        overrides_with = "no_end_diagonal",
        hide_short_help = true,
        help_heading = "Read ends"
    )]
    pub end_diagonal: bool,
    /// Widen read-end bands over the window's slack.
    #[arg(
        long,
        overrides_with = "end_diagonal",
        hide_short_help = true,
        help_heading = "Read ends"
    )]
    pub no_end_diagonal: bool,
    /// Align read ends in raw bases (needs raw transcripts in the index) [default: on except
    /// drna].
    #[arg(
        long,
        overrides_with = "no_raw_ends",
        hide_short_help = true,
        help_heading = "Read ends"
    )]
    pub raw_ends: bool,
    /// Align read ends in homopolymer-compressed bases.
    #[arg(
        long,
        overrides_with = "raw_ends",
        hide_short_help = true,
        help_heading = "Read ends"
    )]
    pub no_raw_ends: bool,
    /// Shortest poly(A) tail recognized, in bases [default: drna 10, others 15].
    #[arg(long, hide_short_help = true, help_heading = "Read ends")]
    pub tail_run: Option<usize>,
    /// Points for isoforms ending where the read's poly(A) tail starts (0: off) [default: 10].
    #[arg(long, hide_short_help = true, help_heading = "Read ends")]
    pub tail_bonus: Option<i32>,
    /// Near-tied read orientations: points for the one whose read ends in a poly(A) tail (0: off)
    /// [default: 30].
    #[arg(long, hide_short_help = true, help_heading = "Read ends")]
    pub orient_tail: Option<i32>,
    /// Tie-break both read orientations when the other's best chain is within this fraction
    /// of the best (or the tie margin) [default: 0.35].
    #[arg(long, value_name = "FRAC", hide_short_help = true, help_heading = "Tie-break")]
    pub orient_margin_rel: Option<f32>,

    // Rescue.
    /// Rescue seed scheme [default: lcp:K=2].
    #[arg(long, hide_short_help = true, help_heading = "Rescue")]
    pub rescue_scheme: Option<String>,
    /// Rescue candidate genes per strand [default: 16].
    #[arg(long, hide_short_help = true, help_heading = "Rescue")]
    pub rescue_genes: Option<usize>,
    /// Reads shorter than this (parsed bases) get proportionally more candidate genes, up to 4×
    /// (0: off) [default: 350].
    #[arg(long, value_name = "LEN", hide_short_help = true, help_heading = "Rescue")]
    pub rescue_genes_len: Option<u32>,
    /// Transcripts re-seeded per gene, those with the most anchors [default: 128].
    #[arg(long, hide_short_help = true, help_heading = "Rescue")]
    pub rescue_tx: Option<usize>,
    /// Candidate genes must score this fraction of the best gene [default: 0.5].
    #[arg(long, value_name = "FRAC", hide_short_help = true, help_heading = "Rescue")]
    pub rescue_alpha: Option<f32>,
    /// What of each transcript to re-seed: PAD,DIV (PAD + read/DIV bases around the read's
    /// projection) or `full` [default: 10,50].
    #[arg(
        long,
        value_name = "PAD,DIV|full",
        hide_short_help = true,
        help_heading = "Rescue"
    )]
    pub rescue_window: Option<RescueWindow>,
    /// Re-seed the other orientation only if its score is within this margin of the chosen
    /// one's (negative: never) [default: 30].
    #[arg(
        long,
        allow_hyphen_values = true,
        conflicts_with = "rescue_other_always",
        hide_short_help = true,
        help_heading = "Rescue"
    )]
    pub rescue_other_margin: Option<f32>,
    /// Re-seed the other orientation regardless of its score.
    #[arg(long, hide_short_help = true, help_heading = "Rescue")]
    pub rescue_other_always: bool,
    /// Skip transcripts longer than this, in compressed bases [default: 50000].
    #[arg(long, hide_short_help = true, help_heading = "Rescue")]
    pub rescue_max_tx_len: Option<u32>,
}

impl MapArgs {
    /// The settings given on the command line.
    pub fn settings(&self) -> MapSettings {
        let chainer = self.chainer.map(|c| match c {
            ChainerArg::Oracle => ChainerKind::Oracle,
            ChainerArg::Refine => ChainerKind::Refine { beam: self.beam },
        });
        MapSettings::builder(self.preset.into())
            .maybe_strands(self.strands.map(Into::into))
            .maybe_report_within(self.report_within)
            .maybe_min_coverage(self.min_coverage)
            .maybe_extend_ends(onoff!(self.extend_ends / no_extend_ends))
            .maybe_extend_within(self.extend_within)
            .maybe_keep_genes(self.keep_genes)
            .maybe_keep_genes_tie(self.keep_genes_tie)
            .maybe_max_genes(self.max_genes)
            .maybe_max_seed_genes(self.max_seed_genes)
            .maybe_min_gene_frac(self.min_gene_frac)
            .maybe_seed_span_cap(self.seed_span_cap)
            .maybe_seed_gene_weight_max(self.seed_gene_weight_max)
            .maybe_chain_alpha(self.chain_alpha)
            .maybe_chain_window(self.chain_window)
            .maybe_chain_max_gap(self.chain_max_gap)
            .maybe_chain_max_overlap(self.chain_max_overlap)
            .maybe_chain_beta(self.chain_beta)
            .maybe_merge_anchors(self.no_merge_anchors.then_some(false))
            .maybe_chainer(chainer)
            .maybe_tiebreak(onoff!(self.tiebreak / no_tiebreak))
            .maybe_tie_margin(self.tie_margin)
            .maybe_tie_margin_rel(self.tie_margin_rel)
            .maybe_tie_cands(self.tie_cands)
            .maybe_tie_kappa(self.tie_kappa)
            .maybe_tie_lambda(self.tie_lambda)
            .maybe_tie_affine_within(self.tie_affine_within)
            .maybe_tie_group(onoff!(self.tie_group / no_tie_group))
            .maybe_tie_cap(self.tie_cap)
            .maybe_tie_unit_cost(self.tie_unit_cost.then_some(true))
            .maybe_tie_whole_read(self.tie_whole_read.then_some(true))
            .maybe_tie_whole_trim(self.no_tie_whole_trim.then_some(false))
            .maybe_tie_orient_edit(self.tie_orient_edit.then_some(true))
            .maybe_tie_flank_min(self.tie_flank_min)
            .maybe_clip_ends(onoff!(self.clip_ends / no_clip_ends))
            .maybe_end_bonus(self.end_bonus)
            .maybe_xdrop(self.xdrop)
            .maybe_band(self.band)
            .maybe_pad(self.pad)
            .maybe_pad_div(self.pad_div)
            .maybe_end_diagonal(onoff!(self.end_diagonal / no_end_diagonal))
            .maybe_raw_ends(onoff!(self.raw_ends / no_raw_ends))
            .maybe_tail_trim(onoff!(self.tail_trim / no_tail_trim))
            .maybe_tail_run(self.tail_run)
            .maybe_tail_bonus(self.tail_bonus)
            .maybe_orient_tail(self.orient_tail)
            .maybe_orient_margin_rel(self.orient_margin_rel)
            .maybe_rescue(onoff!(self.rescue / no_rescue))
            .maybe_rescue_scheme(self.rescue_scheme.clone())
            .maybe_rescue_genes(self.rescue_genes)
            .maybe_rescue_genes_len(self.rescue_genes_len)
            .maybe_rescue_tx(self.rescue_tx)
            .maybe_rescue_alpha(self.rescue_alpha)
            .maybe_rescue_window(self.rescue_window)
            .maybe_rescue_other_margin(self.rescue_other_margin)
            .maybe_rescue_other_always(self.rescue_other_always.then_some(true))
            .maybe_rescue_max_tx_len(self.rescue_max_tx_len)
            .build()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[derive(Parser)]
    struct Cli {
        #[command(flatten)]
        m: MapArgs,
    }

    fn settings(args: &[&str]) -> MapSettings {
        Cli::try_parse_from(std::iter::once("x").chain(args.iter().copied()))
            .unwrap()
            .m
            .settings()
    }

    /// The defaults quoted in the help text are the library's: giving a flag its documented
    /// default (for ont) leaves the options unchanged.
    #[test]
    fn documented_defaults_are_the_defaults() {
        use clap::CommandFactory;
        let ont = spearfish::MapOptions::preset(Preset::Ont);
        let cmd = Cli::command();
        let mut checked = 0;
        for arg in cmd.get_arguments() {
            let Some(long) = arg.get_long() else { continue };
            let help = arg
                .get_long_help()
                .or(arg.get_help())
                .map(|h| h.to_string())
                .unwrap_or_default()
                .replace('\n', " ");
            let Some(d) = help.split("[default: ").nth(1).and_then(|r| r.split(']').next()) else {
                continue;
            };
            let flag = if d == "on" || d.starts_with("on for ont") || d.starts_with("on except") {
                format!("--{long}")
            } else if d == "hifi" {
                format!("--no-{long}")
            } else if d.contains(' ') {
                // Preset-dependent (`drna 4, others all`): checked by preset tests.
                continue;
            } else {
                format!("--{long}={d}")
            };
            let s = Cli::try_parse_from(["x", flag.as_str()])
                .unwrap_or_else(|e| panic!("{flag}: {e}"))
                .m
                .settings();
            assert_eq!(s.resolve().unwrap(), ont, "{flag} is not the default");
            checked += 1;
        }
        assert!(checked > 30, "only {checked} documented defaults found");
    }

    #[test]
    fn unset_flags_change_nothing() {
        assert_eq!(settings(&[]), MapSettings::preset(Preset::Ont));
        assert_eq!(settings(&["--preset", "hifi"]), MapSettings::preset(Preset::Hifi));
    }

    /// Every setting is reachable from the command line: with every flag given, no field of
    /// the settings is left unset.
    #[test]
    fn every_setting_has_a_flag() {
        let s = settings(&[
            "--strands=fwd",
            "--report-within=0.2",
            "--min-coverage=0.4",
            "--no-extend-ends",
            "--extend-within=0.1",
            "--keep-genes=4",
            "--keep-genes-tie=0.02",
            "--max-genes=32",
            "--max-seed-genes=500",
            "--min-gene-frac=0.4",
            "--seed-span-cap=50",
            "--seed-gene-weight-max=100",
            "--chain-alpha=1",
            "--chain-window=40",
            "--chain-max-gap=4000",
            "--chain-max-overlap=50",
            "--chain-beta=0.1",
            "--no-merge-anchors",
            "--chainer=refine",
            "--beam=8",
            "--tiebreak",
            "--tie-margin=20",
            "--tie-margin-rel=0.1",
            "--tie-cands=100",
            "--tie-kappa=5",
            "--tie-lambda=2",
            "--tie-affine-within=3",
            "--no-tie-group",
            "--tie-cap=8",
            "--tie-unit-cost",
            "--tie-whole-read",
            "--no-tie-whole-trim",
            "--tie-orient-edit",
            "--tie-flank-min=10",
            "--no-clip-ends",
            "--end-bonus=3",
            "--xdrop=30",
            "--band=8",
            "--pad=10",
            "--pad-div=5",
            "--end-diagonal",
            "--no-raw-ends",
            "--no-tail-trim",
            "--tail-run=12",
            "--tail-bonus=0",
            "--orient-tail=0",
            "--orient-margin-rel=0.2",
            "--rescue",
            "--rescue-scheme=lcp:K=2",
            "--rescue-genes=8",
            "--rescue-genes-len=0",
            "--rescue-tx=64",
            "--rescue-alpha=0.4",
            "--rescue-window=full",
            "--rescue-other-margin=-1",
            "--rescue-max-tx-len=1000",
        ]);
        let j = serde_json::to_value(&s).unwrap();
        let missing: Vec<_> = MapSettings::FIELDS
            .iter()
            // rescue_other_always conflicts with rescue_other_margin; its flag is tested below.
            .filter(|f| j.get(**f).is_none() && **f != "rescue_other_always")
            .collect();
        assert!(missing.is_empty(), "settings without a flag: {missing:?}");
        assert_eq!(
            settings(&["--rescue-other-always"]).rescue_other_always,
            Some(true)
        );
        assert!(s.resolve().is_ok());
    }
}
