//! Mapping settings: a preset plus the settings the caller chose to change.
//!
//! [`MapSettings`] holds a [`Preset`] and one optional field per setting. A field left unset
//! takes the preset's value, so a preset's choices are never overridden by accident, and
//! [`MapSettings::resolve`] turns the settings into the complete [`MapOptions`] a mapper runs
//! with. Build settings with the builder:
//!
//! ```no_run
//! # use std::sync::Arc;
//! use spearfish::{Index, MapSettings, Preset};
//!
//! let index = Arc::new(Index::open("transcripts.sfi".as_ref())?);
//! let mapper = MapSettings::builder(Preset::Hifi)
//!     .extend_ends(true)
//!     .tie_kappa(5.0)
//!     .mapper(index)?;
//! # Ok::<(), anyhow::Error>(())
//! ```
//!
//! The settings serialize to JSON with unset fields left out; [`MapOptions`] serializes to
//! every value used, which `spearfish map` records next to its output.

use std::sync::Arc;

use anyhow::{Result, bail, ensure};
use serde::{Deserialize, Serialize};

use crate::{ChainerKind, Index, MapOptions, Mapper, Preset};
use map::{RescueParams, Strands};

/// Where the rescue re-seeds a candidate transcript.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum RescueWindow {
    /// The window the read can occupy: `pad` + 1/`div` of the read beyond its projection.
    Window {
        /// Bases beyond the read's projection.
        pad: u32,
        /// Plus 1/`div` of the read's length.
        div: u32,
    },
    /// The whole transcript.
    Full,
}

impl std::str::FromStr for RescueWindow {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self> {
        if s == "full" {
            return Ok(RescueWindow::Full);
        }
        let (pad, div) = s
            .split_once(',')
            .ok_or_else(|| anyhow::anyhow!("rescue window {s:?}: expected PAD,DIV or full"))?;
        let (pad, div) = (pad.trim().parse()?, div.trim().parse()?);
        ensure!(div > 0, "rescue window {s:?}: DIV must be positive");
        Ok(RescueWindow::Window { pad, div })
    }
}

/// Mapping settings: a preset and the settings changed from it (`None`: the preset's value).
///
/// Settings marked *experimental* exist for method development; their defaults are what was
/// evaluated.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, bon::Builder)]
#[serde(default, deny_unknown_fields)]
#[non_exhaustive]
pub struct MapSettings {
    /// Read type; sets the defaults of every other setting.
    #[builder(start_fn)]
    pub preset: Preset,
    /// Read orientations to map.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub strands: Option<Strands>,

    // Reported hits.
    /// Report transcripts scoring within this fraction of the read's best.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub report_within: Option<f32>,
    /// Flag a read `low_coverage` when its best hit spans less than this fraction of it.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_coverage: Option<f32>,
    /// Extend reported spans past the outermost seeds to where the read's alignment ends.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub extend_ends: Option<bool>,
    /// Extend only hits scoring within this fraction of the read's best.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub extend_within: Option<f32>,

    // Candidate genes.
    /// Candidate genes kept per read (best first).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub keep_genes: Option<usize>,
    /// Also keep genes within this relative score of the last kept one...
    #[serde(skip_serializing_if = "Option::is_none")]
    pub keep_genes_tie: Option<f32>,
    /// ...up to this many in total.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_genes: Option<usize>,
    /// Ignore seeds found in more than this many genes.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_seed_genes: Option<usize>,
    /// Keep genes scoring at least this fraction of the best gene.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_gene_frac: Option<f32>,
    /// Cap on a seed's span in the gene score.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub seed_span_cap: Option<u32>,
    /// Experimental: gene count at which a seed's weight reaches zero.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub seed_gene_weight_max: Option<f32>,

    // Chaining.
    /// Chaining cost per base of gap discrepancy |Δq − Δt|.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chain_alpha: Option<f32>,
    /// Chaining predecessor window (anchors).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chain_window: Option<usize>,
    /// Longest gap between chained anchors (bases).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chain_max_gap: Option<i64>,
    /// Largest overlap between chained anchors (bases).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chain_max_overlap: Option<i64>,
    /// Experimental: chaining cost per log2 of gap length.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chain_beta: Option<f32>,
    /// Experimental: merge collinear anchors within a unitig before chaining.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub merge_anchors: Option<bool>,
    /// Experimental: the chainer (both give identical scores without a beam).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chainer: Option<ChainerKind>,

    // Tie-break.
    /// Re-rank near-best transcripts by base-level alignment where they differ.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tiebreak: Option<bool>,
    /// Absolute chain-score margin defining the tie set.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tie_margin: Option<f32>,
    /// Relative margin defining the tie set (the larger margin applies).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tie_margin_rel: Option<f32>,
    /// Align at most this many tied transcripts.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tie_cands: Option<usize>,
    /// Score charged per edit (unit-cost scoring).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tie_kappa: Option<f32>,
    /// Score charged per affine-score point below the best.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tie_lambda: Option<f32>,
    /// Affine-score only texts within this many edits of the best (unset: all).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tie_affine_within: Option<u32>,
    /// When tied hits share no anchor, group them by shared anchors instead of aligning
    /// each whole read.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tie_group: Option<bool>,
    /// Distances beyond the best candidate's plus this are not computed exactly.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tie_cap: Option<u32>,
    /// Experimental: score segments by unit-cost edit distance instead of affine gaps.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tie_unit_cost: Option<bool>,
    /// Experimental: align whole reads instead of only where tied transcripts differ.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tie_whole_read: Option<bool>,
    /// Experimental: whole-read distances cover only the span the tied hits anchor.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tie_whole_trim: Option<bool>,
    /// Experimental: compare near-tied orientations by whole-read edit distance instead of
    /// clipped affine scores.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tie_orient_edit: Option<bool>,
    /// Experimental: points a candidate must gain past the shortest annotated end at a read
    /// end for them to count (0: all count).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tie_flank_min: Option<i64>,

    // Read ends (tie-break and span extension).
    /// Clip read ends that don't align (adapters, primers) instead of forcing them in.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub clip_ends: Option<bool>,
    /// Bonus for aligning through to a read's end.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end_bonus: Option<i32>,
    /// X-drop for read-end alignment.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub xdrop: Option<i32>,
    /// Band of the affine alignments (cells off the diagonals).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub band: Option<usize>,
    /// Transcript bases searched beyond a read end's projection: `pad` + 1/`pad_div` of
    /// the unanchored read end.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pad: Option<u32>,
    /// See `pad`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pad_div: Option<u32>,
    /// Keep read-end bands on the diagonal from the last shared anchor.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end_diagonal: Option<bool>,
    /// Align read ends in raw (uncompressed) bases (needs an HPC index with raw transcripts).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub raw_ends: Option<bool>,
    /// Trim a detected poly(A) tail from the read's 3' end before aligning it.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tail_trim: Option<bool>,
    /// Shortest poly(A) tail recognized: a window this long with at most one other base.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tail_run: Option<usize>,
    /// Points for isoforms whose 3' end lies where the read's poly(A) tail starts (0: off).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tail_bonus: Option<i32>,
    /// Near-tied read orientations: points for the one whose read ends in a poly(A) tail (0: off).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub orient_tail: Option<i32>,
    /// Relative margin within which both read orientations are tie-broken.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub orient_margin_rel: Option<f32>,

    // Rescue.
    /// Re-seed reads that map poorly at a lower LCP level within their candidate genes.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rescue: Option<bool>,
    /// Rescue seed scheme.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rescue_scheme: Option<String>,
    /// Rescue candidate genes per strand.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rescue_genes: Option<usize>,
    /// Rescue: reads shorter than this (parsed bases) get proportionally more candidate genes.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rescue_genes_len: Option<u32>,
    /// Rescue: transcripts re-seeded per gene (those with the most anchors).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rescue_tx: Option<usize>,
    /// Rescue: candidate genes must score this fraction of the best gene.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rescue_alpha: Option<f32>,
    /// Rescue: what of each transcript to re-seed.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rescue_window: Option<RescueWindow>,
    /// Rescue: re-seed the other orientation only if its best score is within this margin
    /// of the chosen one's (negative: never).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rescue_other_margin: Option<f32>,
    /// Rescue: re-seed the other orientation regardless of its score.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rescue_other_always: Option<bool>,
    /// Rescue: skip transcripts longer than this (compressed bases).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rescue_max_tx_len: Option<u32>,
}

/// Set `$dst` from an optional setting.
macro_rules! set {
    ($src:expr => $dst:expr) => {
        if let Some(v) = $src.clone() {
            $dst = v;
        }
    };
}

impl MapSettings {
    /// Every field's name, as serialized.
    pub const FIELDS: &[&str] = &[
        "preset",
        "strands",
        "report_within",
        "min_coverage",
        "extend_ends",
        "extend_within",
        "keep_genes",
        "keep_genes_tie",
        "max_genes",
        "max_seed_genes",
        "min_gene_frac",
        "seed_span_cap",
        "seed_gene_weight_max",
        "chain_alpha",
        "chain_window",
        "chain_max_gap",
        "chain_max_overlap",
        "chain_beta",
        "merge_anchors",
        "chainer",
        "tiebreak",
        "tie_margin",
        "tie_margin_rel",
        "tie_cands",
        "tie_kappa",
        "tie_lambda",
        "tie_affine_within",
        "tie_group",
        "tie_cap",
        "tie_unit_cost",
        "tie_whole_read",
        "tie_whole_trim",
        "tie_orient_edit",
        "tie_flank_min",
        "clip_ends",
        "end_bonus",
        "xdrop",
        "band",
        "pad",
        "pad_div",
        "end_diagonal",
        "raw_ends",
        "tail_trim",
        "tail_run",
        "tail_bonus",
        "orient_tail",
        "orient_margin_rel",
        "rescue",
        "rescue_scheme",
        "rescue_genes",
        "rescue_genes_len",
        "rescue_tx",
        "rescue_alpha",
        "rescue_window",
        "rescue_other_margin",
        "rescue_other_always",
        "rescue_max_tx_len",
    ];

    /// The preset's settings, nothing changed.
    pub fn preset(preset: Preset) -> Self {
        MapSettings {
            preset,
            ..Default::default()
        }
    }

    /// Every setting that only applies when the tie-break runs, by name, if set.
    fn tiebreak_only(&self) -> Vec<&'static str> {
        let mut v = Vec::new();
        macro_rules! chk {
            ($($f:ident),*) => { $( if self.$f.is_some() { v.push(stringify!($f)); } )* };
        }
        chk!(
            tie_margin,
            tie_margin_rel,
            tie_cands,
            tie_kappa,
            tie_lambda,
            tie_affine_within,
            tie_group,
            tie_cap,
            tie_unit_cost,
            tie_whole_read,
            tie_whole_trim,
            tie_orient_edit,
            tie_flank_min,
            raw_ends,
            tail_trim,
            tail_run,
            tail_bonus,
            orient_tail,
            orient_margin_rel,
            end_diagonal
        );
        v
    }

    fn rescue_only(&self) -> Vec<&'static str> {
        let mut v = Vec::new();
        macro_rules! chk {
            ($($f:ident),*) => { $( if self.$f.is_some() { v.push(stringify!($f)); } )* };
        }
        chk!(
            rescue_scheme,
            rescue_genes,
            rescue_genes_len,
            rescue_tx,
            rescue_alpha,
            rescue_window,
            rescue_other_margin,
            rescue_other_always,
            rescue_max_tx_len
        );
        v
    }

    /// The complete options: the preset's, with the settings given applied. Fails on
    /// settings that contradict each other; settings that would have no effect (e.g.
    /// tie-break settings with the tie-break off) are logged as warnings.
    pub fn resolve(&self) -> Result<MapOptions> {
        let mut o = MapOptions::preset(self.preset);
        let p = &mut o.params;
        set!(self.strands => p.strands);
        set!(self.report_within => p.delta_rep);
        set!(self.min_coverage => p.pi_cov);
        set!(self.merge_anchors => p.merge_collinear);

        let a = &mut p.anchor;
        set!(self.keep_genes => a.f_keep);
        set!(self.keep_genes_tie => a.f_tie);
        set!(self.max_genes => a.f_max);
        set!(self.max_seed_genes => a.h_max);
        set!(self.min_gene_frac => a.alpha_fam);
        set!(self.seed_span_cap => a.s_max);
        set!(self.seed_gene_weight_max => a.g_max);
        ensure!(a.f_keep >= 1, "keep_genes must be at least 1");
        ensure!(a.f_max >= a.f_keep, "max_genes must be at least keep_genes");

        let c = &mut p.chain;
        set!(self.chain_alpha => c.alpha);
        set!(self.chain_beta => c.beta);
        set!(self.chain_window => c.h);
        set!(self.chain_max_gap => c.g_max);
        set!(self.chain_max_overlap => c.o_max);
        c.check()?;

        if let Some(k) = self.chainer {
            o.chainer = k;
        }

        match (self.tiebreak, o.params.tiebreak.is_some()) {
            (Some(false), _) => o.params.tiebreak = None,
            (Some(true), false) => o.params.tiebreak = Some(Default::default()),
            _ => {}
        }
        if let Some(tb) = o.params.tiebreak.as_mut() {
            set!(self.tie_margin => tb.margin_abs);
            set!(self.tie_margin_rel => tb.margin_rel);
            set!(self.tie_cands => tb.max_cands);
            set!(self.tie_kappa => tb.kappa);
            set!(self.tie_lambda => tb.lambda);
            if self.tie_affine_within.is_some() {
                tb.affine_within = self.tie_affine_within;
            }
            set!(self.tie_group => tb.group_fallback);
            set!(self.tie_cap => tb.cap);
            if self.tie_unit_cost == Some(true) {
                tb.affine = None;
            }
            if let Some(w) = self.tie_whole_read {
                tb.segmented = !w;
            }
            set!(self.tie_whole_trim => tb.whole_trim);
            if let Some(e) = self.tie_orient_edit {
                tb.orient_clipped = !e;
            }
            set!(self.tie_flank_min => tb.flank_min);
            set!(self.clip_ends => tb.clip_ends);
            set!(self.end_bonus => tb.end_bonus);
            set!(self.xdrop => tb.xdrop);
            set!(self.band => tb.band);
            set!(self.pad => tb.pad);
            set!(self.pad_div => tb.pad_div);
            set!(self.end_diagonal => tb.end_diagonal);
            set!(self.raw_ends => tb.raw_ends);
            set!(self.tail_trim => tb.tail_trim);
            set!(self.tail_run => tb.tail_run);
            set!(self.tail_bonus => tb.tail_bonus);
            set!(self.orient_tail => tb.orient_tail);
            set!(self.orient_margin_rel => tb.orient_margin_rel);
            ensure!(tb.pad_div >= 1, "pad_div must be at least 1");
            ensure!(tb.max_cands >= 1, "tie_cands must be at least 1");
            ensure!(tb.tail_run >= 2, "tail_run must be at least 2");
            ensure!(
                tb.orient_margin_rel >= 0.0,
                "orient_margin_rel must be non-negative"
            );
        } else {
            let unused = self.tiebreak_only();
            if !unused.is_empty() {
                log::warn!("tie-break is off; ignoring {}", unused.join(", "));
            }
        }

        match (self.rescue, o.params.rescue.is_some()) {
            (Some(false), _) => o.params.rescue = None,
            (Some(true), false) => o.params.rescue = Some(RescueParams::default()),
            _ => {}
        }
        if let Some(r) = o.params.rescue.as_mut() {
            set!(self.rescue_scheme => r.spec);
            set!(self.rescue_genes => r.max_genes);
            set!(self.rescue_genes_len => r.genes_len);
            set!(self.rescue_tx => r.max_tx_per_gene);
            set!(self.rescue_alpha => r.alpha);
            set!(self.rescue_max_tx_len => r.max_tx_len);
            match self.rescue_window {
                Some(RescueWindow::Full) => r.window = None,
                Some(RescueWindow::Window { pad, div }) => r.window = Some((pad, div)),
                None => {}
            }
            if self.rescue_other_always == Some(true) {
                if self.rescue_other_margin.is_some() {
                    bail!("rescue_other_margin and rescue_other_always contradict each other");
                }
                r.other_margin = None;
            } else if self.rescue_other_margin.is_some() {
                r.other_margin = self.rescue_other_margin;
            }
        } else {
            let unused = self.rescue_only();
            if !unused.is_empty() {
                log::warn!("rescue is off; ignoring {}", unused.join(", "));
            }
        }

        if let Some(on) = self.extend_ends {
            o.set_extend_ends(on);
        }
        match o.params.ends.as_mut() {
            Some(e) => {
                set!(self.extend_within => e.within);
                // Read-end alignment settings apply to span extension as to the tie-break.
                set!(self.end_bonus => e.end_bonus);
                set!(self.xdrop => e.xdrop);
                set!(self.band => e.band);
                set!(self.pad => e.pad);
                set!(self.pad_div => e.pad_div);
                ensure!(e.pad_div >= 1, "pad_div must be at least 1");
            }
            None => {
                if self.extend_within.is_some() {
                    log::warn!("span extension is off; ignoring extend_within");
                }
            }
        }
        if let ChainerKind::Refine { beam: Some(0) } = o.chainer {
            bail!("the refine chainer's beam must be at least 1");
        }
        Ok(o)
    }

    /// Resolve the settings and build a mapper over `index`.
    pub fn mapper(&self, index: Arc<Index>) -> Result<Mapper> {
        Mapper::new(index, &self.resolve()?)
    }
}

impl<S: map_settings_builder::IsComplete> MapSettingsBuilder<S> {
    /// Finish the settings and build a mapper over `index`.
    pub fn mapper(self, index: Arc<Index>) -> Result<Mapper> {
        self.build().mapper(index)
    }

    /// Finish the settings and resolve them to the complete options.
    pub fn resolve(self) -> Result<MapOptions> {
        self.build().resolve()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unset_settings_are_the_preset() {
        for p in [Preset::Ont, Preset::Drna, Preset::Hifi] {
            assert_eq!(MapSettings::preset(p).resolve().unwrap(), MapOptions::preset(p));
            assert_eq!(MapSettings::builder(p).build(), MapSettings::preset(p));
        }
    }

    #[test]
    fn settings_override_and_validate() {
        let o = MapSettings::builder(Preset::Hifi)
            .extend_ends(true)
            .rescue(true)
            .tie_group(false)
            .xdrop(40)
            .resolve()
            .unwrap();
        assert!(o.params.rescue.is_some());
        let tb = o.params.tiebreak.as_ref().unwrap();
        assert!(!tb.group_fallback);
        assert_eq!(tb.xdrop, 40);
        let e = o.params.ends.as_ref().unwrap();
        assert_eq!(e.xdrop, 40);
        // Extended spans are scored as the preset scores.
        assert_eq!(e.sc, map::align::Affine::HIFI);
        assert!(MapSettings::builder(Preset::Ont).pad_div(0).resolve().is_err());
        // Chain limits the kernels can't represent (i32 gaps, window arithmetic) are errors.
        assert!(
            MapSettings::builder(Preset::Ont)
                .chain_window(0)
                .resolve()
                .is_err()
        );
        assert!(
            MapSettings::builder(Preset::Ont)
                .chain_window(usize::MAX)
                .resolve()
                .is_err()
        );
        assert!(
            MapSettings::builder(Preset::Ont)
                .chain_max_gap(3_000_000_000)
                .resolve()
                .is_err()
        );
        assert!(
            MapSettings::builder(Preset::Ont)
                .chain_max_overlap(1 << 31)
                .resolve()
                .is_err()
        );
        let o = MapSettings::builder(Preset::Ont)
            .chain_max_gap(i32::MAX as i64)
            .chain_max_overlap(-10)
            .chain_window(1000)
            .resolve()
            .unwrap();
        assert_eq!(
            (o.params.chain.g_max, o.params.chain.o_max),
            (i32::MAX as i64, -10)
        );
        assert!(
            MapSettings::builder(Preset::Ont)
                .rescue_other_always(true)
                .rescue_other_margin(3.0)
                .resolve()
                .is_err()
        );
        let w: RescueWindow = "10,50".parse().unwrap();
        assert_eq!(w, RescueWindow::Window { pad: 10, div: 50 });
        assert_eq!("full".parse::<RescueWindow>().unwrap(), RescueWindow::Full);
        assert!("10".parse::<RescueWindow>().is_err());
    }

    #[test]
    fn field_list_is_complete() {
        // An exhaustive pattern: adding a field without listing it in FIELDS fails here.
        let MapSettings {
            preset: _,
            strands: _,
            report_within: _,
            min_coverage: _,
            extend_ends: _,
            extend_within: _,
            keep_genes: _,
            keep_genes_tie: _,
            max_genes: _,
            max_seed_genes: _,
            min_gene_frac: _,
            seed_span_cap: _,
            seed_gene_weight_max: _,
            chain_alpha: _,
            chain_window: _,
            chain_max_gap: _,
            chain_max_overlap: _,
            chain_beta: _,
            merge_anchors: _,
            chainer: _,
            tiebreak: _,
            tie_margin: _,
            tie_margin_rel: _,
            tie_cands: _,
            tie_kappa: _,
            tie_lambda: _,
            tie_affine_within: _,
            tie_group: _,
            tie_cap: _,
            tie_unit_cost: _,
            tie_whole_read: _,
            tie_whole_trim: _,
            tie_orient_edit: _,
            tie_flank_min: _,
            clip_ends: _,
            end_bonus: _,
            xdrop: _,
            band: _,
            pad: _,
            pad_div: _,
            end_diagonal: _,
            raw_ends: _,
            tail_trim: _,
            tail_run: _,
            tail_bonus: _,
            orient_tail: _,
            orient_margin_rel: _,
            rescue: _,
            rescue_scheme: _,
            rescue_genes: _,
            rescue_genes_len: _,
            rescue_tx: _,
            rescue_alpha: _,
            rescue_window: _,
            rescue_other_margin: _,
            rescue_other_always: _,
            rescue_max_tx_len: _,
        } = MapSettings::default();
        assert_eq!(
            MapSettings::FIELDS,
            [
                stringify!(preset),
                stringify!(strands),
                stringify!(report_within),
                stringify!(min_coverage),
                stringify!(extend_ends),
                stringify!(extend_within),
                stringify!(keep_genes),
                stringify!(keep_genes_tie),
                stringify!(max_genes),
                stringify!(max_seed_genes),
                stringify!(min_gene_frac),
                stringify!(seed_span_cap),
                stringify!(seed_gene_weight_max),
                stringify!(chain_alpha),
                stringify!(chain_window),
                stringify!(chain_max_gap),
                stringify!(chain_max_overlap),
                stringify!(chain_beta),
                stringify!(merge_anchors),
                stringify!(chainer),
                stringify!(tiebreak),
                stringify!(tie_margin),
                stringify!(tie_margin_rel),
                stringify!(tie_cands),
                stringify!(tie_kappa),
                stringify!(tie_lambda),
                stringify!(tie_affine_within),
                stringify!(tie_group),
                stringify!(tie_cap),
                stringify!(tie_unit_cost),
                stringify!(tie_whole_read),
                stringify!(tie_whole_trim),
                stringify!(tie_orient_edit),
                stringify!(tie_flank_min),
                stringify!(clip_ends),
                stringify!(end_bonus),
                stringify!(xdrop),
                stringify!(band),
                stringify!(pad),
                stringify!(pad_div),
                stringify!(end_diagonal),
                stringify!(raw_ends),
                stringify!(tail_trim),
                stringify!(tail_run),
                stringify!(tail_bonus),
                stringify!(orient_tail),
                stringify!(orient_margin_rel),
                stringify!(rescue),
                stringify!(rescue_scheme),
                stringify!(rescue_genes),
                stringify!(rescue_genes_len),
                stringify!(rescue_tx),
                stringify!(rescue_alpha),
                stringify!(rescue_window),
                stringify!(rescue_other_margin),
                stringify!(rescue_other_always),
                stringify!(rescue_max_tx_len)
            ]
        );
        let mut j = serde_json::to_value(MapSettings::preset(Preset::Ont)).unwrap();
        for f in MapSettings::FIELDS.iter().filter(|f| **f != "preset") {
            assert!(j.get(*f).is_none());
            j[*f] = serde_json::Value::Null;
        }
        // Every listed name is a real field (unknown fields are rejected).
        assert!(serde_json::from_value::<MapSettings>(j).is_ok());
    }

    #[test]
    fn settings_json_roundtrip() {
        let s = MapSettings::builder(Preset::Drna)
            .tie_kappa(5.0)
            .rescue_window(RescueWindow::Full)
            .build();
        let j = serde_json::to_string(&s).unwrap();
        assert_eq!(serde_json::from_str::<MapSettings>(&j).unwrap(), s);
        assert!(!j.contains("null"), "{j}");
    }
}
