//! Global work counters for profiling (`SPEARFISH_STATS` set in the environment); printed by
//! `spearfish map`. Off by default: one cached flag check per event.

use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

pub fn on() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("SPEARFISH_STATS").is_some())
}

macro_rules! counters {
    ($($name:ident),* $(,)?) => {
        $(pub static $name: AtomicU64 = AtomicU64::new(0);)*
        /// All counters as (name, value).
        pub fn snapshot() -> Vec<(&'static str, u64)> {
            vec![$((stringify!($name), $name.load(Relaxed))),*]
        }
    };
}

counters!(
    READS,
    RESCUE_STRANDS,
    RESCUE_GENES,
    RESCUE_TX,
    RESCUE_TX_BASES,
    TIEBREAK_CALLS,
    TIEBREAK_TIED,
    TIEBREAK_WHOLE_READ,
    TIEBREAK_SEGMENTS,
    AFFINE_CALLS,
    AFFINE_CELLS,
    AFFINE_END_CELLS,
    AFFINE_TEXT_BASES,
    AFFINE_TEXT_NEW_BASES,
    MYERS_CALLS,
    MYERS_CELLS,
    BOTH_ORIENTATIONS,
    SPAN_ENDS,
    SPAN_END_CALLS,
    SPAN_END_CELLS,
    DP_SEG_CELLS,
    DP_SEG_ROWS,
    DP_TBEND_CELLS,
    DP_TBEND_ROWS,
    DP_CLIP_MID_CELLS,
    DP_CLIP_MID_ROWS,
    DP_CLIP_END_CELLS,
    DP_CLIP_END_ROWS,
    DP_SPAN_CELLS,
    DP_SPAN_ROWS,
    ANCHORS,
    FAMILIES,
    KEPT_FAMILIES,
    KEPT_ANCHORS,
    MERGED_ANCHORS,
    CHAIN_TX,
    CHAIN_PLACEMENTS,
    RESCUE_PTS,
    RESCUE_PTS_DP,
);

#[inline]
pub fn add(c: &AtomicU64, v: u64) {
    if on() {
        c.fetch_add(v, Relaxed);
    }
}
