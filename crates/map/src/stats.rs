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
    END_SEGMENTS,
    AFFINE_TEXT_BASES,
    AFFINE_TEXT_NEW_BASES,
    MYERS_CALLS,
    MYERS_CELLS,
    BOTH_ORIENTATIONS,
);

#[inline]
pub fn add(c: &AtomicU64, v: u64) {
    if on() {
        c.fetch_add(v, Relaxed);
    }
}
