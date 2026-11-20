//! Spliced genome mapping (in development; see the `genome-mapping` branch plan).

/// A kernel compiled twice: as written, and again with AVX2 allowed, chosen when the program
/// runs. `$body` is an `#[inline(always)]` function; its loops are plain and the compiler
/// vectorizes them for whichever instruction set the copy is built for.
macro_rules! avx2_variant {
    ($(#[$meta:meta])* $vis:vis fn $name:ident = $body:ident, $avx:ident ($($a:ident : $t:ty),* $(,)?) -> $r:ty) => {
        #[cfg(target_arch = "x86_64")]
        #[target_feature(enable = "avx2")]
        #[allow(clippy::too_many_arguments)]
        unsafe fn $avx($($a: $t),*) -> $r {
            $body($($a),*)
        }
        $(#[$meta])*
        #[inline]
        $vis fn $name($($a: $t),*) -> $r {
            #[cfg(target_arch = "x86_64")]
            if std::is_x86_feature_detected!("avx2") {
                // SAFETY: AVX2 was just detected.
                return unsafe { $avx($($a),*) };
            }
            $body($($a),*)
        }
    };
}

pub mod bam;
mod gap;
pub mod index;
pub mod junctions;
pub mod map;
