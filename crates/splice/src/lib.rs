//! Spliced genome mapping (in development; see the `genome-mapping` branch plan).

/// A kernel compiled twice: as written, and again with AVX2 allowed, chosen when the program
/// runs. `$body` is an `#[inline(always)]` function; its loops are plain and the compiler
/// vectorizes them for whichever instruction set the copy is built for.
macro_rules! avx2_variant {
    ($(#[$meta:meta])* $vis:vis fn $name:ident = $body:ident, $avx:ident ($($a:ident : $t:ty),* $(,)?) -> $r:ty) => {
        #[cfg(target_arch = "x86_64")]
        #[allow(clippy::too_many_arguments)]
        mod $avx {
            #[allow(unused_imports)]
            use super::*;
            #[target_feature(enable = "avx2")]
            pub(super) unsafe fn avx2($($a: $t),*) -> $r {
                $body($($a),*)
            }
            #[target_feature(enable = "avx512f,avx512bw,avx512vl,avx512dq")]
            pub(super) unsafe fn avx512($($a: $t),*) -> $r {
                $body($($a),*)
            }
        }
        $(#[$meta])*
        #[inline]
        $vis fn $name($($a: $t),*) -> $r {
            #[cfg(target_arch = "x86_64")]
            {
                if std::is_x86_feature_detected!("avx512bw") && std::is_x86_feature_detected!("avx512vl") && std::is_x86_feature_detected!("avx512dq") {
                    // SAFETY: the features were just detected.
                    return unsafe { $avx::avx512($($a),*) };
                }
                if std::is_x86_feature_detected!("avx2") {
                    // SAFETY: AVX2 was just detected.
                    return unsafe { $avx::avx2($($a),*) };
                }
            }
            $body($($a),*)
        }
    };
}

pub mod bam;
pub mod base;
mod gap;
pub mod index;
pub mod junctions;
pub mod map;
