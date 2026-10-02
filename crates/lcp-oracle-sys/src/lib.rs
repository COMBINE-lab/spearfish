//! Test-only FFI to lcptools (vendored C, BSD-3-Clause). Exposes core boundaries per level
//! so `lcp-core` can be checked against the reference implementation.

use std::os::raw::{c_char, c_int};
use std::sync::Once;

#[repr(C)]
struct RawCore {
    start: u32,
    end: u32,
    bit_rep: *mut u32,
    label: u64,
    bit_size: u32,
}

#[repr(C)]
struct RawLps {
    level: c_int,
    size: c_int,
    cores: *mut RawCore,
}

unsafe extern "C" {
    fn LCP_INIT();
    fn init_lps(lps: *mut RawLps, s: *const c_char, len: c_int);
    fn lps_deepen(lps: *mut RawLps, level: c_int) -> c_int;
    fn free_lps(lps: *mut RawLps);
}

static INIT: Once = Once::new();

/// Core boundaries `(start, end)` at levels `1..=max_level`, as lcptools reports them.
/// At level > 1 lcptools' `start` includes the DCT context (our `Level::ctx_start`).
/// `seq` must be ACGT only.
pub fn lcptools_boundaries(seq: &[u8], max_level: u8) -> Vec<Vec<(u32, u32)>> {
    assert!(
        seq.iter().all(|b| b"ACGTacgt".contains(b)),
        "lcptools oracle needs ACGT input"
    );
    INIT.call_once(|| unsafe { LCP_INIT() });
    let mut out = Vec::with_capacity(max_level as usize);
    if max_level == 0 {
        return out;
    }
    let mut lps = RawLps {
        level: 0,
        size: 0,
        cores: std::ptr::null_mut(),
    };
    unsafe {
        init_lps(&mut lps, seq.as_ptr() as *const c_char, seq.len() as c_int);
        for k in 1..=max_level {
            if k > 1 {
                lps_deepen(&mut lps, k as c_int);
            }
            let cores = if lps.size > 0 {
                std::slice::from_raw_parts(lps.cores, lps.size as usize)
            } else {
                &[]
            };
            out.push(cores.iter().map(|c| (c.start, c.end)).collect());
        }
        free_lps(&mut lps);
    }
    out
}

#[cfg(test)]
mod tests {
    #[test]
    fn smoke() {
        let b = super::lcptools_boundaries(
            b"CACACAGTGTCAGGGTTTACGATCGATCGGGATCGATCGTAGCTAGCATCGACTGACTAGC",
            2,
        );
        assert_eq!(b.len(), 2);
        assert!(!b[0].is_empty());
    }
}
