//! `Bump` as the global allocator of a program that never calls `enable`: every
//! allocation is the wrapped allocator's.

use std::alloc::System;

use nudge_alloc::{self as alloc, Bump};

#[global_allocator]
static GLOBAL: Bump<System> = Bump::new(System);

#[test]
fn nothing_comes_from_the_arenas_until_they_are_enabled() {
    assert!(!alloc::enabled());
    let blocks: Vec<Vec<u8>> = (0..1000usize).map(|i| vec![i as u8; 1 + i]).collect();
    assert!(blocks.iter().all(|b| !alloc::in_arena(b.as_ptr())));
    assert!(
        blocks
            .iter()
            .enumerate()
            .all(|(i, b)| b.len() == 1 + i && b.iter().all(|&x| x == i as u8))
    );
    assert_eq!(alloc::stats(), alloc::Stats::default());
}
