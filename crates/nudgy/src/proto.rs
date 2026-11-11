//! What the allocator's threads share, and the order they read and write it in.
//!
//! - An arena has one thread at a time (a bit of a `TAKEN` word). What the thread leaves in the
//!   arena's counts when it gives the arena back must be visible to the next thread.
//! - A chunk has two counts, in bytes and modulo 2^32: `LIVE`, what was handed out of it
//!   minus what the arena's thread freed, which only that thread reads and writes, and
//!   `FREED`, what other threads freed. The chunk's bytes are handed out again when the two
//!   are equal, and whatever the freeing threads wrote to the blocks must come before that.
//!   A thread keeps the `LIVE` of the chunk it is allocating from to itself, and writes it
//!   here when it moves to another chunk or gives the arena back.
//! - `enable` publishes the reservation's address before any thread can get a block of it.
//!
//! The functions take the atomics by reference, of the types `crate::sync` names:
//! `tests/loom.rs` includes this file with loom's atomics in their place and runs small models
//! of the cases above through every interleaving. The statics themselves are in `lib.rs`.

use crate::sync::Ordering::{Acquire, Relaxed, Release};
use crate::sync::{AtomicBool, AtomicU32, AtomicU64, AtomicUsize};

/// Takes the lowest arena of `word` that no thread has, leaving out those in `below`.
pub(crate) fn take(word: &AtomicU64, below: u64) -> Option<u32> {
    let mut cur = word.load(Relaxed);
    while cur | below != u64::MAX {
        let bit = (!(cur | below)).trailing_zeros();
        // Acquire: the counters the arena's last thread left in LIVE are visible.
        match word.compare_exchange_weak(cur, cur | 1 << bit, Acquire, Relaxed) {
            Ok(_) => return Some(bit),
            Err(now) => cur = now,
        }
    }
    None
}

/// Gives arena `bit` of `word` back.
pub(crate) fn give_back(word: &AtomicU64, bit: u32) {
    word.fetch_and(!(1 << bit), Release);
}

/// The arena's thread stops allocating from a chunk whose `LIVE` it counted up to `now`.
pub(crate) fn left(live: &AtomicU32, now: u32) {
    live.store(now, Relaxed);
}

/// A chunk's `LIVE`, as its arena's thread left it. For the arena's thread only.
#[inline(always)]
pub(crate) fn live(live: &AtomicU32) -> u32 {
    live.load(Relaxed)
}

/// Counts `n` bytes that the arena's thread freed in a chunk it is not allocating from.
#[inline(always)]
pub(crate) fn freed_by_owner(live: &AtomicU32, n: u32) {
    live.store(live.load(Relaxed).wrapping_sub(n), Relaxed);
}

/// Counts `n` bytes that a thread other than the arena's freed.
#[inline(always)]
pub(crate) fn freed_by_other(freed: &AtomicU32, n: u32) {
    freed.fetch_add(n, Release);
}

/// Whether a chunk whose `LIVE` is `live` holds no block. For the arena's thread only.
#[inline(always)]
pub(crate) fn holds_none(live: u32, freed: &AtomicU32) -> bool {
    // Acquire: what the threads that freed the chunk's blocks wrote to them came first.
    live == freed.load(Acquire)
}

/// `enable`: the reservation is at `at`, and blocks of up to `limit` bytes come from it.
pub(crate) fn turn_on(base: &AtomicUsize, on: &AtomicBool, small: &AtomicUsize, at: usize, limit: usize) {
    base.store(at, Relaxed);
    on.store(true, Release);
    small.store(limit, Release);
}

/// Whether the arenas are on. Acquire: with `on`, the reservation's address is visible, to the
/// thread that asks and (through the block it then gets) to whichever thread frees the block.
pub(crate) fn is_on(on: &AtomicBool) -> bool {
    on.load(Acquire)
}

/// The largest size served from an arena now: zero until `turn_on`.
#[inline(always)]
pub(crate) fn limit(small: &AtomicUsize) -> usize {
    small.load(Relaxed)
}

/// The reservation's address: the placeholder that `base` started with until `turn_on`.
#[inline(always)]
pub(crate) fn base(base: &AtomicUsize) -> usize {
    base.load(Relaxed)
}
