//! Loom models of what the allocator's threads share. `src/proto.rs` is compiled here against
//! loom's atomics, so the models run the allocator's own loads and stores, with the orderings
//! as written, through every interleaving loom allows (weaker orderings included). A block's
//! bytes are a loom `UnsafeCell`: loom reports any two accesses that are not ordered.
//!
//! The models are one arena of one chunk. They do not cover what a single thread does with its
//! own chunk (`tests/alloc.rs`, `tests/stress.rs` and Miri do).

#![cfg(target_os = "linux")]

/// What `proto` takes for `crate::sync`.
mod sync {
    pub(crate) use loom::sync::atomic::{AtomicBool, AtomicU32, AtomicU64, AtomicUsize, Ordering};
}
#[path = "../src/proto.rs"]
mod proto;

use loom::cell::UnsafeCell;
use loom::sync::Arc;
use loom::thread;
use sync::Ordering::{Acquire, Relaxed, Release};
use sync::{AtomicBool, AtomicU32, AtomicU64, AtomicUsize};

/// One arena of one chunk, and the bytes of one block of it.
struct Arena {
    taken: AtomicU64,
    live: AtomicU32,
    freed: AtomicU32,
    block: UnsafeCell<u64>,
    /// The test's own: the block was sent to another thread.
    sent: AtomicBool,
}

impl Arena {
    fn new() -> Arc<Arena> {
        Arc::new(Arena {
            // Every arena but the first is taken.
            taken: AtomicU64::new(!1),
            live: AtomicU32::new(0),
            freed: AtomicU32::new(0),
            block: UnsafeCell::new(0),
            sent: AtomicBool::new(false),
        })
    }
    fn take(&self) -> bool {
        proto::take(&self.taken, 0).is_some()
    }
    fn give_back(&self) {
        proto::give_back(&self.taken, 0);
    }
    /// What the arena's thread asks of a chunk it is not allocating from.
    fn is_empty(&self) -> bool {
        proto::holds_none(proto::live(&self.live), &self.freed)
    }
    fn write(&self, v: u64) {
        // SAFETY: loom checks that no other access is concurrent.
        self.block.with_mut(|p| unsafe { *p = v });
    }
    fn read(&self) -> u64 {
        // SAFETY: as in `write`.
        self.block.with(|p| unsafe { *p })
    }
}

/// A allocates a block, sends it to C, and exits; C reads and frees it; B takes the arena and
/// uses the chunk again once it is empty.
#[test]
fn a_block_freed_by_a_third_thread_after_its_arena_changed_hands() {
    loom::model(|| {
        let s = Arena::new();
        // Before B exists, or B could take the arena first.
        assert!(s.take());
        let a = {
            let s = s.clone();
            thread::spawn(move || {
                // A block of 16 bytes: the thread counts it itself until it leaves the chunk.
                s.write(7);
                s.sent.store(true, Release);
                proto::left(&s.live, 16);
                s.give_back();
            })
        };
        let c = {
            let s = s.clone();
            thread::spawn(move || {
                while !s.sent.load(Acquire) {
                    thread::yield_now();
                }
                assert_eq!(s.read(), 7);
                proto::freed_by_other(&s.freed, 16);
            })
        };
        let b = {
            let s = s.clone();
            thread::spawn(move || {
                if s.take() {
                    // What A counted is visible to the arena's next thread.
                    assert_eq!(s.live.load(Relaxed), 16);
                    if s.is_empty() {
                        s.write(9);
                    }
                }
            })
        };
        a.join().unwrap();
        b.join().unwrap();
        c.join().unwrap();
    });
}

/// A allocates a block, frees it, and exits; B takes the arena and uses the chunk again.
#[test]
fn an_arena_changes_hands_after_its_thread_freed_its_block() {
    loom::model(|| {
        let s = Arena::new();
        assert!(s.take());
        let a = {
            let s = s.clone();
            thread::spawn(move || {
                // A block of 16 bytes, allocated and freed: the thread's own count is 0 again.
                s.write(7);
                assert!(proto::holds_none(0, &s.freed));
                proto::left(&s.live, 0);
                s.give_back();
            })
        };
        let b = {
            let s = s.clone();
            thread::spawn(move || {
                if s.take() {
                    assert!(s.is_empty());
                    s.write(9);
                }
            })
        };
        a.join().unwrap();
        b.join().unwrap();
    });
}

/// The arena's thread has two blocks and gives one to C. C frees it while the thread frees the
/// other; the thread starts the chunk again once it sees it empty, at the free or at its next
/// look.
#[test]
fn the_owner_and_another_thread_free_the_last_two_blocks() {
    loom::model(|| {
        let s = Arena::new();
        assert!(s.take());
        // Two blocks of 16 bytes, of which the thread frees one: its own count is 16.
        s.write(1);
        let c = {
            let s = s.clone();
            thread::spawn(move || {
                assert_eq!(s.read(), 1);
                proto::freed_by_other(&s.freed, 16);
            })
        };
        if proto::holds_none(16, &s.freed) {
            s.write(3);
        } else {
            // Its next look, after it left the chunk for another.
            proto::left(&s.live, 16);
            if s.is_empty() {
                s.write(3);
            }
        }
        c.join().unwrap();
        assert!(proto::holds_none(16, &s.freed));
    });
}

/// The same in a chunk the arena's thread has left for another: its frees go to the chunk's
/// `LIVE`, and it uses the chunk again when a later look finds it empty.
#[test]
fn the_owner_and_another_thread_free_the_last_two_blocks_of_a_chunk_it_left() {
    loom::model(|| {
        let s = Arena::new();
        assert!(s.take());
        s.write(1);
        proto::left(&s.live, 32);
        let c = {
            let s = s.clone();
            thread::spawn(move || {
                assert_eq!(s.read(), 1);
                proto::freed_by_other(&s.freed, 16);
            })
        };
        proto::freed_by_owner(&s.live, 16);
        if s.is_empty() {
            s.write(3);
        }
        c.join().unwrap();
        assert!(s.is_empty());
    });
}

/// Two threads other than the arena's free its last two blocks at once: both frees count,
/// and the arena's thread, which waits for neither, does not use the chunk before both.
#[test]
fn two_other_threads_free_the_last_two_blocks() {
    loom::model(|| {
        let s = Arena::new();
        assert!(s.take());
        s.write(1);
        proto::left(&s.live, 32);
        let free = |s: &Arc<Arena>| {
            let s = s.clone();
            thread::spawn(move || {
                assert_eq!(s.read(), 1);
                proto::freed_by_other(&s.freed, 16);
            })
        };
        let (c, d) = (free(&s), free(&s));
        if s.is_empty() {
            s.write(3);
        }
        c.join().unwrap();
        d.join().unwrap();
        assert!(s.is_empty());
    });
}

/// Two threads ask for the one free arena: one gets it.
#[test]
fn one_arena_goes_to_one_thread() {
    loom::model(|| {
        let s = Arena::new();
        let got = Arc::new(AtomicU32::new(0));
        let hs: Vec<_> = (0..2)
            .map(|_| {
                let (s, got) = (s.clone(), got.clone());
                thread::spawn(move || {
                    if s.take() {
                        got.fetch_add(1, Relaxed);
                        // The arena's thread alone writes LIVE.
                        proto::left(&s.live, 16);
                        s.write(5);
                    }
                })
            })
            .collect();
        for h in hs {
            h.join().unwrap();
        }
        assert_eq!(got.load(Relaxed), 1);
        assert_eq!(s.live.load(Relaxed), 16);
    });
}

/// `enable` on one thread, an allocation on a second, and the block's free on a third: a
/// thread that gets a block, and the thread that frees it, see the reservation's address.
#[test]
fn enable_publishes_the_reservation_before_any_block() {
    const PLACEHOLDER: usize = usize::MAX;
    const AT: usize = 1 << 20;
    loom::model(|| {
        let base = Arc::new(AtomicUsize::new(PLACEHOLDER));
        let on = Arc::new(AtomicBool::new(false));
        let small = Arc::new(AtomicUsize::new(0));
        // The test's own: the address of a block sent to the third thread.
        let block = Arc::new(AtomicUsize::new(0));
        let enable = {
            let (base, on, small) = (base.clone(), on.clone(), small.clone());
            thread::spawn(move || proto::turn_on(&base, &on, &small, AT, 8192))
        };
        let allocate = {
            let (base, on, small, block) = (base.clone(), on.clone(), small.clone(), block.clone());
            thread::spawn(move || {
                // `serves`, then `refill` for a thread without an arena.
                if 64 <= proto::limit(&small) && proto::is_on(&on) {
                    let b = proto::base(&base);
                    assert_eq!(b, AT);
                    block.store(b + 16, Release);
                }
            })
        };
        let free = {
            let (base, block) = (base.clone(), block.clone());
            thread::spawn(move || {
                if block.load(Acquire) != 0 {
                    // `dealloc` finds the block in the reservation.
                    assert_eq!(proto::base(&base), AT);
                }
            })
        };
        enable.join().unwrap();
        allocate.join().unwrap();
        free.join().unwrap();
    });
}
