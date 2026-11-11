//! The paths an allocation or a free takes within one thread's chunk, each reached on purpose
//! and checked by address: growth and shrinking in place, a chunk started again where it is,
//! and a chunk that other threads emptied. `Bump` is called directly here (the test program's
//! own allocator is the system's), so a test knows every block of its thread's arena.

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::Mutex;
use std::thread;

use nudgy::Bump;

static B: Bump<System> = Bump::new(System);

/// The largest block an arena serves.
const M: usize = nudgy::DEFAULT_MAX_BLOCK;
const CHUNK: usize = 1 << 16;

/// One test at a time: each counts on its thread's arena holding its own blocks only.
static ONE: Mutex<()> = Mutex::new(());

fn lay(n: usize) -> Layout {
    Layout::from_size_align(n, 8).unwrap()
}

struct Sent(*mut u8);
// SAFETY: a block of the allocator, which the receiving thread writes and frees.
unsafe impl Send for Sent {}

/// Runs `f` on a thread of its own, which so starts with an arena of empty chunks. Does
/// nothing on a target without the arenas.
fn on_a_new_thread(f: impl FnOnce() + Send + 'static) {
    if !nudgy::enable() {
        return;
    }
    let _one = ONE.lock().unwrap_or_else(|e| e.into_inner());
    thread::spawn(f).join().unwrap();
}

/// Writes the blocks (each of `n` bytes) and frees them, on another thread.
fn free_elsewhere(blocks: Vec<Sent>, n: usize) -> thread::JoinHandle<()> {
    thread::spawn(move || {
        for Sent(p) in blocks {
            // SAFETY: a live block of `n` bytes that this thread now owns, freed once.
            unsafe {
                p.write_bytes(0xAB, n);
                B.dealloc(p, lay(n));
            }
        }
    })
}

#[test]
fn the_last_block_grows_and_shrinks_where_it_is() {
    // SAFETY (here and below): every block is freed or reallocated once, with its layout.
    on_a_new_thread(|| unsafe {
        let keep = B.alloc(lay(48));
        let p = B.alloc(lay(8));
        assert_eq!(p as usize, keep as usize + 48);
        p.write_bytes(5, 8);
        let mut size = 8;
        while size < M {
            assert_eq!(B.realloc(p, lay(size), size * 2), p);
            p.add(size).write_bytes(5, size);
            size *= 2;
        }
        assert!((0..M).all(|i| *p.add(i) == 5));
        // Shrunk, it gives its end back: the next block follows it.
        assert_eq!(B.realloc(p, lay(M), 100), p);
        let q = B.alloc(lay(16));
        assert_eq!(q as usize, p as usize + 112);
        // No longer the last block: growing copies it.
        let moved = B.realloc(p, lay(100), 200);
        assert_ne!(moved, p);
        assert!(nudgy::in_arena(moved));
        assert!((0..100).all(|i| *moved.add(i) == 5));
        B.dealloc(moved, lay(200));
        B.dealloc(q, lay(16));
        B.dealloc(keep, lay(48));
        // All freed: the chunk starts again where `keep` was.
        let again = B.alloc(lay(16));
        assert_eq!(again, keep);
        B.dealloc(again, lay(16));
    });
}

#[test]
fn growth_at_the_end_of_a_chunk_moves() {
    on_a_new_thread(|| unsafe {
        let first = B.alloc(lay(M));
        assert!(nudgy::in_arena(first));
        let mut held = vec![first];
        for _ in 0..CHUNK / M - 2 {
            held.push(B.alloc(lay(M)));
        }
        // The chunk's last 8 KB: a block that ends 16 bytes short grows to the end in place,
        // and one byte more takes it out of the arena.
        let last = B.alloc(lay(M - 16));
        assert_eq!(last as usize + M, first as usize + CHUNK);
        assert_eq!(B.realloc(last, lay(M - 16), M), last);
        last.write_bytes(2, M);
        let out = B.realloc(last, lay(M), M + 1);
        assert!(!nudgy::in_arena(out));
        assert!((0..M).all(|i| *out.add(i) == 2));
        B.dealloc(out, lay(M + 1));
        // The 8 KB are free again. A block that ends with the chunk cannot grow there.
        let a = B.alloc(lay(M - 16));
        let b = B.alloc(lay(16));
        assert_eq!(b as usize + 16, first as usize + CHUNK);
        b.write_bytes(3, 16);
        let b2 = B.realloc(b, lay(16), 32);
        assert_ne!(b2, b);
        assert!((0..16).all(|i| *b2.add(i) == 3));
        B.dealloc(b2, lay(32));
        B.dealloc(a, lay(M - 16));
        for p in held {
            B.dealloc(p, lay(M));
        }
    });
}

/// Other threads free every block of the chunk but the last, which its thread frees: the
/// chunk is empty then, and starts again where it is.
#[test]
fn a_chunk_starts_again_after_frees_by_another_thread() {
    on_a_new_thread(|| unsafe {
        let first = B.alloc(lay(M));
        let mut sent = vec![Sent(first)];
        for _ in 0..6 {
            sent.push(Sent(B.alloc(lay(M))));
        }
        let freeing = free_elsewhere(sent, M);
        // Not joined first: the allocator's own ordering is what lets this thread write where
        // the other one wrote.
        loop {
            let q = B.alloc(lay(16));
            q.write_bytes(1, 16);
            B.dealloc(q, lay(16));
            if q == first {
                break;
            }
            thread::yield_now();
        }
        freeing.join().unwrap();
    });
}

/// A chunk filled to its end and then emptied by another thread: the thread that filled it,
/// looking for an empty chunk, finds the one it is in.
#[test]
fn a_full_chunk_emptied_by_another_thread_is_taken_again() {
    on_a_new_thread(|| unsafe {
        let first = B.alloc(lay(M));
        let mut sent = vec![Sent(first)];
        for _ in 0..CHUNK / M - 1 {
            sent.push(Sent(B.alloc(lay(M))));
        }
        free_elsewhere(sent, M).join().unwrap();
        // Twice: the chunk keeps its count of what the other thread freed when it starts again.
        for _ in 0..2 {
            let q = B.alloc(lay(48));
            assert_eq!(q, first);
            let r = B.alloc(lay(48));
            assert_eq!(r as usize, first as usize + 48);
            // The first block (counted), then the last (given back): empty, and started
            // again.
            B.dealloc(q, lay(48));
            B.dealloc(r, lay(48));
        }
        let again = B.alloc(lay(16));
        assert_eq!(again, first);
        B.dealloc(again, lay(16));
    });
}

/// A thread leaves its blocks to another and exits; the next thread has its arena, and uses
/// the chunk again once the blocks are freed, which may be while it is looking.
#[test]
fn the_next_thread_reuses_a_chunk_freed_after_the_hand_over() {
    if !nudgy::enable() {
        return;
    }
    let _one = ONE.lock().unwrap_or_else(|e| e.into_inner());
    let sent = thread::spawn(|| {
        // SAFETY: as above; the blocks are freed by `free_elsewhere`.
        let first = unsafe { B.alloc(lay(M)) };
        let mut sent = vec![Sent(first)];
        for _ in 0..5 {
            sent.push(Sent(unsafe { B.alloc(lay(M)) }));
        }
        sent
    })
    .join()
    .unwrap();
    let first = sent[0].0 as usize;
    let freeing = free_elsewhere(sent, M);
    thread::spawn(move || unsafe {
        let mut held = Vec::new();
        loop {
            let q = B.alloc(lay(M));
            q.write_bytes(1, M);
            if q as usize == first {
                B.dealloc(q, lay(M));
                break;
            }
            held.push(q);
            if held.len() > 12 {
                for q in held.drain(..) {
                    B.dealloc(q, lay(M));
                }
            }
            thread::yield_now();
        }
        for q in held {
            B.dealloc(q, lay(M));
        }
    })
    .join()
    .unwrap();
    freeing.join().unwrap();
}
