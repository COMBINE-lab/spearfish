//! `alloc::Bump` as this test program's global allocator, with the arenas on: blocks keep
//! their contents whichever thread frees them, through thread exits and when arenas or chunks
//! run out.

use std::alloc::System;
use std::sync::{Arc, Barrier, Mutex, mpsc};

use spearfish::alloc::{self, Bump};

#[global_allocator]
static GLOBAL: Bump<System> = Bump::new(System);

/// Whether this target has the arenas; the tests pass (and test little) where it does not.
fn on() -> bool {
    alloc::enable()
}

/// Held by the tests that count on a free arena, and by the one that takes them all.
static ARENAS: Mutex<()> = Mutex::new(());

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as usize
    }
}

/// A block of `len` bytes that says who made it.
fn block(tag: usize, len: usize) -> Vec<u8> {
    let mut v = Vec::with_capacity(len);
    v.extend((0..len).map(|i| (tag + i) as u8));
    v
}

fn check(tag: usize, v: &[u8]) {
    assert!(
        v.iter().enumerate().all(|(i, &b)| b == (tag + i) as u8),
        "block {tag} was overwritten"
    );
}

#[test]
fn small_blocks_come_from_the_arena() {
    let on = on();
    assert_eq!(alloc::enabled(), on);
    let small = Box::new([7u8; 64]);
    let edge = vec![1u8; 1024];
    let large = vec![1u8; 1025];
    #[repr(align(64))]
    struct Wide(#[allow(dead_code)] u8);
    let wide = Box::new(Wide(1));
    assert_eq!(alloc::in_arena(small.as_ptr()), on);
    assert_eq!(alloc::in_arena(edge.as_ptr()), on);
    assert!(!alloc::in_arena(large.as_ptr()));
    assert!(!alloc::in_arena(&raw const *wide as *const u8));
    assert_eq!(small.as_ptr() as usize % 16, 0);
    // Zeroed blocks are zero in a chunk that was used before.
    for _ in 0..1000 {
        drop(block(3, 512));
        assert!(vec![0u8; 512].iter().all(|&b| b == 0));
    }
}

#[test]
fn blocks_grow_and_shrink_with_their_contents() {
    on();
    // Growth in place (the last block), by moving (another block follows), and out of the arena.
    let mut a: Vec<u8> = Vec::new();
    let mut b: Vec<u8> = Vec::new();
    for i in 0..5000usize {
        a.push(i as u8);
        if i % 3 == 0 {
            b.push((i / 3) as u8);
        }
    }
    check(0, &a);
    check(0, &b);
    assert!(!alloc::in_arena(a.as_ptr()));
    a.truncate(100);
    a.shrink_to_fit();
    b.truncate(10);
    b.shrink_to_fit();
    check(0, &a);
    check(0, &b);
    let mut s = String::new();
    for i in 0..300 {
        s.push_str("ab");
        let t = s.clone();
        assert_eq!(t.len(), 2 * (i + 1));
    }
    assert!(s.bytes().step_by(2).all(|c| c == b'a'));
}

#[test]
fn blocks_survive_churn_and_frees_on_other_threads() {
    on();
    let (tx, rx) = mpsc::sync_channel::<(usize, Vec<u8>)>(1024);
    let consumer = std::thread::spawn(move || {
        let mut n = 0usize;
        // Hold some blocks for a while, so that chunks fill up behind them.
        let mut held: Vec<(usize, Vec<u8>)> = Vec::new();
        for (tag, v) in rx {
            check(tag, &v);
            n += 1;
            if tag % 7 == 0 {
                held.push((tag, v));
                if held.len() > 500 {
                    for (tag, v) in held.drain(..250) {
                        check(tag, &v);
                    }
                }
            }
        }
        for (tag, v) in &held {
            check(*tag, v);
        }
        n
    });
    let workers: Vec<_> = (0..8usize)
        .map(|w| {
            let tx = tx.clone();
            std::thread::spawn(move || {
                let mut rng = Rng(0x9e37_79b9_7f4a_7c15 ^ (w as u64 + 1));
                let mut pool: Vec<(usize, Vec<u8>)> = Vec::new();
                for i in 0..200_000usize {
                    let tag = w * 1_000_003 + i;
                    let len = 1 + rng.next() % if i % 50 == 0 { 3000 } else { 300 };
                    let v = block(tag, len);
                    match rng.next() % 8 {
                        0 => tx.send((tag, v)).unwrap(),
                        1 | 2 => pool.push((tag, v)),
                        _ => check(tag, &v),
                    }
                    if pool.len() > 64 {
                        let k = rng.next() % pool.len();
                        let (tag, mut v) = pool.swap_remove(k);
                        check(tag, &v);
                        // Grow one of the survivors, too.
                        let grown = v.len() + rng.next() % 200;
                        v.extend((v.len()..grown).map(|i| (tag + i) as u8));
                        check(tag, &v);
                    }
                }
                for (tag, v) in &pool {
                    check(*tag, v);
                }
            })
        })
        .collect();
    drop(tx);
    for w in workers {
        w.join().unwrap();
    }
    assert!(consumer.join().unwrap() > 0);
}

#[test]
fn a_thread_that_exits_passes_its_arena_on() {
    let on = on();
    let _arenas = ARENAS.lock().unwrap();
    // More threads, one after the other, than there are arenas: each still gets one, and the
    // blocks that outlive their thread keep their contents while later threads use its arena.
    let mut kept: Vec<(usize, Vec<u8>)> = Vec::new();
    for t in 0..600usize {
        let (in_arena, v) = std::thread::spawn(move || {
            let scratch: Vec<Vec<u8>> = (0..100).map(|i| block(t + i, 48)).collect();
            let v = block(t, 200);
            (
                alloc::in_arena(scratch[99].as_ptr()) && alloc::in_arena(v.as_ptr()),
                v,
            )
        })
        .join()
        .unwrap();
        assert_eq!(in_arena, on, "thread {t}");
        kept.push((t, v));
    }
    for (t, v) in &kept {
        check(*t, v);
    }
}

#[test]
fn more_threads_than_arenas_use_the_wrapped_allocator() {
    on();
    let _arenas = ARENAS.lock().unwrap();
    const THREADS: usize = 300;
    let barrier = Arc::new(Barrier::new(THREADS));
    let handles: Vec<_> = (0..THREADS)
        .map(|t| {
            let barrier = barrier.clone();
            std::thread::spawn(move || {
                let v = block(t, 100);
                // Every thread is alive, and has asked for an arena, at once.
                barrier.wait();
                let more: Vec<Vec<u8>> = (0..50).map(|i| block(t + i, 1 + i * 20)).collect();
                check(t, &v);
                for (i, m) in more.iter().enumerate() {
                    check(t + i, m);
                }
                alloc::in_arena(v.as_ptr())
            })
        })
        .collect();
    let in_arena = handles
        .into_iter()
        .map(|h| h.join().unwrap())
        .filter(|&a| a)
        .count();
    assert!(in_arena <= 256);
    assert!(alloc::arenas_in_use() <= 256);
}

#[test]
fn a_full_arena_uses_the_wrapped_allocator() {
    let on = on();
    let _arenas = ARENAS.lock().unwrap();
    // 16 MB of live small blocks from one thread, whose arena holds 4 MB.
    std::thread::spawn(move || {
        let blocks: Vec<Vec<u8>> = (0..64_000usize).map(|i| block(i, 250)).collect();
        let n = blocks.iter().filter(|b| alloc::in_arena(b.as_ptr())).count();
        if on {
            assert!(n > 10_000 && n < 64_000, "{n} blocks in the arena");
        }
        for (i, b) in blocks.iter().enumerate() {
            check(i, b);
        }
        // Freed in another order than allocated. The thread then looks for a free chunk again
        // within a few thousand allocations, and finds its arena empty.
        let (even, odd): (Vec<_>, Vec<_>) = blocks.into_iter().enumerate().partition(|(i, _)| i % 2 == 0);
        drop(even);
        for (i, b) in &odd {
            check(*i, b);
        }
        drop(odd);
        let again: Vec<Vec<u8>> = (0..20_000usize).map(|i| block(i, 100)).collect();
        assert_eq!(again[5000..].iter().all(|b| alloc::in_arena(b.as_ptr())), on);
        for (i, b) in again.iter().enumerate() {
            check(i, b);
        }
    })
    .join()
    .unwrap();
}
