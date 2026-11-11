//! `Bump` as this test program's global allocator, with the arenas on: blocks keep
//! their contents whichever thread frees them, through thread exits and when arenas or chunks
//! run out.

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::RefCell;
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
use std::sync::{Arc, Barrier, Mutex, mpsc};

use nudgy::{self as alloc, Bump};

#[global_allocator]
static GLOBAL: Bump<System> = Bump::new(System);

/// The largest block an arena serves.
const M: usize = alloc::DEFAULT_MAX_BLOCK;

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
    let edge = vec![1u8; M];
    let large = vec![1u8; M + 1];
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
    for i in 0..2 * M {
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
                    let len = 1 + rng.next() % if i % 50 == 0 { M + M / 2 } else { 300 };
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
    let stats = alloc::stats();
    assert!(stats.arenas_in_use <= 256);
    if stats.enabled {
        // Every arena was taken, and the threads that found none were counted.
        assert_eq!(stats.arenas_peak, 256);
        assert!(stats.no_arena >= (THREADS - 256) as u64, "{stats:?}");
        assert!(stats.fallback_allocations() >= stats.no_arena);
    }
}

#[test]
fn a_full_arena_uses_the_wrapped_allocator() {
    let on = on();
    let _arenas = ARENAS.lock().unwrap();
    // 16 MB of live small blocks from one thread, whose arena holds 4 MB.
    std::thread::spawn(move || {
        let before = alloc::stats();
        let blocks: Vec<Vec<u8>> = (0..64_000usize).map(|i| block(i, 250)).collect();
        let n = blocks.iter().filter(|b| alloc::in_arena(b.as_ptr())).count();
        if on {
            assert!(n > 10_000 && n < 64_000, "{n} blocks in the arena");
            // The thread went through its 64 chunks, then found its arena full.
            let after = alloc::stats();
            assert!(
                after.chunk_refills >= before.chunk_refills + 64,
                "{before:?} {after:?}"
            );
            assert!(after.arena_full > before.arena_full, "{before:?} {after:?}");
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

#[test]
fn the_next_owner_and_other_threads_free_what_an_arena_held() {
    on();
    let _arenas = ARENAS.lock().unwrap();
    for round in 0..50usize {
        // A thread leaves part-used chunks behind. The next thread takes its arena, frees
        // half of the blocks itself and passes the rest to a third thread, allocating all
        // the while.
        let left: Vec<(usize, Vec<u8>)> = std::thread::spawn(move || {
            let mut v: Vec<_> = (0..2000usize)
                .map(|i| (i, block(i, 1 + (i * 37) % 600)))
                .collect();
            for i in (0..v.len()).rev().step_by(3) {
                v.swap_remove(i);
            }
            v
        })
        .join()
        .unwrap();
        let (tx, rx) = mpsc::channel::<(usize, Vec<u8>)>();
        let third = std::thread::spawn(move || {
            for (tag, v) in rx {
                check(tag, &v);
            }
        });
        let kept = std::thread::spawn(move || {
            let mut rng = Rng(0x9e37_79b9_7f4a_7c15 ^ (round as u64 + 1));
            let mut kept: Vec<(usize, Vec<u8>)> = Vec::new();
            for (n, (tag, v)) in left.into_iter().enumerate() {
                if n % 2 == 0 {
                    check(tag, &v);
                } else {
                    tx.send((tag, v)).unwrap();
                }
                for i in 0..rng.next() % 4 {
                    let tag = 1_000_000 + n * 4 + i;
                    // One in eight just under the largest size served, to grow past it.
                    let near = rng.next().is_multiple_of(8);
                    let len = if near {
                        M - rng.next() % 200
                    } else {
                        1 + rng.next() % 1024
                    };
                    let mut v = block(tag, len);
                    let grown = v.len() + rng.next() % 300;
                    v.extend((v.len()..grown).map(|i| (tag + i) as u8));
                    if rng.next().is_multiple_of(4) {
                        kept.push((tag, v));
                    }
                }
            }
            kept
        })
        .join()
        .unwrap();
        third.join().unwrap();
        for (tag, v) in &kept {
            check(*tag, v);
        }
    }
}

/// Blocks that a thread's destructors free, and one that allocates as it runs.
struct AtExit(Vec<(usize, Vec<u8>)>, mpsc::Sender<(usize, Vec<u8>)>);
impl Drop for AtExit {
    fn drop(&mut self) {
        let mut more: Vec<_> = (0..50usize).map(|i| (i, block(i, 1 + i * 13))).collect();
        for (tag, v) in more.iter_mut().step_by(3) {
            v.extend((v.len()..v.len() + 17).map(|i| (*tag + i) as u8));
        }
        for (tag, v) in self.0.iter().chain(&more) {
            check(*tag, v);
        }
        // Some outlive the thread.
        for b in more.drain(..10) {
            self.1.send(b).unwrap();
        }
    }
}
thread_local! {
    static BEFORE: RefCell<Option<AtExit>> = const { RefCell::new(None) };
    static AFTER: RefCell<Option<AtExit>> = const { RefCell::new(None) };
}

#[test]
fn thread_local_destructors_allocate_and_free() {
    on();
    let _arenas = ARENAS.lock().unwrap();
    let (tx, rx) = mpsc::channel();
    for t in 0..200usize {
        let tx = tx.clone();
        std::thread::spawn(move || {
            // One destructor is registered before the thread has an arena, and so runs after
            // the arena was given back; the other one before.
            BEFORE.with(|_| ());
            let held = |k: usize| {
                (0..100usize)
                    .map(|i| (t + k + i, block(t + k + i, 1 + i * 7)))
                    .collect()
            };
            let (a, b): (Vec<_>, Vec<_>) = (held(0), held(1000));
            BEFORE.with(|c| *c.borrow_mut() = Some(AtExit(a, tx.clone())));
            AFTER.with(|c| *c.borrow_mut() = Some(AtExit(b, tx)));
        })
        .join()
        .unwrap();
        // The next thread has the arena while the last one's blocks are still here.
        if t % 10 == 9 {
            for (tag, v) in rx.try_iter() {
                check(tag, &v);
            }
        }
    }
}

/// An allocator that counts what reaches it, and answers a request for no bytes (which
/// `System` need not) with an address no block has.
struct Probe;
static PROBED: AtomicUsize = AtomicUsize::new(0);
unsafe impl GlobalAlloc for Probe {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        PROBED.fetch_add(1, Relaxed);
        if l.size() == 0 {
            return l.align() as *mut u8;
        }
        unsafe { System.alloc(l) }
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        PROBED.fetch_add(1, Relaxed);
        if l.size() != 0 {
            unsafe { System.dealloc(p, l) }
        }
    }
}

#[test]
fn direct_calls_keep_to_the_contract() {
    let on = on();
    // One thread, so that the blocks follow one another in its chunk.
    std::thread::spawn(move || unsafe {
        let a = Bump::new(Probe);
        let fill = |p: *mut u8, n: usize| (0..n).for_each(|i| *p.add(i) = (i * 7) as u8);
        let filled = |p: *mut u8, n: usize| (0..n).all(|i| *p.add(i) == (i * 7) as u8);
        // Every alignment, sizes around the largest served.
        for align in [1usize, 2, 4, 8, 16, 32, 64] {
            for size in [1usize, 15, 16, 17, M - 1, M, M + 1, 5 * M] {
                let l = Layout::from_size_align(size, align).unwrap();
                let mut p = a.alloc_zeroed(l);
                assert_eq!(p as usize % align, 0);
                assert!((0..size).all(|i| *p.add(i) == 0));
                assert_eq!(alloc::in_arena(p), on && size <= M && align <= 16);
                fill(p, size);
                // To the other side of the largest size and back, growing and shrinking. A
                // block that left the arenas stays the wrapped allocator's.
                for new in [size + 1, M, M + 1, 3 * M, 700, 1, size] {
                    let was = alloc::in_arena(p);
                    let q = a.realloc(p, l, new);
                    assert_eq!(q as usize % align, 0);
                    assert!(filled(q, size.min(new)));
                    assert_eq!(alloc::in_arena(q), was && new <= M);
                    p = a.realloc(q, Layout::from_size_align(new, align).unwrap(), size);
                    assert_eq!(alloc::in_arena(p), was && new <= M);
                    assert!(filled(p, size.min(new)));
                    fill(p, size);
                }
                a.dealloc(p, l);
            }
        }
        // A request for no bytes is the wrapped allocator's, also when the thread's chunk is
        // full to its last byte, as it is three times here.
        let (word, none) = (Layout::new::<[u8; 16]>(), Layout::from_size_align(0, 1).unwrap());
        let mut blocks: Vec<*mut u8> = Vec::with_capacity(3 * 4096);
        for _ in 0..3 * 4096 {
            blocks.push(a.alloc(word));
            let before = PROBED.load(Relaxed);
            let z = a.alloc(none);
            assert!(!alloc::in_arena(z));
            a.dealloc(z, none);
            assert_eq!(PROBED.load(Relaxed), before + 2);
        }
        assert_eq!(blocks.iter().all(|&p| alloc::in_arena(p)), on);
        for p in blocks {
            a.dealloc(p, word);
        }
    })
    .join()
    .unwrap();
}
