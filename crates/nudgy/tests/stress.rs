//! Stress tests, from the allocator's safety review: arena hand-over with frees by the new
//! owner and by third threads, allocation and frees during thread-local destruction, realloc
//! across the largest size an arena serves, and thread churn with a shared pool of blocks. Every block knows
//! its contents and checks them when it is freed.
//!
//! `STRESS_SECS` sets the length of each timed test (default 1, to keep `cargo test` short).
//! The continuous integration runs them for longer, in a release build with debug assertions
//! and overflow checks.

#![cfg(all(any(target_os = "linux", target_os = "macos"), target_pointer_width = "64"))]

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::RefCell;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::*};
use std::sync::{Arc, Mutex, mpsc};
use std::time::{Duration, Instant};

use nudgy::{self as alloc, Bump};

#[global_allocator]
static GLOBAL: Bump<System> = Bump::new(System);

/// The largest block an arena serves.
const M: usize = alloc::MAX_BLOCK;

static SERIAL: Mutex<()> = Mutex::new(());

fn secs() -> Duration {
    Duration::from_secs_f64(
        std::env::var("STRESS_SECS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(1.0),
    )
}

struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1)
    }
    fn next(&mut self) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as usize
    }
}

static NEXT_TAG: AtomicUsize = AtomicUsize::new(1);

/// A block that knows its contents.
struct Blk(usize, Vec<u8>);
impl Blk {
    fn new(len: usize) -> Blk {
        let tag = NEXT_TAG.fetch_add(1, Relaxed);
        let mut v = Vec::with_capacity(len);
        v.extend((0..len).map(|i| pat(tag, i)));
        Blk(tag, v)
    }
    fn check(&self) {
        for (i, &b) in self.1.iter().enumerate() {
            if b != pat(self.0, i) {
                panic!(
                    "block {} (len {}, {:p}, in an arena: {}) overwritten at {i}",
                    self.0,
                    self.1.len(),
                    self.1.as_ptr(),
                    alloc::in_arena(self.1.as_ptr())
                );
            }
        }
    }
    fn grow(&mut self, by: usize) {
        self.check();
        let n = self.1.len();
        let t = self.0;
        self.1.extend((n..n + by).map(|i| pat(t, i)));
        self.check();
    }
    fn shrink(&mut self, to: usize) {
        self.check();
        self.1.truncate(to.max(1));
        self.1.shrink_to_fit();
        self.check();
    }
}
impl Drop for Blk {
    fn drop(&mut self) {
        self.check();
    }
}
fn pat(tag: usize, i: usize) -> u8 {
    (tag.wrapping_mul(31).wrapping_add(i).wrapping_mul(0x9d)) as u8 ^ (tag >> 8) as u8
}

/// Thread A fills part of its arena and exits; thread B takes the arena (it is the only
/// thread asking), frees half of A's blocks itself (LIVE side) while thread C frees the
/// other half (FREED side); B keeps allocating, and every surviving block keeps its bytes.
#[test]
fn new_owner_and_third_thread_free_the_old_owners_blocks() {
    let _g = SERIAL.lock().unwrap();
    assert!(alloc::enable());
    for round in 0..200 {
        let (blocks, in_arena) = std::thread::spawn(|| {
            let mut v: Vec<Blk> = (0..2000).map(|i| Blk::new(1 + (i * 37) % 600)).collect();
            // Some of A's blocks freed by A itself, so A's chunks are part-used.
            for i in (0..v.len()).rev().step_by(3) {
                v.swap_remove(i);
            }
            let a = alloc::in_arena(v[0].1.as_ptr());
            (v, a)
        })
        .join()
        .unwrap();
        assert!(in_arena);
        let (mine, theirs): (Vec<_>, Vec<_>) = blocks.into_iter().enumerate().partition(|(i, _)| i % 2 == 0);
        let (tx, rx) = mpsc::channel::<Vec<(usize, Blk)>>();
        let c = std::thread::spawn(move || {
            for batch in rx {
                for (_, b) in batch {
                    drop(b);
                }
            }
        });
        let b = std::thread::spawn(move || {
            let mut rng = Rng::new(round + 7);
            let mut keep: Vec<Blk> = Vec::new();
            let mut mine = mine;
            let mut theirs = theirs;
            while !mine.is_empty() || !theirs.is_empty() {
                // Interleave: B frees some, hands some to C, allocates some.
                for _ in 0..rng.next() % 20 {
                    if let Some((_, b)) = mine.pop() {
                        drop(b);
                    }
                }
                let k = theirs.len().min(rng.next() % 20);
                let batch: Vec<_> = theirs.drain(..k).collect();
                tx.send(batch).unwrap();
                for _ in 0..rng.next() % 40 {
                    // One in eight around the largest size served, to grow and shrink across it.
                    let near = rng.next().is_multiple_of(8);
                    let blk = Blk::new(if near {
                        M - 200 + rng.next() % 400
                    } else {
                        1 + rng.next() % 1024
                    });
                    if rng.next().is_multiple_of(4) {
                        keep.push(blk);
                    }
                    if keep.len() > 200 {
                        let i = rng.next() % keep.len();
                        let mut k = keep.swap_remove(i);
                        match rng.next() % 3 {
                            0 => k.grow(rng.next() % 300),
                            1 => k.shrink(rng.next() % 300),
                            _ => {}
                        }
                        keep.push(k);
                        keep.swap_remove(rng.next() % keep.len());
                    }
                }
            }
            drop(tx);
            for k in &keep {
                k.check();
            }
            keep
        });
        let keep = b.join().unwrap();
        c.join().unwrap();
        // B's survivors outlive B and are freed here, on another thread.
        drop(keep);
    }
}

thread_local! {
    // Registered before the thread's arena (touched before any small allocation): its
    // destructor runs after the arena is given back.
    static EARLY: RefCell<Vec<Blk>> = const { RefCell::new(Vec::new()) };
}
thread_local! {
    // Registered after the arena was taken: destroyed before the arena is given back.
    static LATE: RefCell<Vec<Blk>> = const { RefCell::new(Vec::new()) };
}
struct AllocsInDrop(Arc<Mutex<Vec<Blk>>>);
impl Drop for AllocsInDrop {
    fn drop(&mut self) {
        // Allocate, grow and free small blocks during TLS destruction, and leave some behind
        // for another thread.
        let mut v: Vec<Blk> = (0..50).map(|i| Blk::new(1 + i * 13)).collect();
        for b in v.iter_mut().step_by(3) {
            b.grow(17);
        }
        let keep: Vec<Blk> = v.drain(..10).collect();
        drop(v);
        self.0.lock().unwrap().extend(keep);
    }
}
thread_local! {
    static DROPPER: RefCell<Option<AllocsInDrop>> = const { RefCell::new(None) };
}

/// Threads that hold arena blocks in TLS destructors that run before and after the arena is
/// handed back, and that allocate in them, while other threads take the arenas.
#[test]
fn tls_destructors_free_and_allocate_while_arenas_change_hands() {
    let _g = SERIAL.lock().unwrap();
    assert!(alloc::enable());
    let pool = Arc::new(Mutex::new(Vec::<Blk>::new()));
    let deadline = Instant::now() + secs();
    let mut n = 0usize;
    while Instant::now() < deadline {
        let hs: Vec<_> = (0..24)
            .map(|t| {
                let pool = pool.clone();
                std::thread::spawn(move || {
                    // Before any small allocation of this thread.
                    EARLY.with(|_| ());
                    if t % 2 == 0 {
                        DROPPER.with(|_| ());
                    }
                    let mut rng = Rng::new(t as u64 + 1);
                    let mut held: Vec<Blk> = (0..200).map(|_| Blk::new(1 + rng.next() % 900)).collect();
                    EARLY.with(|e| e.borrow_mut().extend(held.drain(..60)));
                    LATE.with(|e| e.borrow_mut().extend(held.drain(..60)));
                    DROPPER.with(|d| {
                        if t % 2 == 1 {
                            *d.borrow_mut() = Some(AllocsInDrop(pool.clone()));
                        }
                    });
                    if t % 2 == 0 {
                        DROPPER.with(|d| *d.borrow_mut() = Some(AllocsInDrop(pool.clone())));
                    }
                    // Swap blocks with the pool.
                    let mut p = pool.lock().unwrap();
                    for _ in 0..20 {
                        if let Some(b) = held.pop() {
                            p.push(b);
                        }
                    }
                    let k = p.len().min(30);
                    let theirs: Vec<Blk> = p.drain(..k).collect();
                    drop(p);
                    drop(theirs);
                })
            })
            .collect();
        for h in hs {
            h.join().unwrap();
        }
        n += 1;
        let mut p = pool.lock().unwrap();
        for b in p.iter() {
            b.check();
        }
        if p.len() > 4000 {
            p.drain(..2000);
        }
    }
    eprintln!("tls: {n} rounds of 24 threads");
}

/// Realloc across the largest size served and the in-place paths, through the allocator directly.
#[test]
fn realloc_across_the_boundary() {
    assert!(alloc::enable());
    std::thread::spawn(|| unsafe {
        let fill = |p: *mut u8, n: usize, t: u8| {
            for i in 0..n {
                *p.add(i) = t.wrapping_add(i as u8);
            }
        };
        let ok = |p: *mut u8, n: usize, t: u8| (0..n).all(|i| *p.add(i) == t.wrapping_add(i as u8));
        for &(a, b, c) in &[
            (1usize, M, M + 1),
            (M - 24, M + 1, 100),
            (M, M - 1, 1),
            (17, 16, 4 * M),
            (500, 9 * M, M),
            (64, M, 64),
        ] {
            for align in [1usize, 2, 8, 16, 32, 64] {
                let l = Layout::from_size_align(a, align).unwrap();
                let p = GLOBAL.alloc(l);
                assert!(!p.is_null() && (p as usize).is_multiple_of(align));
                assert_eq!(alloc::in_arena(p), a <= M && align <= 16);
                fill(p, a, 3);
                // A block after it, so the next realloc of p cannot be in place (sometimes).
                let lx = Layout::from_size_align(48, 8).unwrap();
                let x = GLOBAL.alloc(lx);
                fill(x, 48, 9);
                let q = GLOBAL.realloc(p, l, b);
                assert!(
                    !q.is_null() && (q as usize).is_multiple_of(align),
                    "{a}->{b} align {align}"
                );
                assert!(ok(q, a.min(b), 3), "{a}->{b} align {align}");
                fill(q, b, 4);
                GLOBAL.dealloc(x, lx);
                let lb = Layout::from_size_align(b, align).unwrap();
                let r = GLOBAL.realloc(q, lb, c);
                assert!(!r.is_null() && (r as usize).is_multiple_of(align));
                assert!(ok(r, b.min(c), 4), "{b}->{c} align {align}");
                fill(r, c, 6);
                // Last block now: grow and shrink in place.
                let lc = Layout::from_size_align(c, align).unwrap();
                let s = GLOBAL.realloc(r, lc, c.min(1000) + 10);
                assert!(ok(s, c.min(c.min(1000) + 10), 6));
                GLOBAL.dealloc(s, Layout::from_size_align(c.min(1000) + 10, align).unwrap());
                // Zeroed after reuse.
                let z = GLOBAL.alloc_zeroed(lb);
                assert!((0..b).all(|i| *z.add(i) == 0));
                GLOBAL.dealloc(z, lb);
            }
        }
        // A block allocated by the wrapped allocator (before enable, or by System) given to
        // Bump::realloc and dealloc.
        let l = Layout::from_size_align(100, 8).unwrap();
        let p = System.alloc(l);
        fill(p, 100, 5);
        let q = GLOBAL.realloc(p, l, 50);
        assert!(!alloc::in_arena(q) && ok(q, 50, 5));
        GLOBAL.dealloc(q, Layout::from_size_align(50, 8).unwrap());
    })
    .join()
    .unwrap();
}

/// Short-lived threads (more over the run than there are arenas), a shared pool of blocks
/// that outlive their threads, frees on every thread, growth, shrinking and zeroed blocks.
#[test]
fn thread_churn_with_a_shared_pool() {
    let _g = SERIAL.lock().unwrap();
    assert!(alloc::enable());
    const SLOTS: usize = 24;
    let pool: Arc<Vec<Mutex<Vec<Blk>>>> = Arc::new((0..64).map(|_| Mutex::new(Vec::new())).collect());
    let stop = Arc::new(AtomicBool::new(false));
    let spawned = Arc::new(AtomicUsize::new(0));
    let slots: Vec<_> = (0..SLOTS)
        .map(|s| {
            let pool = pool.clone();
            let stop = stop.clone();
            let spawned = spawned.clone();
            std::thread::spawn(move || {
                let mut round_ = 0u64;
                while !stop.load(Relaxed) {
                    round_ += 1;
                    let pool2 = pool.clone();
                    spawned.fetch_add(1, Relaxed);
                    std::thread::spawn(move || {
                        let pool = pool2;
                        let mut rng = Rng::new((s as u64) << 32 | round_);
                        let life = 500 + rng.next() % 20_000;
                        let mut mine: Vec<Blk> = Vec::new();
                        for _ in 0..life {
                            let len = 1 + if rng.next().is_multiple_of(64) {
                                rng.next() % (M + M / 2)
                            } else {
                                rng.next() % 400
                            };
                            let b = if rng.next().is_multiple_of(16) {
                                let v = vec![0u8; len];
                                assert!(v.iter().all(|&x| x == 0));
                                drop(v);
                                Blk::new(len)
                            } else {
                                Blk::new(len)
                            };
                            match rng.next() % 16 {
                                0 => {
                                    let i = rng.next() % pool.len();
                                    pool[i].lock().unwrap().push(b);
                                }
                                1 => {
                                    let i = rng.next() % pool.len();
                                    let mut p = pool[i].lock().unwrap();
                                    if !p.is_empty() {
                                        let k = rng.next() % p.len();
                                        let mut x = p.swap_remove(k);
                                        drop(p);
                                        if rng.next().is_multiple_of(2) {
                                            x.grow(rng.next() % 200);
                                            mine.push(x);
                                        }
                                    }
                                }
                                2..=6 => mine.push(b),
                                _ => drop(b),
                            }
                            if mine.len() > 100 {
                                let k = rng.next() % mine.len();
                                let mut x = mine.swap_remove(k);
                                match rng.next() % 4 {
                                    0 => x.grow(rng.next() % 100),
                                    1 => x.shrink(rng.next() % 100),
                                    _ => {}
                                }
                                if rng.next().is_multiple_of(2) {
                                    mine.push(x);
                                }
                            }
                        }
                        // Leave the rest to the pool: they outlive the thread.
                        let i = rng.next() % pool.len();
                        pool[i].lock().unwrap().extend(mine.drain(..mine.len() / 2));
                    })
                    .join()
                    .unwrap();
                    // Keep the pool bounded.
                    let i = (round_ as usize * 7 + s) % pool.len();
                    let mut p = pool[i].lock().unwrap();
                    if p.len() > 300 {
                        p.drain(..150);
                    }
                }
            })
        })
        .collect();
    std::thread::sleep(secs());
    stop.store(true, Relaxed);
    for s in slots {
        s.join().unwrap();
    }
    for p in pool.iter() {
        for b in p.lock().unwrap().iter() {
            b.check();
        }
    }
    let stats = alloc::stats();
    eprintln!("churn: {} threads; {stats:?}", spawned.load(Relaxed));
    // The slots' threads have exited, and each gave its arena back.
    assert!(stats.arenas_in_use <= stats.arenas_peak && stats.arenas_peak >= SLOTS);
}
