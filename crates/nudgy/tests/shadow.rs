//! Random allocations, reallocations and frees through a `Bump` that is called directly (the
//! test program's own allocator is the system's), checked against a record of every live
//! block: a block that overlaps a live one fails the test, and so does one whose contents
//! changed by the time it is freed or reallocated. Blocks pass between threads, outlive the
//! thread that made them, and are made and freed in thread-local destructors that run before
//! and after the thread's arena is given back.
//!
//! `GENS` sets how many generations of threads run (8 in `cargo test`).

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

use nudgy::Bump;

static B: Bump<System> = Bump::new(System);
const M: usize = nudgy::MAX_BLOCK;

/// Address to size and tag, of every live block, from an arena or not.
static LIVE: Mutex<BTreeMap<usize, (usize, usize)>> = Mutex::new(BTreeMap::new());
/// Blocks any thread may free or reallocate.
static POOL: Mutex<Vec<Blk>> = Mutex::new(Vec::new());
static TAG: AtomicUsize = AtomicUsize::new(1);
static IN_ARENA: AtomicUsize = AtomicUsize::new(0);
static TOTAL: AtomicUsize = AtomicUsize::new(0);

/// The bytes a block of `size` at `p` takes: an arena rounds a size up to 16.
fn ext(p: usize, size: usize) -> usize {
    if nudgy::in_arena(p as *const u8) {
        (size + 15) & !15
    } else {
        size
    }
}

fn insert(p: usize, size: usize, tag: usize) {
    let e = ext(p, size);
    let mut m = LIVE.lock().unwrap();
    if let Some((&q, &(qs, qt))) = m.range(..p + e).next_back() {
        let qe = ext(q, qs);
        assert!(
            q + qe <= p,
            "overlap: new {p:#x}+{e} (tag {tag}) with live {q:#x}+{qe} (tag {qt})"
        );
    }
    m.insert(p, (size, tag));
}
fn remove(p: usize) {
    LIVE.lock().unwrap().remove(&p).expect("freeing a block not live");
}

/// Byte `i` of the block with this tag.
fn pat(tag: usize, i: usize) -> u8 {
    (tag.wrapping_mul(0x9e37).wrapping_add(i * 7)) as u8 ^ (tag >> 9) as u8
}

struct Blk {
    p: *mut u8,
    size: usize,
    align: usize,
    tag: usize,
}
// SAFETY: a block of the allocator, which any one thread at a time may use and free.
unsafe impl Send for Blk {}

impl Blk {
    fn new(size: usize, align: usize, zeroed: bool) -> Blk {
        let l = Layout::from_size_align(size, align).unwrap();
        // SAFETY: `size` is not zero.
        let p = unsafe { if zeroed { B.alloc_zeroed(l) } else { B.alloc(l) } };
        assert!(!p.is_null());
        assert_eq!(p as usize % align, 0);
        if zeroed {
            // SAFETY: within the block.
            assert!(
                (0..size).all(|i| unsafe { *p.add(i) } == 0),
                "zeroed block not zero"
            );
        }
        let tag = TAG.fetch_add(1, Relaxed);
        TOTAL.fetch_add(1, Relaxed);
        if nudgy::in_arena(p) {
            IN_ARENA.fetch_add(1, Relaxed);
        }
        insert(p as usize, size, tag);
        for i in 0..size {
            // SAFETY: within the block.
            unsafe { *p.add(i) = pat(tag, i) };
        }
        Blk { p, size, align, tag }
    }
    fn check(&self, n: usize) {
        for i in 0..n {
            // SAFETY: within the block.
            let b = unsafe { *self.p.add(i) };
            assert_eq!(
                b,
                pat(self.tag, i),
                "block tag {} {:p} size {} corrupted at {i}",
                self.tag,
                self.p,
                self.size
            );
        }
    }
    fn realloc(&mut self, new: usize) {
        self.check(self.size);
        remove(self.p as usize);
        let l = Layout::from_size_align(self.size, self.align).unwrap();
        // SAFETY: the block is live with layout `l`, and `new` is not zero.
        let q = unsafe { B.realloc(self.p, l, new) };
        assert!(!q.is_null());
        assert_eq!(q as usize % self.align, 0);
        let keep = self.size.min(new);
        self.p = q;
        self.check(keep);
        insert(q as usize, new, self.tag);
        for i in keep..new {
            // SAFETY: within the block.
            unsafe { *q.add(i) = pat(self.tag, i) };
        }
        self.size = new;
    }
}
impl Drop for Blk {
    fn drop(&mut self) {
        self.check(self.size);
        remove(self.p as usize);
        // SAFETY: the block is live with this layout, and freed once.
        unsafe { B.dealloc(self.p, Layout::from_size_align(self.size, self.align).unwrap()) };
    }
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as usize
    }
}

/// Mostly small, some near the largest an arena serves, some past it.
fn size(r: &mut Rng) -> usize {
    match r.next() % 10 {
        0 => 1 + r.next() % (M + 600),
        1 => M - r.next() % 64,
        2 => 1 + r.next() % 16,
        _ => 1 + r.next() % 500,
    }
}

thread_local! {
    static EARLY: RefCell<Vec<Blk>> = const { RefCell::new(Vec::new()) };
    static LATE: RefCell<Vec<Blk>> = const { RefCell::new(Vec::new()) };
}

/// Allocates, reallocates and frees when its thread's destructors run.
struct Churn(u64);
impl Drop for Churn {
    fn drop(&mut self) {
        let mut r = Rng(self.0 | 1);
        let mut v: Vec<Blk> = (0..20).map(|_| Blk::new(size(&mut r), 8, false)).collect();
        for b in v.iter_mut().step_by(2) {
            let n = size(&mut r);
            b.realloc(n);
        }
        POOL.lock().unwrap().extend(v.drain(..5));
    }
}
thread_local! {
    static CHURN: RefCell<Option<Churn>> = const { RefCell::new(None) };
}

fn worker(seed: u64, ops: usize) {
    let mut r = Rng(seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1);
    if seed.is_multiple_of(2) {
        // Registered before the allocator's own destructor: these run after the arena is
        // given back.
        EARLY.with(|_| ());
        CHURN.with(|c| *c.borrow_mut() = Some(Churn(seed)));
    }
    let mut mine: Vec<Blk> = Vec::new();
    for _ in 0..ops {
        match r.next() % 16 {
            0..=5 => {
                let al = [1, 2, 8, 16, 16, 32][r.next() % 6];
                let z = r.next().is_multiple_of(8);
                mine.push(Blk::new(size(&mut r), al, z));
            }
            6..=8 => {
                if !mine.is_empty() {
                    // Often the last one, which gives its bytes back.
                    let k = if r.next().is_multiple_of(2) {
                        mine.len() - 1
                    } else {
                        r.next() % mine.len()
                    };
                    drop(mine.swap_remove(k));
                }
            }
            9..=11 => {
                if !mine.is_empty() {
                    let k = if r.next().is_multiple_of(2) {
                        mine.len() - 1
                    } else {
                        r.next() % mine.len()
                    };
                    let b = &mut mine[k];
                    let n = match r.next() % 5 {
                        0 => b.size + 1 + r.next() % 64,
                        1 => (b.size / 2).max(1),
                        2 => ((b.size + 15) & !15).max(1),
                        3 => M + r.next() % 3,
                        _ => size(&mut r),
                    };
                    b.realloc(n);
                }
            }
            12 => {
                if let Some(b) = mine.pop() {
                    POOL.lock().unwrap().push(b);
                }
            }
            13 => {
                let got = {
                    let mut p = POOL.lock().unwrap();
                    if p.is_empty() {
                        None
                    } else {
                        let k = r.next() % p.len();
                        Some(p.swap_remove(k))
                    }
                };
                if let Some(mut b) = got {
                    match r.next() % 3 {
                        0 => drop(b),
                        1 => {
                            let n = size(&mut r);
                            b.realloc(n);
                            mine.push(b)
                        }
                        _ => mine.push(b),
                    }
                }
            }
            14 => {
                if let Some(b) = mine.pop() {
                    if r.next().is_multiple_of(2) {
                        EARLY.with(|e| e.borrow_mut().push(b))
                    } else {
                        LATE.with(|e| e.borrow_mut().push(b))
                    }
                }
            }
            _ => {
                // A run of frees, newest first, to empty chunks.
                let k = mine.len().min(r.next() % 40);
                let tail = mine.len() - k;
                for b in mine.drain(tail..).rev() {
                    drop(b);
                }
            }
        }
    }
    if seed % 2 == 1 {
        // Registered after the allocator's destructor: runs while the thread has its arena.
        CHURN.with(|c| *c.borrow_mut() = Some(Churn(seed)));
    }
    let half = mine.len() / 2;
    POOL.lock().unwrap().extend(mine.drain(..half));
}

#[test]
fn no_block_overlaps_a_live_one() {
    let on = nudgy::enable();
    let (gens, threads, ops) = if cfg!(miri) {
        (3, 3, 150)
    } else {
        let gens = std::env::var("GENS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(8);
        (gens, 6, 20_000)
    };
    for g in 0..gens {
        let hs: Vec<_> = (0..threads)
            .map(|t| std::thread::spawn(move || worker((g * 100 + t + 1) as u64, ops)))
            .collect();
        for h in hs {
            h.join().unwrap();
        }
        let mut p = POOL.lock().unwrap();
        if p.len() > 3000 {
            let k = p.len() - 1500;
            let old: Vec<Blk> = p.drain(..k).collect();
            drop(p);
            drop(old);
        }
    }
    drop(std::mem::take(&mut *POOL.lock().unwrap()));
    assert!(LIVE.lock().unwrap().is_empty());
    // Or the test checked the system allocator.
    assert!(!on || IN_ARENA.load(Relaxed) * 2 > TOTAL.load(Relaxed));
}
