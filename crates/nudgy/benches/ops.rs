//! Time per allocation of a few patterns, for the wrapped allocator alone, for `Bump` around
//! it before `enable`, and after: `cargo bench -p nudgy`. The patterns call the allocators
//! directly, so the program's own allocations stay out of the arenas. No statistics: each
//! figure is the best of several rounds, in nanoseconds for one allocation and its free.

use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::time::Instant;

use nudgy::Bump;

static BUMP: Bump<System> = Bump::new(System);

fn layout(size: usize) -> Layout {
    Layout::from_size_align(size, 8).unwrap()
}

/// Block sizes of 16 to 527 bytes, the same on every run.
fn sizes(n: usize) -> Vec<usize> {
    let mut x = 0x9e37_79b9_7f4a_7c15u64;
    (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            16 + (x % 512) as usize
        })
        .collect()
}

/// One block at a time: allocated, written, freed.
fn pair<A: GlobalAlloc>(a: &A, rounds: usize) -> usize {
    for _ in 0..rounds {
        // SAFETY: a nonzero size; the block is freed with its layout.
        unsafe {
            let p = a.alloc(layout(64));
            p.write(1);
            a.dealloc(black_box(p), layout(64));
        }
    }
    rounds
}

/// What handling one item of work looks like: `sizes.len()` blocks allocated and written,
/// then freed, in the order `order` gives.
fn batch<A: GlobalAlloc>(a: &A, rounds: usize, sizes: &[usize], order: &[usize]) -> usize {
    let mut held = vec![std::ptr::null_mut::<u8>(); sizes.len()];
    for _ in 0..rounds {
        for (slot, &size) in held.iter_mut().zip(sizes) {
            // SAFETY: a nonzero size.
            unsafe {
                *slot = a.alloc(layout(size));
                slot.write(1);
            }
        }
        for &i in order {
            // SAFETY: allocated above with this layout, freed once.
            unsafe { a.dealloc(black_box(held[i]), layout(sizes[i])) };
        }
    }
    rounds * sizes.len()
}

/// A vector that grows from 8 bytes to 4 KB by doubling, as `Vec::push` does, next to a
/// block that stays: the last block of a chunk grows where it is.
fn grow<A: GlobalAlloc>(a: &A, rounds: usize) -> usize {
    for _ in 0..rounds {
        // SAFETY: nonzero sizes; each block is freed or reallocated with its current layout.
        unsafe {
            let keep = a.alloc(layout(48));
            let (mut p, mut size) = (a.alloc(layout(8)), 8);
            while size < 4096 {
                p.add(size - 1).write(1);
                p = a.realloc(p, layout(size), size * 2);
                size *= 2;
            }
            a.dealloc(black_box(p), layout(size));
            a.dealloc(keep, layout(48));
        }
    }
    rounds * 10
}

/// Blocks allocated on one thread and freed on another.
fn handed_over<A: GlobalAlloc + Sync>(a: &'static A, rounds: usize, sizes: &'static [usize]) -> usize {
    let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<usize>>(4);
    let free = std::thread::spawn(move || {
        for batch in rx {
            for (p, &size) in batch.into_iter().zip(sizes) {
                // SAFETY: allocated by the sender with this layout, freed once.
                unsafe { a.dealloc(std::ptr::with_exposed_provenance_mut(p), layout(size)) };
            }
        }
    });
    for _ in 0..rounds {
        // SAFETY: nonzero sizes.
        let batch = sizes
            .iter()
            .map(|&size| unsafe { a.alloc(layout(size)).expose_provenance() });
        tx.send(batch.collect()).unwrap();
    }
    drop(tx);
    free.join().unwrap();
    rounds * sizes.len()
}

/// Nanoseconds per allocation: the best of five runs of `f`, which returns how many it made.
fn time(mut f: impl FnMut() -> usize) -> f64 {
    (0..5)
        .map(|_| {
            let start = Instant::now();
            let n = f();
            start.elapsed().as_nanos() as f64 / n as f64
        })
        .fold(f64::INFINITY, f64::min)
}

fn run<A: GlobalAlloc + Sync>(a: &'static A) -> Vec<f64> {
    let sizes: &'static [usize] = sizes(300).leak();
    let reverse: Vec<usize> = (0..sizes.len()).rev().collect();
    let forward: Vec<usize> = (0..sizes.len()).collect();
    // Mostly the reverse of the allocations, as nested scopes free, with every fifth block
    // kept to the end.
    let mut nested: Vec<usize> = reverse.iter().copied().filter(|i| i % 5 != 0).collect();
    nested.extend(forward.iter().filter(|&i| i % 5 == 0));
    let threads = std::thread::available_parallelism().map_or(4, |n| n.get().min(8));
    vec![
        time(|| pair(a, 2_000_000)),
        time(|| batch(a, 10_000, sizes, &reverse)),
        time(|| batch(a, 10_000, sizes, &forward)),
        time(|| batch(a, 10_000, sizes, &nested)),
        time(|| grow(a, 200_000)),
        time(|| handed_over(a, 5_000, sizes)),
        time(|| {
            let (nested, n) = (&nested, std::sync::atomic::AtomicUsize::new(0));
            std::thread::scope(|s| {
                for _ in 0..threads {
                    s.spawn(|| {
                        n.fetch_add(
                            batch(a, 10_000, sizes, nested),
                            std::sync::atomic::Ordering::Relaxed,
                        )
                    });
                }
            });
            // Per thread, so that it compares with the row above it.
            n.into_inner() / threads
        }),
    ]
}

fn main() {
    let rows = [
        "one block at a time",
        "300 blocks, freed in reverse",
        "300 blocks, freed in order",
        "300 blocks, freed as nested scopes do",
        "a vector growing to 4 KB",
        "300 blocks, freed on another thread",
        "nested, on several threads at once",
    ];
    // Once unmeasured: the system allocator is fastest on a heap nothing has used yet.
    run(&System);
    let system = run(&System);
    let off = run(&BUMP);
    if !nudgy::enable() {
        println!("no arenas on this target");
        return;
    }
    let on = run(&BUMP);
    println!(
        "{:<40}{:>10}{:>12}{:>10}",
        "ns per allocation", "System", "Bump, off", "Bump, on"
    );
    for (i, row) in rows.iter().enumerate() {
        println!("{row:<40}{:>10.1}{:>12.1}{:>10.1}", system[i], off[i], on[i]);
    }
    println!("{:?}", nudgy::stats());
}
