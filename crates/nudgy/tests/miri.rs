//! Small enough for Miri (`cargo +nightly miri test --test miri`), which checks every access
//! and looks for data races: arena hand-over, frees on other threads, thread-local
//! destructors, and blocks that grow and shrink. Under Miri the crate has fewer and smaller
//! arenas, since Miri backs the whole reservation. It also runs as an ordinary test.

#![cfg(all(any(target_os = "linux", target_os = "macos"), target_pointer_width = "64"))]

use std::alloc::System;
use std::cell::RefCell;
use std::sync::mpsc;

use nudgy::{self as alloc, Bump};

#[global_allocator]
static GLOBAL: Bump<System> = Bump::new(System);

thread_local! {
    static EARLY: RefCell<Vec<Vec<u8>>> = const { RefCell::new(Vec::new()) };
}

fn blk(t: u8, n: usize) -> Vec<u8> {
    let mut v = Vec::with_capacity(n);
    v.extend((0..n).map(|i| t.wrapping_add(i as u8)));
    v
}
fn ok(t: u8, v: &[u8]) {
    assert!(v.iter().enumerate().all(|(i, &b)| b == t.wrapping_add(i as u8)));
}

#[test]
fn miri_protocol() {
    assert!(alloc::enable());
    let (tx, rx) = mpsc::channel::<(u8, Vec<u8>)>();
    let c = std::thread::spawn(move || {
        let mut n = 0;
        for (t, v) in rx {
            ok(t, &v);
            n += 1;
        }
        n
    });
    let mut kept = Vec::new();
    for r in 0..3u8 {
        let tx = tx.clone();
        let v = std::thread::spawn(move || {
            EARLY.with(|_| ());
            let mut mine: Vec<(u8, Vec<u8>)> = (0..20u8)
                .map(|i| (i ^ r, blk(i ^ r, 1 + i as usize * 9)))
                .collect();
            for (t, v) in mine.drain(..5) {
                tx.send((t, v)).unwrap();
            }
            EARLY.with(|e| e.borrow_mut().push(blk(r, 33)));
            let (t, mut g) = mine.pop().unwrap();
            g.extend((g.len()..g.len() + 40).map(|i| t.wrapping_add(i as u8)));
            ok(t, &g);
            g.truncate(3);
            g.shrink_to_fit();
            ok(t, &g);
            // Out of the arena as it grows.
            let n = alloc::MAX_BLOCK - 24;
            let mut big: Vec<u8> = blk(1, n);
            big.extend((n..n + 100).map(|i| 1u8.wrapping_add(i as u8)));
            ok(1, &big);
            assert!(alloc::in_arena(mine[0].1.as_ptr()));
            mine
        })
        .join()
        .unwrap();
        kept.push(v);
    }
    drop(tx);
    assert_eq!(c.join().unwrap(), 15);
    for v in kept.iter().flatten() {
        ok(v.0, &v.1);
    }
}
