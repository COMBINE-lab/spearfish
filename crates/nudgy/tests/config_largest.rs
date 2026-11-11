//! The arenas turned on with more than the largest limit a `Config` takes: the limit is the
//! largest one, half a chunk, and blocks of that size fill a chunk two at a time. `Bump` is
//! called directly, as in `tests/paths.rs`.

use std::alloc::{GlobalAlloc, Layout, System};
use std::thread;

use nudgy::{Bump, Config};

static B: Bump<System> = Bump::new(System);

const L: usize = nudgy::LARGEST_MAX_BLOCK;
const CHUNK: usize = 1 << 16;

fn lay(n: usize) -> Layout {
    Layout::from_size_align(n, 8).unwrap()
}

#[test]
fn blocks_of_the_largest_limit_fill_a_chunk_two_at_a_time() {
    assert_eq!(2 * L, CHUNK);
    if !nudgy::enable_with(Config::new().with_max_block(usize::MAX)) {
        return;
    }
    assert_eq!(nudgy::stats().max_block, L);
    // A thread of its own, which so starts with an arena of empty chunks.
    // SAFETY: every block is freed or reallocated once, with its layout.
    thread::spawn(|| unsafe {
        let a = B.alloc(lay(L));
        let b = B.alloc(lay(L));
        let c = B.alloc_zeroed(lay(L));
        assert!(nudgy::in_arena(a));
        // The first chunk is full after two; the third block starts the next chunk.
        assert_eq!(b as usize, a as usize + L);
        assert_eq!(c as usize, a as usize + CHUNK);
        a.write_bytes(1, L);
        b.write_bytes(2, L);
        assert!((0..L).all(|i| *c.add(i) == 0));
        c.write_bytes(3, L);
        assert!((0..L).all(|i| *a.add(i) == 1 && *b.add(i) == 2));
        let over = B.alloc(lay(L + 1));
        assert!(!nudgy::in_arena(over));
        B.dealloc(over, lay(L + 1));
        B.dealloc(a, lay(L));
        B.dealloc(b, lay(L));
        B.dealloc(c, lay(L));

        // The second chunk, empty again, is filled from its start: a block grows to the limit
        // where it is, and past the limit it moves to the wrapped allocator.
        let p = B.alloc(lay(16));
        assert_eq!(p, c);
        p.write_bytes(9, 16);
        let q = B.realloc(p, lay(16), L);
        assert_eq!(q, p);
        let r = B.realloc(q, lay(L), L + 1);
        assert!(!nudgy::in_arena(r));
        assert!((0..16).all(|i| *r.add(i) == 9));
        B.dealloc(r, lay(L + 1));
    })
    .join()
    .unwrap();
}
