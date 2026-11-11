//! The arenas turned on with a limit of the program's choosing, here one that was asked for
//! as a number that is not a multiple of 16: which blocks an arena serves, how far the last block of a
//! chunk grows where it is, and that the first configuration is the one that holds. `Bump` is
//! called directly, as in `tests/paths.rs`.

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::Mutex;
use std::thread;

use nudgy::{Bump, Config};

static B: Bump<System> = Bump::new(System);

/// The limit this program asks for, and the one that makes: the next multiple of 16.
const ASKED: usize = 100;
const LIMIT: usize = 112;

/// One test at a time: each counts on its thread's arena holding its own blocks only.
static ONE: Mutex<()> = Mutex::new(());

fn lay(n: usize) -> Layout {
    Layout::from_size_align(n, 8).unwrap()
}

/// Whether this target has the arenas, turned on with `LIMIT` by whichever test asks first.
fn on() -> bool {
    nudgy::enable_with(Config::new().with_max_block(ASKED))
}

/// Runs `f` on a thread of its own, which so starts with an arena of empty chunks. Does
/// nothing on a target without the arenas.
fn on_a_new_thread(f: impl FnOnce() + Send + 'static) {
    if !on() {
        return;
    }
    let _one = ONE.lock().unwrap_or_else(|e| e.into_inner());
    thread::spawn(f).join().unwrap();
}

#[test]
fn a_config_keeps_what_it_is_given_up_to_the_largest_limit() {
    assert_eq!(Config::new().max_block(), nudgy::DEFAULT_MAX_BLOCK);
    assert_eq!(Config::default(), Config::new());
    assert_eq!(Config::new().with_max_block(ASKED).max_block(), LIMIT);
    assert_eq!(Config::new().with_max_block(LIMIT).max_block(), LIMIT);
    assert_eq!(Config::new().with_max_block(1).max_block(), 16);
    assert_eq!(Config::new().with_max_block(0).max_block(), 0);
    let largest = nudgy::LARGEST_MAX_BLOCK;
    assert_eq!(Config::new().with_max_block(largest).max_block(), largest);
    assert_eq!(Config::new().with_max_block(largest + 1).max_block(), largest);
    assert_eq!(Config::new().with_max_block(usize::MAX).max_block(), largest);
}

#[test]
fn a_config_is_read_from_a_size_as_a_user_writes_it() {
    let limit = |s: &str| s.parse::<Config>().map(|c| c.max_block());
    for (s, bytes) in [
        ("4096", 4096),
        ("512B", 512),
        ("512b", 512),
        ("8KB", 8192),
        ("8kb", 8192),
        ("8K", 8192),
        ("8KiB", 8192),
        (" 1 KB ", 1024),
        ("100", LIMIT),
        ("0", 0),
        ("0KB", 0),
        ("32KB", nudgy::LARGEST_MAX_BLOCK),
        ("32768", nudgy::LARGEST_MAX_BLOCK),
    ] {
        assert_eq!(limit(s), Ok(bytes), "{s:?}");
    }
    // A size above the largest limit is an error, where `with_max_block` gives that limit.
    for s in [
        "32769",
        "33KB",
        "64KB",
        "99999999999999999999999999",
        "18446744073709551615KB",
    ] {
        let e = limit(s).expect_err(s);
        assert_eq!(e.to_string(), "the largest block is at most 32 KB", "{s:?}");
    }
    for s in [
        "", "KB", "8MB", "8GB", "-1", "+8KB", "1.5KB", "8 K B", "0x10", "off", "8KB8",
    ] {
        let e = limit(s).expect_err(s);
        assert!(e.to_string().starts_with("expected a size"), "{s:?}: {e}");
    }
}

#[test]
fn the_first_configuration_holds() {
    if !on() {
        assert_eq!(nudgy::stats().max_block, 0);
        return;
    }
    assert_eq!(nudgy::stats().max_block, LIMIT);
    assert!(nudgy::enable());
    assert!(nudgy::enable_with(Config::new().with_max_block(4096)));
    assert_eq!(nudgy::stats().max_block, LIMIT);
    // SAFETY: a block freed once, with its layout.
    unsafe {
        let p = B.alloc(lay(4096));
        assert!(!nudgy::in_arena(p));
        B.dealloc(p, lay(4096));
    }
}

#[test]
fn an_arena_serves_blocks_up_to_the_limit() {
    // SAFETY (here and below): every block is freed or reallocated once, with its layout.
    on_a_new_thread(|| unsafe {
        let first = B.alloc(lay(1));
        assert!(nudgy::in_arena(first));
        let at_limit = B.alloc(lay(LIMIT));
        assert_eq!(at_limit as usize, first as usize + 16);
        at_limit.write_bytes(7, LIMIT);
        let next = B.alloc_zeroed(lay(ASKED));
        assert_eq!(next as usize, at_limit as usize + LIMIT);
        assert!((0..LIMIT).all(|i| *at_limit.add(i) == 7));
        for n in [LIMIT + 1, 2 * LIMIT, nudgy::DEFAULT_MAX_BLOCK] {
            let p = B.alloc(lay(n));
            assert!(!nudgy::in_arena(p), "{n} bytes");
            B.dealloc(p, lay(n));
            let p = B.alloc_zeroed(lay(n));
            assert!(!nudgy::in_arena(p), "{n} bytes, zeroed");
            B.dealloc(p, lay(n));
        }
        B.dealloc(next, lay(ASKED));
        B.dealloc(at_limit, lay(LIMIT));
        B.dealloc(first, lay(1));
    });
}

#[test]
fn the_last_block_grows_where_it_is_up_to_the_limit() {
    on_a_new_thread(|| unsafe {
        let keep = B.alloc(lay(16));
        let p = B.alloc(lay(8));
        assert_eq!(p as usize, keep as usize + 16);
        p.write_bytes(3, 8);
        let q = B.realloc(p, lay(8), LIMIT);
        assert_eq!(q, p);
        q.add(8).write_bytes(4, LIMIT - 8);
        // One byte more than the limit: the block moves to the wrapped allocator, and its
        // bytes with it.
        let r = B.realloc(q, lay(LIMIT), LIMIT + 1);
        assert!(!nudgy::in_arena(r));
        assert!((0..8).all(|i| *r.add(i) == 3) && (8..LIMIT).all(|i| *r.add(i) == 4));
        // The arena has the block's place back.
        let again = B.alloc(lay(LIMIT));
        assert_eq!(again, p);
        B.dealloc(again, lay(LIMIT));
        B.dealloc(r, lay(LIMIT + 1));
        B.dealloc(keep, lay(16));
    });
}
