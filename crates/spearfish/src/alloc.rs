//! An optional allocator for mapping: small blocks from per-thread arenas.
//!
//! Mapping a read makes a few hundred small allocations and frees nearly all of them before
//! the next read. [`Bump`] wraps a program's global allocator and, once [`enable`] has been
//! called, serves blocks of up to 1 KB from memory that belongs to the allocating thread: a
//! block is the next free bytes of the thread's current 64 KB chunk, and a chunk whose blocks
//! have all been freed is used again from its start. A thread mapping reads then stays in the
//! same few pages, which is where the gain comes from: mapping takes about 2% less CPU time,
//! for memory of at most 4 MB a thread. Every other allocation, and every allocation before
//! `enable`, goes to the wrapped allocator, at the cost of one comparison. What is mapped
//! does not change.
//!
//! A program opts in twice: it declares the allocator, and it turns the arenas on.
//!
//! ```no_run
//! use spearfish::alloc::Bump;
//!
//! #[global_allocator]
//! static GLOBAL: Bump<std::alloc::System> = Bump::new(std::alloc::System);
//!
//! fn main() {
//!     spearfish::alloc::enable();
//!     // ... build a Mapper and map reads on any threads ...
//! }
//! ```
//!
//! The wrapped allocator can be any [`GlobalAlloc`] (`spearfish map` wraps mimalloc). A block
//! may be freed by a thread other than the one that allocated it, and a thread's arena passes
//! to a later thread when it exits. Up to 256 threads have an arena at a time, and a thread
//! without one, or whose 64 chunks all hold live blocks, uses the wrapped allocator. Blocks
//! that stay allocated keep their chunk from being reused, so the arenas suit programs whose
//! small allocations are mostly short-lived, as a mapper's are.
//!
//! The arenas are one reservation of 1 GB of address space, of which only the chunks in use
//! are backed by memory. They are there on 64-bit Linux and macOS; on other targets `enable`
//! returns false and `Bump` only forwards.

use std::alloc::{GlobalAlloc, Layout};
use std::cell::Cell;
use std::sync::atomic::Ordering::{Acquire, Relaxed, Release};
use std::sync::atomic::{AtomicBool, AtomicU32, AtomicU64, AtomicUsize};

/// Largest block served from an arena, and the alignment of every such block.
const SMALL: usize = 1024;
const ALIGN: usize = 16;
const CHUNK_BITS: usize = 16;
const CHUNK: usize = 1 << CHUNK_BITS;
/// Chunks in a thread's arena.
const ARENA_BITS: usize = 6;
const ARENA_CHUNKS: usize = 1 << ARENA_BITS;
/// Arenas, so threads that can have one at a time.
const ARENAS: usize = 256;
const CHUNKS: usize = ARENAS * ARENA_CHUNKS;
const REGION: usize = CHUNKS * CHUNK;
/// Small allocations a thread passes to the wrapped allocator after it found no free chunk (or
/// no free arena) before it looks again.
const BACKOFF: u32 = 4096;
/// `Local::skip` of a thread that never looks again.
const NEVER: u32 = u32::MAX;
const NO_ARENA: usize = usize::MAX;

/// Targets with the arenas: their thread-local storage is the linker's, so reading `LOCAL`
/// inside the allocator never allocates, and no address is within `REGION` of the top of
/// their address space.
const SUPPORTED: bool = cfg!(all(
    any(target_os = "linux", target_os = "macos"),
    target_pointer_width = "64"
));

static ON: AtomicBool = AtomicBool::new(false);
/// The largest size served from an arena now: zero (no allocation has it) until `enable`.
/// The one load an allocation costs while the arenas are off.
static SMALL_NOW: AtomicUsize = AtomicUsize::new(0);
/// The reservation is `BASE..BASE + REGION`. Until `enable` reserved it, the top of the
/// address space, where no block is: the one load a free costs while the arenas are off.
static BASE: AtomicUsize = AtomicUsize::new(0usize.wrapping_sub(REGION));

/// Bytes from the reservation's start to `p`, when `p` is in it.
#[inline(always)]
fn offset(base: usize, p: *const u8) -> Option<usize> {
    let off = (p as usize).wrapping_sub(base);
    (SUPPORTED && off < REGION).then_some(off)
}

/// One bit per arena: a thread has it.
static TAKEN: [AtomicU64; ARENAS / 64] = [const { AtomicU64::new(0) }; ARENAS / 64];

#[repr(align(64))]
struct PerChunk([AtomicU32; CHUNKS]);
/// Blocks handed out of each chunk minus those its arena's thread freed. Only that thread
/// reads or writes it.
static LIVE: PerChunk = PerChunk([const { AtomicU32::new(0) }; CHUNKS]);
/// Blocks of each chunk that other threads freed. `LIVE == FREED`: the chunk holds no block.
static FREED: PerChunk = PerChunk([const { AtomicU32::new(0) }; CHUNKS]);

/// A thread's arena and its place in the current chunk.
#[derive(Clone, Copy)]
struct Local {
    /// The next free byte of the current chunk, and the chunk's end.
    ptr: usize,
    end: usize,
    chunk: usize,
    arena: usize,
    /// Small allocations still to pass on before `refill` looks for a chunk again.
    skip: u32,
}

/// Gives the thread's arena back when the thread exits.
struct Exit;

thread_local! {
    // No destructor, so it stays readable while the thread's other destructors run.
    static LOCAL: Cell<Local> = const {
        Cell::new(Local { ptr: 0, end: 0, chunk: 0, arena: NO_ARENA, skip: 0 })
    };
    static EXIT: Exit = const { Exit };
}

impl Drop for Exit {
    fn drop(&mut self) {
        let _ = LOCAL.try_with(|cell| {
            let s = cell.get();
            // Destructors that run after this one allocate from the wrapped allocator, and
            // free this arena's blocks as any other thread would.
            cell.set(Local {
                ptr: 0,
                end: 0,
                chunk: 0,
                arena: NO_ARENA,
                skip: NEVER,
            });
            if s.arena != NO_ARENA {
                give_back(s.arena);
            }
        });
    }
}

/// The first arena from `from` on that no thread has.
fn take_arena(from: usize) -> Option<usize> {
    for (w, word) in TAKEN.iter().enumerate().skip(from / 64) {
        let below = if w == from / 64 {
            (1u64 << (from % 64)) - 1
        } else {
            0
        };
        let mut cur = word.load(Relaxed);
        while cur | below != u64::MAX {
            let bit = (!(cur | below)).trailing_zeros();
            // Acquire: the counters the arena's last thread left in LIVE are visible.
            match word.compare_exchange_weak(cur, cur | 1 << bit, Acquire, Relaxed) {
                Ok(_) => return Some(w * 64 + bit as usize),
                Err(now) => cur = now,
            }
        }
    }
    None
}

fn give_back(arena: usize) {
    TAKEN[arena / 64].fetch_and(!(1 << (arena % 64)), Release);
}

/// The first chunk of `arena` that holds no block. For the arena's thread only.
fn empty_chunk(arena: usize) -> Option<usize> {
    let first = arena << ARENA_BITS;
    // Acquire: what the threads that freed the chunk's blocks wrote to them came first.
    (first..first + ARENA_CHUNKS).find(|&c| LIVE.0[c].load(Relaxed) == FREED.0[c].load(Acquire))
}

/// Points `cell` at an empty chunk of the thread's arena, the first one, so that a thread
/// keeps to the start of its arena. False: the allocation goes to the wrapped allocator.
#[cold]
#[inline(never)]
fn refill(cell: &Cell<Local>) -> bool {
    let mut s = cell.get();
    if s.skip > 0 {
        if s.skip != NEVER {
            s.skip -= 1;
            cell.set(s);
        }
        return false;
    }
    if s.arena == NO_ARENA {
        // Acquire: with ON, BASE is visible, here and (through the block) to whichever thread
        // frees it.
        if !ON.load(Acquire) {
            return false;
        }
        // Registering the exit hook may allocate, and a thread that is exiting has none.
        s.skip = NEVER;
        cell.set(s);
        if EXIT.try_with(|_| ()).is_err() {
            return false;
        }
        // An arena whose chunks all hold blocks that outlived its earlier threads is passed
        // over.
        let mut from = 0;
        s.skip = BACKOFF;
        while let Some(arena) = take_arena(from) {
            if empty_chunk(arena).is_some() {
                s.arena = arena;
                s.skip = 0;
                break;
            }
            give_back(arena);
            from = arena + 1;
        }
        if s.arena == NO_ARENA {
            cell.set(s);
            return false;
        }
    }
    match empty_chunk(s.arena) {
        Some(c) => {
            s.chunk = c;
            s.ptr = BASE.load(Relaxed) + (c << CHUNK_BITS);
            s.end = s.ptr + CHUNK;
        }
        None => s.skip = BACKOFF,
    }
    cell.set(s);
    s.skip == 0
}

#[inline(always)]
fn round(n: usize) -> usize {
    (n + (ALIGN - 1)) & !(ALIGN - 1)
}

/// `n` bytes (a multiple of `ALIGN`) of the thread's current chunk, or 0.
#[inline(always)]
fn bump(n: usize) -> usize {
    LOCAL.with(|cell| {
        let mut s = cell.get();
        if s.ptr + n > s.end {
            if !refill(cell) {
                return 0;
            }
            s = cell.get();
        }
        let p = s.ptr;
        s.ptr += n;
        cell.set(s);
        let live = &LIVE.0[s.chunk];
        live.store(live.load(Relaxed).wrapping_add(1), Relaxed);
        p
    })
}

/// Frees the arena block at `p` (`off` bytes into the reservation at `base`).
#[inline(always)]
fn release(base: usize, off: usize, p: usize, size: usize) {
    let c = off >> CHUNK_BITS;
    LOCAL.with(|cell| {
        let mut s = cell.get();
        if c >> ARENA_BITS == s.arena {
            let live = LIVE.0[c].load(Relaxed).wrapping_sub(1);
            LIVE.0[c].store(live, Relaxed);
            if c == s.chunk {
                if live == FREED.0[c].load(Acquire) {
                    // The current chunk is empty: start it again.
                    s.ptr = base + (c << CHUNK_BITS);
                    cell.set(s);
                } else if p + round(size) == s.ptr {
                    // The last block handed out: take it back.
                    s.ptr = p;
                    cell.set(s);
                }
            }
        } else {
            FREED.0[c].fetch_add(1, Release);
        }
    })
}

/// Turns the arenas on for every [`Bump`] of the process, from now on; there is no turning
/// them off. Returns whether they are on: false where the address space could not be
/// reserved, or on a target other than 64-bit Linux and macOS. Without a `Bump` as the
/// global allocator this reserves the address space and changes nothing else.
pub fn enable() -> bool {
    static RESERVE: std::sync::Once = std::sync::Once::new();
    RESERVE.call_once(|| {
        if let Some(base) = reserve() {
            BASE.store(base, Relaxed);
            ON.store(true, Release);
            SMALL_NOW.store(SMALL, Release);
        }
    });
    enabled()
}

/// Whether [`enable`] turned the arenas on.
pub fn enabled() -> bool {
    ON.load(Acquire)
}

/// Whether `p` points into the arenas (for tests and diagnostics).
pub fn in_arena(p: *const u8) -> bool {
    offset(BASE.load(Relaxed), p).is_some()
}

/// Threads that have an arena now (for tests and diagnostics).
pub fn arenas_in_use() -> usize {
    TAKEN.iter().map(|w| w.load(Relaxed).count_ones() as usize).sum()
}

#[cfg(all(any(target_os = "linux", target_os = "macos"), target_pointer_width = "64"))]
fn reserve() -> Option<usize> {
    // Address space only: a page is backed when it is first written.
    #[cfg(target_os = "linux")]
    const NORESERVE: libc::c_int = libc::MAP_NORESERVE;
    #[cfg(not(target_os = "linux"))]
    const NORESERVE: libc::c_int = 0;
    // SAFETY: an anonymous private mapping at an address the kernel chooses; nothing else
    // refers to it.
    let p = unsafe {
        libc::mmap(
            std::ptr::null_mut(),
            REGION,
            libc::PROT_READ | libc::PROT_WRITE,
            libc::MAP_PRIVATE | libc::MAP_ANON | NORESERVE,
            -1,
            0,
        )
    };
    (p != libc::MAP_FAILED).then_some(p as usize)
}

#[cfg(not(all(any(target_os = "linux", target_os = "macos"), target_pointer_width = "64")))]
fn reserve() -> Option<usize> {
    None
}

/// A global allocator that serves small blocks from per-thread arenas once [`enable`] has been
/// called, and forwards everything else to `A` (see the [module](self)).
///
/// The arenas are the process's, not the value's: declare one `Bump` as the global allocator.
pub struct Bump<A>(A);

impl<A> Bump<A> {
    /// Wraps the allocator that serves whatever the arenas do not.
    pub const fn new(inner: A) -> Self {
        Bump(inner)
    }
}

// SAFETY: a block is either the wrapped allocator's, which gets it back with the layout it
// was given, or `round(size)` bytes of the reservation, aligned to ALIGN, that no other live
// block overlaps: a chunk's bytes are handed out again only below the current thread's own
// last block or once the chunk's counts say that every block in it was freed. The two kinds
// are told apart by address, so a block allocated before `enable` is never taken for an
// arena's.
unsafe impl<A: GlobalAlloc> GlobalAlloc for Bump<A> {
    #[inline]
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        if l.size() <= SMALL_NOW.load(Relaxed) && l.align() <= ALIGN {
            let p = bump(round(l.size()));
            if p != 0 {
                return p as *mut u8;
            }
        }
        // SAFETY: the caller's contract, passed on.
        unsafe { self.0.alloc(l) }
    }

    #[inline]
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        let base = BASE.load(Relaxed);
        match offset(base, p) {
            Some(off) => release(base, off, p as usize, l.size()),
            // SAFETY: not an arena block, so `p` came from the wrapped allocator with `l`.
            None => unsafe { self.0.dealloc(p, l) },
        }
    }

    #[inline]
    unsafe fn alloc_zeroed(&self, l: Layout) -> *mut u8 {
        if l.size() <= SMALL_NOW.load(Relaxed) && l.align() <= ALIGN {
            let p = bump(round(l.size()));
            if p != 0 {
                // SAFETY: `p` is a fresh block of at least `l.size()` bytes. A reused chunk
                // holds old bytes.
                unsafe { std::ptr::write_bytes(p as *mut u8, 0, l.size()) };
                return p as *mut u8;
            }
        }
        // SAFETY: the caller's contract, passed on.
        unsafe { self.0.alloc_zeroed(l) }
    }

    #[inline]
    unsafe fn realloc(&self, p: *mut u8, l: Layout, new: usize) -> *mut u8 {
        let base = BASE.load(Relaxed);
        let Some(off) = offset(base, p) else {
            // SAFETY: not an arena block, so `p` came from the wrapped allocator with `l`.
            return unsafe { self.0.realloc(p, l, new) };
        };
        // The last block of this thread's current chunk grows or shrinks where it is.
        let (old, n) = (round(l.size()), round(new));
        let c = off >> CHUNK_BITS;
        let in_place = new <= SMALL
            && LOCAL.with(|cell| {
                let mut s = cell.get();
                let last = c == s.chunk && c >> ARENA_BITS == s.arena && p as usize + old == s.ptr;
                let fits = last && p as usize + n <= s.end;
                if fits {
                    s.ptr = p as usize + n;
                    cell.set(s);
                }
                fits
            });
        if in_place {
            return p;
        }
        // SAFETY: `new` is nonzero and valid for `l.align()` (the caller's contract); the two
        // blocks are distinct, and `p` is freed with the layout it was allocated with.
        unsafe {
            let q = self.alloc(Layout::from_size_align_unchecked(new, l.align()));
            if !q.is_null() {
                std::ptr::copy_nonoverlapping(p, q, l.size().min(new));
                self.dealloc(p, l);
            }
            q
        }
    }
}
