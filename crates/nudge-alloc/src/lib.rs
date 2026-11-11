//! Small blocks from per-thread arenas, in front of any global allocator.
//!
//! Some programs make many small allocations and free nearly all of them soon after: a read
//! mapper makes a few hundred for a read and frees them before the next one. [`Bump`] wraps a
//! program's global allocator and, once [`enable`] has been called, serves blocks of up to
//! 1 KB from memory that belongs to the allocating thread: a block is the next free bytes of
//! the thread's current 64 KB chunk, and a chunk whose blocks have all been freed is used
//! again from its start. Such a thread then stays in the same few pages, which is where the
//! gain comes from, for memory of at most 4 MB a thread. Every other allocation, and every
//! allocation before `enable`, goes to the wrapped allocator, at the cost of one comparison.
//! (It was written for [spearfish](https://github.com/COMBINE-lab/spearfish), which maps reads
//! in about 2% less CPU time with it.)
//!
//! A program opts in twice: it declares the allocator, and it turns the arenas on.
//!
//! ```no_run
//! use nudge_alloc::Bump;
//!
//! #[global_allocator]
//! static GLOBAL: Bump<std::alloc::System> = Bump::new(std::alloc::System);
//!
//! fn main() {
//!     nudge_alloc::enable();
//!     // ... the work, on any threads ...
//!     eprintln!("{:?}", nudge_alloc::stats());
//! }
//! ```
//!
//! The wrapped allocator can be any [`GlobalAlloc`]. A block may be freed by a thread other
//! than the one that allocated it, and a thread's arena passes to a later thread when it
//! exits. Up to 256 threads have an arena at a time, and a thread without one, or whose 64
//! chunks all hold live blocks, uses the wrapped allocator. Blocks that stay allocated keep
//! their chunk from being reused, so the arenas suit programs whose small allocations are
//! mostly short-lived. Nothing says so when a thread goes without: [`stats`] counts how often
//! it happened.
//!
//! The arenas are one reservation of 1 GB of address space, of which only the chunks in use
//! are backed by memory. Where transparent huge pages are always on, a chunk's first write
//! backs the 2 MB around it, half its thread's arena, which is a little faster and stays
//! within the 4 MB a thread. Where the kernel does not overcommit
//! (`vm.overcommit_memory=2`), the whole reservation counts against the commit limit. The
//! arenas are there on 64-bit Linux and macOS; on other targets `enable` returns false and
//! `Bump` only forwards.
//!
//! As with most allocators, a signal handler must not allocate or free through it. A thread
//! that never exits normally (after `fork`, or a detached thread when the process ends) keeps
//! its arena.
//!
//! # Changing it
//!
//! What threads share (who has an arena, and whether a chunk still holds a block) is in
//! `src/proto.rs`, a file that `tests/loom.rs` also compiles, against
//! [loom](https://docs.rs/loom)'s atomics, to check every interleaving of small models. A
//! change to an ordering there is checked by `cargo test`. `tests/miri.rs` is small enough for
//! Miri (`cargo +nightly miri test --test miri`), and `tests/stress.rs` runs for as long as
//! `STRESS_SECS` says.

#![warn(missing_docs)]

use std::alloc::{GlobalAlloc, Layout};
use std::cell::Cell;
use std::ptr::with_exposed_provenance_mut;
use std::sync::atomic::Ordering::Relaxed;

mod proto;
/// The atomics `proto` is written against (`tests/loom.rs` gives it loom's).
mod sync {
    pub(crate) use std::sync::atomic::{AtomicBool, AtomicU32, AtomicU64, AtomicUsize, Ordering};
}
use sync::{AtomicBool, AtomicU32, AtomicU64, AtomicUsize};

/// Miri backs the whole reservation and runs the tests slowly: under it, fewer and smaller
/// arenas and a shorter back-off.
const TINY: bool = cfg!(miri);

/// Largest block served from an arena, and the alignment of every such block.
const SMALL: usize = 1024;
const ALIGN: usize = 16;
const CHUNK_BITS: usize = 16;
const CHUNK: usize = 1 << CHUNK_BITS;
/// Chunks in a thread's arena.
const ARENA_BITS: usize = if TINY { 2 } else { 6 };
const ARENA_CHUNKS: usize = 1 << ARENA_BITS;
/// Arenas, so threads that can have one at a time.
const ARENAS: usize = if TINY { 64 } else { 256 };
const CHUNKS: usize = ARENAS * ARENA_CHUNKS;
const REGION: usize = CHUNKS * CHUNK;
/// Small allocations a thread passes to the wrapped allocator after it found no free chunk (or
/// no free arena) before it looks again.
const BACKOFF: u32 = if TINY { 64 } else { 4096 };
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
/// The largest size served from an arena now: zero, so that none is, until `enable`. The one
/// load an allocation costs while the arenas are off.
static SMALL_NOW: AtomicUsize = AtomicUsize::new(0);
/// The reservation is `BASE..BASE + REGION`. Until `enable` reserved it, the top of the
/// address space, where no block is: the one load a free costs while the arenas are off.
static BASE: AtomicUsize = AtomicUsize::new(0usize.wrapping_sub(REGION));

/// Bytes from the reservation's start to `p`, when `p` is in it.
#[inline(always)]
fn offset(base: usize, p: *const u8) -> Option<usize> {
    let off = p.addr().wrapping_sub(base);
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

/// What an arena's threads did, for [`stats`]. Only the arena's thread writes it, so each
/// arena has a cache line of its own.
#[repr(align(64))]
struct Tally {
    refills: AtomicU64,
    full: AtomicU64,
}
static TALLY: [Tally; ARENAS] = [const {
    Tally {
        refills: AtomicU64::new(0),
        full: AtomicU64::new(0),
    }
}; ARENAS];
/// Times a thread found every arena taken, and the most arenas taken at once.
static NO_ARENA_FOUND: AtomicU64 = AtomicU64::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);

/// Adds one to a count that only this thread writes.
fn tick(count: &AtomicU64) {
    count.store(count.load(Relaxed) + 1, Relaxed);
}

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
    // EXIT has one, and a thread's first small allocation registers it from inside the
    // allocator. That is sound where registering does not come back here for memory in a way
    // std forbids: glibc and macOS register through libc, and std's own list (musl) is on the
    // system allocator from Rust 1.89, the crate's `rust-version`; an earlier std aborts.
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
        if let Some(bit) = proto::take(word, below) {
            return Some(w * 64 + bit as usize);
        }
    }
    None
}

fn give_back(arena: usize) {
    proto::give_back(&TAKEN[arena / 64], (arena % 64) as u32);
}

fn taken() -> usize {
    TAKEN.iter().map(|w| w.load(Relaxed).count_ones() as usize).sum()
}

/// The first chunk of `arena` that holds no block. For the arena's thread only.
fn empty_chunk(arena: usize) -> Option<usize> {
    let first = arena << ARENA_BITS;
    (first..first + ARENA_CHUNKS).find(|&c| proto::is_empty(&LIVE.0[c], &FREED.0[c]))
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
        if !proto::is_on(&ON) {
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
            NO_ARENA_FOUND.fetch_add(1, Relaxed);
            cell.set(s);
            return false;
        }
        PEAK.fetch_max(taken(), Relaxed);
    }
    match empty_chunk(s.arena) {
        Some(c) => {
            s.chunk = c;
            s.ptr = proto::base(&BASE) + (c << CHUNK_BITS);
            s.end = s.ptr + CHUNK;
            tick(&TALLY[s.arena].refills);
        }
        None => {
            s.skip = BACKOFF;
            tick(&TALLY[s.arena].full);
        }
    }
    cell.set(s);
    s.skip == 0
}

#[inline(always)]
fn round(n: usize) -> usize {
    (n + (ALIGN - 1)) & !(ALIGN - 1)
}

/// Whether an arena serves `l` now: 1 to `SMALL_NOW` bytes. A size of zero, which a caller
/// must not ask for, goes to the wrapped allocator like any size while the arenas are off: a
/// block of no bytes at the end of a full chunk would be counted in one chunk and freed in
/// the next.
#[inline(always)]
fn serves(l: Layout) -> bool {
    l.size().wrapping_sub(1) < proto::limit(&SMALL_NOW) && l.align() <= ALIGN
}

/// The address of `n` bytes (a nonzero multiple of `ALIGN`) of the thread's current chunk,
/// or 0.
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
        proto::handed_out(&LIVE.0[s.chunk]);
        p
    })
}

/// Frees the arena block at address `p` (`off` bytes into the reservation at `base`).
#[inline(always)]
fn release(base: usize, off: usize, p: usize, size: usize) {
    let c = off >> CHUNK_BITS;
    LOCAL.with(|cell| {
        let mut s = cell.get();
        if c >> ARENA_BITS == s.arena {
            let live = proto::freed_by_owner(&LIVE.0[c]);
            if c == s.chunk {
                if proto::holds_none(live, &FREED.0[c]) {
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
            proto::freed_by_other(&FREED.0[c]);
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
            proto::turn_on(&BASE, &ON, &SMALL_NOW, base, SMALL);
        }
    });
    enabled()
}

/// Whether [`enable`] turned the arenas on.
pub fn enabled() -> bool {
    proto::is_on(&ON)
}

/// What the arenas have done since [`enable`], from [`stats`].
///
/// A thread that finds no arena free, or no empty chunk in its arena, passes its next 4096
/// small allocations to the wrapped allocator before it looks again. A program that sees
/// many of either has more threads allocating than there are arenas, or keeps more small
/// blocks alive than an arena holds, and gains less from the arenas than it could.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct Stats {
    /// Whether the arenas are on.
    pub enabled: bool,
    /// Threads that have an arena now.
    pub arenas_in_use: usize,
    /// The most threads that had an arena at once.
    pub arenas_peak: usize,
    /// Times a thread moved on to an empty chunk of its arena, its first one included. A
    /// chunk that empties while it is the current one is used again without being counted.
    pub chunk_refills: u64,
    /// Times a thread looked for an arena and every one was taken.
    pub no_arena: u64,
    /// Times a thread's arena had no empty chunk.
    pub arena_full: u64,
}

impl Stats {
    /// At most this many small allocations went to the wrapped allocator because a thread
    /// had no arena or a full one: 4097 for each time counted in [`no_arena`](Self::no_arena)
    /// and [`arena_full`](Self::arena_full), fewer where the thread stopped allocating first.
    pub fn fallback_allocations(&self) -> u64 {
        (self.no_arena + self.arena_full) * (u64::from(BACKOFF) + 1)
    }
}

/// What the arenas have done so far. The counts are read while threads go on allocating, so
/// they are a moment's view, and all zero before [`enable`].
pub fn stats() -> Stats {
    Stats {
        enabled: enabled(),
        arenas_in_use: taken(),
        arenas_peak: PEAK.load(Relaxed),
        chunk_refills: TALLY.iter().map(|t| t.refills.load(Relaxed)).sum(),
        no_arena: NO_ARENA_FOUND.load(Relaxed),
        arena_full: TALLY.iter().map(|t| t.full.load(Relaxed)).sum(),
    }
}

/// Whether `p` points into the arenas (for the crate's tests).
#[doc(hidden)]
pub fn in_arena(p: *const u8) -> bool {
    offset(proto::base(&BASE), p).is_some()
}

#[cfg(all(any(target_os = "linux", target_os = "macos"), target_pointer_width = "64"))]
fn reserve() -> Option<usize> {
    // Address space only: a page is backed when it is first written.
    #[cfg(all(target_os = "linux", not(miri)))]
    const NORESERVE: libc::c_int = libc::MAP_NORESERVE;
    #[cfg(not(all(target_os = "linux", not(miri))))]
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
    // Exposed: every block is made from an address in the mapping.
    (p != libc::MAP_FAILED).then(|| p.expose_provenance())
}

#[cfg(not(all(any(target_os = "linux", target_os = "macos"), target_pointer_width = "64")))]
fn reserve() -> Option<usize> {
    None
}

/// A global allocator that serves small blocks from per-thread arenas once [`enable`] has been
/// called, and forwards everything else to `A` (see the [crate] documentation).
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
        if serves(l) {
            let p = bump(round(l.size()));
            if p != 0 {
                return with_exposed_provenance_mut(p);
            }
        }
        // SAFETY: the caller's contract, passed on.
        unsafe { self.0.alloc(l) }
    }

    #[inline]
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        let base = proto::base(&BASE);
        match offset(base, p) {
            Some(off) => release(base, off, p.addr(), l.size()),
            // SAFETY: not an arena block, so `p` came from the wrapped allocator with `l`.
            None => unsafe { self.0.dealloc(p, l) },
        }
    }

    #[inline]
    unsafe fn alloc_zeroed(&self, l: Layout) -> *mut u8 {
        if serves(l) {
            let p = bump(round(l.size()));
            if p != 0 {
                let p = with_exposed_provenance_mut::<u8>(p);
                // SAFETY: `p` is a fresh block of at least `l.size()` bytes. A reused chunk
                // holds old bytes.
                unsafe { std::ptr::write_bytes(p, 0, l.size()) };
                return p;
            }
        }
        // SAFETY: the caller's contract, passed on.
        unsafe { self.0.alloc_zeroed(l) }
    }

    #[inline]
    unsafe fn realloc(&self, p: *mut u8, l: Layout, new: usize) -> *mut u8 {
        let base = proto::base(&BASE);
        let Some(off) = offset(base, p) else {
            // SAFETY: not an arena block, so `p` came from the wrapped allocator with `l`.
            return unsafe { self.0.realloc(p, l, new) };
        };
        // The last block of this thread's current chunk grows or shrinks where it is.
        let (old, n) = (round(l.size()), round(new));
        let c = off >> CHUNK_BITS;
        let in_place = new.wrapping_sub(1) < SMALL
            && LOCAL.with(|cell| {
                let mut s = cell.get();
                let last = c == s.chunk && c >> ARENA_BITS == s.arena && p.addr() + old == s.ptr;
                let fits = last && p.addr() + n <= s.end;
                if fits {
                    s.ptr = p.addr() + n;
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
