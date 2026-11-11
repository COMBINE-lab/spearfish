//! Small blocks from per-thread arenas, in front of any global allocator.
//!
//! Some programs make many small allocations and free nearly all of them soon after: a read
//! mapper makes a few hundred for a read and frees them before the next one. [`Bump`] wraps a
//! program's global allocator and, once [`enable`] has been called, serves blocks of up to
//! 8 KB ([`MAX_BLOCK`]) from memory that belongs to the allocating thread: a block is the next
//! free bytes of the thread's current 64 KB chunk, and a chunk whose blocks have all been
//! freed is used again from its start. Such a thread then stays in the same few pages, which
//! is where the gain comes from, for memory of at most 4 MB a thread. Every other allocation,
//! and every allocation before `enable`, goes to the wrapped allocator, at the cost of one
//! comparison. (It was written for [spearfish](https://github.com/COMBINE-lab/spearfish),
//! which maps reads in about 3% less CPU time with it.)
//!
//! A program opts in twice: it declares the allocator, and it turns the arenas on.
//!
//! ```no_run
//! use nudgy::Bump;
//!
//! #[global_allocator]
//! static GLOBAL: Bump<std::alloc::System> = Bump::new(std::alloc::System);
//!
//! fn main() {
//!     nudgy::enable();
//!     // ... the work, on any threads ...
//!     eprintln!("{:?}", nudgy::stats());
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
//! An allocation from an arena moves a pointer in the thread's own storage, and a free moves
//! it back or adds to a count there; neither writes anything that another thread reads. The
//! last block of a chunk grows and shrinks where it is, so a vector that is being filled is
//! not copied until another block follows it. `cargo bench` times the common patterns with
//! and without the arenas.
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
//! Miri (`cargo +nightly miri test --test miri`), `tests/stress.rs` runs for as long as
//! `STRESS_SECS` says, and `benches/ops.rs` (`cargo bench`) times an allocation and its free.

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

/// The largest block, in bytes, that an arena serves. A larger block, or one aligned to more
/// than 16 bytes, is the wrapped allocator's.
///
/// The limit is a trade. A vector that grows past it is copied to the wrapped allocator, so
/// a low limit sends vectors there that would have been filled and freed within their chunk;
/// larger blocks fill a thread's 4 MB sooner, and leave more of a chunk unused at its end
/// (at most this many bytes of its 64 KB). For spearfish, 8 KB was the best of the limits
/// from 512 bytes to 32 KB.
pub const MAX_BLOCK: usize = 8192;
/// The alignment of every block of an arena.
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
const NO_CHUNK: usize = usize::MAX;

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
/// Bytes handed out of each chunk minus those its arena's thread freed, modulo 2^32. Only
/// that thread reads or writes it, and for the chunk it is allocating from it keeps the count
/// in `Local` instead: this one is behind until the thread leaves the chunk (`leave`).
static LIVE: PerChunk = PerChunk([const { AtomicU32::new(0) }; CHUNKS]);
/// Bytes of each chunk that other threads freed, modulo 2^32. `LIVE == FREED`: the chunk
/// holds no block (the two differ by the bytes of its live blocks, at most a chunk's).
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
struct Local {
    /// The next free byte of the current chunk, and the chunk's end: both zero without one.
    ptr: Cell<usize>,
    end: Cell<usize>,
    /// `ptr - bias`, in its low 32 bits, is the current chunk's `LIVE`. An allocation moves
    /// `ptr` and so counts itself; a free moves `ptr` back or `bias` up.
    bias: Cell<usize>,
    /// The current chunk (one of the arena's), or `NO_CHUNK`.
    chunk: Cell<usize>,
    arena: Cell<usize>,
    /// Small allocations still to pass on before `refill` looks for a chunk again.
    skip: Cell<u32>,
}

impl Local {
    const fn new() -> Local {
        Local {
            ptr: Cell::new(0),
            end: Cell::new(0),
            bias: Cell::new(0),
            chunk: Cell::new(NO_CHUNK),
            arena: Cell::new(NO_ARENA),
            skip: Cell::new(0),
        }
    }

    /// The current chunk's `LIVE`.
    #[inline(always)]
    fn live(&self) -> usize {
        self.ptr.get().wrapping_sub(self.bias.get())
    }

    /// Writes the current chunk's `LIVE` where whoever looks at the chunk next reads it:
    /// before this thread looks for an empty chunk, and before it gives the arena back.
    fn leave(&self) {
        let c = self.chunk.get();
        if c != NO_CHUNK {
            proto::left(&LIVE.0[c], self.live() as u32);
        }
    }

    /// Makes `c`, which holds no block and whose `LIVE` is `live`, the current chunk.
    fn start(&self, base: usize, c: usize, live: usize) {
        let first = base + (c << CHUNK_BITS);
        self.chunk.set(c);
        self.ptr.set(first);
        self.end.set(first + CHUNK);
        self.bias.set(first.wrapping_sub(live));
    }
}

/// Gives the thread's arena back when the thread exits.
struct Exit;

thread_local! {
    // No destructor, so it stays readable while the thread's other destructors run.
    // EXIT has one, and a thread's first small allocation registers it from inside the
    // allocator. That is sound where registering does not come back here for memory in a way
    // std forbids: glibc and macOS register through libc, and std's own list (musl) is on the
    // system allocator from Rust 1.89, the crate's `rust-version`; an earlier std aborts.
    static LOCAL: Local = const { Local::new() };
    static EXIT: Exit = const { Exit };
}

impl Drop for Exit {
    fn drop(&mut self) {
        let _ = LOCAL.try_with(|l| {
            let arena = l.arena.get();
            l.leave();
            // Destructors that run after this one allocate from the wrapped allocator, and
            // free this arena's blocks as any other thread would.
            l.ptr.set(0);
            l.end.set(0);
            l.chunk.set(NO_CHUNK);
            l.arena.set(NO_ARENA);
            l.skip.set(NEVER);
            if arena != NO_ARENA {
                give_back(arena);
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

/// The first chunk of `arena` that holds no block, and its `LIVE`. For the arena's thread
/// only, once it has left its current chunk.
fn empty_chunk(arena: usize) -> Option<(usize, u32)> {
    let first = arena << ARENA_BITS;
    (first..first + ARENA_CHUNKS).find_map(|c| {
        let live = proto::live(&LIVE.0[c]);
        proto::holds_none(live, &FREED.0[c]).then_some((c, live))
    })
}

/// Points the thread at an empty chunk of its arena, the first one, so that a thread keeps
/// to the start of its arena. False: the allocation goes to the wrapped allocator.
fn refill(l: &Local) -> bool {
    let skip = l.skip.get();
    if skip > 0 {
        if skip != NEVER {
            l.skip.set(skip - 1);
        }
        return false;
    }
    let mut arena = l.arena.get();
    if arena == NO_ARENA {
        if !proto::is_on(&ON) {
            return false;
        }
        // Registering the exit hook may allocate, and a thread that is exiting has none.
        l.skip.set(NEVER);
        if EXIT.try_with(|_| ()).is_err() {
            return false;
        }
        // An arena whose chunks all hold blocks that outlived its earlier threads is passed
        // over.
        let mut from = 0;
        while let Some(a) = take_arena(from) {
            if empty_chunk(a).is_some() {
                arena = a;
                break;
            }
            give_back(a);
            from = a + 1;
        }
        if arena == NO_ARENA {
            l.skip.set(BACKOFF);
            NO_ARENA_FOUND.fetch_add(1, Relaxed);
            return false;
        }
        l.skip.set(0);
        l.arena.set(arena);
        PEAK.fetch_max(taken(), Relaxed);
    }
    l.leave();
    match empty_chunk(arena) {
        Some((c, live)) => {
            l.start(proto::base(&BASE), c, live as usize);
            tick(&TALLY[arena].refills);
            true
        }
        None => {
            l.skip.set(BACKOFF);
            tick(&TALLY[arena].full);
            false
        }
    }
}

#[inline(always)]
fn round(n: usize) -> usize {
    (n + (ALIGN - 1)) & !(ALIGN - 1)
}

/// Whether an arena serves `l` now: 1 to `SMALL_NOW` bytes. A size of zero, which a caller
/// must not ask for, goes to the wrapped allocator like any size while the arenas are off.
#[inline(always)]
fn serves(l: Layout) -> bool {
    l.size().wrapping_sub(1) < proto::limit(&SMALL_NOW) && l.align() <= ALIGN
}

/// The address of `n` bytes (a nonzero multiple of `ALIGN`) of the thread's current chunk,
/// when it has that many left. All an allocation costs when it does.
#[inline(always)]
fn bump(n: usize) -> Option<usize> {
    LOCAL.with(|l| {
        let p = l.ptr.get();
        (n <= l.end.get() - p).then(|| {
            l.ptr.set(p + n);
            p
        })
    })
}

/// `bump` for a thread whose current chunk has fewer than `n` bytes left, or that has none:
/// from an empty chunk of its arena, or 0.
#[cold]
#[inline(never)]
fn bump_slow(n: usize) -> usize {
    LOCAL.with(|l| {
        if !refill(l) {
            return 0;
        }
        let p = l.ptr.get();
        l.ptr.set(p + n);
        p
    })
}

/// Frees the `n` bytes (a block's size, rounded) at address `p`, `off` bytes into the
/// reservation at `base`.
#[inline(always)]
fn release(base: usize, off: usize, p: usize, n: usize) {
    let c = off >> CHUNK_BITS;
    LOCAL.with(|l| {
        if c == l.chunk.get() {
            if p + n == l.ptr.get() {
                // The last block handed out: take it back.
                l.ptr.set(p);
            } else {
                l.bias.set(l.bias.get().wrapping_add(n));
            }
            let live = l.live();
            if proto::holds_none(live as u32, &FREED.0[c]) {
                // The current chunk is empty: start it again, with the same `LIVE`.
                let first = base + (c << CHUNK_BITS);
                l.ptr.set(first);
                l.bias.set(first.wrapping_sub(live));
            }
        } else if c >> ARENA_BITS == l.arena.get() {
            proto::freed_by_owner(&LIVE.0[c], n as u32);
        } else {
            proto::freed_by_other(&FREED.0[c], n as u32);
        }
    })
}

/// Makes the block of `old` bytes at `p`, in chunk `c`, one of `n` bytes where it is, when
/// it is the last block of this thread's current chunk and the chunk has the room.
#[inline(always)]
fn resize_last(c: usize, p: usize, old: usize, n: usize) -> bool {
    LOCAL.with(|l| {
        let fits = c == l.chunk.get() && p + old == l.ptr.get() && p + n <= l.end.get();
        if fits {
            l.ptr.set(p + n);
        }
        fits
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
            proto::turn_on(&BASE, &ON, &SMALL_NOW, base, MAX_BLOCK);
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

impl<A: GlobalAlloc> Bump<A> {
    /// `alloc` and `alloc_zeroed` for a small block that the thread's current chunk has no
    /// room for.
    #[cold]
    #[inline(never)]
    unsafe fn alloc_slow(&self, l: Layout, zeroed: bool) -> *mut u8 {
        match bump_slow(round(l.size())) {
            // SAFETY: the caller's contract, passed on.
            0 if zeroed => unsafe { self.0.alloc_zeroed(l) },
            // SAFETY: as above.
            0 => unsafe { self.0.alloc(l) },
            p => {
                let p = with_exposed_provenance_mut::<u8>(p);
                if zeroed {
                    // SAFETY: `p` is a fresh block of at least `l.size()` bytes.
                    unsafe { std::ptr::write_bytes(p, 0, l.size()) };
                }
                p
            }
        }
    }
}

impl<A: GlobalAlloc> Bump<A> {
    /// `realloc` for an arena block that cannot stay where it is: a new block, the bytes
    /// copied, the old one freed.
    #[inline(never)]
    unsafe fn realloc_moving(&self, p: *mut u8, l: Layout, new: usize) -> *mut u8 {
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

// SAFETY: a block is either the wrapped allocator's, which gets it back with the layout it
// was given, or `round(size)` bytes of the reservation, aligned to ALIGN, that no other live
// block overlaps: a chunk's bytes are handed out again only above the current thread's own
// last block or once the chunk's counts say that every byte handed out of it was freed (the
// counts are of the same rounded sizes on both sides, because a block is freed with the size
// it was allocated or last reallocated with). The two kinds are told apart by address, so a
// block allocated before `enable` is never taken for an arena's.
unsafe impl<A: GlobalAlloc> GlobalAlloc for Bump<A> {
    #[inline]
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        if serves(l) {
            return match bump(round(l.size())) {
                Some(p) => with_exposed_provenance_mut(p),
                // SAFETY: the caller's contract, passed on.
                None => unsafe { self.alloc_slow(l, false) },
            };
        }
        // SAFETY: the caller's contract, passed on.
        unsafe { self.0.alloc(l) }
    }

    #[inline]
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        let base = proto::base(&BASE);
        match offset(base, p) {
            Some(off) => release(base, off, p.addr(), round(l.size())),
            // SAFETY: not an arena block, so `p` came from the wrapped allocator with `l`.
            None => unsafe { self.0.dealloc(p, l) },
        }
    }

    #[inline]
    unsafe fn alloc_zeroed(&self, l: Layout) -> *mut u8 {
        if serves(l) {
            return match bump(round(l.size())) {
                Some(p) => {
                    let p = with_exposed_provenance_mut::<u8>(p);
                    // SAFETY: `p` is a fresh block of at least `l.size()` bytes. A reused
                    // chunk holds old bytes.
                    unsafe { std::ptr::write_bytes(p, 0, l.size()) };
                    p
                }
                // SAFETY: the caller's contract, passed on.
                None => unsafe { self.alloc_slow(l, true) },
            };
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
        // A block keeps its place when its rounded size does not change, and the last block
        // of this thread's current chunk grows (to MAX_BLOCK at most) or shrinks where it is.
        let (old, n) = (round(l.size()), round(new));
        if n == old || (new.wrapping_sub(1) < MAX_BLOCK && resize_last(off >> CHUNK_BITS, p.addr(), old, n)) {
            return p;
        }
        // SAFETY: the caller's contract, passed on.
        unsafe { self.realloc_moving(p, l, new) }
    }
}
