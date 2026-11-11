# nudge-alloc

Small blocks from per-thread arenas, in front of any global allocator.

Some programs make many small allocations and free nearly all of them soon after. `Bump<A>`
wraps a global allocator `A` and, once `nudge_alloc::enable()` has been called, serves blocks
of up to 1 KB from memory that belongs to the allocating thread: a block is the next free
bytes of the thread's current 64 KB chunk, and a chunk whose blocks have all been freed is
used again from its start. Such a thread stays in the same few pages. Everything else, and
everything before `enable`, goes to `A`, at the cost of one comparison.

```rust
use nudge_alloc::Bump;

#[global_allocator]
static GLOBAL: Bump<std::alloc::System> = Bump::new(std::alloc::System);

fn main() {
    nudge_alloc::enable();
    // ... the work, on any threads ...
    eprintln!("{:?}", nudge_alloc::stats());
}
```

- A block may be freed on any thread, and a thread's arena passes to a later thread when it
  exits.
- Up to 256 threads have an arena (4 MB) at a time. A thread without one, or whose arena is
  full of live blocks, uses the wrapped allocator; `stats()` counts how often that happened.
- The arenas are there on 64-bit Linux and macOS. Elsewhere `enable()` returns false and
  `Bump` only forwards.
- It suits programs whose small allocations are mostly short-lived. Blocks that stay
  allocated keep their chunk from being reused.

It was written for [spearfish](https://github.com/COMBINE-lab/spearfish), a long-read mapper
(`spearfish map --bump-alloc`), and lives in its repository.

## Checking a change

The allocator is unsafe code that every allocation of a program goes through, so the crate
carries the means to check it:

- `cargo test -p nudge-alloc` runs the tests, among them [loom](https://docs.rs/loom) models
  of what threads share (`tests/loom.rs`, on Linux). The models call the allocator's own
  functions (`src/proto.rs`), so a changed memory ordering is checked as written.
- `cargo +nightly miri test -p nudge-alloc --test miri` runs a small program under Miri, with
  its data-race detection.
- `STRESS_SECS=20 cargo test -p nudge-alloc --release --test stress` runs the stress tests for
  longer: arenas changing hands, frees on other threads, and allocation inside thread-local
  destructors.

The repository's continuous integration runs all three.
