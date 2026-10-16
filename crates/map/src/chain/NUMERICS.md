# Bidirectional chaining and exact fallback

The default library/CLI `oracle` selection uses `BidirectionalChainer` through
`AnyChainer::Oracle`. Its public spelling and default are unchanged.
`OracleChainer` remains the original forward normalized-prefix implementation,
usable directly as a reference. `RefineChainer` remains available unchanged.
No scoring parameters or seeding, alignment, or rescue rules are changed.

The default first checks whether bidirectional work is worthwhile. Single-transcript,
tiny, zero-window, and numerically ineligible candidates go directly to the
forward reference. For other candidates a planner examines the original placement
buckets before constructing per-transcript DP arrays. It counts forward predecessor
lanes after full-prefix reuse, prefix/suffix lanes after sharing, and crossing pairs.
Bidirectional work is selected only if it saves enough lanes to cover a setup
allowance of four edge evaluations per placement. This is a work-count heuristic,
not a measured runtime prediction. Identical complete lists and candidates without
enough suffix reuse keep the forward implementation's existing sharing.

The planner retains the full normalized forward ordering for both prefix processing
and any fallback, avoiding a redundant sort. It computes a separate suffix ordering
and caches adjacent common-prefix/suffix lengths. Sliding identity-window cursors
replace per-endpoint binary searches, and same-identity placements are excluded
before edge arithmetic. Prefix-only DP arrays contain only prefix placements.
Numerical certificate checks reduce each batch once rather than applying an f64
bound calculation to every edge.

The split is `candidate.anchors.len() / 2`, in candidate anchor identities,
not transcript placement indices. Every occurrence of one identity stays on
one side. Missing identities do not shrink the identity window.

Prefix placement lists are sorted lexicographically by `(identity, target start
relative to first prefix placement)`. Common prefixes copy F scores, predecessor
indices, and incoming discrepancies. A separate ordering uses reversed suffix
lists, normalized to their last placement. Common reversed prefixes copy H
continuations, next-placement distances, and outgoing discrepancies. Identities,
relative target starts, order, and multiplicity are all part of the signatures;
query spans and weights are fixed by the candidate identity. Scores and pointers
are retained in per-transcript scratch arrays, so a second ordering cannot
invalidate them. Scratch capacities are reused across candidates and reads.

H excludes its starting anchor. Its recurrence iterates right to left but calls
the original directed `transition(i, k)` arithmetic, using k's query-overlap gain
and weight.
The batched kernel uses structure-of-arrays inputs and architecture-dispatched
independent lanes. It does not reverse coordinates or interchange weights.

The result is the best of F[j], w[i]+H[i], and F[j]+e(j,i)+H[i]. Crossing edges
only involve prefix identities at least b-h and suffix identities below b+h;
each pair must also satisfy the exact identity difference and original validity
conditions. Missing identities and repeated placements are handled explicitly.
There are at most h(h+1)/2 crossing identity pairs, with all corresponding
placement combinations. The join uses compact batched arrays, not a tree, and
never propagates joined scores forward through the suffix.

Every valid edge strictly increases identity, so a path crosses the split at
most once. A crossing path therefore decomposes into a prefix ending at j,
the edge j -> i, and a suffix continuing after i. Endpoint-specific F and H
values cover all such paths, including paths skipping the split's neighbors.
Traceback reverses the prefix predecessor path and then follows suffix next
pointers; the maximum discrepancy includes the crossing edge.

## Numerical certificate and ties

The real-arithmetic decomposition alone does not preserve f32 accumulation.
This first implementation deliberately has a narrow certified fast path:

- All candidate weights and every evaluated valid edge increment must be finite
  integral f32 values. Negative zero also falls back to preserve score bits.
- With n candidate identities and M the largest absolute weight or increment,
  `(n+1)*M <= 2^24` must hold. The bound is checked in f64 using integer-valued
  quantities; the configured score calculations themselves remain f32.
- Query/normalized target coordinates must fit the original forward SIMD
  oracle's signed coordinates without overflow, spans must be nonempty, and
  gap bounds must be nonnegative i32 values. Parameters must be finite.
- Prefix ties retain the original forward strict-greater/earliest-predecessor
  behavior directly. A suffix continuation optimum tie (including stopping) or
  final winner tie causes fallback for that transcript. A tie later exceeded by
  a better alternative is harmless.
- Uncertifiable edge increments or sums also fall back per transcript. A failed
  transcript is excluded as a sharing source; if a planned neighbor fails, common
  lengths are recomputed against the last successful transcript. Certificate
  bounds are restored on failure so one transcript cannot poison the others.

A path uses at most n identities. Every forward partial sum, backward partial
sum, crossing sum, and losing alternative considered by the DP is thus an
integer of magnitude at most 2^24 and is exactly representable in f32. This
certifies every alternative, rather than just recomputing the selected path.
The edge increments are evaluated with the existing `transition` arithmetic;
its zero-overlap gain agrees with the original forward SIMD expression for
nonempty spans. Normalized target translations leave all coordinate differences
unchanged. Reused signatures have the same edges and therefore inherit their
certificate. Nonfinite valid increments reject the certificate.

The original forward tie behavior in F, uniqueness of retained suffix optima,
and uniqueness of the final decomposition suffice to preserve traceback. All
prefix paths have already been resolved exactly as in the forward reference;
no different prefix path can make a strictly worse endpoint sum win. Otherwise
the original forward chainer reproduces its earliest-predecessor/endpoint tie
behavior for the affected transcript. Failed transcripts are processed together
in their already sorted forward order, retaining prefix sharing within that subset.
Certified results for other transcripts are kept. Fallback is never a forward
recomputation of a joined suffix accepted as a bidirectional result.

Fractional weights (including ordinary specificity-weighted anchors), fractional
overlap/penalty increments, large sums, and ambiguous ties commonly fall back.
No approximation, scoring-precision change, or relaxed test tolerance is used.
This implementation makes no claim of a measured speedup. Broader numerical
certificates can be added later without changing the correctness contract.
