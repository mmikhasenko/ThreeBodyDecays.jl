# [Performance notes](@id performance_notes)

This page gives a performance model for `ThreeBodyDecays.jl`: where a typical amplitude evaluation spends time, which blocks are plausible optimization targets, and which attractive changes have already been checked and found not to move the full workload.

The guiding principle is to optimize at the level where the computation is actually expensive. Some kinematic helpers run in only a few nanoseconds. Improving them is useful when the change is local and readable, but the effect is quickly hidden inside the microsecond-scale amplitude path. Larger gains are more likely to come from reducing repeated Wigner matrix work across chains, or from evaluating many phase-space points together in fitting workloads.

## Performance model

A pointwise model evaluation has a few distinct layers:

- Kinematic primitives convert between invariant masses, decay angles, and phase-space limits.
- `aligned_amplitude` builds the chain amplitude in the natural frame of a single topology.
- `amplitude` rotates that aligned amplitude into the chosen helicity reference frame using Wigner small-`d` matrices.
- `unpolarized_intensity` evaluates one or more chains and sums over helicities.

The first layer is cheap. The middle layers do more work: they construct vertex matrices, evaluate Wigner functions, apply Clebsch-Gordan coefficients through the recoupling objects, and contract small tensors. In multi-chain models, the same reference-frame rotations can also be recomputed many times for the same phase-space point.

## Kinematic primitives

Functions such as `cosθij`, `σjofk`, and `lims` sit on hot call paths, so they should stay type-stable and easy for Julia to inline. Inlining matters because it lets the compiler propagate a literal topology index through wrappers such as:

```julia
cosθ23(σs, msq) = cosθij(σs, msq; k = 1)
```

Benchmarks of these helpers show nanosecond-scale runtimes. With inlining, the plain integer keyword path and a `Val`-based path reached the same performance for the tested angle and invariant conversions. `@constprop :aggressive` alone was not enough in the original shape, because the topology lookup still hid the constant from the tuple indexing sites.

The practical conclusion is modest: keep the primitives inlineable and type-stable, but do not expect this layer to dominate complete amplitude evaluations.

## Type inference

JET and `@inferred` checks are useful here as guardrails. They can catch type-instability before it turns into a performance regression, and they are especially helpful around generic tuple and keyword-heavy code.

The current findings do not indicate that inference is the limiting factor for the main amplitude path. On Julia 1.11.5, inference succeeds for the tested uniform-tuple indexing cases. The same pattern would deserve renewed attention if the code started to use heterogeneous tuples, if indexed values became type parameters for later dispatch, or if a caller introduced an outer loop that could not specialize.

## Aligned amplitude

`aligned_amplitude` allocates small arrays for the two vertex matrices and the output tensor. It also uses `permutedims` and a Tullio contraction to express the tensor layout and summation compactly.

This block is a plausible target for allocation reduction, but not automatically a good target for speed. Measurements indicate that the small heap allocations are not the main cost in typical scalar-amplitude use. The expensive work is the physics computation around Wigner functions and recouplings. Removing the arrays without changing that computation can improve allocation counts while leaving runtime nearly unchanged.

The best reason to revisit this block would be a profile from a real fitting workload that shows allocation pressure or garbage collection in this function. Without that signal, preserving the readable tensor expression is a reasonable tradeoff.

## Helicity rotations

The full `amplitude` call builds Wigner small-`d` matrices for the parent and final-state particles, then contracts them with the aligned amplitude. For a single chain this is already a meaningful part of the computation. For a model with many chains, the same alignment rotations can be recomputed for every chain at the same phase-space point.

This makes caching a stronger optimization candidate than rewriting the local tensor containers. A model-level path could compute the alignment matrices once for a given `(σs, refζs, two_js)` and pass them into the chain evaluations. The expected benefit depends on the number of chains and helicity dimensions, but the direction is aligned with the observed cost structure: avoid repeated Wigner work rather than only changing where small arrays live.

## Static spin dimensions

One tempting idea is to encode spin dimensions as compile-time values, build Wigner matrices as `SMatrix` objects, and replace Tullio contractions with explicit static loops. This was prototyped and benchmarked.

The result was mixed:

| Case | Result |
| :-- | :-- |
| small-spin aligned amplitude | about neutral |
| small-spin scalar amplitude | slightly slower |
| small-spin array amplitude | faster, roughly 1.6x in the tested setup |
| larger-spin scalar and array paths | slower |
| 20-chain unpolarized intensity | slower in the tested setup |

This is an important negative result. The static version reduced allocations, but it also added dispatch barriers, more type parameters, less readable contractions, and more specialization per spin configuration. Since the typical scalar and multi-chain workloads did not improve, static spin dimensions are not a recommended general refactor.

## Batch evaluation

Fitting workflows usually evaluate a model over many phase-space points. Optimizing only the single-point call can miss opportunities that appear at the batch level: shared setup, better memory locality, and possible vectorization of repeated kinematic and Wigner computations.

A batched API is therefore one of the more promising directions for large speedups. It should be designed around the data layout used by downstream fits, because the best representation depends on whether the caller needs all helicities, scalar intensities, chain-level amplitudes, or derivatives.

## Practical guidance

When asking whether `ThreeBodyDecays.jl` can be faster, start by identifying the workload:

- For single-chain scalar amplitudes, kinematic helper improvements are unlikely to change the total time substantially.
- For multi-chain intensities, look first for repeated Wigner matrix construction and shared reference-frame work.
- For large fits, consider batched evaluation before rewriting the pointwise tensor code.
- For allocation-heavy downstream profiles, consider a workspace API, but only after the profile shows that allocation pressure is a real bottleneck.
- Avoid broad `Val` or `StaticArrays` rewrites unless a benchmark for the target workload shows a clear win.

The current performance picture is encouraging: the simple paths are already type-stable, the cheap helpers are in the nanosecond range, and the main opportunities are at well-defined higher levels of the computation rather than in hidden type-instability.
