# [Performance notes](@id performance_notes)

This page collects the main performance findings from the JET and benchmarking pass around the amplitude code. It is written as a map for the next time someone asks whether `ThreeBodyDecays.jl` can run faster: yes, but not every attractive micro-optimization moves the end-to-end workload.

## The short version

The hot physics path is not dominated by small kinematic helper functions. A few nanosecond-level improvements are real, but they disappear inside a microsecond-level amplitude evaluation. For realistic multi-chain models, the more interesting work is to avoid repeated Wigner matrix construction across chains, or to evaluate many phase-space points as a batch.

The most useful ranking from the investigation is:

- Cache the alignment Wigner matrices across chains when evaluating a model at one phase-space point.
- Add a batched phase-space evaluation API for fitting loops.
- Consider a reusable workspace if allocation pressure shows up in a real profile.
- Do not rewrite the amplitude stack around `Val`-wrapped spin dimensions unless a new benchmark changes the conclusion.

## Kinematic primitives

Functions such as `cosθij`, `σjofk`, and `lims` are small enough that inference and inlining matter. Adding `@inline` to the hot kinematic functions lets Julia constant-fold through the topology index when it is passed as a literal keyword argument:

```julia
cosθ23(σs, msq) = cosθij(σs, msq; k = 1)
```

In the benchmark used for this pass, the plain integer path and a `Val` path both reached the same performance once `cosθij` was inlined. The important result was not that `Val` is bad; it was that `Val` did not add more once inlining gave the compiler enough context.

This also explains an initially surprising result: `@constprop :aggressive` alone did not help the original code, because the topology lookup still hid the constant from the tuple indexing sites.

## Inference was not the bottleneck

One of the original concerns was that `@inferred` might fail for tuple indexing with a runtime integer. On Julia 1.11.5, the baseline inference checks passed for the variants tested in this pass. The abstract interpreter can resolve indexing into uniform tuples well enough here.

That does not make inference checks useless. It means this package should treat them as a guardrail, not as proof of a performance problem. The same pattern could become real trouble if:

- tuples become heterogeneous,
- indexed values are used later as type parameters,
- the call sits inside an outer loop that cannot specialize.

Those conditions did not apply to the tested amplitude path.

## The amplitude path

`aligned_amplitude` still allocates small arrays: vertex matrices, an output tensor, a permuted view of that tensor, and Tullio-generated temporaries. The array-returning `amplitude` path allocates more because it builds Wigner matrices and contracts them with the aligned amplitude tensor.

That sounds like an obvious target for `StaticArrays`, but the measurements were less dramatic. In the tested spin-1/2 setup, replacing dynamic arrays with static arrays and manual contractions improved the array amplitude case, but the scalar amplitude case was neutral to slower. For a larger spin setup and a 20-chain model, the static version was slower.

The reason is simple: Wigner `d` evaluations and Clebsch-Gordan coefficient work dominated the runtime. Container choice reduced allocation counts, but did not remove the main computation.

## The rejected optimization

A full prototype encoded spin dimensions as compile-time values, built Wigner matrices as `SMatrix` objects, and replaced Tullio contractions with explicit static loops. The result reduced allocations, but did not meet the bar for carrying the extra complexity:

| Case | Result |
| :-- | :-- |
| small-spin aligned amplitude | about neutral |
| small-spin scalar amplitude | slightly slower |
| small-spin array amplitude | faster, roughly 1.6x in the tested setup |
| larger-spin scalar and array paths | slower |
| 20-chain unpolarized intensity | slower in the tested setup |

The complexity cost was high: an extra static spin representation, more dispatch barriers, more type parameters, less readable contractions, and likely more compilation per spin configuration. Since the typical model path did not improve, this direction was rejected for now.

## Better next targets

The promising ideas operate at the level where the time is actually spent.

First, cache alignment Wigner matrices across chains. In `unpolarized_intensity(model, σs)`, many chains share the same `σs` and reference-frame choice, so the same four alignment matrices can be recomputed repeatedly. Building them once per model evaluation and passing them down could be a small API change with a meaningful gain for multi-chain models.

Second, batch the evaluation over many phase-space points. Fits usually call the model over large samples, not isolated points. A batch API could let the Wigner and kinematic computations use vectorization or better memory locality, and it avoids optimizing only the per-point overhead.

Third, add a workspace API only if a profile shows allocation pressure in a downstream fitting loop. It is lower risk than a `Val` rewrite, but it should still be guided by a real caller.

## How to answer the speed question

If someone asks whether the project can be faster, the best answer is:

> Probably, but the easy-looking static-array rewrite was tested and did not pay off for the typical workload. The next serious attempts should cache repeated Wigner matrices across chains and support batched phase-space evaluation.

Micro-optimizing kinematic helpers is still worthwhile when it is local and readable. It should not be expected to change the runtime of full amplitude evaluations by itself.
