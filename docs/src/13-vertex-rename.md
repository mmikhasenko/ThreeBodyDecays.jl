# [Renaming `VertexFunction` → `Vertex`](@id vertex_rename)

ThreeBodyDecays.jl v0.15 renames the public type [`VertexFunction`](@ref) to [`Vertex`](@ref).

## Why

A decay chain is built from local two-body vertices (`HRk`, `Hij` in [`DecayChain`](@ref)), each carrying a
[`Recoupling`](@ref) scheme and an optional form-factor functor. The old name `VertexFunction` did not add
useful distinction—the type is already a vertex payload, not a callable amplitude function.

The shorter name aligns with the wider decay-chain ecosystem (e.g. `Propagator` in
[CascadeDecays.jl](https://github.com/RUB-EP1/CascadeDecays.jl)) and reads naturally in model code:

```julia
Vertex(RecouplingLS((2, 0)), my_form_factor)
```

## What changed

| Before | After |
| ------ | ----- |
| `VertexFunction(h)` | `Vertex(h)` |
| `VertexFunction(h, ff)` | `Vertex(h, ff)` |
| `HRk::VertexFunction` in docstrings/types | `HRk::Vertex` |

The struct fields are unchanged: `h` (recoupling) and `ff` (form factor).

## Backward compatibility

For one release cycle, the old name remains available as a deprecated alias.
On Julia 1.11 and later, using `VertexFunction` emits a deprecation warning; on Julia 1.9–1.10 it is a silent alias.

```julia
using ThreeBodyDecays

VertexFunction === Vertex  # true

dc = DecayChain(;
    k = 1,
    two_j = 2,
    Xlineshape = σ -> 1.0,
    HRk = VertexFunction(RecouplingLS((2, 0))),  # still works
    Hij = Vertex(RecouplingLS((2, 0))),          # preferred
    tbs = my_system,
)
```

`VertexFunction` is still exported but documented as deprecated. It will be removed in a future release.

## Migration checklist

1. Replace `VertexFunction` with `Vertex` in your model definitions and tests.
2. Update type annotations (`<: VertexFunction` → `<: Vertex`) if you extend [`DecayChain`](@ref).
3. In downstream packages that re-export the old name (e.g.
   `import ThreeBodyDecays: VertexFunction as Vertex`), switch to `import ThreeBodyDecays: Vertex`.

No changes are required to recoupling or form-factor implementations—only the wrapper type name changed.

## Related docs

- [Energy dependence of the amplitude](@ref energy_dependence) — constructing vertices with form factors
- [`Vertex`](@ref) — API reference
- [`DecayChain`](@ref) — where `HRk` and `Hij` vertices are used
