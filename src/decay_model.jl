"""
    ThreeBodyDecay

Model for a three-body decay amplitude defined as a coherent sum of decay chains.

The model stores:
- `chains`: a concretely typed `Tuple` of [`AbstractDecayChain`](@ref) objects,
- `couplings`: complex (or real) coefficients multiplying each chain,
- `names`: labels for the chains (typically resonance names).

Evaluate with [`amplitude`](@ref) and summarize over spin projections with
[`unpolarized_intensity`](@ref).
"""
struct ThreeBodyDecay{N,Ch<:Tuple{Vararg{AbstractDecayChain}},L,S}
    chains::Ch
    couplings::NTuple{N,L}
    names::NTuple{N,S}
end

"""
ThreeBodyDecay(chains, couplings, names)

Constructs a `ThreeBodyDecay` object with the given parameters.

# Arguments
- `chains`: Chains involved in the decay, as a tuple or vector. The length should match `couplings` and `names`.
- `couplings`: Coupling constants for each chain, as a tuple or vector.
- `names`: Names for each chain, as a tuple or vector.

# Returns
- A `ThreeBodyDecay` object with the specified chains, couplings, and names.

# Examples
```julia
ThreeBodyDecay(
    (chain1, chain2, chain3),
    (1.0, -1.0, 0.2im),
    ("L1405", "L1405", "K892"),
)
```
"""
function ThreeBodyDecay(
    chains::Tuple{Vararg{AbstractDecayChain}},
    couplings::Tuple{Vararg{Number}},
    names::Tuple{Vararg{AbstractString}},
)
    N = length(chains)
    length(couplings) == N == length(names) ||
        throw(ArgumentError("The lengths of chains, couplings, and names must be equal"))
    L = promote_type(typeof.(couplings)...)
    coupling_tuple = ntuple(i -> convert(L, couplings[i]), N)
    name_tuple = ntuple(i -> String(names[i]), N)
    return ThreeBodyDecay{N,typeof(chains),L,eltype(name_tuple)}(
        chains,
        coupling_tuple,
        name_tuple,
    )
end

function ThreeBodyDecay(
    chains::AbstractVector{<:AbstractDecayChain},
    couplings::AbstractVector{<:Number},
    names::AbstractVector{<:AbstractString},
)
    N = length(chains)
    length(couplings) == N == length(names) ||
        throw(ArgumentError("The lengths of chains, couplings, and names must be equal"))
    return ThreeBodyDecay((chains...,), (couplings...,), (names...,))
end

"""
ThreeBodyDecay(descriptor)

Constructs a `ThreeBodyDecay` object using one argument, a descriptor.
The `descriptor` is a list of pairs, `names .=> zip(couplings, chains)`.

# Examples
```julia
ThreeBodyDecay("K892" .=> zip([1.0, -1.0, 0.2im], [chain1, chain2, chain3]))
```
"""
ThreeBodyDecay(descriptor::Pair) = ThreeBodyDecay([descriptor])

function ThreeBodyDecay(descriptor)
    names = first.(descriptor)
    cd = getindex.(descriptor, 2)
    couplings = first.(cd)
    chains = getindex.(cd, 2)
    ThreeBodyDecay((chains...,), (couplings...,), (names...,))
end

amplitude(model::ThreeBodyDecay, p...; kw...) =
    sum(c * amplitude(d, p...; kw...) for (c, d) in zip(model.couplings, model.chains))

import Base: getindex, length
import LinearAlgebra: adjoint

adjoint(c::NTuple{N,L}) where {N,L<:Number} = adjoint(SVector(c))
function Base.:*(a::Adjoint{T,<:AbstractVector}, c::NTuple{N,L}) where {T,N,L<:Number}
    return a * SVector(c)
end
function Base.:*(a::Transpose{T,<:AbstractVector}, c::NTuple{N,L}) where {T,N,L<:Number}
    return a * SVector(c)
end

function _submodel(model::ThreeBodyDecay, idx::Tuple{Vararg{Int}})
    N = length(idx)
    chains = ntuple(i -> model.chains[idx[i]], N)
    couplings = ntuple(i -> model.couplings[idx[i]], N)
    names = ntuple(i -> model.names[idx[i]], N)
    ThreeBodyDecay(chains, couplings, names)
end

getindex(model::ThreeBodyDecay, i::Integer) = _submodel(model, (i,))
getindex(model::ThreeBodyDecay, r::AbstractRange{<:Integer}) = _submodel(model, Tuple(r))
getindex(model::ThreeBodyDecay, idx::AbstractVector{<:Integer}) = _submodel(model, Tuple(idx))
getindex(model::ThreeBodyDecay, mask::Union{NTuple{N,Bool},AbstractVector{Bool}} where {N}) =
    getindex(model, findall(mask))

length(model::ThreeBodyDecay{N}) where {N} = N

"""
    system(model::ThreeBodyDecay) -> ThreeBodySystem

Return the [`ThreeBodySystem`](@ref) associated with `model` (taken from the first chain).
"""
system(model::ThreeBodyDecay) = first(model.chains).tbs

"""
    masses(model::ThreeBodyDecay) -> MassTuple

Convenience wrapper for `masses(system(model))`.
"""
masses(model::ThreeBodyDecay) = masses(system(model))

"""
    spins(model::ThreeBodyDecay) -> SpinTuple

Convenience wrapper for `spins(system(model))`.
"""
spins(model::ThreeBodyDecay) = spins(system(model))

"""
unpolarized_intensity(model::ThreeBodyDecay, σs; kw...)

Computes squared amplitude summed over spin projections.
"""
unpolarized_intensity(model, σs; kw...) = sum(abs2, amplitude(model, σs; kw...))

"""
    Base.vcat(models::ThreeBodyDecay...)

Concatenates multiple `ThreeBodyDecay` objects into a single `ThreeBodyDecay`.
Argument is variable number of `ThreeBodyDecay` objects.

# An example
```julia
extended_model = vcat(model[2], model[2:3], model)
```
"""
function Base.vcat(models::ThreeBodyDecay...)
    chains = mapreduce(m -> m.chains, (a, b) -> (a..., b...), models)
    couplings = mapreduce(m -> m.couplings, (a, b) -> (a..., b...), models)
    names = mapreduce(m -> m.names, (a, b) -> (a..., b...), models)
    ThreeBodyDecay(chains, couplings, names)
end
