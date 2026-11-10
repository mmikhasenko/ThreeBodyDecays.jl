"""
    ThreeBodyDecay

Model for a three-body decay amplitude defined as a coherent sum of decay chains.

The model stores:
- `chains`: a concretely typed `Tuple` of [`AbstractDecayChain`](@ref) objects,
- `couplings`: complex (or real) coefficients multiplying each chain,
- `names`: labels for the chains (typically resonance names),
- `reference_k`: the common reference topology for all helicities and couplings.

The reference is required at construction and cannot be overridden at evaluation.

Evaluate with [`amplitude`](@ref) and summarize over spin projections with
[`unpolarized_intensity`](@ref).
"""
struct ThreeBodyDecay{N,Ch<:Tuple{Vararg{AbstractDecayChain}},L,S}
    chains::Ch
    couplings::NTuple{N,L}
    names::NTuple{N,S}
    reference_k::Int
    function ThreeBodyDecay{N,Ch,L,S}(
        chains::Ch,
        couplings::NTuple{N,L},
        names::NTuple{N,S},
        reference_k::Integer,
    ) where {N,Ch<:Tuple{Vararg{AbstractDecayChain}},L,S}
        N > 0 || throw(ArgumentError("a model must contain at least one chain"))
        length(chains) == N ||
            throw(ArgumentError("chains, couplings, and names must have equal lengths"))
        # Custom AbstractDecayChain implementations need not expose a system.
        all(chain -> chain isa DecayChain, chains) &&
            !all(
                chain ->
                    masses(chain) == masses(first(chains)) &&
                    spins(chain) == spins(first(chains)),
                chains,
            ) &&
            throw(ArgumentError("all chains must share external masses and spins"))
        new{N,Ch,L,S}(chains, couplings, names, _validate_reference_k(reference_k))
    end
end

function _validate_reference_k(k::Integer)
    (k isa Bool || !(1 <= k <= 3)) && throw(ArgumentError("reference_k must be 1, 2, or 3"))
    return Int(k)
end

"""
    reference_k(model::ThreeBodyDecay) -> Int

Return the common reference topology stored with the model's coupling coefficients.
The cyclic topologies are `1 => ((2,3),1)`, `2 => ((3,1),2)`, and `3 => ((1,2),3)`.
"""
reference_k(model::ThreeBodyDecay) = model.reference_k

"""
ThreeBodyDecay(chains, couplings, names; reference_k)

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
    ("L1405", "L1405", "K892");
    reference_k = 1,
)
```
"""
function ThreeBodyDecay(
    chains::Tuple{Vararg{AbstractDecayChain}},
    couplings::Tuple{Vararg{Number}},
    names::Tuple{Vararg{AbstractString}};
    reference_k::Integer,
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
        reference_k,
    )
end

function ThreeBodyDecay(
    chains::AbstractVector{<:AbstractDecayChain},
    couplings::AbstractVector{<:Number},
    names::AbstractVector{<:AbstractString};
    reference_k::Integer,
)
    N = length(chains)
    length(couplings) == N == length(names) ||
        throw(ArgumentError("The lengths of chains, couplings, and names must be equal"))
    return ThreeBodyDecay((chains...,), (couplings...,), (names...,); reference_k)
end

"""
ThreeBodyDecay(descriptor; reference_k)

Constructs a `ThreeBodyDecay` object using one argument, a descriptor.
The `descriptor` is a list of pairs, `names .=> zip(couplings, chains)`.
The keyword `reference_k` is required and must be 1, 2, or 3.

# Examples
```julia
ThreeBodyDecay("K892" .=> zip([1.0, -1.0, 0.2im], [chain1, chain2, chain3]); reference_k = 1)
```
"""
ThreeBodyDecay(descriptor::Pair; reference_k::Integer) =
    ThreeBodyDecay([descriptor]; reference_k)

function ThreeBodyDecay(descriptor; reference_k::Integer)
    names = first.(descriptor)
    cd = getindex.(descriptor, 2)
    couplings = first.(cd)
    chains = getindex.(cd, 2)
    ThreeBodyDecay((chains...,), (couplings...,), (names...,); reference_k)
end

amplitude(model::ThreeBodyDecay, p...) = sum(
    c * amplitude(d, p...; reference_k = reference_k(model)) for
    (c, d) in zip(model.couplings, model.chains)
)

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
    ThreeBodyDecay(chains, couplings, names; reference_k = reference_k(model))
end

getindex(model::ThreeBodyDecay, i::Integer) = _submodel(model, (i,))
getindex(model::ThreeBodyDecay, r::AbstractRange{<:Integer}) = _submodel(model, Tuple(r))
getindex(model::ThreeBodyDecay, idx::AbstractVector{<:Integer}) =
    _submodel(model, Tuple(idx))
getindex(
    model::ThreeBodyDecay,
    mask::Union{NTuple{N,Bool},AbstractVector{Bool}} where {N},
) = getindex(model, findall(mask))

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
    unpolarized_intensity(model::ThreeBodyDecay, σs)
    unpolarized_intensity(model::ThreeBodyDecay, angles, σs)

Computes squared amplitude summed over spin projections. Model evaluation uses the
stored reference topology and accepts no reference override.
"""
unpolarized_intensity(model, p...; kw...) = sum(abs2, amplitude(model, p...; kw...))

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
    isempty(models) && throw(ArgumentError("provide at least one model"))
    k = reference_k(first(models))
    all(model -> reference_k(model) == k, models) || throw(
        ArgumentError("cannot concatenate models with different reference topologies"),
    )
    chains = mapreduce(m -> m.chains, (a, b) -> (a..., b...), models)
    couplings = mapreduce(m -> m.couplings, (a, b) -> (a..., b...), models)
    names = mapreduce(m -> m.names, (a, b) -> (a..., b...), models)
    ThreeBodyDecay(chains, couplings, names; reference_k = k)
end

ThreeBodyDecay(; chains, couplings, names, reference_k::Integer) =
    ThreeBodyDecay(chains, couplings, names; reference_k)
