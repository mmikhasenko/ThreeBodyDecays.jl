@with_kw struct ThreeBodyDecay{N,T<:AbstractDecayChain,L<:Number,S<:AbstractString}
    chains::SVector{N,T}
    couplings::SVector{N,L}
    names::SVector{N,S}
end

"""
	ThreeBodyDecay(; chains, couplings, names)

Constructs a `ThreeBodyDecay` object with the given parameters.

# Arguments
- `chains`: An array of chains involved in the decay. The length of this array should match the lengths of `couplings` and `names`.
- `couplings`: An array of coupling constants for each chain in the decay. The length of this array should match the lengths of `chains` and `names`.
- `names`: An array of names for each chain, or names of resonances in the decay. The length of this array should match the lengths of `chains` and `couplings`.

# Returns
- A `ThreeBodyDecay` object with the specified chains, couplings, and names.

# Examples
```julia
ThreeBodyDecay(
	chains=[chain1, chain2, chain3],
	couplings=[1.0, -1.0, 0.2im],
	names=["L1405", "L1405", "K892"])
```
"""
function ThreeBodyDecay(
    chains::Vector{T}, couplings::Vector{L}, names::Vector{S}) where {T<:AbstractDecayChain,L<:Number,S<:AbstractString}
    N = length(chains)
    @assert length(couplings) == N && length(names) == N "The lengths of chains, couplings, and names must be equal"
    return ThreeBodyDecay(SVector{N,T}(chains), SVector{N,L}(couplings), SVector{N,S}(names))
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
function ThreeBodyDecay(descriptor)
    N = length(descriptor)
    #
    names = first.(descriptor)
    cd = getindex.(descriptor, 2)
    couplings = first.(cd)
    chains = getindex.(cd, 2)
    #
    N = length(chains)
    ThreeBodyDecay(SVector{N}(chains), SVector{N}(couplings), SVector{N}(names))
end

amplitude(model::ThreeBodyDecay, p...; kw...) =
    sum(c * amplitude(d, p...; kw...)
        for (c, d) in zip(model.couplings, model.chains))

import Base: getindex, length

function getindex(model::ThreeBodyDecay, key...)
    description = getindex(model.names, key...) .=> zip(
        getindex(model.couplings, key...),
        getindex(model.chains, key...))
    ThreeBodyDecay(description)
end

length(model::ThreeBodyDecay{N}) where {N} = length(model.chains)
system(model::ThreeBodyDecay) = first(model.chains).tbs
masses(model::ThreeBodyDecay) = masses(system(model))
spins(model::ThreeBodyDecay) = spins(system(model))

"""
	unpolarized_intensity(model::ThreeBodyDecay, σs; kw...)

Computes squared amplitude summed over spin projections.
"""
unpolarized_intensity(model, σs; kw...) =
    sum(abs2, amplitude(model, σs, two_λs; kw...)
              for two_λs in itr(spins(model)))


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
    names = vcat(getproperty.(models, :names)...)
    couplings = vcat(getproperty.(models, :couplings)...)
    chains = vcat(getproperty.(models, :chains)...)
    ThreeBodyDecay(names .=> zip(couplings, chains))
end
