@with_kw struct ThreeBodyDecay{N, T, L<:Number}
    chains::SVector{N,T}
    couplings::SVector{N,L}
    names::SVector{N,String}
end

const VectPairStringChain = Vector{Pair{String, Tuple{F, DecayChain{X,T}}}} where {F<:Number, X,T}
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
function ThreeBodyDecay(descriptor::VectPairStringChain)
    N = length(descriptor)

    #
    names = first.(descriptor)
    cd = getindex.(descriptor, 2)
    couplings = first.(cd)
    chains = getindex.(cd, 2)
    # 
    sv_chains = (SVector{N})(chains)
    sv_couplings = SVector{N}(couplings)
    sv_names = SVector{N}(names)
    #
    ThreeBodyDecay(sv_chains, sv_couplings, sv_names)
end

amplitude(model::ThreeBodyDecay, p...; kw...) =
    sum(c * amplitude(d, p...; kw...)
        for (c, d) in zip(model.couplings, model.chains))

import Base: getindex, length

getindex(model::ThreeBodyDecay, key...) = 
    ThreeBodyDecay(
        getindex(model.names, key...) .=> zip(
            getindex(model.couplings, key...),
            getindex(model.chains, key...)))

length(model::ThreeBodyDecay{N}) where N = length(model.chains)
masses(model::ThreeBodyDecay) = masses(first(model.chains).tbs)
spins(model::ThreeBodyDecay) = spins(first(model.chains).tbs)

"""
    unpolarized_intensity(model::ThreeBodyDecay, σs; kw...)

Computes squared amplitude summed over spin projections.
"""
unpolarized_intensity(model::ThreeBodyDecay, σs; kw...) =
    sum(abs2, amplitude(model, σs, two_λs; kw...)
              for two_λs in itr(spins(model)))