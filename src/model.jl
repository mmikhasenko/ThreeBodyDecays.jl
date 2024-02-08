@with_kw struct ThreeBodyDecay{N, T, L<:Number}
    names::SVector{N,String}
    chains::SVector{N,T}
    couplings::SVector{N,L}
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
ThreeBodyDecay(chains=[chain1, chain2, chain3], couplings=[1.0, -1.0, 0.2im],
    names=["L1405", "L1405", "K892"])
```
"""
function ThreeBodyDecay(; chains, couplings, names)
    N = length(chains)
    N != length(couplings) && error("Length of couplings does not match the length of the chains")
    N != length(names) && error("Length of names does not match the length of the chains")
    #
    sv_chains = (SVector{N,ChainΛb2Λγ{T,R1,R2} where {T,R1,R2}})(chains)
    sv_couplings = SVector{N}(couplings)
    sv_names = SVector{N}(names)
    #
    ThreeBodyDecay(sv_chains, sv_couplings, sv_names)
end

amplitude(model::ThreeBodyDecay, σs, two_λs) =
    sum(c * amplitude(d, σs, two_λs)
        for (c, d) in zip(model.couplings, model.chains))

import Base: getindex, length

getindex(model::ThreeBodyDecay, key...) = 
    ThreeBodyDecay(;
        chains=getindex(model.chains, key...),
        couplings=getindex(model.couplings, key...),
        names=getindex(model.names, key...))

length(model::ThreeBodyDecay{N}) where N = length(model.chains)
masses(model::ThreeBodyDecay) = masses(first(model.chains).tbs)
spins(model::ThreeBodyDecay) = spins(first(model.chains).tbs)

masses(tbs::ThreeBodySystem) = tbs.ms
spins(tbs::ThreeBodySystem) = tbs.two_js

export masses, spins
export ThreeBodyDecay