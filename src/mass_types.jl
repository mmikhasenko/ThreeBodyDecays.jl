"""
    MassTuple{T}

A named tuple representing the masses of a three-body system.
Contains masses m₁, m₂, m₃ of the decay products and m₀ of the parent particle.
"""
const MassTuple{T} = NamedTuple{(:m1, :m2, :m3, :m0), NTuple{4, T}}

"""
    ThreeBodyMasses(m1, m2, m3; m0)
    ThreeBodyMasses(; m1, m2, m3, m0)

Construct a MassTuple for a three-body system.

# Arguments
- `m1`, `m2`, `m3`: Masses of the decay products
- `m0`: Mass of the parent particle

# Returns
- `MassTuple{T}`: A named tuple containing the masses

# Throws
- `ErrorException` if m₀ is less than the sum of m₁, m₂, and m₃
"""
function ThreeBodyMasses(m1, m2, m3; m0)
    tm0 = typeof(m0)
    tm0 <: Number && (m0 < m1 + m2 + m3) && error("m₀ should be bigger than m₁+m₂+m₃")
    MassTuple{tm0}((m1, m2, m3, m0))
end

ThreeBodyMasses(; m1, m2, m3, m0) = ThreeBodyMasses(m1, m2, m3; m0)

"""
    lims(ms::MassTuple; k::Int)
    lims(k::Int, ms::MassTuple)
    lims1(ms::MassTuple)
    lims2(ms::MassTuple)
    lims3(ms::MassTuple)

Calculate the kinematic limits (boundaries) for the Mandelstam variables σᵢ in a three-body decay.
For each pair of particles (i,j), the invariant mass squared σₖ = (pᵢ + pⱼ)² must lie within physical limits:
- Lower bound: (mᵢ + mⱼ)² (threshold for producing particles i and j)
- Upper bound: (M - mₖ)² (maximum energy available when particle k recoils)

# Arguments
- `ms`: Tuple of masses (m₁,m₂,m₃,M)
- `k`: Index specifying which pair of particles (1,2,3)

# Returns
- Tuple (min,max) of the allowed range for σₖ

# Example
```julia
ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0=4.0)
lims1(ms)  # limits for σ₁ = (p₂ + p₃)²
lims2(ms)  # limits for σ₂ = (p₃ + p₁)²
lims3(ms)  # limits for σ₃ = (p₁ + p₂)²
```

See also [`isphysical`](@ref), [`Kibble`](@ref).
"""
function lims(ms::MassTuple; k::Int)
    i, j = ij_from_k(k)
    ((ms[i] + ms[j])^2, (ms[4] - ms[k])^2)
end
lims(k::Int, ms::MassTuple) = lims(ms; k)
lims1(ms::MassTuple) = lims(ms; k = 1)
lims2(ms::MassTuple) = lims(ms; k = 2)
lims3(ms::MassTuple) = lims(ms; k = 3)

import Base: getindex, ^, length, iterate
^(ms::MassTuple, i::Int) = Tuple(ms) .^ i
