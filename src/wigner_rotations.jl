"""
    AbstractWignerRotation

Abstract type for representing Wigner rotations in three-body decays.
Subtypes include `TrivialWignerRotation` and `WignerRotation{N}` for N = 0,2,3.
"""
abstract type AbstractWignerRotation end

"""
    TrivialWignerRotation <: AbstractWignerRotation

Represents a trivial Wigner rotation (identity transformation).
Used when the reference frame and system frame are the same.

# Fields
- `k::Int`: Index of the spectator particle
"""
struct TrivialWignerRotation <: AbstractWignerRotation
    k::Int
end

"""
    WignerRotation{N} <: AbstractWignerRotation

Represents a Wigner rotation in three-body decays.
The type parameter N determines the kind of rotation:
- N=0: Rotation between total system frames
- N=2: Rotation involving two distinct indices
- N=3: Rotation involving three distinct indices

# Fields
- `k::Int`: Index of the spectator particle
- `ispositive::Bool`: Direction of rotation
- `iseven::Bool`: Parity of the rotation
"""
struct WignerRotation{N} <: AbstractWignerRotation
    k::Int
    ispositive::Bool
    iseven::Bool
end

const Arg0WignerRotation = WignerRotation{0}
const Arg2WignerRotation = WignerRotation{2}
const Arg3WignerRotation = WignerRotation{3}

WignerRotation{N}(k::Int, ispositive::Bool) where {N} =
    WignerRotation{N}(k, ispositive, true)

"""
    ispositive(wr::AbstractWignerRotation)

Determine if a Wigner rotation is in the positive direction.

# Arguments
- `wr`: A Wigner rotation object

# Returns
- `true` for positive rotations, `false` otherwise
"""
ispositive(wr::TrivialWignerRotation) = true
ispositive(wr::WignerRotation) = wr.ispositive

import Base: iseven
iseven(wr::TrivialWignerRotation) = true
iseven(wr::WignerRotation{0}) = true
iseven(wr::WignerRotation{3}) = true
iseven(wr::WignerRotation{2}) = wr.iseven

ijk(wr::AbstractWignerRotation) = ijk(wr.k)
issequential(i, j) = (j - i) ∈ (1, -2)

"""
    wr(system_a, reference_b, particle_c=0)

Create a WignerRotation object for transforming between different reference frames in a three-body system.

# Arguments
- `system_a`: Index of the isobar system being considered (1,2,3)
- `reference_b`: Index of the reference system (1,2,3)
- `particle_c`: Index of the particle (0 for parent particle, 1,2,3 for daughters)

# Returns
- A WignerRotation object of the appropriate type

# Example
```julia
# Rotation from system 2 to system 1 for particle 1 => zeta_21_for1
w = wr(2, 1, 1)
# Rotation between total system frames => zeta_12_for0
w0 = wr(1, 2, 0)
```
"""
function wr(system_a, reference_b, particle_c = 0)
    system_a == reference_b && return TrivialWignerRotation(particle_c)
    S = issequential(system_a, reference_b)
    A, B = S ? (system_a, reference_b) : (reference_b, system_a)
    #
    particle_c == 0 && return Arg0WignerRotation(A, S)
    #
    particle_c ∉ (system_a, reference_b) && return Arg3WignerRotation(particle_c, S)
    #
    T = (particle_c == A)
    return Arg2WignerRotation(particle_c, !(S), T)
end

"""
    cosζ(wr::AbstractWignerRotation, σs, msq)

Calculate the cosine of the Wigner rotation angle ζ for a given kinematic configuration.
Different methods are implemented for different types of Wigner rotations.

# Arguments
- `wr`: A Wigner rotation object
- `σs`: Tuple of Mandelstam variables
- `msq`: Tuple of squared masses

# Returns
- Cosine of the Wigner rotation angle

# Example
```julia
ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0=4.0)
σs = x2σs([0.5, 0.5], ms; k=1)
w = wr(2, 1, 1)
cosζ(w, σs, ms^2)  # Get rotation angle cosine
```
"""
cosζ(wr::TrivialWignerRotation, σs, msq) = one(σs[1])

function cosζ(wr::Arg0WignerRotation, σs, msq)
    i, j, k = ijk(wr)
    #
    s = msq[4]
    EE4s = (s + msq[i] - σs[i]) * (s + msq[k] - σs[k])
    pp4s = sqrt(Kallen(s, msq[i], σs[i]) * Kallen(s, msq[k], σs[k]))
    rest = σs[j] - msq[i] - msq[k]
    return (EE4s - 2s * rest) / pp4s
end

function cosζ(wr::Arg2WignerRotation, σs, msq)
    i, j, k = ijk(wr)
    #
    if !(iseven(wr))
        i, j = j, i
    end
    #
    msq[k] ≈ 0 && return one(σs[1])
    #
    s = msq[4]
    EE4mksq = (s + msq[k] - σs[k]) * (σs[i] - msq[k] - msq[j])
    pp4mksq = sqrt(Kallen(s, msq[k], σs[k]) * Kallen(msq[k], msq[j], σs[i]))
    rest = σs[j] - s - msq[j]
    return (2msq[k] * rest + EE4mksq) / pp4mksq
end

#
function cosζ(wr::Arg3WignerRotation, σs, msq)
    i, j, k = ijk(wr)
    #
    msq[k] ≈ 0 && return one(σs[1])
    #
    s = msq[4]
    EE4m1sq = (σs[i] - msq[j] - msq[k]) * (σs[j] - msq[k] - msq[i])
    pp4m1sq = sqrt(Kallen(σs[i], msq[j], msq[k]) * Kallen(σs[j], msq[k], msq[i]))
    rest = msq[i] + msq[j] - σs[k]
    return (2msq[k] * rest + EE4m1sq) / pp4m1sq
end

# explicit
cosζ21_for1(σs, ms²) = cosζ(wr(2, 1, 1), σs, ms²)
cosζ21_for2(σs, ms²) = cosζ(wr(2, 1, 2), σs, ms²)
cosζ13_for1(σs, ms²) = cosζ(wr(1, 3, 1), σs, ms²)
cosζ13_for3(σs, ms²) = cosζ(wr(1, 3, 3), σs, ms²)
cosζ32_for3(σs, ms²) = cosζ(wr(3, 2, 3), σs, ms²)
cosζ32_for2(σs, ms²) = cosζ(wr(3, 2, 2), σs, ms²)

cosζ12_for3(σs, ms²) = cosζ(wr(1, 2, 3), σs, ms²)
cosζ23_for1(σs, ms²) = cosζ(wr(2, 3, 1), σs, ms²)
cosζ31_for2(σs, ms²) = cosζ(wr(3, 1, 2), σs, ms²)

cosζ12_for0(σs, ms²) = cosζ(wr(1, 2, 0), σs, ms²)
cosζ23_for0(σs, ms²) = cosζ(wr(2, 3, 0), σs, ms²)
cosζ31_for0(σs, ms²) = cosζ(wr(3, 1, 0), σs, ms²)

#
cosζk1_for1(k, σs, ms²) = cosζ(wr(k, 1, 1), σs, ms²)
cosζk2_for2(k, σs, ms²) = cosζ(wr(k, 2, 2), σs, ms²)
cosζk3_for3(k, σs, ms²) = cosζ(wr(k, 3, 3), σs, ms²)

"""
Phase for wigner d-functions for clockwise rotations
"""
phase(two_λ1_minus_λ2) =
    (abs(two_λ1_minus_λ2) % 4 == 2 ? -one(two_λ1_minus_λ2) : one(two_λ1_minus_λ2))
phase(two_λ1, two_λ2) = phase(two_λ1 - two_λ2)
