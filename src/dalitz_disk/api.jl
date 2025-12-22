# High-level API for Dalitz to Unit Circle Mapping
# This file defines the type-stable API structure

# Include boundary module (defines BoundaryFunction)
include("boundary.jl")

# Include domain module (defines DemocraticMap)
include("domain.jl")

# Include zipper module (defines ZipperMap)
include("zipper.jl")

# Include normalize module (defines normalization functions)
include("normalize.jl")

# DalitzMap is defined below, so diagnostics.jl needs to be included after
# But we'll include it at the end after DalitzMap is defined

# Export normalize_landmarks for use in dalitzmap
# (Other functions are internal)

# DemocraticMap struct and methods are defined in domain.jl

# ZipperMap struct and methods are defined in zipper.jl

"""
    MobiusTransform{T}

A callable struct representing a Möbius transformation of the unit disk.

Maps three points to specified targets, used for channel-democratic normalization.

# Fields
- `a::Complex{T}`, `b::Complex{T}`, `c::Complex{T}`, `d::Complex{T}`: Möbius coefficients

# Formula
A(w) = (aw + b) / (cw + d)
"""
struct MobiusTransform{T}
    a::Complex{T}
    b::Complex{T}
    c::Complex{T}
    d::Complex{T}
end

"""
    (mt::MobiusTransform)(w::Complex{T}) -> Complex{T}

Apply Möbius transformation.

# Formula
A(w) = (aw + b) / (cw + d)
"""
function (mt::MobiusTransform)(w::Complex{T}) where {T}
    return (mt.a * w + mt.b) / (mt.c * w + mt.d)
end

"""
    DalitzMap{T}

Top-level struct representing the complete conformal mapping from Dalitz region
to normalized unit disk.

# Fields
- `ms::MassTuple{T}`: Masses of the three-body system
- `boundary::BoundaryFunction{T}`: Boundary parameterization
- `democratic::DemocraticMap{T}`: Democratic coordinate transformation
- `zipper::ZipperMap{T}`: Zipper conformal mapping
- `normalization::MobiusTransform{T}`: Channel-democratic normalization

# Example
```julia
ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0=4.0)
map = dalitzmap(ms)
w = map(σs)  # Map Mandelstam invariants to unit disk
```
"""
struct DalitzMap{T}
    ms::MassTuple{T}
    boundary::BoundaryFunction{T}
    democratic::DemocraticMap{T}
    zipper::ZipperMap{T}
    normalization::MobiusTransform{T}
end

"""
    (map::DalitzMap)(σs::MandelstamTuple) -> Complex{T}

Apply the complete mapping: Mandelstam invariants → normalized unit disk.

# Formula
w = normalization(zipper(democratic(σs)))
"""
function (map::DalitzMap)(σs::MandelstamTuple)
    z = map.democratic(σs)
    w_raw = map.zipper(z)
    return map.normalization(w_raw)
end

"""
    (map::DalitzMap)(s12::Real, s23::Real) -> Complex{T}

Convenience method: map from (s₁₂, s₂₃) coordinates.

# Arguments
- `s12::Real`: Invariant mass squared of particles 1 and 2 (σ₃)
- `s23::Real`: Invariant mass squared of particles 2 and 3 (σ₁)

# Returns
- `Complex{T}`: Mapped point in normalized unit disk

# Note
This computes σ₂ from the sum rule: σ₁ + σ₂ + σ₃ = Σ
"""
function (map::DalitzMap)(s12::Real, s23::Real)
    # s12 = σ₃, s23 = σ₁
    # Compute σ₂ from sum rule
    Σ = sum(map.ms^2)
    σ2 = Σ - s12 - s23
    σs = Invariants(map.ms; σ1 = s23, σ2 = σ2, σ3 = s12)
    return map(σs)
end

"""
    dalitzmap(ms::MassTuple; method=:zipper) -> DalitzMap

Construct a conformal mapping from the Dalitz region to the unit disk.

# Arguments
- `ms::MassTuple`: Masses of the three-body system
- `method::Symbol`: Mapping method (currently only `:zipper` supported)

# Returns
- `DalitzMap{T}`: Complete mapping structure

# Example
```julia
ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0=4.0)
map = dalitzmap(ms)
w = map(σs)  # |w| < 1 for physical points
```
"""
function dalitzmap(ms::MassTuple{T}; method::Symbol = :zipper) where {T}
    # Construct all components
    bf = BoundaryFunction(ms)
    dm = DemocraticMap(ms)
    zm = ZipperMap(bf, dm; n = 200)
    mt = normalize_landmarks(bf, dm, zm)

    return DalitzMap(ms, bf, dm, zm, mt)
end

# Include diagnostics module after DalitzMap is defined
include("diagnostics.jl")
