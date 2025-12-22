# Democratic coordinate transformation for Dalitz region

"""
    DemocraticMap{T}

A callable struct for the democratic coordinate transformation.

Democratic coordinates map Mandelstam invariants to a complex plane using
cube-root-weighted linear combination, preserving D₃ symmetry.

# Fields
- `ms::MassTuple{T}`: Masses of the three-body system
- `Λ::T`: Normalization parameter
- `Σ::T`: Sum of squared masses (m₁² + m₂² + m₃² + M²)
"""
struct DemocraticMap{T}
    ms::MassTuple{T}
    Λ::T
    Σ::T
end

"""
    DemocraticMap(ms::MassTuple, Λ::Real=nothing)

Construct a `DemocraticMap` for the given mass configuration and normalization.

# Arguments
- `ms::MassTuple`: Masses of the three-body system
- `Λ::Real`: Normalization parameter (optional, will be computed if not provided)

# Formula
The democratic coordinate is defined as:
z = (1/Λ) Σᵢ (σᵢ - σ̄) ω^(i-1)
where σ̄ = Σ/3, ω = exp(2πi/3), and Σ = m₁² + m₂² + m₃² + M²

# Returns
- `DemocraticMap{T}`: Callable struct for democratic coordinate transformation

# Example
```julia
ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0=4.0)
dm = DemocraticMap(ms)  # Λ computed automatically
z = dm(σs)  # Forward: Mandelstam → complex
```
"""
function DemocraticMap(ms::MassTuple{T}, Λ::Union{T,Nothing} = nothing) where {T}
    Σ = sum(ms^2)
    if Λ === nothing
        Λ = compute_Λ(ms)
    end
    return DemocraticMap(ms, Λ, Σ)
end

"""
    compute_Λ(ms::MassTuple, bf=nothing; n::Int=200)

Compute normalization parameter Λ by finding maximum distance from center
along the boundary.

# Formula
Λ = max_{θ ∈ [0, 2π)} |Σᵢ (σᵢ(θ) - σ̄) ω^(i-1)|

# Arguments
- `ms::MassTuple`: Masses of the system
- `bf`: Boundary function (optional, will be created if not provided)
- `n::Int`: Number of sample points

# Returns
- `T`: Normalization parameter Λ
"""
function compute_Λ(ms::MassTuple{T}, bf = nothing; n::Int = 200) where {T}
    if bf === nothing
        # BoundaryFunction will be available when domain.jl is included after boundary.jl
        bf = BoundaryFunction(ms)
    end
    Σ = sum(ms^2)
    σ̄ = Σ / 3
    ω = exp(2π * im / 3)

    # Sample boundary and find maximum |z|
    θs = range(0, 2π, length = n+1)[1:(end-1)]
    max_abs_z = zero(T)

    for θ in θs
        σs = bf(θ)
        # Compute z without normalization
        z_raw = (σs.σ1 - σ̄) + (σs.σ2 - σ̄) * ω + (σs.σ3 - σ̄) * ω^2
        max_abs_z = max(max_abs_z, abs(z_raw))
    end

    return max_abs_z
end

"""
    (dm::DemocraticMap{T})(σs::MandelstamTuple) -> Complex{T}

Forward transformation: Mandelstam invariants → democratic coordinates.

# Formula
z = (1/Λ) Σᵢ (σᵢ - σ̄) ω^(i-1)
where:
- σ̄ = Σ/3 is the mean invariant
- ω = exp(2πi/3) is the cube root of unity
- Λ is the normalization parameter

# Arguments
- `σs::MandelstamTuple`: Mandelstam invariants

# Returns
- `Complex{T}`: Democratic coordinate z
"""
function (dm::DemocraticMap{T})(σs::MandelstamTuple) where {T}
    σ̄ = dm.Σ / 3
    ω = exp(2π * im / 3)

    # Compute z = (1/Λ) Σᵢ (σᵢ - σ̄) ω^(i-1)
    z = ((σs.σ1 - σ̄) + (σs.σ2 - σ̄) * ω + (σs.σ3 - σ̄) * ω^2) / dm.Λ

    return z
end

"""
    (dm::DemocraticMap)(z::Complex{T}) -> MandelstamTuple{T}

Inverse transformation: democratic coordinates → Mandelstam invariants.

# Formula
σ₁ = σ̄ + Λ·Re(z)
σ₂ = σ̄ + Λ·Re(ω⁻¹z)
σ₃ = Σ - σ₁ - σ₂

where ω⁻¹ = ω² = exp(-2πi/3)

# Arguments
- `z::Complex{T}`: Democratic coordinate

# Returns
- `MandelstamTuple{T}`: Mandelstam invariants
"""
function (dm::DemocraticMap{T})(z::Complex{T}) where {T}
    σ̄ = dm.Σ / 3

    # Recover invariants from z = (1/Λ)[(σ₁ - σ̄) + (σ₂ - σ̄)ω + (σ₃ - σ̄)ω²]
    # Expanding: Λz = (σ₁ - σ̄) + (σ₂ - σ̄)ω + (σ₃ - σ̄)ω²
    # Real part: Re(Λz) = (σ₁ - σ̄) + (σ₂ - σ̄)Re(ω) + (σ₃ - σ̄)Re(ω²)
    #           = (σ₁ - σ̄) - (σ₂ - σ̄)/2 - (σ₃ - σ̄)/2
    #           = (σ₁ - σ̄) - (σ₂ + σ₃ - 2σ̄)/2
    #           = (σ₁ - σ̄) - (Σ - σ₁ - 2σ̄)/2
    #           = (3σ₁ - Σ)/2
    # So: σ₁ = (2Re(Λz) + Σ)/3

    # Imaginary part: Im(Λz) = (σ₂ - σ̄)Im(ω) + (σ₃ - σ̄)Im(ω²)
    #                = (σ₂ - σ̄)√3/2 - (σ₃ - σ̄)√3/2
    #                = √3/2(σ₂ - σ₃)
    # So: σ₂ - σ₃ = 2Im(Λz)/√3

    # Combined with σ₁ + σ₂ + σ₃ = Σ:
    # σ₂ = (Σ - σ₁ + (σ₂ - σ₃))/2
    # σ₃ = (Σ - σ₁ - (σ₂ - σ₃))/2

    Λz = dm.Λ * z
    σ1 = (2 * real(Λz) + dm.Σ) / 3
    σ2_minus_σ3 = 2 * imag(Λz) / sqrt(3)
    σ2 = (dm.Σ - σ1 + σ2_minus_σ3) / 2
    σ3 = (dm.Σ - σ1 - σ2_minus_σ3) / 2

    return MandelstamTuple{T}((σ1, σ2, σ3))
end
