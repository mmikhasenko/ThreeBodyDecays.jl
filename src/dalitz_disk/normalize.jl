# Channel-democratic normalization via Möbius transformation

"""
    find_landmarks(bf::BoundaryFunction, n::Int=500)

Find three boundary landmarks by maximizing each σᵢ along the boundary.

# Arguments
- `bf::BoundaryFunction`: Boundary function
- `n::Int`: Number of sample points for optimization

# Returns
- `Tuple{Real, Real, Real}`: Angles θ₁, θ₂, θ₃ where each σᵢ is maximized
"""
function find_landmarks(bf::BoundaryFunction{T}, n::Int = 500) where {T}
    θs = range(0, 2π, length = n+1)[1:(end-1)]

    # Find maximum of each σᵢ
    σ1_vals = [bf(θ).σ1 for θ in θs]
    σ2_vals = [bf(θ).σ2 for θ in θs]
    σ3_vals = [bf(θ).σ3 for θ in θs]

    idx1 = argmax(σ1_vals)
    idx2 = argmax(σ2_vals)
    idx3 = argmax(σ3_vals)

    θ1 = θs[idx1]
    θ2 = θs[idx2]
    θ3 = θs[idx3]

    return (θ1, θ2, θ3)
end

"""
    mobius_transform_three_points(w1, w2, w3, target1, target2, target3)

Compute Möbius transformation mapping three points to three targets.

# Formula
A(w) = (aw + b) / (cw + d) where A(wᵢ) = targetᵢ

# Arguments
- `w1, w2, w3`: Source points (complex)
- `target1, target2, target3`: Target points (complex)

# Returns
- `MobiusTransform`: Möbius transformation
"""
function mobius_transform_three_points(
    w1::Complex{T},
    w2::Complex{T},
    w3::Complex{T},
    target1::Complex{T},
    target2::Complex{T},
    target3::Complex{T},
) where {T}
    # Solve for Möbius transformation A(w) = (aw + b) / (cw + d)
    # that maps w1→target1, w2→target2, w3→target3

    # With normalization d=1, we have:
    # target_i = (a*w_i + b) / (c*w_i + 1)
    # Rearranging: target_i * (c*w_i + 1) = a*w_i + b
    # This gives: a*w_i + b - c*target_i*w_i = target_i

    # Linear system: [w1  1  -target1*w1] [a]   [target1]
    #                [w2  1  -target2*w2] [b] = [target2]
    #                [w3  1  -target3*w3] [c]   [target3]

    M = [
        w1 1 -target1*w1;
        w2 1 -target2*w2;
        w3 1 -target3*w3
    ]
    rhs = [target1; target2; target3]
    sol = M \ rhs
    a, b, c = sol[1], sol[2], sol[3]
    d = one(Complex{T})

    # MobiusTransform is defined in api.jl, but we need to reference it
    # It should be available when normalize.jl is included after api.jl
    return MobiusTransform(a, b, c, d)
end

"""
    normalize_landmarks(bf::BoundaryFunction, dm::DemocraticMap, zm::ZipperMap)

Compute channel-democratic normalization by finding landmarks and applying Möbius transformation.

# Arguments
- `bf::BoundaryFunction`: Boundary function
- `dm::DemocraticMap`: Democratic coordinate map
- `zm::ZipperMap`: Zipper map

# Returns
- `MobiusTransform`: Möbius transformation mapping landmarks to 1, ω, ω²
"""
function normalize_landmarks(
    bf::BoundaryFunction{T},
    dm::DemocraticMap{T},
    zm::ZipperMap{T},
) where {T}
    # Find landmarks
    θ1, θ2, θ3 = find_landmarks(bf)

    # Get landmark points in Mandelstam coordinates
    σs1 = bf(θ1)
    σs2 = bf(θ2)
    σs3 = bf(θ3)

    # Map to democratic coordinates
    z1 = dm(σs1)
    z2 = dm(σs2)
    z3 = dm(σs3)

    # Map through zipper
    w1 = zm(z1)
    w2 = zm(z2)
    w3 = zm(z3)

    # Target points: 1, ω, ω²
    ω = exp(2π * im / 3)
    target1 = one(Complex{T})
    target2 = ω
    target3 = ω^2

    # Compute Möbius transformation
    return mobius_transform_three_points(w1, w2, w3, target1, target2, target3)
end
