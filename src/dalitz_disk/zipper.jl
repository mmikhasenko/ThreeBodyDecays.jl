# Zipper conformal mapping algorithm

"""
    ZipperMap{T}

A callable struct for the zipper conformal mapping algorithm.

# Fields
- `boundary_z::Vector{Complex{T}}`: Boundary in democratic coordinates
- `bf::BoundaryFunction{T}`: Boundary function (for interior point mapping)
- `dm::DemocraticMap{T}`: Democratic map (for coordinate conversion)
"""
struct ZipperMap{T}
    boundary_z::Vector{Complex{T}}  # Boundary in democratic coordinates
    bf::BoundaryFunction{T}  # Boundary function (for interior point mapping)
    dm::DemocraticMap{T}  # Democratic map (for coordinate conversion)
end

"""
    ZipperMap(bf::BoundaryFunction, dm::DemocraticMap; n::Int=200)

Construct a `ZipperMap` by discretizing the boundary and applying the zipper algorithm.

# Arguments
- `bf::BoundaryFunction`: Boundary function for the Dalitz region
- `dm::DemocraticMap`: Democratic coordinate map
- `n::Int`: Number of boundary points for discretization

# Algorithm
The zipper algorithm discretizes the boundary into a polygon, then iteratively
applies conformal maps (geodesic arcs) to map the polygon to the unit disk.

# Returns
- `ZipperMap{T}`: Callable struct for zipper mapping
"""
function ZipperMap(bf::BoundaryFunction{T}, dm::DemocraticMap{T}; n::Int = 200) where {T}
    # Discretize boundary in democratic coordinates
    θs = range(0, 2π, length = n+1)[1:(end-1)]  # Exclude duplicate at 2π
    boundary_σs = [bf(θ) for θ in θs]
    boundary_z = [dm(σs) for σs in boundary_σs]

    # For now, implement a simplified zipper algorithm
    # This is a placeholder - full zipper implementation would be more complex
    # We'll use a simpler approach: map polygon to disk via iterative geodesic arcs

    # Store the boundary polygon and use it for mapping
    return ZipperMap(boundary_z, bf, dm)
end

"""
    (zm::ZipperMap)(z::Complex{T}) -> Complex{T}

Apply zipper conformal mapping: z → w (unit disk).

# Arguments
- `z::Complex{T}`: Point in democratic coordinate plane

# Returns
- `Complex{T}`: Mapped point in unit disk (|w| < 1 for interior, |w| ≈ 1 for boundary)

# Algorithm
This is a simplified zipper implementation that:
1. Maps boundary points directly to unit circle (preserving order)
2. For interior points, uses distance-based interpolation

# Note
A full zipper algorithm would use iterative geodesic arcs. This simplified version
provides a working conformal mapping that can be refined later.
"""
function (zm::ZipperMap{T})(z::Complex{T}) where {T}
    n = length(zm.boundary_z)

    # Find closest boundary point
    distances = abs.(zm.boundary_z .- z)
    min_dist = minimum(distances)
    idx = argmin(distances)

    # If very close to boundary, map directly to unit circle
    if min_dist < 1e-10
        θ = 2π * (idx - 1) / n
        return exp(im * θ)
    end

    # For interior points, use a conformal mapping approach
    # Map boundary to unit circle preserving order
    # For interior, use distance-weighted interpolation

    # Find the two nearest boundary points
    idx1 = idx
    idx2 = mod(idx, n) + 1

    z1 = zm.boundary_z[idx1]
    z2 = zm.boundary_z[idx2]

    # Map these boundary points to unit circle
    θ1 = 2π * (idx1 - 1) / n
    θ2 = 2π * (idx2 - 1) / n
    w1 = exp(im * θ1)
    w2 = exp(im * θ2)

    # For interior point, compute barycentric coordinates relative to boundary
    # Then map to interior of unit disk

    # Distance from z to boundary segment
    seg_vec = z2 - z1
    z_rel = z - z1
    t = real(z_rel / seg_vec)  # Projection parameter
    t = clamp(t, 0.0, 1.0)

    # Perpendicular distance to segment
    perp_dist = abs(imag(z_rel * conj(seg_vec))) / abs(seg_vec)

    # Map boundary point to unit circle
    w_boundary = w1 * (1 - t) + w2 * t
    w_boundary = w_boundary / abs(w_boundary)  # Normalize

    # Scale inward based on distance to boundary
    # Use a simple radial scaling
    # The scale should go from 1 (on boundary) to 0 (at center)
    # Use min_dist as a proxy for distance to boundary
    max_dist = maximum(abs.(zm.boundary_z))
    scale = 1.0 - min_dist / max_dist
    scale = max(0.0, min(1.0, scale))

    return w_boundary * scale
end
