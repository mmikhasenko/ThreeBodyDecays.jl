# Zipper conformal mapping algorithm

using ConformalMaps

"""
    ZipperMap{T}

A callable struct for the zipper conformal mapping algorithm using ConformalMaps.jl.

# Fields
- `boundary_z::Vector{Complex{T}}`: Boundary in democratic coordinates
- `bf::BoundaryFunction{T}`: Boundary function (for interior point mapping)
- `dm::DemocraticMap{T}`: Democratic map (for coordinate conversion)
- `conformal_map::ConformalMap`: The conformal map from polygon to unit disk
"""
struct ZipperMap{T}
    boundary_z::Vector{Complex{T}}  # Boundary in democratic coordinates
    bf::BoundaryFunction{T}  # Boundary function (for interior point mapping)
    dm::DemocraticMap{T}  # Democratic map (for coordinate conversion)
    conformal_map::ConformalMaps.ConformalMap  # Conformal map from polygon to unit disk
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
function ZipperMap(
    bf::BoundaryFunction{T},
    dm::DemocraticMap{T};
    n::Int = 200,
    point_spacing = nothing,
) where {T}
    # Discretize boundary in democratic coordinates
    θs = range(0, 2π, length = n + 1)[1:(end-1)]  # Exclude duplicate at 2π
    boundary_σs = [bf(θ) for θ in θs]
    boundary_z = [dm(σs) for σs in boundary_σs]

    # Create conformal map using ConformalMaps.jl
    # The center in democratic coordinates is always at the origin (0.0 + 0.0im)
    center = 0.0 + 0.0im

    # Convert boundary_z to the format ConformalMaps expects (Vector{Complex})
    # It's already in the right format, but ensure it's Complex{T}
    boundary_complex = Complex{T}.(boundary_z)

    # Create the conformal map
    # point_spacing controls accuracy: smaller = more accurate but slower
    if point_spacing === nothing
        # Default: 1% of domain diameter
        domain_diameter = maximum(abs.(boundary_complex)) - minimum(abs.(boundary_complex))
        point_spacing = 0.01 * domain_diameter
    end

    conformal_map =
        ConformalMaps.ConformalMap(boundary_complex, center; point_spacing = point_spacing)

    return ZipperMap(boundary_z, bf, dm, conformal_map)
end

"""
    (zm::ZipperMap)(z::Complex{T}) -> Complex{T}

Apply zipper conformal mapping: z → w (unit disk) using ConformalMaps.jl.

# Arguments
- `z::Complex{T}`: Point in democratic coordinate plane

# Returns
- `Complex{T}`: Mapped point in unit disk (|w| < 1 for interior, |w| ≈ 1 for boundary)

# Algorithm
Uses the zipper algorithm from ConformalMaps.jl to perform a proper conformal mapping
from the polygon (boundary in democratic coordinates) to the unit disk.
"""
function (zm::ZipperMap{T})(z::Complex{T}) where {T}
    # Use ConformalMaps.jl to perform the mapping
    # Convert z to the appropriate type if needed
    z_mapped = zm.conformal_map(z)
    return z_mapped
end
