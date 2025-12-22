# Boundary parameterization and validation for Dalitz region

"""
    BoundaryFunction{T}

A callable struct representing the parametric boundary of the Dalitz region.

The boundary is parameterized by angle θ ∈ [0, 2π) and uses polynomial root finding
to determine the boundary point at each angle.

# Fields
- `ms::MassTuple{T}`: Masses of the three-body system
- `expansion_point::MandelstamTuple{T}`: Expansion point for polynomial parameterization
- `msq::NTuple{4,T}`: Cached squared masses for efficiency
"""
struct BoundaryFunction{T}
    ms::MassTuple{T}
    expansion_point::MandelstamTuple{T}
    msq::NTuple{4,T}
end

"""
    BoundaryFunction(ms::MassTuple)

Construct a `BoundaryFunction` for the given mass configuration.

The boundary is parameterized using polynomial root finding, following the
approach in `border()` from `src/dalitz.jl`.

# Formula
The boundary at angle θ is found by solving φ(σ₁(θ,r), σ₂(θ,r), σ₃(θ,r)) = 0
for the smallest positive real root r(θ), where φ is the Kibble function.

# Arguments
- `ms::MassTuple`: Masses of the three-body system

# Returns
- `BoundaryFunction{T}`: Callable struct for boundary evaluation

# Example
```julia
ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0=4.0)
bf = BoundaryFunction(ms)
σs = bf(π/4)  # Get boundary point at angle π/4
```
"""
function BoundaryFunction(ms::MassTuple{T}) where {T}
    # Calculate expansion point (same as in border())
    expansion_point = let
        f = 0.5
        z = 0.0
        σ1 = (ms[2] + ms[3])^2 + f * ((ms[4] - ms[1])^2 - (ms[2] + ms[3])^2)
        σ3 = σ3of1(z, σ1, ms^2)
        Invariants(ms; σ1, σ3)
    end

    return BoundaryFunction(ms, expansion_point, ms^2)
end

"""
    (bf::BoundaryFunction)(θ::Real) -> MandelstamTuple

Evaluate the boundary at angle θ.

# Arguments
- `θ::Real`: Polar angle in [0, 2π)

# Returns
- `MandelstamTuple{T}`: Boundary point in Mandelstam invariants

# Formula
The boundary is found by solving φ(σ₁(θ,r), σ₂(θ,r), σ₃(θ,r)) = 0 for the smallest
positive real root r(θ), where:
- σᵢ(θ,r) = σᵢ⁽⁰⁾ + r·pᵢ(θ)
- p₁(θ) = -cos(θ), p₂(θ) = cos(θ + π/3), p₃(θ) = cos(θ - π/3)
- φ is the Kibble function
"""
function (bf::BoundaryFunction{T})(θ::Real) where {T}
    # Get direction polynomials
    σs_poly(θ) = polardalitz2invariants(θ, bf.expansion_point |> Tuple)

    # Find boundary radius by solving Kibble = 0
    function rborder(θ)
        σs_θ = σs_poly(θ)
        ϕ_poly = Kibble(σs_θ, bf.msq)
        _roots = PolynomialRoots.roots(coeffs(ϕ_poly))
        # Filter for real, positive roots and convert to real
        real_roots = Float64[]
        for r in _roots
            if abs(imag(r)) < 1e-10 && real(r) > 0.0
                push!(real_roots, real(r))
            end
        end
        isempty(real_roots) && error("No valid boundary root found at θ=$θ")
        return minimum(real_roots)
    end

    r = rborder(θ)
    σs_θ = σs_poly(θ)
    # Evaluate polynomials at r
    σs_values = map(P -> P(r), σs_θ)
    return MandelstamTuple{T}(σs_values)
end

"""
    boundary_closure_error(bf::BoundaryFunction; tolerance=1e-10)

Check boundary closure: |boundary(0) - boundary(2π)| < tolerance.

# Returns
- `Float64`: Maximum absolute difference between boundary(0) and boundary(2π)
"""
function boundary_closure_error(bf::BoundaryFunction; tolerance = 1e-10)
    σs_0 = bf(0.0)
    σs_2π = bf(2π)
    max_diff = maximum(abs.(Tuple(σs_0) .- Tuple(σs_2π)))
    return max_diff
end

"""
    boundary_orientation(bf::BoundaryFunction, n::Int=100)

Check boundary orientation: compute signed area to verify CCW orientation.

# Arguments
- `bf::BoundaryFunction`: Boundary function
- `n::Int`: Number of sample points

# Returns
- `Float64`: Signed area (should be > 0 for CCW)
"""
function boundary_orientation(bf::BoundaryFunction, n::Int = 100)
    θs = range(0, 2π, length = n+1)[1:(end-1)]  # Exclude duplicate at 2π
    points = [bf(θ) for θ in θs]

    # Compute signed area using shoelace formula in (σ₁, σ₂) plane
    area = 0.0
    for i = 1:length(points)
        j = mod(i, length(points)) + 1
        area += points[i].σ1 * points[j].σ2 - points[j].σ1 * points[i].σ2
    end
    return area / 2.0
end

"""
    inside_D(σs::MandelstamTuple, ms::MassTuple) -> Bool

Check if a point is inside the Dalitz region using the Kibble function.

# Formula
A point is inside if φ(σ₁, σ₂, σ₃) < 0 (Kibble condition) and all
invariants are within their physical ranges.

# Arguments
- `σs::MandelstamTuple`: Point to check
- `ms::MassTuple`: Masses of the system

# Returns
- `Bool`: `true` if inside, `false` otherwise
"""
inside_D(σs::MandelstamTuple, ms::MassTuple) = isphysical(σs, ms)

"""
    boundary_self_intersections(bf::BoundaryFunction, n::Int=300)

Check for self-intersections in the boundary polygon.

# Arguments
- `bf::BoundaryFunction`: Boundary function
- `n::Int`: Number of sample points

# Returns
- `Bool`: `true` if no intersections found, `false` if intersections detected
"""
function boundary_self_intersections(bf::BoundaryFunction, n::Int = 300)
    θs = range(0, 2π, length = n+1)[1:(end-1)]
    points = [bf(θ) for θ in θs]

    # Check all pairs of non-adjacent segments for intersection
    for i = 1:length(points)
        p1 = (points[i].σ1, points[i].σ2)
        p2 = (points[mod(i, length(points))+1].σ1, points[mod(i, length(points))+1].σ2)

        for j = (i+2):length(points)
            if j == length(points) && i == 1
                continue  # Skip last-first segment pair
            end
            p3 = (points[j].σ1, points[j].σ2)
            p4 = (points[mod(j, length(points))+1].σ1, points[mod(j, length(points))+1].σ2)

            if segments_intersect(p1, p2, p3, p4)
                return false
            end
        end
    end
    return true
end

"""
    segments_intersect(p1, p2, p3, p4) -> Bool

Check if two line segments intersect using cross product method.

# Arguments
- `p1, p2`: Endpoints of first segment
- `p3, p4`: Endpoints of second segment

# Returns
- `Bool`: `true` if segments intersect
"""
function segments_intersect(p1, p2, p3, p4)
    # Helper function for orientation
    orientation(p, q, r) = (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])

    o1 = orientation(p1, p2, p3)
    o2 = orientation(p1, p2, p4)
    o3 = orientation(p3, p4, p1)
    o4 = orientation(p3, p4, p2)

    # General case: segments intersect if orientations differ
    if o1 * o2 < 0 && o3 * o4 < 0
        return true
    end

    # Special cases: collinear points
    if o1 == 0 && on_segment(p1, p3, p2)
        return true
    end
    if o2 == 0 && on_segment(p1, p4, p2)
        return true
    end
    if o3 == 0 && on_segment(p3, p1, p4)
        return true
    end
    if o4 == 0 && on_segment(p3, p2, p4)
        return true
    end

    return false
end

"""
    on_segment(p, q, r) -> Bool

Check if point q lies on segment pr.
"""
function on_segment(p, q, r)
    return (
        q[1] <= max(p[1], r[1]) &&
        q[1] >= min(p[1], r[1]) &&
        q[2] <= max(p[2], r[2]) &&
        q[2] >= min(p[2], r[2])
    )
end

"""
    boundary_cosine_scan(ms::MassTuple, k::Int=1; n::Int=300)

Alternative boundary construction using cosine scan in sequential frame.

This provides cross-validation for the polynomial-based boundary.
Uses cosθ = ±1 to get the two branches of the boundary.

# Arguments
- `ms::MassTuple`: Masses of the system
- `k::Int`: Spectator index (default 1)
- `n::Int`: Number of points

# Returns
- `Vector{MandelstamTuple}`: Boundary points
"""
function boundary_cosine_scan(ms::MassTuple{T}, k::Int = 1; n::Int = 300) where {T}
    i, j = ij_from_k(k)
    σk_min, σk_max = lims(ms; k = k)
    σk_range = range(σk_min, σk_max, length = n)

    boundary_points = MandelstamTuple{T}[]
    for σk in σk_range
        # Get limits for σj given σk using cosθ = ±1
        σj_at_cos1 = σjofk(1.0, σk, ms^2; k = k)   # cosθ = 1
        σj_at_cosm1 = σjofk(-1.0, σk, ms^2; k = k) # cosθ = -1

        # Compute corresponding σi values
        σi_at_cos1 = sum(ms^2) - σk - σj_at_cos1
        σi_at_cosm1 = sum(ms^2) - σk - σj_at_cosm1

        # Create invariants tuples with proper ordering
        σt1 = circleorigin(-k, (σi_at_cos1, σj_at_cos1, σk))
        σt2 = circleorigin(-k, (σi_at_cosm1, σj_at_cosm1, σk))

        push!(boundary_points, MandelstamTuple{T}(σt1))
        push!(boundary_points, MandelstamTuple{T}(σt2))
    end

    # Filter to physical points only
    physical_points = filter(σs -> isphysical(σs, ms), boundary_points)
    return physical_points
end
