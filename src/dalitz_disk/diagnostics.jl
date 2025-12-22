# Diagnostics functions for Dalitz mapping
# Note: DalitzMap is defined in api.jl, which includes this file after DalitzMap is defined

"""
    boundary_error(map::DalitzMap; n::Int=100)

Compute maximum error in boundary mapping: how far from unit circle are boundary points?

# Returns
- `Float64`: Maximum |w| - 1| for boundary points
"""
function boundary_error(map::DalitzMap{T}; n::Int = 100) where {T}
    θs = range(0, 2π, length = n+1)[1:(end-1)]
    boundary_σs = [map.boundary(θ) for θ in θs]
    boundary_w = [map(σs) for σs in boundary_σs]
    radii = abs.(boundary_w)
    return maximum(abs.(radii .- 1.0))
end

"""
    conformality_error(map::DalitzMap; n::Int=20)

Compute Cauchy-Riemann residual error on a grid of interior points.

# Returns
- `Float64`: Maximum Cauchy-Riemann residual
"""
function conformality_error(map::DalitzMap{T}; n::Int = 20) where {T}
    ms = map.ms
    max_error = 0.0
    h = 1e-6

    # Sample interior points
    for _ = 1:n
        x = rand(2)
        σs = y2σs(x, ms)
        if isphysical(σs, ms)
            z = map.democratic(σs)
            w = map.zipper(z)

            # Check Cauchy-Riemann equations
            w_dx = map.zipper(z + h)
            w_dy = map.zipper(z + im * h)

            du_dx = real(w_dx - w) / h
            dv_dx = imag(w_dx - w) / h
            du_dy = real(w_dy - w) / h
            dv_dy = imag(w_dy - w) / h

            cr1 = abs(du_dx - dv_dy)
            cr2 = abs(du_dy + dv_dx)
            max_error = max(max_error, cr1, cr2)
        end
    end

    return max_error
end
