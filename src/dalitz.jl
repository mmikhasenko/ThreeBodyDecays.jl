const DEFAULT_BORDER_POINTS = 300

"""
    polardalitz2invariants(θ, expansion_point)

For given polar angle θ, returns an (σ1,σ2,σ3) Tuple of polynomials of radius r(θ) around the expansion point.
The polynomial works as a function of the r coordinate.

# Arguments
- `θ`: Polar angle
- `expansion_point`: Tuple of expansion point coordinates

# Returns
- `Tuple{Polynomial,Polynomial,Polynomial}`: Polynomials for each invariant
"""
polardalitz2invariants(θ, expansion_point::Tuple) =
    (
        Polynomial([0, -cos(θ)]),
        Polynomial([0, cos(θ + π / 3)]),
        Polynomial([0, cos(θ - π / 3)]),
    ) .+ expansion_point

function invariants(σx, σy; iσx, iσy, ms)
    iσx == 1 && iσy == 2 && return Invariants(ms; σ1 = σx, σ2 = σy)
    iσx == 2 && iσy == 1 && return Invariants(ms; σ2 = σx, σ1 = σy)
    iσx == 1 && iσy == 3 && return Invariants(ms; σ1 = σx, σ3 = σy)
    iσx == 3 && iσy == 1 && return Invariants(ms; σ3 = σx, σ1 = σy)
    iσx == 3 && iσy == 2 && return Invariants(ms; σ3 = σx, σ2 = σy)
    # remaining case σx = 2 && σy = 3
    return Invariants(ms; σ2 = σx, σ3 = σy)
end

"""
    border(ms::MassTuple{T}; Nx::Int=300) where T

Calculate the border of the Dalitz plot.

# Arguments
- `ms::MassTuple{T}`: Masses of the system
- `Nx::Int`: Number of points to generate

# Returns
- `Vector{MandelstamTuple{T}}`: Points on the border
"""
function border(ms::MassTuple{T}; Nx::Int = DEFAULT_BORDER_POINTS) where {T}
    # Calculate a physically valid expansion point
    expansion_point = let
        f = 0.5
        z = 0.0
        σ1 = (ms[2] + ms[3])^2 + f * ((ms[4] - ms[1])^2 - (ms[2] + ms[3])^2)
        σ3 = σ3of1(z, σ1, ms^2)
        Invariants(ms; σ1, σ3)
    end

    σs(θ) = polardalitz2invariants(θ, expansion_point |> Tuple)
    ϕ0 = Kibble(expansion_point, ms^2)
    ϕ(σs) = Kibble(σs, ms^2)

    function rborder(θ)
        _roots = PolynomialRoots.roots(coeffs(ϕ(σs(θ))))
        filter(_roots) do r
            (abs(imag(r)) < 1e-10) && real(r) > 0.0
        end |> real |> minimum
    end

    function σs_border(θ)
        r = rborder(θ)
        return map(P -> P(r), σs(θ))
    end

    θs = range(-π / 9, 2π - π / 9, length = Nx)
    σs_tuple = σs_border.(θs)
    return MandelstamTuple{T}.(σs_tuple)
end

# border13, border12, border21, border23, border32
for (i, j) in ((1, 2), (2, 1), (2, 3), (3, 2), (3, 1), (1, 3))
    eval(
        quote
            $(Symbol(:border, i, j))(ms; Nx::Int = DEFAULT_BORDER_POINTS) =
                NamedTuple{$(Symbol(:σ, i), Symbol(:σ, j))}.(border(ms; Nx))
        end,
    )
end

@recipe function f(ms::MassTuple, intensity::Function)
    return (intensity, ms)
end

@recipe function f(
    intensity::Function,
    ms::MassTuple;
    iσx = 1,
    iσy = 2,
    grid_size = 100,
    xlims = lims(ms; k = iσx),
    ylims = lims(ms; k = iσy),
)
    #
    σxv = range(xlims..., length = grid_size + 1) |> shift_by_half
    σyv = range(ylims..., length = grid_size + 1) |> shift_by_half
    #
    values = [
        (#
            _σs = invariants(σx, σy; iσx = iσx, iσy = iσy, ms = ms);
            Kibble(_σs, ms^2) > 0 ? NaN : intensity(_σs)
        ) for σy in σyv, σx in σxv
    ]

    seriestype := :heatmap
    colorbar --> false
    #
    σxv, σyv, values
end


"""
    DalitzPlot # only for documentation, see example below
    plot(intensity, ms)
    plot(ms, intensity)

A plotting recipe for `DalitzPlot`.

This recipe generates a Dalitz plot as a heatmap, visualizing the intensity of a function over
a specified range of invariants. The plot provides insights into the kinematic regions of a three-body decay or similar processes.

# Parameters:
- `intensity::Function`: A real function of the invariants, `(m23², m31², m12²)`, returning the intensity at a given kinematic point.
- `ms::MassTuple`: A tuple representing the masses of the particles involved in the system.

# Keyword Arguments:
- `iσx`: Index of the first invariant to use for the x-axis, `1->m23²`, `2->m31²`, and `3->m12²`. Defaults to 1.
- `iσy`: Index of the second invariant to use for the y-axis. Defaults to 2.
- `grid_size`: The resolution of the plot grid. Higher values result in finer detail. Defaults to 100.
- `xlims`: Limits for the x-axis in terms of the invariant range. Defaults to `lims(iσx, ms)` (calculated automatically).
- `ylims`: Limits for the y-axis in terms of the invariant range. Defaults to `lims(iσy, ms)` (calculated automatically).

# Output:
The recipe generates:
1. A grid of invariant values for `σx` and `σy` axes.
2. A 2D array of `values`, where each element is either:
   - `NaN` if the corresponding kinematic region is forbidden by the Kibble condition.
   - The intensity calculated from the provided `intensity` function.

# Usage:
Just call a plot command,
```julia
plot(intensity, ms)
plot(ms, intensity)
```
"""
struct DalitzPlot end
