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

# Helper function to process xlims/ylims which can be tuples with :auto
function process_lims(lim_tuple, default_lims)
    if lim_tuple isa Tuple && length(lim_tuple) == 2
        lower = lim_tuple[1] === :auto ? default_lims[1] : lim_tuple[1]
        upper = lim_tuple[2] === :auto ? default_lims[2] : lim_tuple[2]
        return (lower, upper)
    else
        return lim_tuple
    end
end

@recipe function f(
    intensity::Function,
    ms::MassTuple,
    args...;
    xbins = 100,
    ybins = 100,
    iσx = 1,
    iσy = 2,
    xlims = (:auto, :auto),
    ylims = (:auto, :auto),
)
    # Process xlims and ylims to handle :auto
    default_xlims = lims(ms; k = iσx)
    default_ylims = lims(ms; k = iσy)
    processed_xlims = process_lims(xlims, default_xlims)
    processed_ylims = process_lims(ylims, default_ylims)
    # +1 and shift by half gets the right picture:
    #   values are computed in the middle of the bin
    #   the bin edge extends exactly to the specified limits
    σxv = range(processed_xlims..., length = xbins+1) |> shift_by_half
    σyv = range(processed_ylims..., length = ybins+1) |> shift_by_half
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


# User plot for dalitzplot function (this also creates the DalitzPlot struct)
@userplot DalitzPlot

"""
    plot(intensity, ms;
        xbins = 100,
        ybins = 100,
        iσx = 1, iσy = 2,
        xlims = (:auto, :auto),
        ylims = (:auto, :auto)
    )
    # also
    plot(ms, intensity)
    dalitzplot(ms, intensity)
    dalitzplot(intensity, ms)

A plotting recipe for a Dalitz plot.

This recipe generates a Dalitz plot as a heatmap, visualizing the intensity of a function over a specified range of invariants.

# Parameters:
- `intensity::Function`: A real function of the invariants, `(m23², m31², m12²)`, returning a value at a given kinematic point.
- `ms::MassTuple`: A tuple representing the masses of the particles involved in the system, used for determination of the borders and function dispatch.

# Keyword Arguments:
- `iσx`: Index of the first invariant to use for the x-axis, `1->m23²`, `2->m31²`, and `3->m12²`. Defaults to 1.
- `iσy`: Index of the second invariant to use for the y-axis. Defaults to 2.
- `xbins`: Number of points for the x-axis grid.
- `ybins`: Number of points for the y-axis grid.
- `xlims`: Limits for the x-axis in terms of the invariant range. Defaults to `lims(iσx, ms)` (calculated automatically). Can be a tuple with `:auto` for automatic limits (e.g., `(:auto, 4.4)`).
- `ylims`: Limits for the y-axis in terms of the invariant range. Defaults to `lims(iσy, ms)` (calculated automatically). Can be a tuple with `:auto` for automatic limits (e.g., `(2.2, :auto)`).

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
dalitzplot(ms, intensity)
dalitzplot(intensity, ms)
```
"""
dalitzplot


@recipe function f(dp::DalitzPlot, args...)
    # Extract arguments - can be (ms, intensity) or (intensity, ms)
    if length(dp.args) == 2
        arg1, arg2 = dp.args
        if arg1 isa MassTuple
            ms, function_of_σs = arg1, arg2
        else
            function_of_σs, ms = arg1, arg2
        end
        # Use the existing recipe by returning the function_of_σs and ms
        return function_of_σs, ms
    else
        error(
            "dalitzplot requires exactly 2 arguments: (ms, function_of_σs) or (function_of_σs, ms)",
        )
    end
end


# User plot for dalitzprojection function (this also creates the DalitzProjection struct)
@userplot DalitzProjection

"""
    dalitzprojection(function_of_σs, ms, integrator;
        k, # must be specified
        bins = 100,
        xlims = (:auto, :auto)
    )
    # also
    dalitzprojection(ms, function_of_σs, integrator; k = 1)

A plotting recipe for Dalitz plot projections.

This recipe generates a 1D projection of the Dalitz plot function_of_σs by integrating
over one of the invariant mass coordinates, visualizing the function_of_σs as a function
of a single kinematic variable.

# Parameters:
- `ms::MassTuple`: A tuple representing the masses of the particles involved in the system.
- `function_of_σs::Function`: A real function of the invariants, `(m23², m31², m12²)`,
  returning a value at a given kinematic point.
- `integrator`: A numerical integration function (e.g., `quadgk` from QuadGK.jl).
  The integrator should accept `(function, lower_limit, upper_limit)` and return
  a tuple where the first element is the integral value.

# Keyword Arguments:
- `k`: Index of the invariant to project onto (1, 2, or 3). The projection integrates
  over the other two coordinates.
- `bins`: Number of bins for the projection axis. Defaults to 100.
- `xlims`: Minimal and maximal values of the bin edges for the projection axis in terms of the invariant range.
  Defaults to `lims(ms; k)` (calculated automatically). Can be a tuple with `:auto`
  for automatic limits (e.g., `(:auto, 4.4)`).

# Usage:
```julia
using QuadGK
dalitzprojection(ms, function_of_σs, quadgk; k = 1)
dalitzprojection(function_of_σs, ms, quadgk; k = 3, bins = 150)
```

# Example:
```julia
using ThreeBodyDecays
using Plots
using QuadGK

ms = ThreeBodyMasses(0.938, 0.493, 0.0; m0 = 5.62)
function_of_σs = σs -> 1.0  # constant function_of_σs

# Project onto σ₁
dalitzprojection(ms, function_of_σs, quadgk; k = 1, bins = 100)

# Project onto σ₃ with custom limits
dalitzprojection(ms, function_of_σs, quadgk; k = 3, xlims = (20.0, 26.0))
```
"""
dalitzprojection

@recipe function f(
    dalitz_projection::DalitzProjection,
    args...;
    k::Int,
    bins::Int = 100,
    xlims = (:auto, :auto),
)
    # Extract arguments - can be (ms, function_of_σs, integrator) or (function_of_σs, ms, integrator)
    if length(dalitz_projection.args) == 3
        integrator = dalitz_projection.args[3]  # integrator is now the last argument
        arg1, arg2 = dalitz_projection.args[1], dalitz_projection.args[2]

        if arg1 isa MassTuple
            ms, function_of_σs = arg1, arg2
        else
            function_of_σs, ms = arg1, arg2
        end
    else
        error(
            "dalitzprojection requires exactly 3 arguments: (ms, function_of_σs, integrator) or (function_of_σs, ms, integrator)",
        )
    end

    # Process xlims to handle :auto
    default_xlims = lims(ms; k)
    processed_xlims = process_lims(xlims, default_xlims)

    # Generate the range of σk values
    σk_range = range(processed_xlims..., length = bins+1) |> shift_by_half

    # Compute the projection for each σk
    projection_values = map(σk_range) do σk
        integrand = projection_integrand(function_of_σs, ms, σk; k)
        result = integrator(integrand, 0, 1)
        result[1]  # Extract the integral value from the result tuple
    end

    # Set up the plot
    seriestype := :path
    xguide --> "σ$k"

    σk_range, projection_values
end
