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
            ms, intensity = arg1, arg2
        else
            intensity, ms = arg1, arg2
        end
        # Use the existing recipe by returning the intensity and ms
        return intensity, ms
    else
        error("dalitzplot requires exactly 2 arguments: (ms, intensity) or (intensity, ms)")
    end
end


# User plot for dalitzprojection function (this also creates the DalitzProjection struct)
@userplot DalitzProjection

"""
    dalitzprojection(quadgk, ms, intensity;
        k = 1,
        bins = 100,
        xlims = (:auto, :auto)
    )
    # also
    dalitzprojection(quadgk, intensity, ms; k = 1)

A plotting recipe for Dalitz plot projections.

This recipe generates a 1D projection of the Dalitz plot intensity by integrating 
over one of the invariant mass coordinates, visualizing the intensity as a function 
of a single kinematic variable.

# Parameters:
- `integrator`: A numerical integration function (e.g., `quadgk` from QuadGK.jl). 
  The integrator should accept `(function, lower_limit, upper_limit)` and return 
  a tuple where the first element is the integral value.
- `ms::MassTuple`: A tuple representing the masses of the particles involved in the system.
- `intensity::Function`: A real function of the invariants, `(m23², m31², m12²)`, 
  returning a value at a given kinematic point.

# Keyword Arguments:
- `k`: Index of the invariant to project onto (1, 2, or 3). The projection integrates 
  over the other two coordinates. Defaults to 1.
- `bins`: Number of points for the projection axis. Defaults to 100.
- `xlims`: Limits for the projection axis in terms of the invariant range. 
  Defaults to `lims(ms; k)` (calculated automatically). Can be a tuple with `:auto` 
  for automatic limits (e.g., `(:auto, 4.4)`).

# Output:
The recipe generates:
1. A range of invariant values for the projection axis.
2. A 1D array of integrated intensity values.

# Usage:
```julia
using QuadGK
dalitzprojection(quadgk, ms, intensity; k = 1)
dalitzprojection(quadgk, intensity, ms; k = 3, bins = 150)
```

# Example:
```julia
using ThreeBodyDecays
using Plots
using QuadGK

ms = ThreeBodyMasses(0.938, 0.493, 0.0; m0 = 5.62)
intensity = σs -> 1.0  # constant intensity

# Project onto σ₁
dalitzprojection(quadgk, ms, intensity; k = 1, bins = 100)

# Project onto σ₃ with custom limits
dalitzprojection(quadgk, ms, intensity; k = 3, xlims = (20.0, 26.0))
```
"""
dalitzprojection

@recipe function f(
    dproj::DalitzProjection,
    args...;
    k::Int = 1,
    bins::Int = 100,
    xlims = (:auto, :auto),
)
    # Extract arguments - can be (integrator, ms, intensity) or (integrator, intensity, ms)
    if length(dproj.args) == 3
        integrator = dproj.args[1]
        arg2, arg3 = dproj.args[2], dproj.args[3]
        
        if arg2 isa MassTuple
            ms, intensity = arg2, arg3
        else
            intensity, ms = arg2, arg3
        end
    else
        error("dalitzprojection requires exactly 3 arguments: (integrator, ms, intensity) or (integrator, intensity, ms)")
    end
    
    # Process xlims to handle :auto
    default_xlims = lims(ms; k = k)
    processed_xlims = process_lims(xlims, default_xlims)
    
    # Generate the range of σk values
    σk_range = range(processed_xlims..., length = bins)
    
    # Compute the projection for each σk
    projection_values = map(σk_range) do σk
        integrand = projection_integrand(intensity, ms, σk; k = k)
        result = integrator(integrand, 0, 1)
        result[1]  # Extract the integral value from the result tuple
    end
    
    # Set up the plot
    seriestype := :path
    xguide --> "σ$k"
    yguide --> "Integrated Intensity"
    
    σk_range, projection_values
end
