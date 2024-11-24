function invariants(σx, σy; iσx, iσy, ms)
    iσx == 1 && iσy == 2 && return Invariants(ms, σ1 = σx, σ2 = σy)
    iσx == 2 && iσy == 1 && return Invariants(ms, σ1 = σy, σ2 = σx)
    iσx == 1 && iσy == 3 && return Invariants(ms, σ1 = σx, σ3 = σy)
    iσx == 3 && iσy == 1 && return Invariants(ms, σ1 = σy, σ3 = σx)
    iσx == 3 && iσy == 2 && return Invariants(ms, σ2 = σy, σ3 = σx)
    return Invariants(ms, σ2 = σx, σ3 = σy) # σx = 2 && σy = 3
    # error()
end

@recipe function f(ms::MassTuple, intensity::Function)
    return (intensity, ms)
end

@recipe function f(intensity::Function, ms::MassTuple;
    iσx = 1, iσy = 2, grid_size = 100,
    xlims = lims(iσx, ms),
    ylims = lims(iσy, ms))
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
````
"""
struct DalitzPlot end
