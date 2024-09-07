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

@recipe function f(intensity::Function, ms::MassTuple; iσx = 1, iσy = 2, grid_size = 100)
    #
    σxv = range(lims(iσx, ms)..., length = grid_size + 1) |> shift_by_half
    σyv = range(lims(iσy, ms)..., length = grid_size + 1) |> shift_by_half
    #
    #
    values = [
        (
            Kibble(invariants(σx, σy; iσx = iσx, iσy = iσy, ms = ms), ms^2) > 0 ? NaN :
            intensity(invariants(σx, σy; iσx = iσx, iσy = iσy, ms = ms))
        ) for σy in σyv, σx in σxv
    ]

    seriestype := :heatmap
    colorbar --> false
    #
    σxv, σyv, values
end
