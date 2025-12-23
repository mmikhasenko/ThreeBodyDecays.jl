using ThreeBodyDecays
using Test
using Random
using Plots
theme(:boxed)


ms = ThreeBodyMasses(0.1, 1.2, 1.3; m0 = 4.0)
my_map = dalitzmap(ms)

σs_rand = x2σs.(rand(2, 3000) |> eachcol, ms |> Ref; k = 3)
wv = my_map.(σs_rand)  # Map Mandelstam invariants to unit disk
let
    plot(border12(ms))
    xv = getproperty.(σs_rand, :σ1)
    yv = getproperty.(σs_rand, :σ2)
    scatter!(xv, yv)
    #
    scale = sum(ms^2) / 10
    shift = lims(ms; k = 1)[end] + 1im * lims(ms; k = 2)[end]
    t(w) = w .* scale .+ shift
    #
    scatter!(t(wv))
    plot!(t(cis.(range(0, 2π, 100))))

    # map(zip(xv .+ 1im .* yv, t.(wv))) do (w1, w2)
    #     plot!([w1, w2], arrow = true, lw = 0.3, alpha = 0.5)
    # end
    # plot!()
end
