using Test
using ThreeBodyDecays.Lineshapes
# using Plots

# theme(:wong2, frame=:box, grid=false, minorticks=true,
#     guidefontvalign=:top, guidefonthalign=:right,
#     xlim=(:auto, :auto), ylim=(:auto, :auto),
#     lw=1.2, lab="", colorbar=false)

Γ0 = 0.2
bw = BreitWigner(1.6, Γ0)
@testset "BW call" begin
    @test bw(bw.m^2) ≈ 1im / bw.m / Γ0
    refA = -1.2039660056657226 + 0.5665722379603404im
    @test bw((bw.m + Γ0)^2) ≈ refA
end

bw0 = BlattWeisskopf{0}(1.5)
bw1 = BlattWeisskopf{1}(1.5)
bw2 = BlattWeisskopf{2}(1.5)
bw3 = BlattWeisskopf{3}(1.5)

@testset "BlattWeisskopf" begin
    @test bw0(0.0) == 1.0
    @test bw1(0.0) == bw2(0.0) == bw3(0.0) == 0.0
    @test bw0(1.1) > bw1(1.1) > bw2(1.1) > bw3(1.1)
    refs = (1, 0.9761870601839528, 0.924462392487166, 0.8353277954487898)
    @test all((bw0(3), bw1(3), bw2(3), bw3(3)) .≈ refs)
end

# let
#     plot()
#     map((1, 2, 3, 4, 5)) do l
#         ff = BlattWeisskopf{l}(1.5)
#         plot!(x -> ff(x), 0, 5.0)
#     end
#     plot!()
# end

bw1_of_sq = bw1(z -> z^2)
bw0_scale = bw0 * 5.5
X = bw * bw3
kill_dep = bw(bw0)

@testset "Operations" begin
    @test bw1(2.5^2) == bw1_of_sq(2.5)
    @test bw0_scale(2.5^2) == 5.5
    ref = 0.002582023399904232 + 0.00036559623361475844im
    @test X(0.3) ≈ ref
    @test kill_dep(4.4) == bw(1.0)
end

# let m = 0.77, Γ = 0.15
#     ρ(e) = sqrt(e - 0.28)
#     plot()
#     bw = BreitWigner(m, Γ)
#     plot!(0.28, 0.7, lab="const width") do e
#         e * angle(bw(e^2))
#     end
#     bw = BreitWigner(m, Γ, 0.14, 0.14, 1, 5.0)
#     plot!(0.28, 0.7, lab="P-wave (5/Gev)") do e
#         e * angle(bw(e^2))
#     end
#     bw = BreitWigner(m, Γ, 0.14, 0.14, 1, 1.5)
#     plot!(0.28, 0.7, lab="P-wave (1.5/Gev)") do e
#         e * angle(bw(e^2))
#     end
#     vline!([m], lab="")
# end
