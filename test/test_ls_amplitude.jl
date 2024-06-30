using ThreeBodyDecays
using Test

@testset "LS amplitude scalars" begin
    tbs = ThreeBodySystem(2.0, 1.0, 1.5; m0 = 6.0)
    dpp = randomPoint(tbs)
    #
    dc = DecayChainLS(;
        k = 1,
        Xlineshape = σ -> 1 / (4.1^2 - σ - 0.1im),
        jp = "2+",
        Ps = ThreeBodyParities('+', '+', '+', P0 = '+'),
        tbs = tbs)
    @test sum(reim(amplitude(dc, dpp)) .≈ 0.0) == 0
end

two_js, Ps = ThreeBodySpinParities("1/2-", "0-", "0-"; jp0 = "1/2+")
tbs = ThreeBodySystem(2.0, 1.0, 1.5; m0 = 6.0, two_js)  # 1/2+ 0- 0- 1/2+
σs = x2σs([0.5, 0.3], tbs.ms; k = 1)
#
bw(σ) = 1 / (4.1^2 - σ - 0.1im)
dc = DecayChainLS(; k = 3, Xlineshape = bw, jp = "3-", Ps, tbs)

@testset "LS amplitude half-integer spin" begin
    dpp = DalitzPlotPoint(; σs = σs, two_λs = [1, 0, 0, 1])
    @test sum(reim(amplitude(dc, dpp)) .≈ 0.0) == 0
    @test length(dc) == 1
end


# Manual construction of the decay chains
ms = ThreeBodyMasses(1.1, 2.2, 3.3; m0 = 7.7)
(two_js, PC), (_, PV) = map(["1/2+", "1/2-"]) do jp0
    ThreeBodySpinParities("1-", "1/2+", "0-"; jp0)
end
tbs = ThreeBodySystem(ms, two_js)

chains = let
    k = 3
    jp = jp"3/2-"
    ls_all = possible_ls_ij(jp, two_js, PC; k)
    LS_all = vcat(
        possible_ls_Rk(jp, two_js, PC; k),
        possible_ls_Rk(jp, two_js, PV; k))
    #
    map(Iterators.product(LS_all, ls_all)) do (two_LS, two_ls)
        L = div(two_LS[1], 2)
        l = div(two_ls[1], 2)
        Xlineshape = identity
        DecayChain(; k, jp.two_j, tbs, Xlineshape,
            HRk = RecouplingLS(two_LS), Hij = RecouplingLS(two_ls))
    end
end

sort_by_Ll(chains) = sort(chains,
    by = x -> x.Hij.two_ls[1] * 10 + x.HRk.two_ls[1])

@testset "DecayChainsLS: all decay chains" begin
    @test sort_by_Ll(vec(chains)) ==
          sort_by_Ll(vec(
        [DecayChainsLS(; k = 3, Xlineshape = identity, jp = "3/2-", Ps = PC, tbs);
            DecayChainsLS(; k = 3, Xlineshape = identity, jp = "3/2-", Ps = PV, tbs)]))
end

σs = x2σs([0.5, 0.3], tbs.ms; k = 1)
@testset "Unpolarized intensity values for chains" begin
    refs = [    183.0603173468474 100.20583627496791 183.0603173468474
        183.0603173468474 100.20583627496791 183.0603173468474]
    @test all(unpolarized_intensity.(chains, Ref(σs)) .≈ refs)
end
