using ThreeBodyDecays
using Test

@testset "LS amplitude integer spin" begin
    tbs = ThreeBodySystem(2.0, 1.0, 1.5; m0=6.0)
    dpp = randomPoint(tbs)
    #
    dc = DecayChainLS(;
        k=1,
        Xlineshape=σ -> 1 / (4.1^2 - σ - 0.1im),
        jp="2+",
        Ps=ThreeBodyParities('+', '+', '+', P0='+'),
        tbs=tbs)
    @test sum(reim(amplitude(dc, dpp)) .≈ 0.0) == 0
end

@testset "LS amplitude half-integer spin" begin
    tbs = ThreeBodySystem(2.0, 1.0, 1.5; m0=6.0,
        two_js=ThreeBodySpins(1, 0, 0; two_h0=1))  # 1/2+ 0- 0- 1/2+
    σs = randomPoint(tbs.ms)
    dpp = DalitzPlotPoint(; σs=σs, two_λs=[1, 0, 0, 1])
    #
    dc = DecayChainLS(; k=3,
        Xlineshape=σ -> 1 / (4.1^2 - σ - 0.1im),
        jp="3-", Ps=ThreeBodyParities('+', '-', '-', P0='+'), tbs=tbs)
    # @show amplitude(dpp, dc)
    @test sum(reim(amplitude(dc, dpp)) .≈ 0.0) == 0
    # testing something else?
end
