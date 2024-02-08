using Test
using ThreeBodyDecays

# test model constraction
@testset "ThreeBodyDecay" begin
    tbs = ThreeBodySystem(
        ms = ThreeBodyMasses(2.0, 1.0, 1.5; m0=6.0),
        two_js=ThreeBodySpins(1, 0, 0; two_h0=1))
    dpp = randomPoint(tbs)
    #
    dc = DecayChainLS(
        k=1,
        two_j=2,
        Xlineshape = σ -> 1 / (4.1^2 - σ - 0.1im),
        # Hij = 
        tbs=tbs)
    @test sum(reim(amplitude(dc, dpp)) .≈ 0.0) == 0
end
