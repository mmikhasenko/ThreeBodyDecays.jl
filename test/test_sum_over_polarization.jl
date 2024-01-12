using Test
using ThreeBodyDecays

@testset "Sum over polarization" begin
    mΛb = 5.61960
    mπ = 0.14
    # intermediate resonances
    mΣb = 5.81065
    ΓΣb = 0.005
    mΣb_x = 5.83032
    ΓΣb_x = 0.0094

    mΛb2S = 6.3 # just a peak of the plot

    tbs = ThreeBodySystem(mπ, mΛb, mπ; m0=mΛb2S,
        two_js=ThreeBodySpins(0, 1, 0; two_h0=1))  # 1/2+ 0- 0- 1/2+
    dpp = randomPoint(tbs)
    σs = dpp.σs

    bw_Σb(σ) = 1 / (mΣb^2 - σ - 1im * mΣb * ΓΣb)
    bwΣb_x(σ) = 1 / (mΣb_x^2 - σ - 1im * mΣb_x * ΓΣb_x)
    # lineshape
    dc_Σb1 = DecayChainLS(1, bw_Σb; tbs=tbs, two_j=1, parity='+', Ps=['-', '+', '-', '+'])
    dc_Σb3 = DecayChainLS(3, bw_Σb; tbs=tbs, two_j=1, parity='+', Ps=['-', '+', '-', '+'])

    dc_Σb_x1 = DecayChainLS(1, bwΣb_x; tbs=tbs, two_j=3, parity='+', Ps=['-', '+', '-', '+'])
    dc_Σb_x3 = DecayChainLS(3, bwΣb_x; tbs=tbs, two_j=3, parity='+', Ps=['-', '+', '-', '+'])

    full_a(σs, two_λs) = sum(amplitude(ch, σs, two_λs) for ch in [dc_Σb1, dc_Σb_x1, dc_Σb3, dc_Σb_x3])
    full_a(dpp) = full_a(dpp.σs, dpp.two_λs)
    #
    total_I(σs) = sum(abs2(full_a(σs, two_λs)) for two_λs in itr(tbs.two_js))

    test_I = summed_over_polarization((σs, two_λs) -> abs2(full_a(σs, two_λs)), tbs.two_js)
    @test test_I(dpp.σs) == total_I(dpp.σs)
end
