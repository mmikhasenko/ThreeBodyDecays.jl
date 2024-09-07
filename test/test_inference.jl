using ThreeBodyDecays
using Test
using ThreeBodyDecays.Parameters

@testset "Inferred type" begin
    tbs = ThreeBodySystem(
        2.0,
        1.0,
        1.5;
        m0 = 6.0,
        two_js = ThreeBodySpins(1, 0, 0; two_h0 = 1),
    )  # 1/2+ 0- 0- 1/2+

    dc = DecayChainLS(;
        k = 3,
        Xlineshape = σ -> 1 / (4.1^2 - σ - 0.1im),
        jp = jp"3-",
        Ps = ThreeBodyParities('+', '-', '-'; P0 = '+'),
        tbs = tbs,
    )

    @unpack σs, two_λs = randomPoint(tbs)

    @inferred Complex{Float64} amplitude(dc, σs, two_λs)

    # using InteractiveUtils
    # @code_warntype amplitude(dc, σs, two_λs)
end
