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
        jp = jp"1/2-",
        Ps = ThreeBodyParities('+', '-', '-'; P0 = '+'),
        tbs = tbs,
    )

    @unpack σs, two_λs = randomPoint(tbs)

    @inferred Complex{Float64} amplitude(dc, σs, two_λs)

    # heterogeneous chains should stay concretely typed in the model tuple
    ch_bw = DecayChain(;
        k = 2,
        two_j = 2,
        Xlineshape = σ -> 1 / (5.0^2 - σ - 0.5im),
        Hij = RecouplingLS((2, 0)) |> VertexFunction,
        HRk = RecouplingLS((2, 2)) |> VertexFunction,
        tbs,
    )
    model = ThreeBodyDecay(
        (dc, ch_bw, DecayChain(dc; k = 3)),
        (1.0, -0.5, 0.2im),
        ("a", "b", "c"),
    )
    @test typeof(model.chains) <: Tuple
    @test length(model.chains) == 3
    @inferred Complex{Float64} amplitude(model, σs, two_λs)
    @inferred Float64 unpolarized_intensity(model, σs)

    # using InteractiveUtils
    # @code_warntype amplitude(dc, σs, two_λs)
end
