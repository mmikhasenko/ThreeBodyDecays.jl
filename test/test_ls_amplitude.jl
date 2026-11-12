using ThreeBodyDecays
using Test

using ThreeBodyDecays.Parameters
using BenchmarkTools
using Random

@testset "LS amplitude scalars" begin
    tbs = ThreeBodySystem(2.0, 1.0, 1.5; m0 = 6.0)
    dpp = randomPoint(tbs)
    #
    dc = DecayChainLS(;
        k = 1,
        Xlineshape = σ -> 1 / (4.1^2 - σ - 0.1im),
        jp = "2+",
        Ps = ThreeBodyParities('+', '+', '+', P0 = '+'),
        tbs = tbs,
    )
    @test sum(reim(amplitude(dc, dpp)) .≈ 0.0) == 0
end

two_js, Ps = ThreeBodySpinParities("1/2-", "0-", "0-"; jp0 = "1/2+")
tbs = ThreeBodySystem(2.0, 1.0, 1.5; m0 = 6.0, two_js)  # 1/2+ 0- 0- 1/2+
σs = x2σs([0.5, 0.3], tbs.ms; k = 1)
#
bw(σ) = 1 / (4.1^2 - σ - 0.1im)
dc = DecayChainLS(; k = 3, Xlineshape = bw, jp = "1/2-", Ps, tbs)

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
    LS_all = vcat(possible_ls_Rk(jp, two_js, PC; k), possible_ls_Rk(jp, two_js, PV; k))
    #
    map(Iterators.product(LS_all, ls_all)) do (two_LS, two_ls)
        Xlineshape = identity
        DecayChain(;
            k,
            jp.two_j,
            tbs,
            Xlineshape,
            HRk = RecouplingLS(two_LS) |> Vertex,
            Hij = RecouplingLS(two_ls) |> Vertex,
        )
    end
end

sort_by_Ll(chains) = sort(chains, by = x -> x.Hij.h.two_ls[1] * 10 + x.HRk.h.two_ls[1])

@testset "DecayChainsLS: all decay chains" begin
    @test sort_by_Ll(vec(chains)) == sort_by_Ll(
        vec(
            [
                DecayChainsLS(; k = 3, Xlineshape = identity, jp = "3/2-", Ps = PC, tbs)
                DecayChainsLS(; k = 3, Xlineshape = identity, jp = "3/2-", Ps = PV, tbs)
            ],
        ),
    )
end

σs = x2σs([0.5, 0.3], tbs.ms; k = 1)
@testset "Unpolarized intensity values for chains" begin
    refs = [
        183.0603173468474 100.20583627496791 183.0603173468474
        183.0603173468474 100.20583627496791 183.0603173468474
    ]
    @test all(unpolarized_intensity.(chains, Ref(σs)) .≈ refs)
end



const model = let
    two_js, Ps = ThreeBodySpinParities("1-", "1/2+", "0-"; jp0 = "3/2+")
    tbs = ThreeBodySystem(1.1, 2.2, 3.3; m0 = 7.7, two_js)
    #
    models = map([
        (name = "R3_3h-", k = 3, two_jp = jp"3/2-"),
        (name = "R1_3h-", k = 1, two_jp = jp"3/2-"),
        (name = "R3_1h-", k = 3, two_jp = jp"1/2-"),
        (name = "R1_1h-", k = 1, two_jp = jp"1/2-"),
        (name = "R2", k = 2, two_jp = jp"3-"),
    ]) do (; k, two_jp, name)
        #
        chains = map(possible_lsLS(two_jp, two_js, Ps; k)) do conf
            @unpack two_LS, two_ls = conf
            DecayChain(;
                k,
                two_jp.two_j,
                tbs,
                Xlineshape = identity,
                HRk = RecouplingLS(two_LS) |> Vertex,
                Hij = RecouplingLS(two_ls) |> Vertex,
            )
        end
        ci = ones(Float64, length(chains))
        ThreeBodyDecay(name .=> zip(ci, chains); reference_k = 1)
    end
    vcat(models...)
end

@testset "Elapsed time" begin
    @test length(model) == 20
    σs = x2σs([0.5, 0.3], masses(model); k = 1)
    ref_I = 17646.88022101494
    @test unpolarized_intensity(model, σs) ≈ ref_I
    @test unpolarized_intensity(
        ThreeBodyDecay(model.chains, model.couplings, model.names; reference_k = 1),
        σs,
    ) ≈ ref_I
    # Convert the reference-1 couplings, rather than relabeling the same numbers.
    sign_for(k, r) = k == r ? 1 : (-1)^spins(model)[6-k-r]
    for r = 2:3
        converted = map(model.chains, model.couplings) do chain, c
            sign_for(chain.k, r) * sign_for(chain.k, 1) * c
        end
        referenced = ThreeBodyDecay(model.chains, converted, model.names; reference_k = r)
        @test unpolarized_intensity(referenced, σs) ≈ ref_I
        unchanged =
            ThreeBodyDecay(model.chains, model.couplings, model.names; reference_k = r)
        @test !isapprox(unpolarized_intensity(unchanged, σs), ref_I)
    end
    # @btime unpolarized_intensity($model, $σs)
    eval_time = @elapsed unpolarized_intensity(model, σs)
    @info """Unpolarized_intensity for a model (3/2->1,1/2,0) with 20 chains is computed in $(round(1000*eval_time, digits=2)) ms
Compare to my usual evaluation time of about 5ms
"""
end

@testset "Full amplitude components" begin
    σs = x2σs([0.5, 0.3], masses(model); k = 1)
    # Fixed polarized oracle from the pre-migration API (main),
    # using the old common-reference tuples (k,k,k,k).
    expected = (
        (24.313600920343276, -16.53265693046843),
        (40.37973824558429, -14.059543287122107),
        (27.740984239852448, -4.141879753609057),
    )
    for k = 1:3
        referenced =
            ThreeBodyDecay(model.chains, model.couplings, model.names; reference_k = k)
        tensor = amplitude(referenced, σs)
        @test tensor[end, end, end, end] ≈ expected[k][1]
        @test tensor[1, 1, end, end] ≈ expected[k][2]
    end
    A_2103 = amplitude(model, σs, [2, 1, 0, 3])
    A = amplitude(model, σs)
    @test A_2103 ≈ A[end, end, end, end] ≈ 24.313600920343276
    A_m2m103 = amplitude(model, σs, [-2, -1, 0, 3])
    @test A_m2m103 ≈ A[1, 1, end, end] ≈ -16.53265693046843
end





struct FF
    R::Float64
end
function (ff::FF)(m0sq, m1sq, m2sq)
    m0, m1, m2 = sqrt.((m0sq, m1sq, m2sq))
    p = breakup(m0, m1, m2)
    z = p * ff.R
    z = p / sqrt(1 + z^2)
end

@testset "Vertex constructor" begin
    @test Vertex(RecouplingLS((2, 1))) isa Vertex
    @test Vertex(RecouplingLS((2, 1)), FF(1.5)) isa Vertex
end

@testset "VertexFunction deprecation alias" begin
    @test VertexFunction === Vertex
    @test VertexFunction(RecouplingLS((2, 1))) isa Vertex
    @test VertexFunction(RecouplingLS((2, 1)), FF(1.5)) isa Vertex
end


Random.seed!(1234)
@testset "Vertex representation property" begin
    ms = ThreeBodyMasses(0.938, 0.49367, 0.13957; m0 = 2.46867)
    tbs = ThreeBodySystem(ms, ThreeBodySpins(1, 1, 1; two_h0 = 1))

    # Test that Vertex preserves representation property
    dc = DecayChain(;
        k = 1,
        two_j = 2,
        Xlineshape = σ -> 1.0,
        HRk = Vertex(RecouplingLS((2, 1)), FF(1.5)),
        Hij = Vertex(RecouplingLS((2, 1)), FF(1.5)),
        tbs,
    )

    σs = randomPoint(tbs.ms)
    @test amplitude(dc, σs) ≈ aligned_amplitude(dc, σs)
    @test sum(abs2, amplitude(dc, σs)) > 0
end
