using ThreeBodyDecays
using Test


ms = ThreeBodyMasses(1.1, 3.3, 5.5; m0 = 20.0)
σs = Invariants(ms; σ3 = 4.5^2, σ2 = 10.1^2)

@testset "Invariants structure" begin
    #
    @test σs.σ1 == σs[1]
    @test σs.σ2 == σs[2]
    @test σs.σ3 == σs[3]
    #
    @test_throws BoundsError σs[5]
    @test σs == Invariants(; σs.σ1, σs.σ2, σs.σ3)
    @test σs == Invariants(σs.σ1, σs.σ2, σs.σ3)
end

@testset "operations and integrate" begin
    @test sum(σs) == sum(ms^2)
    @test length(σs) == 3
end

@testset "Creation from two given" begin
    @test σs == Invariants(ms; σ1 = σs.σ1, σ2 = σs.σ2)
    @test σs == Invariants(ms; σ2 = σs.σ2, σ3 = σs.σ3)
    @test σs == Invariants(ms; σ3 = σs.σ3, σ1 = σs.σ1)
    #
    @test_throws ErrorException Invariants(ms; σ1 = 1.0, σ2 = 1.0, σ3 = 1.0)
    #
    @test Kibble(σs, ms^2) < 0
end


@testset "consistency of cosθij and σiofk, σjofk functions" begin
    z23 = cosθ23(σs, ms^2)
    @test σs.σ3 ≈ σ3of1(z23, σs.σ1, ms^2)
    @test σs.σ2 ≈ σ2of1(z23, σs.σ1, ms^2)
    #
    z12 = cosθ12(σs, ms^2)
    @test σs.σ1 ≈ σ1of3(z12, σs.σ3, ms^2)
    @test σs.σ2 ≈ σ2of3(z12, σs.σ3, ms^2)
    #
    z31 = cosθ31(σs, ms^2)
    @test σs.σ3 ≈ σ3of2(z31, σs.σ2, ms^2)
    @test σs.σ1 ≈ σ1of2(z31, σs.σ2, ms^2)
end

@testset "Circle to Origin" begin
    @test circleorigin(-3, (1, 2, 3)) == (1, 2, 3)
    @test circleorigin(-2, (3, 1, 2)) == (1, 2, 3)
    @test circleorigin(-1, (2, 3, 1)) == (1, 2, 3)
end

@testset "Mapping of [0,1]x[0,1] to invariants" begin
    @test Kibble(x2σs([0.01, 0.99], ms; k = 3), ms^2) < 0
    @test Kibble(x2σs([0.99, 0.01], ms; k = 3), ms^2) < 0
    #
    @test Kibble(x2σs([0.01, 0.99], ms; k = 3), ms^2) < 0
    @test Kibble(x2σs([0.99, 0.01], ms; k = 3), ms^2) < 0
    #
    @test Kibble(x2σs([0.01, 0.99], ms; k = 3), ms^2) < 0
    @test Kibble(x2σs([0.99, 0.01], ms; k = 3), ms^2) < 0
end

@testset "Random points" begin
    σs = randomPoint(ms)
    @test Kibble(σs, ms^2) < 0
end

tbs = ThreeBodySystem(ms)

@testset "Three-body system" begin
    @test tbs == ThreeBodySystem(ms = ms, two_js = ThreeBodySpins(0, 0, 0; two_h0 = 0))
    @test tbs == ThreeBodySystem(ms, ThreeBodySpins(0, 0, 0; two_h0 = 0))
    @test tbs == ThreeBodySystem(
        ms.m1,
        ms.m2,
        ms.m3;
        m0 = ms.m0,
        two_js = ThreeBodySpins(0, 0, 0; two_h0 = 0),
    )
end

tbsS = ThreeBodySystem(ms, ThreeBodySpins(1, 3, 5; two_h0 = 1))
rp = randomPoint(tbsS)

@testset "Random points" begin
    @test Kibble(rp.σs, tbsS.ms^2) < 0
    @test Tuple(rp.two_λs) ∈ collect(itr(tbsS.two_js))
end


@testset "Aligned four vectors" begin
    # Test 1: Check four-momentum conservation
    p1, p2, p3 = aligned_four_vectors(σs, ms; k = 1)
    @test all(p1 .+ p2 .+ p3 .≈ (0, 0, 0, ms.m0))
    #
    # aligned vector k, px=0, is on kth place
    @test aligned_four_vectors(σs, ms; k = 1)[1][1] == 0.0
    @test aligned_four_vectors(σs, ms; k = 2)[2][1] == 0.0
    @test aligned_four_vectors(σs, ms; k = 3)[3][1] == 0.0

    metric = (-1, -1, -1, 1)
    # Test 2: Check mass shell conditions
    @test sum(p1 .^ 2 .* metric) ≈ ms.m1^2
    @test sum(p2 .^ 2 .* metric) ≈ ms.m2^2
    @test sum(p3 .^ 2 .* metric) ≈ ms.m3^2

    # Test 3: Verify invariant masses
    @test sum((p2 .+ p3) .^ 2 .* metric) ≈ σs.σ1
    @test sum((p3 .+ p1) .^ 2 .* metric) ≈ σs.σ2
    @test sum((p1 .+ p2) .^ 2 .* metric) ≈ σs.σ3
end

@testset "Aligned four vectors are cyclic" begin
    p1_1, p2_1, p3_1 = aligned_four_vectors(σs, ms; k = 1)
    #
    p1_2, p2_2, p3_2 = aligned_four_vectors(
        MandelstamTuple{Float64}((σs[3], σs[1], σs[2])),
        ThreeBodyMasses(ms[3], ms[1], ms[2]; ms.m0);
        k = 2,
    )
    #
    p1_3, p2_3, p3_3 = aligned_four_vectors(
        MandelstamTuple{Float64}((σs[2], σs[3], σs[1])),
        ThreeBodyMasses(ms[2], ms[3], ms[1]; ms.m0);
        k = 3,
    )
    #
    @test all(p1_1 .≈ p2_2 .≈ p3_3)
    @test all(p2_1 .≈ p3_2 .≈ p1_3)
    @test all(p3_1 .≈ p1_2 .≈ p2_3)
end
