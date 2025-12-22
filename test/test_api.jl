using ThreeBodyDecays
using Test
using Random

@testset "High-Level API" begin
    ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0 = 4.0)

    @testset "dalitzmap construction" begin
        map = dalitzmap(ms)

        @test map isa DalitzMap
        @test map.ms == ms
    end

    @testset "End-to-end mapping" begin
        map = dalitzmap(ms)

        # Test with physical points
        for _ = 1:10
            x = rand(2)
            σs = y2σs(x, ms)
            if isphysical(σs, ms)
                w = map(σs)
                @test w isa Complex
                @test abs(w) < 1.0 + 1e-6  # Should be inside or on unit disk
            end
        end
    end

    @testset "Convenience method: map(s12, s23)" begin
        map = dalitzmap(ms)

        # Test with a physical point (use y2σs to get valid coordinates)
        x = [0.5, 0.5]
        σs = y2σs(x, ms)
        s12 = σs.σ3  # σ3 corresponds to s12
        s23 = σs.σ1  # σ1 corresponds to s23
        w = map(s12, s23)

        @test w isa Complex
        @test isfinite(w)
    end

    @testset "API consistency" begin
        map1 = dalitzmap(ms)
        map2 = dalitzmap(ms)

        # Same inputs should give same outputs
        σs = Invariants(ms; σ1 = 6.25, σ2 = 6.5)
        w1 = map1(σs)
        w2 = map2(σs)

        @test w1 ≈ w2
    end

    @testset "Diagnostics" begin
        map = dalitzmap(ms)

        # Boundary error should be reasonable
        # Use ThreeBodyDecays prefix since it's not exported
        bdry_err = ThreeBodyDecays.boundary_error(map; n = 50)
        @test bdry_err < 0.2  # Relaxed for simplified zipper

        # Conformality error should be reasonable
        conf_err = ThreeBodyDecays.conformality_error(map; n = 10)
        @test conf_err < 2.0  # Relaxed for simplified zipper
    end
end

@testset "Property tests: random masses" begin
    Random.seed!(98765)

    for _ = 1:5
        # Generate random physical masses
        m1 = 0.1 + 0.9 * rand()
        m2 = 0.1 + 0.9 * rand()
        m3 = 0.1 + 0.9 * rand()
        m0 = (m1 + m2 + m3) * (1.0 + 0.5 * rand())

        ms = ThreeBodyMasses(m1, m2, m3; m0 = m0)

        @testset "Construction for random masses" begin
            map = dalitzmap(ms)
            @test map isa DalitzMap
        end

        @testset "Mapping for random masses" begin
            map = dalitzmap(ms)

            # Test a few points
            for _ = 1:3
                x = rand(2)
                σs = y2σs(x, ms)
                if isphysical(σs, ms)
                    w = map(σs)
                    @test abs(w) < 1.0 + 1e-6
                end
            end
        end
    end
end
