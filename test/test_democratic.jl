using ThreeBodyDecays
using Test
using Random

# Include boundary.jl first (domain.jl depends on it)
ThreeBodyDecays.include(joinpath(@__DIR__, "..", "src", "dalitz_disk", "boundary.jl"))
# Then include domain.jl
ThreeBodyDecays.include(joinpath(@__DIR__, "..", "src", "dalitz_disk", "domain.jl"))

# Import functions
import ThreeBodyDecays: DemocraticMap, compute_Λ, BoundaryFunction

@testset "Democratic Coordinates" begin
    ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0 = 4.0)
    dm = DemocraticMap(ms)

    @testset "Basic forward transformation" begin
        # Test with a known point
        σs = Invariants(ms; σ1 = 6.25, σ2 = 6.25)
        z = dm(σs)

        @test z isa Complex
        @test isfinite(z)
    end

    @testset "Round-trip: forward then inverse" begin
        # Test multiple points
        for _ = 1:10
            # Generate random physical point
            x = rand(2)
            σs = y2σs(x, ms)
            if isphysical(σs, ms)
                z = dm(σs)
                σs_recovered = dm(z)

                # Check round-trip accuracy
                @test σs_recovered.σ1 ≈ σs.σ1 atol = 1e-10
                @test σs_recovered.σ2 ≈ σs.σ2 atol = 1e-10
                @test σs_recovered.σ3 ≈ σs.σ3 atol = 1e-10
            end
        end
    end

    @testset "Permutation covariance: cyclic permutation" begin
        σs = Invariants(ms; σ1 = 6.25, σ2 = 6.5)
        z_original = dm(σs)

        # Cyclic permutation: (σ₁, σ₂, σ₃) → (σ₂, σ₃, σ₁)
        σs_perm = MandelstamTuple{Float64}((σs.σ2, σs.σ3, σs.σ1))
        z_perm = dm(σs_perm)

        # Should rotate by ω² (one forward cyclic permutation)
        ω = exp(2π * im / 3)
        @test z_perm ≈ z_original * ω^2 atol = 1e-10
    end

    @testset "Permutation covariance: transposition" begin
        σs = Invariants(ms; σ1 = 6.25, σ2 = 6.5)
        z_original = dm(σs)

        # Transposition: (σ₁, σ₂, σ₃) → (σ₂, σ₁, σ₃)
        σs_trans = MandelstamTuple{Float64}((σs.σ2, σs.σ1, σs.σ3))
        z_trans = dm(σs_trans)

        # Should reflect (up to rotation)
        # For transposition, z should be related by complex conjugation and rotation
        # This is a weaker test - just check it's different but related
        @test abs(z_trans) ≈ abs(z_original) atol = 1e-10
    end

    @testset "D₃ symmetry for equal masses" begin
        ms_equal = ThreeBodyMasses(1.0, 1.0, 1.0; m0 = 4.0)
        dm_equal = DemocraticMap(ms_equal)
        bf = BoundaryFunction(ms_equal)

        # Sample boundary points
        θs = [0.0, 2π / 3, 4π / 3]
        zs = [dm_equal(bf(θ)) for θ in θs]

        # Rotate first point by 120° - should match second point
        ω = exp(2π * im / 3)
        @test abs(zs[2] - zs[1] * ω) < 0.1 || abs(zs[3] - zs[1] * ω) < 0.1
    end
end

@testset "Property tests: random masses" begin
    Random.seed!(54321)

    for _ = 1:10
        # Generate random physical masses
        m1 = 0.1 + 0.9 * rand()
        m2 = 0.1 + 0.9 * rand()
        m3 = 0.1 + 0.9 * rand()
        m0 = (m1 + m2 + m3) * (1.0 + 0.5 * rand())

        ms = ThreeBodyMasses(m1, m2, m3; m0 = m0)
        dm = DemocraticMap(ms)

        @testset "Round-trip for random masses" begin
            # Test a few random points
            for _ = 1:5
                x = rand(2)
                σs = y2σs(x, ms)
                if isphysical(σs, ms)
                    z = dm(σs)
                    σs_recovered = dm(z)

                    @test σs_recovered.σ1 ≈ σs.σ1 atol = 1e-8
                    @test σs_recovered.σ2 ≈ σs.σ2 atol = 1e-8
                    @test σs_recovered.σ3 ≈ σs.σ3 atol = 1e-8
                end
            end
        end

        @testset "Permutation covariance for random masses" begin
            x = rand(2)
            σs = y2σs(x, ms)
            if isphysical(σs, ms)
                z_original = dm(σs)

                # Cyclic permutation
                σs_perm = MandelstamTuple{typeof(ms.m0)}((σs.σ2, σs.σ3, σs.σ1))
                z_perm = dm(σs_perm)
                ω = exp(2π * im / 3)

                @test z_perm ≈ z_original * ω^2 atol = 1e-8
            end
        end
    end
end
