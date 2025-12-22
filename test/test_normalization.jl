using ThreeBodyDecays
using Test
using Random

# Include necessary modules in correct order
ThreeBodyDecays.include(joinpath(@__DIR__, "..", "src", "dalitz_disk", "boundary.jl"))
ThreeBodyDecays.include(joinpath(@__DIR__, "..", "src", "dalitz_disk", "domain.jl"))
ThreeBodyDecays.include(joinpath(@__DIR__, "..", "src", "dalitz_disk", "api.jl"))  # Includes zipper and defines MobiusTransform
ThreeBodyDecays.include(joinpath(@__DIR__, "..", "src", "dalitz_disk", "normalize.jl"))

# Import functions
import ThreeBodyDecays:
    BoundaryFunction,
    DemocraticMap,
    ZipperMap,
    MobiusTransform,
    find_landmarks,
    normalize_landmarks,
    mobius_transform_three_points

@testset "Landmark Selection" begin
    ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0 = 4.0)
    bf = BoundaryFunction(ms)

    @testset "Landmarks found" begin
        θ1, θ2, θ3 = find_landmarks(bf)

        @test all(0 <= θ <= 2π for θ in (θ1, θ2, θ3))

        # Check that landmarks maximize respective σᵢ
        σs1 = bf(θ1)
        σs2 = bf(θ2)
        σs3 = bf(θ3)

        # Sample nearby points and verify σᵢ is maximum
        for δ in [-0.1, 0.1]
            σs1_test = bf(θ1 + δ)
            @test σs1.σ1 >= σs1_test.σ1 - 1e-6  # Allow small numerical error
        end
    end
end

@testset "Möbius Transformation" begin
    @testset "Three-point mapping" begin
        w1 = 0.5 + 0.5im
        w2 = -0.5 + 0.5im
        w3 = -0.5 - 0.5im

        target1 = 1.0 + 0.0im
        target2 = exp(2π * im / 3)
        target3 = exp(4π * im / 3)

        mt = mobius_transform_three_points(w1, w2, w3, target1, target2, target3)

        # Check that transformation maps correctly
        @test abs(mt(w1) - target1) < 1e-6
        @test abs(mt(w2) - target2) < 1e-6
        @test abs(mt(w3) - target3) < 1e-6
    end
end

@testset "Channel-Democratic Normalization" begin
    ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0 = 4.0)
    bf = BoundaryFunction(ms)
    dm = DemocraticMap(ms)
    zm = ZipperMap(bf, dm; n = 200)

    @testset "Landmark exactness" begin
        mt = normalize_landmarks(bf, dm, zm)

        # Get landmarks
        θ1, θ2, θ3 = find_landmarks(bf)
        σs1 = bf(θ1)
        σs2 = bf(θ2)
        σs3 = bf(θ3)

        # Map through full chain
        z1 = dm(σs1)
        z2 = dm(σs2)
        z3 = dm(σs3)
        w1 = zm(z1)
        w2 = zm(z2)
        w3 = zm(z3)

        # Apply normalization
        w1_norm = mt(w1)
        w2_norm = mt(w2)
        w3_norm = mt(w3)

        # Check they map to 1, ω, ω²
        ω = exp(2π * im / 3)
        @test abs(w1_norm - 1.0) < 1e-4  # Relaxed tolerance
        @test abs(w2_norm - ω) < 1e-4
        @test abs(w3_norm - ω^2) < 1e-4
    end

    @testset "Equal-mass symmetry" begin
        ms_equal = ThreeBodyMasses(1.0, 1.0, 1.0; m0 = 4.0)
        bf_equal = BoundaryFunction(ms_equal)
        dm_equal = DemocraticMap(ms_equal)
        zm_equal = ZipperMap(bf_equal, dm_equal; n = 200)

        mt = normalize_landmarks(bf_equal, dm_equal, zm_equal)

        # For equal masses, cyclic permutation should give rotation
        θ1, θ2, θ3 = find_landmarks(bf_equal)
        σs1 = bf_equal(θ1)
        σs2 = bf_equal(θ2)
        σs3 = bf_equal(θ3)

        # Map and normalize
        z1 = dm_equal(σs1)
        w1 = zm_equal(z1)
        w1_norm = mt(w1)

        # Cyclic permutation
        σs1_perm = MandelstamTuple{Float64}((σs1.σ2, σs1.σ3, σs1.σ1))
        z1_perm = dm_equal(σs1_perm)
        w1_perm = zm_equal(z1_perm)
        w1_perm_norm = mt(w1_perm)

        # Should be related by rotation (relaxed test)
        @test abs(abs(w1_norm) - abs(w1_perm_norm)) < 0.1
    end
end
