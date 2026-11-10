using Test
using LinearAlgebra
using ThreeBodyDecays
using ThreeBodyDecays.PartialWaveFunctions: wignerD_doublearg

@testset "Reference topology belongs to model" begin
    tbs = ThreeBodySystem(
        1.0,
        2.0,
        3.0;
        m0 = 7.9,
        two_js = ThreeBodySpins(1, 2, 0; two_h0 = 1),
    )
    chain1 = DecayChain(;
        k = 1,
        Xlineshape = x -> 1.0,
        two_j = 4,
        HRk = Vertex(RecouplingLS((4, 3))),
        Hij = Vertex(RecouplingLS((4, 2))),
        tbs,
    )
    chain3 = DecayChain(;
        k = 3,
        Xlineshape = x -> 1.0,
        two_j = 1,
        HRk = Vertex(RecouplingLS((2, 1))),
        Hij = Vertex(RecouplingLS((2, 3))),
        tbs,
    )
    chains = (chain1, chain3)
    coefficients = (1.0 + 0.3im, -0.4 + 0.2im)
    names = ("R1", "R3")
    σs = Invariants(tbs.ms; σ1 = 32.5, σ2 = 21.6)
    orientation = (α = 0.3, cosβ = 0.4, γ = 0.5)

    @test_throws UndefKeywordError ThreeBodyDecay(chains, coefficients, names)
    @test_throws ArgumentError ThreeBodyDecay((), (), (); reference_k = 1)
    for invalid in (0, 4, -1, true)
        @test_throws ArgumentError ThreeBodyDecay(
            chains,
            coefficients,
            names;
            reference_k = invalid,
        )
        @test_throws ArgumentError amplitude(chain1, σs; reference_k = invalid)
    end

    for other_tbs in (
        ThreeBodySystem(1.0, 2.0, 3.0; m0 = 8.0, two_js = tbs.two_js),
        ThreeBodySystem(tbs.ms, ThreeBodySpins(1, 0, 0; two_h0 = 1)),
    )
        incompatible = DecayChain(;
            k = chain1.k,
            Xlineshape = chain1.Xlineshape,
            two_j = chain1.two_j,
            HRk = chain1.HRk,
            Hij = chain1.Hij,
            tbs = other_tbs,
        )
        @test_throws ArgumentError ThreeBodyDecay(
            (chain1, incompatible),
            coefficients,
            names;
            reference_k = 1,
        )
    end

    reference_intensity = nothing
    for k = 1:3
        model = ThreeBodyDecay(chains, coefficients, names; reference_k = k)
        @test reference_k(model) == model.reference_k == k
        @test reference_k(
            ThreeBodyDecay(; chains, couplings = coefficients, names, reference_k = k),
        ) == k
        @test reference_k(
            ThreeBodyDecay("R1" => (coefficients[1], chain1); reference_k = k),
        ) == k
        @test reference_k(model[1]) == k
        @test reference_k(model[1:2]) == k
        @test reference_k(model[[2, 1]]) == k
        @test reference_k(model[(true, false)]) == k
        @test reference_k(vcat(model[1], model[2])) == k
        other = ThreeBodyDecay(chains, coefficients, names; reference_k = mod1(k + 1, 3))
        @test_throws ArgumentError vcat(model, other)
        @test_throws MethodError amplitude(model, σs; reference_k = mod1(k + 1, 3))
        @test_throws MethodError unpolarized_intensity(
            model,
            σs;
            reference_k = mod1(k + 1, 3),
        )
        @test_throws MethodError unpolarized_intensity(
            model,
            orientation,
            σs;
            reference_k = mod1(k + 1, 3),
        )
        @test_throws MethodError amplitude(model, σs; refζs = (k, k, k, k))
        @test_throws MethodError amplitude(chain1, σs; refζs = (k, k, k, k))

        expected = sum(
            c * amplitude(chain, σs; reference_k = k) for
            (c, chain) in zip(coefficients, chains)
        )
        @test amplitude(model, σs) ≈ expected
        for helicities in itr(spins(model))
            indices = map((j, λ) -> (j + λ) ÷ 2 + 1, spins(model), helicities)
            @test amplitude(model, σs, helicities) ≈ expected[indices...]
        end
        @test eltype(amplitude(chain1, σs; reference_k = k)) <: Real
        @test any(!iszero, imag.(amplitude(chain1, orientation, σs; reference_k = k)))
        rotated = amplitude(model, orientation, σs)
        D = conj.(wignerD_doublearg(spins(model).two_h0, orientation...))
        expected_rotation = similar(expected)
        for I in CartesianIndices(expected)
            expected_rotation[I] =
                sum(D[I[4], n] * expected[I[1], I[2], I[3], n] for n in axes(expected, 4))
        end
        @test rotated ≈ expected_rotation
        @test unpolarized_intensity(model, orientation, σs) ≈ sum(abs2, expected)
        intensity = unpolarized_intensity(model, σs)
        if isnothing(reference_intensity)
            reference_intensity = intensity
        else
            @test intensity ≈ reference_intensity
        end

        # Checking the polarized cache catches dropping model.reference_k even when
        # a spin-summed overlap happens to be invariant under a basis change.
        cache = chain_amplitudes(model, [σs])
        for (i, chain) in enumerate(chains)
            @test cache.amplitudes[i, :, 1] ≈ vec(amplitude(chain, σs; reference_k = k))
        end
        weighted = sum(c * cache.amplitudes[i, :, 1] for (i, c) in enumerate(coefficients))
        @test weighted ≈ vec(expected)
        @test total_intensity(physical_overlap(chain_overlap_matrix(cache), coefficients)) ≈
              intensity
    end
    for chain in chains
        @test amplitude(chain, σs; reference_k = chain.k) ≈
              ThreeBodyDecays.aligned_amplitude(chain, σs)
    end
end
