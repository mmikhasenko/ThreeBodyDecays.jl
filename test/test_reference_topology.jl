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

    sign_for(k, r) = k == r ? 1 : (-1)^tbs.two_js[6-k-r]
    reference_model = ThreeBodyDecay(chains, coefficients, names; reference_k = 1)
    reference_intensity = unpolarized_intensity(reference_model, σs)
    for k = 1:3
        model = ThreeBodyDecay(chains, coefficients, names; reference_k = k)
        @test ThreeBodyDecays.reference_k(model) == model.reference_k == k
        @test ThreeBodyDecays.reference_k(
            ThreeBodyDecay(; chains, couplings = coefficients, names, reference_k = k),
        ) == k
        @test ThreeBodyDecays.reference_k(
            ThreeBodyDecay("R1" => (coefficients[1], chain1); reference_k = k),
        ) == k
        @test ThreeBodyDecays.reference_k(model[1]) == k
        @test ThreeBodyDecays.reference_k(model[1:2]) == k
        @test ThreeBodyDecays.reference_k(model[[2, 1]]) == k
        @test ThreeBodyDecays.reference_k(model[(true, false)]) == k
        @test ThreeBodyDecays.reference_k(vcat(model[1], model[2])) == k
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
        converted = map(chains, coefficients) do chain, c
            sign_for(chain.k, k) * sign_for(chain.k, 1) * c
        end
        converted_model = ThreeBodyDecay(chains, converted, names; reference_k = k)
        @test unpolarized_intensity(converted_model, σs) ≈ reference_intensity
        if k == 2
            @test !isapprox(intensity, reference_intensity)
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

@testset "Coupling conversion preserves coherent intensities" begin
    # All eight integer/half-integer final-spin sectors, including several fermions.
    # Maximal spin couplings with L=l=0 give allowed, nonzero vertices in every sector.
    for bits in Iterators.product(0:1, 0:1, 0:1)
        two_js = ThreeBodySpins(bits...; two_h0 = sum(bits))
        tbs = ThreeBodySystem(1.1, 1.3, 1.5; m0 = 5.5, two_js)
        chains = map(1:3) do k
            two_j = sum(bits) - bits[k]
            DecayChain(;
                k,
                two_j,
                tbs,
                Xlineshape = σ -> 1 / (8.0 - σ - 0.7im),
                HRk = Vertex(RecouplingLS((0, sum(bits)))),
                Hij = Vertex(RecouplingLS((0, two_j))),
            )
        end
        names = ["chain $k" for k = 1:3]
        c = [1.0 + 0.3im, -0.4 + 0.7im, 0.6 - 0.2im]
        sign_for(k, r) = k == r ? 1 : (-1)^two_js[6-k-r]
        for r = 1:3
            source = ThreeBodyDecay(chains, c, names; reference_k = r)
            points = [x2σs(x, tbs.ms; k = 1) for x in ([0.2, 0.3], [0.5, 0.6], [0.8, 0.4])]
            expected = [unpolarized_intensity(source, σs) for σs in points]
            @test all(>(1e-8), expected)
            for s = 1:3
                r == s && continue # Relabeling a reference to itself adds no coverage.
                converted = [sign_for(k, s) * sign_for(k, r) * c[k] for k = 1:3]
                target = ThreeBodyDecay(chains, converted, names; reference_k = s)
                for (σs, value) in zip(points, expected)
                    @test unpolarized_intensity(target, σs) ≈ value
                end
            end
        end
    end
end
