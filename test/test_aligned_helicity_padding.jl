using Test
using LinearAlgebra
using ThreeBodyDecays

# Independent helicity formula: no shifted array indices or Tullio padding.
function explicit_aligned_helicities(chain, σs)
    k = chain.k
    i, j, _ = ThreeBodyDecays.ij_from_k(k)
    js = chain.tbs.two_js
    J = chain.two_j
    ms² = chain.tbs.ms^2
    cosθ = cosθij(σs, ms²; k)
    lineshape =
        chain.Xlineshape(σs[k]) *
        chain.HRk.ff(ms²[4], σs[k], ms²[k]) *
        chain.Hij.ff(σs[k], ms²[i], ms²[j])
    result = zeros(ComplexF64, Tuple(js) .+ 1)
    for I in CartesianIndices(result)
        λ = ntuple(a -> -js[a] + 2 * (I[a] - 1), 4)
        λR, λR′ = λ[4] + λ[k], λ[i] - λ[j]
        (abs(λR) > J || abs(λR′) > J || isodd(J - λR) || isodd(J - λR′)) && continue
        production =
            amplitude(chain.HRk.h, (λR, λ[k]), (js[4], J, js[k])) *
            (-1)^((js[k] - λ[k]) ÷ 2)
        decay =
            amplitude(chain.Hij.h, (λ[i], λ[j]), (J, js[i], js[j])) *
            (-1)^((js[j] - λ[j]) ÷ 2)
        rotation = ThreeBodyDecays.wignerd_doublearg_sign(J, λR, λR′, cosθ, true)
        result[I] = sqrt(J + 1) * lineshape * production * rotation * decay
    end
    return result
end

@testset "Aligned amplitude covers all external helicities" begin
    # Local spin order is daughter i, daughter j, spectator, parent.
    # First case needs extra production padding, second extra decay padding;
    # the scalar isobar requires both, despite its one-element Wigner matrix.
    cases = (
        (spins = (1, 2, 3, 4), J = 1, LS = (4, 2), ls = (4, 3)),
        (spins = (3, 4, 1, 4), J = 1, LS = (4, 0), ls = (0, 1)),
        (spins = (2, 2, 2, 2), J = 0, LS = (0, 2), ls = (0, 0)),
    )
    ms = ThreeBodyMasses(0.7, 1.1, 0.9; m0 = 4.2)
    for case in cases, swap in (false, true), k = 1:3
        i, j, _ = ThreeBodyDecays.ij_from_k(k)
        local_spins = swap ? case.spins[[2, 1, 3, 4]] : case.spins
        js = ntuple(
            a -> a == i ? local_spins[1] : a == j ? local_spins[2] : local_spins[3],
            3,
        )
        tbs = ThreeBodySystem(ms, ThreeBodySpins(js...; two_h0 = local_spins[4]))
        chain = DecayChain(;
            k,
            two_j = case.J,
            Xlineshape = σ -> 1.0 + 0.2im,
            HRk = Vertex(RecouplingLS(case.LS)),
            Hij = Vertex(RecouplingLS(case.ls)),
            tbs,
        )
        for point in ([0.3, 0.4], [0.65, 0.7])
            σs = x2σs(point, ms; k)
            expected = explicit_aligned_helicities(chain, σs)
            @test norm(expected) > 0.1
            @test aligned_amplitude(chain, σs) ≈ expected
        end
    end
end
