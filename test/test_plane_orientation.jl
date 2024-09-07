using ThreeBodyDecays
using Test

reference_model = let
    ms = ThreeBodyMasses(1, 2, 3; m0 = 7.9)
    two_js = ThreeBodySpins(1, 2, 0; two_h0 = 1)
    tbs = ThreeBodySystem(ms, two_js)

    chain1 = DecayChain(;
        k = 1,
        Xlineshape = x -> 1.0 + 0.0im,
        two_j = 4,
        HRk = RecouplingLS((4, 3)),
        Hij = RecouplingLS((4, 2)),
        tbs,
    )

    chain3 = DecayChain(;
        k = 3,
        Xlineshape = x -> 1.0 + 0.0im,
        two_j = 1,
        HRk = RecouplingLS((2, 1)),
        Hij = RecouplingLS((2, 3)),
        tbs,
    )

    ThreeBodyDecay(["ξ(23)", "ξ(12)"] .=> zip([1.0 + 0.0im, 1.0 + 0.0im], [chain1, chain3]))
end

# kinematics
angles_test = (α = 0.3, cosβ = 0.4, γ = 0.5)
σs_test = Invariants(masses(reference_model); σ1 = 32.5, σ2 = 21.6)

# test 1: same unpolarized intensity
I_dpd = sum(abs2, amplitude(reference_model, σs_test))
I_ori = sum(abs2, amplitude(reference_model, angles_test, σs_test))

@testset "Same unpolarized intensity" begin
    @test I_dpd ≈ I_ori
end

@testset "non-zero interference" begin
    incoherent =
        sum(abs2, amplitude(reference_model[1], angles_test, σs_test)) +
        sum(abs2, amplitude(reference_model[2], angles_test, σs_test))

    coherent = sum(abs2, amplitude(reference_model, angles_test, σs_test))
    #
    interference_fraction = (incoherent - coherent) / coherent
    @test round(interference_fraction * 100; digits = 2) == 2.87 # 2.87%
end



# Other conventions
using ThreeBodyDecays.Parameters
using ThreeBodyDecays.Tullio

function amplitude_helicity(model, angles, pars...; reference_chain, kw...)
    #
    @unpack γ = angles
    two_js = spins(model)
    all_two_js_ijk = itr(collect(two_js)[1:3])
    #
    phase = map(all_two_js_ijk) do two_λs
        cis(γ * two_λs[reference_chain] / 2)
    end
    #
    refζs = fill(reference_chain, 4)
    F0 = amplitude(model, angles, pars...; refζs, kw...)
    #
    F = similar(F0)
    @tullio F[_i, _j, _k, _z] = phase[_i, _j, _k] * F0[_i, _j, _k, _z]
    return F
end

function amplitude_minusphi(model, angles, pars...; reference_chain, kw...)
    refζi, refζj = ij_from_k(reference_chain)
    #
    @unpack α, cosβ, γ = angles
    γ′ = α + γ
    #
    two_js = spins(model)
    all_two_js_ijk = itr(collect(two_js)[1:3])
    phase = map(all_two_js_ijk) do two_λs
        cis(γ′ * (two_λs[reference_chain] - two_λs[refζi] + two_λs[refζj]) / 2)
    end
    refζs = fill(reference_chain, 4)
    F0 = amplitude(model, angles, pars...; refζs, kw...)
    #
    F = similar(F0)
    @tullio F[_i, _j, _k, _z] = phase[_i, _j, _k] * F0[_i, _j, _k, _z]
    return F
end

A_helicity = amplitude_helicity(reference_model, angles_test, σs_test; reference_chain = 1)
A_minusphi = amplitude_minusphi(reference_model, angles_test, σs_test; reference_chain = 1)

@testset "Minus phi and helicity conventions" begin

    @test A_helicity[2, 3, 1, 1] ≈ -0.3203152615528447 - 0.14539779706573902im
    @test A_minusphi[2, 3, 1, 1] ≈ -0.34299056035837566 + 0.07810161125280723im

    @test sum(abs2, A_helicity) ≈ sum(abs2, A_minusphi)
end
