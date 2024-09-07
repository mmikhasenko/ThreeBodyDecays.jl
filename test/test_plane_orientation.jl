using ThreeBodyDecays
using Test

reference_model = let
    ms = ThreeBodyMasses(1, 2, 3; m0=7.9)
    two_js = ThreeBodySpins(1,2,0; two_h0=1)
    tbs = ThreeBodySystem(ms, two_js);

    chain1 = DecayChain(; k=1,
        Xlineshape=x->1.0+0.0im,
        two_j=4,
        HRk = RecouplingLS((4, 3)),
        Hij = RecouplingLS((4, 2)),
        tbs);

    chain3 = DecayChain(; k=3,
        Xlineshape=x->1.0+0.0im,
        two_j=1,
        HRk = RecouplingLS((2, 1)),
        Hij = RecouplingLS((2, 3)),
        tbs);

    ThreeBodyDecay(["ξ(23)", "ξ(12)"] .=>
        zip([1.0+0.0im, 1.0+0.0im], [chain1, chain3]));
end

# kinematics
angles_test = (α = 0.3, cosβ=0.4, γ=0.5)
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
	interference_fraction = (incoherent-coherent)/coherent
    @test round(interference_fraction * 100; digits=2) == 2.87 # 2.87%
end
