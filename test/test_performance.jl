using Test
using ThreeBodyDecays
using StaticArrays

using BenchmarkTools


# test model constraction
model = let
	tbs = ThreeBodySystem(
		ms = ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09),
		two_js = ThreeBodySpins(0, 0, 0; two_h0 = 2))
	dpp = randomPoint(tbs)
	#
	two_j = 2
	ch1 = DecayChain(;
		k = 1,
		two_j,
		Xlineshape = σ -> 1 / (4.1^2 - σ - 0.1im),
		Hij = RecouplingLS(two_j, (two_j, 0), 0, 0),
		HRk = RecouplingLS(tbs.two_js[4], (two_j, two_j), two_j, 0),
		tbs)
	ch2 = DecayChain(ch1; k = 2)
	ch3 = DecayChain(ch1; k = 3)
	# 
	ThreeBodyDecay(
		"K(892)" .=> [(4.0, ch1), (2.0, ch2), (3.0, ch3)])
end


@testset "Amplitude with ThreeBodyDecay" begin
	ms = masses(model)
	σs = x2σs([0.1, 0.9], ms; k = 2)
	refA_xxx1 = 0.5360686514001554 + 0.005104210009658071im
	@test amplitude(model, σs, spins(model)) ≈ refA_xxx1
	@test amplitude(model, σs, spins(model); refζs = (1, 2, 3, 1)) ≈ refA_xxx1
	@test amplitude(model, σs, spins(model); refζs = (2, 3, 2, 1)) ≈ refA_xxx1
	@test amplitude(model, σs, spins(model); refζs = (2, 2, 3, 1)) ≈ refA_xxx1
end

σs = x2σs([0.1, 0.9], masses(model); k = 2)
refI = 0.574791303947608
unpolarized_intensity(model, σs) ≈ refI


@btime unpolarized_intensity($model, $σs) # 14 us
@btime aligned_amplitude($model, $σs) # 14 us


ch1 = model.chains[1]
two_λs = spins(model)
@btime amplitude($ch1, $σs, $two_λs) # 1.3 us

@btime aligned_amplitude($ch1, $σs) # 1.19 us

@btime sum(abs2, aligned_amplitude($ch1, $σs)) # 1.2 us

unpolarized_intensity(model[1], σs) / abs2(model.couplings[1])

