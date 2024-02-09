using Test
using ThreeBodyDecays

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
		"K(892)" .=> [(1.0, ch1), (1.0, ch2), (1.0, ch3)])
end

@testset "masses with ThreeBodySystem" begin
	@test masses(model)[1] == 0.141
	@test masses(model)[2] == 0.142
	@test masses(model)[3] == 0.143
	@test masses(model).m0 == 3.09
end

@testset "spins with ThreeBodySystem" begin
	@test spins(model)[1] == 0
	@test spins(model)[2] == 0
	@test spins(model)[3] == 0
	@test spins(model).two_h0 == 2
end

@testset "Amplitude with ThreeBodyDecay" begin
    ms = masses(model)
	σs = x2σs_ki(2, [0.1, 0.9], ms)
	@test amplitude(model, σs, spins(model)) == 0.1582599633375949 + 0.0014291623024238807im
end
