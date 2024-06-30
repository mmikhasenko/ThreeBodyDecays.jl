using Test
using ThreeBodyDecays
using StaticArrays


# test model constraction
model = let
    tbs = ThreeBodySystem(
        ms=ThreeBodyMasses(0.141, 0.142, 0.143; m0=3.09),
        two_js=ThreeBodySpins(0, 0, 0; two_h0=2))
    dpp = randomPoint(tbs)
    #
    two_j = 2
    ch1 = DecayChain(;
        k=1,
        two_j,
        Xlineshape=σ -> 1 / (4.1^2 - σ - 0.1im),
        Hij=RecouplingLS((two_j, 0)),
        HRk=RecouplingLS((two_j, two_j)),
        tbs)
    ch2 = DecayChain(ch1; k=2)
    ch3 = DecayChain(ch1; k=3)
    #
    ThreeBodyDecay(
        "K(892)" .=> [(4.0, ch1), (2.0, ch2), (3.0, ch3)])
end

@testset "Construction with SVectors" begin
    _model = ThreeBodyDecay(
        ("K(892)", "K(892)", "K(892)") .=> zip(model.couplings, model.chains))
    @test _model == model
    _model = ThreeBodyDecay(
        SubString("K(892))", 1, 6) .=> zip(model.couplings, model.chains))
    @test _model isa ThreeBodyDecay
end

@testset "Properties of the model" begin
    @test length(model) == 3
end

@testset "system, masses, spins with ThreeBodySystem" begin
    @test system(model) == model.chains[1].tbs
    #
    @test masses(model)[1] == 0.141
    @test masses(model)[2] == 0.142
    @test masses(model)[3] == 0.143
    @test masses(model).m0 == 3.09
    #
    @test spins(model)[1] == 0
    @test spins(model)[2] == 0
    @test spins(model)[3] == 0
    @test spins(model).two_h0 == 2
end


@testset "Amplitude with ThreeBodyDecay" begin
    ms = masses(model)
    σs = x2σs([0.1, 0.9], ms; k=2)
    refA_xxx1 = 0.5360686514001554 + 0.005104210009658071im
    @test amplitude(model, σs, spins(model)) ≈ refA_xxx1
    @test amplitude(model, σs, spins(model); refζs=(1, 2, 3, 1)) ≈ refA_xxx1
    @test amplitude(model, σs, spins(model); refζs=(2, 3, 2, 1)) ≈ refA_xxx1
    @test amplitude(model, σs, spins(model); refζs=(2, 2, 3, 1)) ≈ refA_xxx1
end

@testset "Intensity with ThreeBodyDecay" begin
    ms = masses(model)
    σs = x2σs([0.1, 0.9], ms; k=2)
    refI = 0.574791303947608
    @test unpolarized_intensity(model, σs) ≈ refI
    @test unpolarized_intensity(model, σs; refζs=(1, 2, 3, 1)) ≈ refI
end

@testset "Range indexing model" begin
    _model = model[2:3]
    @test _model.chains[1] == model.chains[2]
    @test length(_model) == 2
end

@testset "Scalar indexing model" begin
    _model = model[2]
    @test _model.chains[1] == model.chains[2]
    @test length(_model) == 1
end


@testset "Vertical Concat model" begin
    _model = vcat(model[2], model[2:3], model)
    @test _model isa ThreeBodyDecay
    @test length(_model) == 1 + 2 + 3
end


struct MyDecayChain{T} <: AbstractDecayChain
    x::T
end
ThreeBodyDecays.amplitude(::MyDecayChain, σs, spins) = 1.5

@testset "amplutude for AbstractDecayChain" begin
    @test amplitude(MyDecayChain(1.0), [1, 2, 3], [1, 2, 3, 4]) == 1.5
    @test amplitude(MyDecayChain(1.0), randomPoint(system(model))) == 1.5
end

model3 = ThreeBodyDecay(
    "MyResonance" .=>
        zip([3.3, 4.4], [MyDecayChain(1.0), MyDecayChain(2.2)]))

@testset "Customary model" begin
    @test model3 isa ThreeBodyDecay
    @test length(model3) == 2
end
