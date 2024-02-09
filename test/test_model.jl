using Test
using ThreeBodyDecays

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

# @testset "amplitude with ThreeBodyDecay" begin
let
    dpp = randomPoint(first(model.chains).tbs)
    amplitude(model, dpp) == 0.1415775268065487 + 0.0012649796556426023im
end