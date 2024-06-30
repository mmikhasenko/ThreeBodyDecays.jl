using ThreeBodyDecays
using Test

two_js = ThreeBodySpins(1, 0, 2; two_h0=3)

@testset "Three Body Masses structure" begin
    #
    @test two_js.two_h1 == two_js[1] == 1
    @test two_js.two_h2 == two_js[2] == 0
    @test two_js.two_h3 == two_js[3] == 2
    @test two_js.two_h0 == two_js[4] == 3
    #
    @test_throws ErrorException ThreeBodySpins(0, 1, 1; two_h0=1)
    @test_throws ErrorException ThreeBodySpins(0, 1, 1)
    @test_throws BoundsError two_js[5]
end

@testset "Spin with floats and strings" begin
    #
    @test ThreeBodySpins("1", "1/2", "0"; h0="1/2") ==
          ThreeBodySpins(2, 1, 0; two_h0=1)
    #
    @test ThreeBodySpins(1, 1 / 2, 0; h0=1 / 2) ==
          ThreeBodySpins(2, 1, 0; two_h0=1)
    #
    @test_throws ErrorException ThreeBodySpins("0", "1/2", "1/2"; h0="1/2")
end


two_js, parities = ThreeBodySpinParities("1-", "1/2+", "0-"; jp0="1/2+")
@testset "Spin-Parity call" begin
    @test two_js == ThreeBodySpins(2, 1, 0; two_h0=1)
    @test parities == ThreeBodyParities('-', '+', '-'; P0='+')
end

(two_js, PC), (_, PV) = [ThreeBodySpinParities("1-", "1/2+", "0-"; jp0) for jp0 in ("1/2+", "1/2-")]
@testset "Spin-Parity call with ±" begin
    @test two_js == ThreeBodySpins(2, 1, 0; two_h0=1)
    @test [PC, PV] == [
        ThreeBodyParities('-', '+', '-'; P0='+'),
        ThreeBodyParities('-', '+', '-'; P0='-')]
end

@testset "operations and integrate" begin
    @test iseven(sum(two_js))
    @test div(sum(Tuple(two_js) .|> x2), 2) == sum(two_js)
end
