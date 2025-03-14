using Test
using ThreeBodyDecays

@testset "jp_str macro" begin
    @test jp"3+" == str2jp("3+")
    @test jp"3^+" == jp"3+"
    @test jp"1/2+" == SpinParity(1, '+')
    @test jp"3/2-" == SpinParity(3, '-')
    @test jp"3+" == SpinParity(6, '+')
    @test SpinParity(6, '+').two_j == 6
    @test SpinParity(6, '+').p == '+'
end


@testset "jp ⊗ jp" begin
    S_wave = SpinParity("1/2+") ⊗ SpinParity("1-")
    P_wave = [sw ⊗ SpinParity("1-") for sw in S_wave]
    D_wave = [sw ⊗ SpinParity("2+") for sw in S_wave]
    #
    @test length(S_wave) == 2
    @test length(vcat(P_wave...)) == 5
    @test length(Set(vcat(P_wave...))) == 3
end

let
    @test length(possible_ls(jp"3/2-", jp"3-"; jp = jp"1/2+")) == 4
    @test possible_ls("3/2-", "3-"; jp = "1/2+") ==
          possible_ls(jp"3/2-", jp"3-"; jp = jp"1/2+")
    #
    two_js, Ps = ThreeBodySpinParities("1+", "1/2+", "0-"; jp0 = "1/2-")
    lsLSv = possible_lsLS(jp"1/2+", two_js, Ps; k = 1)
    @test size(lsLSv) == (1, 2)
    #
    @test length(possible_ls_ij(jp"1+", two_js, Ps; k = 1)) == 1
    @test length(possible_ls_Rk(jp"1+", two_js, Ps; k = 1)) == 2
    lsLSv = possible_lsLS(jp"1+", two_js, Ps; k = 1)
    @test size(lsLSv) == (1, 2)
end

@testset "x2, times 2" begin
    @test x2(1 // 2) == 1
    @test x2(1.5) == 3
    @test x2("5/2") == 5
    @test x2(3) == 6
    @test x2((2, 3 // 2)) == (4, 3)
end

@testset "d2, div 2" begin
    @test d2(3) == "3/2"
    @test d2(4) == "2"

    @test d2((2, 1)) == ("1", "1/2")
    @test d2([3, 2]) == ["3/2", "1"]
end

@testset "letterL" begin
    @test letterL(0) == 'S'
    @test letterL(3) == 'F'
    @test letterL("2") == 'D'
    @test letterL(4 |> d2) == 'D'
end


(two_js, pc), (_, pv) = map(["1/2+", "1/2-"]) do jp0
    ThreeBodySpinParities("1/2+", "0-", "0-"; jp0)
end;

@testset "Completing l, s, L, S" begin
    qn = (jp = jp"1/2-", k = 2)
    updated_qn = complete_l_s_L_S(qn.jp, two_js, pc, qn; qn.k)
    @test updated_qn.l == "0"
    @test updated_qn.s == "1/2"
    @test updated_qn.L == "0"
    @test updated_qn.S == "1/2"
    #
    qn = (jp = jp"2+", k = 1)
    @test_throws ErrorException complete_l_s_L_S(qn.jp, two_js, pc, qn; qn.k)
    qn = (jp = jp"2+", S = 3 / 2, k = 1)
    @test complete_l_s_L_S(qn.jp, two_js, pc, qn; qn.k).L == "2"
    #
    qn = (jp = jp"0+", k = 1)
    @test_throws ErrorException complete_l_s_L_S(qn.jp, two_js, [pc, pv], qn; qn.k)
    #
    qn = (jp = jp"0+", k = 1, L = 0)
    @test complete_l_s_L_S(qn.jp, two_js, [pc, pv], qn; qn.k).S == "1/2"
end
