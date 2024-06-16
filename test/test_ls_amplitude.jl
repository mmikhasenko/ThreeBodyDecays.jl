using ThreeBodyDecays
using Test

using Parameters
using BenchmarkTools

@testset "LS amplitude scalars" begin
    tbs = ThreeBodySystem(2.0, 1.0, 1.5; m0=6.0)
    dpp = randomPoint(tbs)
    #
    dc = DecayChainLS(;
        k=1,
        Xlineshape=σ -> 1 / (4.1^2 - σ - 0.1im),
        jp="2+",
        Ps=ThreeBodyParities('+', '+', '+', P0='+'),
        tbs=tbs)
    @test sum(reim(amplitude(dc, dpp)) .≈ 0.0) == 0
end

two_js, Ps = ThreeBodySpinParities("1/2-", "0-", "0-"; jp0="1/2+")
tbs = ThreeBodySystem(2.0, 1.0, 1.5; m0=6.0, two_js)  # 1/2+ 0- 0- 1/2+
σs = x2σs([0.5, 0.3], tbs.ms; k=1)
#
bw(σ) = 1 / (4.1^2 - σ - 0.1im)
dc = DecayChainLS(; k=3, Xlineshape=bw, jp="3-", Ps, tbs)

@testset "LS amplitude half-integer spin" begin
    dpp = DalitzPlotPoint(; σs=σs, two_λs=[1, 0, 0, 1])
    @test sum(reim(amplitude(dc, dpp)) .≈ 0.0) == 0
    @test length(dc) == 1
end


# Manual construction of the decay chains
ms = ThreeBodyMasses(1.1, 2.2, 3.3; m0=7.7)
(two_js, PC), (_, PV) = map(["1/2+", "1/2-"]) do jp0
    ThreeBodySpinParities("1-", "1/2+", "0-"; jp0)
end
tbs = ThreeBodySystem(ms, two_js)

chains = let
    k = 3
    jp = jp"3/2-"
    ls_all = possible_ls_ij(jp, two_js, PC; k)
    LS_all = vcat(
        possible_ls_Rk(jp, two_js, PC; k),
        possible_ls_Rk(jp, two_js, PV; k))
    # 
    map(Iterators.product(LS_all, ls_all)) do (two_LS, two_ls)
        Xlineshape = identity
        DecayChain(; k, jp.two_j, tbs, Xlineshape,
            HRk=RecouplingLS(two_LS), Hij=RecouplingLS(two_ls))
    end
end

sortbyLl(chains) = sort(chains,
    by=x -> x.Hij.two_ls[1] * 10 + x.HRk.two_ls[1])

@testset "DecayChainsLS: all decay chains" begin
    @test sortbyLl(vec(chains)) ==
          sortbyLl(vec(
        [DecayChainsLS(; k=3, Xlineshape=identity, jp="3/2-", Ps=PC, tbs);
            DecayChainsLS(; k=3, Xlineshape=identity, jp="3/2-", Ps=PV, tbs)]))
end

σs = x2σs([0.5, 0.3], tbs.ms; k=1)
@testset "Unpol. intensity values for chains" begin
    refs = [183.0603173468474 100.20583627496791 183.0603173468474
        183.0603173468474 100.20583627496791 183.0603173468474]
    @test all(unpolarized_intensity.(chains, Ref(σs)) .≈ refs)
end



const model = let
    two_js, Ps = ThreeBodySpinParities("1-", "1/2+", "0-"; jp0="3/2+")
    tbs = ThreeBodySystem(1.1, 2.2, 3.3; m0=7.7, two_js)
    #
    models = map([
        (name="R3_3h-", k=3, two_jp=jp"3/2-"),
        (name="R1_3h-", k=1, two_jp=jp"3/2-"),
        (name="R3_1h-", k=3, two_jp=jp"1/2-"),
        (name="R1_1h-", k=1, two_jp=jp"1/2-"),
        (name="R2", k=2, two_jp=jp"3-"),]) do (; k, two_jp, name)
        #
        chains = map(possible_lsLS(two_jp, two_js, Ps; k)) do conf
            @unpack two_LS, two_ls = conf
            DecayChain(; k, two_jp.two_j, tbs, Xlineshape=identity,
                HRk=RecouplingLS(two_LS), Hij=RecouplingLS(two_ls))
        end
        ci = ones(Float64, length(chains))
        ThreeBodyDecay(name .=> zip(ci, chains))
    end
    vcat(models...)
end

@testset "Elapsed time" begin
    @test length(model) == 20
    σs = x2σs([0.5, 0.3], masses(model); k=1)
    ref_I = 17646.88022101494
    @test unpolarized_intensity(model, σs) ≈ ref_I
    @test unpolarized_intensity(model, σs; refζs=(1, 1, 1, 1)) ≈ ref_I
    @test !(unpolarized_intensity(model, σs; refζs=(2, 2, 2, 2)) ≈ ref_I)
    @test !(unpolarized_intensity(model, σs; refζs=(3, 3, 3, 3)) ≈ ref_I)
    # @btime unpolarized_intensity($model, $σs)
    evaltime = @elapsed unpolarized_intensity(model, σs)
    @info """Unpolarized_intensity for a model (3/2->1,1/2,0) with 20 chains is computed in $(round(1000*evaltime, digits=2)) ms
Compare to my usual evaluation time of about 5ms
"""
end

@testset "Full amplitude components" begin
    σs = x2σs([0.5, 0.3], masses(model); k=1)
    A_2103 = amplitude(model, σs, [2, 1, 0, 3])
    A = amplitude(model, σs)
    @test A_2103 ≈ A[end, end, end, end] ≈ 25.650736877020776
    A_m2m103 = amplitude(model, σs, [-2, -1, 0, 3])
    @test A_m2m103 ≈ A[1, 1, end, end] ≈ -14.353968898544204
end
