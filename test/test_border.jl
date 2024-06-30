using ThreeBodyDecays
using Test

ms = let
    mpi = 0.13957
    meta = 0.547862
    mp = 0.938
    #
    ThreeBodyMasses(mpi, meta, mp; m0 = 19.0)
end

σs_b = border(ms, Nx = 10)

Kibble(σs_b[1], ms^2)
Kibble.(σs_b, Ref(ms^2))

@testset "Border is an array of MandelstamTuple, vanishes Kibble" begin
    @test σs_b isa Vector{T} where {T <: ThreeBodyDecays.MandelstamTuple{Float64}}
    @test all(sum.(σs_b) .≈ sum(ms^2))
    @test all(Kibble.(σs_b, Ref(ms^2)) .< 1e-3)
end

# Check if the border function works for the case of unbalanced masses.

const m_small = 0.1
const m_large = 10.0

@testset "Border for unbalanced masses" begin
    for ms in [
        ThreeBodyMasses(m_small, m_small, m_large; m0 = 11.0),
        ThreeBodyMasses(m_small, m_large, m_small; m0 = 11.0),
        ThreeBodyMasses(m_large, m_small, m_small; m0 = 11.0),
        #
        ThreeBodyMasses(m_small, m_large, m_large; m0 = 21.0),
        ThreeBodyMasses(m_large, m_small, m_large; m0 = 21.0),
        ThreeBodyMasses(m_large, m_large, m_small; m0 = 21.0),
    ]
        #
        σs_b = border(ms)
        extremaKibble = Kibble.(σs_b, Ref(ms^2)) |> extrema
        @test all(abs.(extremaKibble) .< 1e-5)
    end
end
