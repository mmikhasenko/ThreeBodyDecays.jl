using ThreeBodyDecays
using Test
using Random

ms = ThreeBodyMasses(1.1, 1.3, 1.5; m0=6.0)

Random.seed!(1234)
N = 100_000
randsample = eachslice(rand(2, N); dims=2)
presample = y2σs.(randsample, Ref(ms))
physical_sample = filter(Base.Fix2(inphrange, ms), presample)
@testset "Generated flat sample" begin
    @test length(physical_sample) / N == 0.70207
end
# histogram2d(phsample, bins=50)

Nx, Ny = 200, 100
grid2d_xv = Iterators.product(range(0, 1, 200), range(0, 1, 100))
grid2d_σsv = y2σs.(grid2d_xv, Ref(ms))
isphys_grid = map(Base.Fix2(inphrange, ms), grid2d_σsv)
@testset "Dalitz grid values" begin
    @test sum(isphys_grid) / (Nx * Ny) == 0.6901
end


b = border(ms)
borderKibblemax = maximum(abs.(Kibble.(b, Ref(ms^2))))
sampleKibblemax = maximum(abs.(Kibble.(physical_sample, Ref(ms^2))))

@testset "Border of Dalitz" begin
    @test borderKibblemax / sampleKibblemax < 1e-10
end
