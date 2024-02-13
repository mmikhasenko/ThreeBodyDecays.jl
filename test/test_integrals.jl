using ThreeBodyDecays
using Test
using Cuba


@testset "phase space mapping" begin
    ms = ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09)
    f = phase_space_integrand(σs -> 1, ms; k=2)
    @test f([0.1, 0.9]) ≈ 7.399872133477282
end

@testset "three-body phase-space integral value" begin
    ms = ThreeBodyMasses(0.938, 0.493, 0.0; m0 = 5.62)
    # 
    integrand = phase_space_integrand(σs -> 1, ms; k=1)
    value_1 = cuhre((x, f)->(f[1] = integrand(x)), 2, 1)[1][1]
    integrand = phase_space_integrand(σs -> 1, ms; k=2)
    value_2 = cuhre((x, f)->(f[1] = integrand(x)), 2, 1)[1][1]
    integrand = phase_space_integrand(σs -> 1, ms; k=3)
    value_3 = cuhre((x, f)->(f[1] = integrand(x)), 2, 1)[1][1]
    # 
    ref_value = 11.516096823773797

    @test isapprox(value_1, ref_value; rtol=1e-6)
    @test isapprox(value_2, ref_value; rtol=1e-6)
    @test isapprox(value_3, ref_value; rtol=1e-6)
end