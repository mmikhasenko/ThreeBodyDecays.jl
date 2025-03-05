using ThreeBodyDecays
using QuadGK
using Test
using Cuba

@testset "phase space mapping" begin
    ms = ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09)
    f = phase_space_integrand(σs -> 1, ms; k = 2)
    @test f([0.1, 0.9]) ≈ 7.399872133477282
end

@testset "three-body phase-space integral value" begin
    ms = ThreeBodyMasses(0.938, 0.493, 0.0; m0 = 5.62)
    #
    integrand = phase_space_integrand(σs -> 1, ms; k = 1)
    value_1 = cuhre((x, f) -> (f[1] = integrand(x)), 2, 1)[1][1]
    integrand = phase_space_integrand(σs -> 1, ms; k = 2)
    value_2 = cuhre((x, f) -> (f[1] = integrand(x)), 2, 1)[1][1]
    integrand = phase_space_integrand(σs -> 1, ms; k = 3)
    value_3 = cuhre((x, f) -> (f[1] = integrand(x)), 2, 1)[1][1]
    #
    ref_value = 11.516096823773797

    @test isapprox(value_1, ref_value; rtol = 1e-6)
    @test isapprox(value_2, ref_value; rtol = 1e-6)
    @test isapprox(value_3, ref_value; rtol = 1e-6)
end

@testset "projection integral" begin
    ms = ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09)
    f = projection_integrand(σs -> 1, ms, 0.7^2; k = 1)
    @test quadgk(f, 0, 1)[1] ≈ 8.253210733598506
    #
    @test quadgk(projection_integrand(σs -> 1, ms, 0.7^2; k = 2), 0, 1)[1] ≈
          8.258640762296254
    @test quadgk(projection_integrand(σs -> 1, ms, 0.7^2; k = 3), 0, 1)[1] ≈
          8.26409477583501
end

@testset "cosθ projection integral" begin
    ms = ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09)
    I = σs -> σs[1] * σs[2] - σs[3]  # Unit amplitude for testing

    # Test integration over cosθ matches phase space integral
    Nb = 500
    # Integration over cosθ
    cosθ_integral = [
        sum(range(-1, 1, Nb)) do z
            integrand = project_cosθij_intergand(I, ms, z; k)
            quadgk(integrand, 0, 1)[1] * 2 / Nb
        end for k = 1:3
    ]

    # Integration over σk
    σk_integral = [
        sum(range(lims(ms; k)..., Nb)) do σk
            integrand = projection_integrand(I, ms, σk; k)
            quadgk(integrand, 0, 1)[1] * diff(collect(lims(ms; k)))[1] / Nb
        end for k = 1:3
    ]

    # Both methods should give approximately the same result
    @test isapprox(cosθ_integral[1], σk_integral[1], rtol = 1e-2)
    @test isapprox(cosθ_integral[2], σk_integral[2], rtol = 1e-2)
    @test isapprox(cosθ_integral[3], σk_integral[3], rtol = 1e-2)
end
