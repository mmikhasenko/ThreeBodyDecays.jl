using ThreeBodyDecays
using Plots
using QuadGK


const mB = 5.2
const mJpsi = 3.09
const mK = 0.489
const mρ = 0.77
const ms = ThreeBodyMasses(mJpsi, mρ, mK; m0=mB)

(two_js, PC), (_, PV) = (ThreeBodySpinParities("1-", "1-", "0-"; jp0) for jp0 = ("0+", "0-"))
tbs = ThreeBodySystem(ms, two_js)

BW(m, Γ) = σ -> 1 / (m^2 - σ - 1im * m * Γ)
const mK1, ΓK1 = 1.453, 0.3

chains_PC = DecayChainsLS(; k=1, jp="1+", Xlineshape=BW(mK1, ΓK1), Ps=PC, tbs)
chains_PV = DecayChainsLS(; k=1, jp="1+", Xlineshape=BW(mK1, ΓK1), Ps=PV, tbs)

let
    chain = chains_PV[4]
    plot(masses(chain), xlab="m²(J/ψρ)", ylab="m²(ρK)", title="Dalitz plot for K1(1270) decay",
        Base.Fix1(unpolarized_intensity, chain); iσx=3, iσy=1)
end

let
    chain = chains_PV[4]
    plot(3.8, 4.0) do e1
        I = Base.Fix1(unpolarized_intensity, chain)
        e1 * quadgk(projection_integrand(I, masses(chain), e1^2; k=3), 0, 1)[1]
    end
    X = BW(3.872, 10e-3)
    plot!(e1 -> 15 * abs2(X(e1^2) / X(3.872^2)), 3.7, 4.5, l=(1, :red))
end