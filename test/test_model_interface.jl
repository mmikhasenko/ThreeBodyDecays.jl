using Test
using ThreeBodyDecays
using Parameters


# II) possible_LS with spin and parity
ms = ThreeBodyMasses(1.1, 2.2, 3.3; m0=7.7)
(two_js, PC), (_, PV) = map(["1/2+", "1/2-"]) do jp0
    ThreeBodySpinParities("1-", "1/2+", "0-"; jp0)
end
tbs = ThreeBodySystem(masses, two_js)

chains = let
    k = 3
    jp = jp"3/2-"
    ls_all = possible_ls_ij(jp, two_js, PC; k)
    LS_all = vcat(
        possible_ls_Rk(jp, two_js, PC; k),
        possible_ls_Rk(jp, two_js, PV; k))
    # 
    map(Iterators.product(LS_all, ls_all)) do (two_LS, two_ls)
        L = div(two_LS[1], 2)
        l = div(two_ls[1], 2)
        Xlineshape =
            BW(m, Γ, ma, mb) *
            BlattWeisskopf{l}(dR)(breakup_ij(ms; k)) *
            BlattWeisskopf{L}(dX)(breakup_Rk(ms; k))
        DecayChain(; k, jp.two_j, tbs, Xlineshape,
            HRk=RecouplingLS(two_LS), Hij=RecouplingLS(two_ls))
    end
end

breakup_Rk(ms; k) = σ -> breakup(ms[4], sqrt(σ), ms[k])
breakup_ij(ms; k) = ((i, j) = ij_from_k(k); σ -> breakup(sqrt(σ), ms[i], ms[j]))

# LS + parity basis
two_λs = possible_rehelicities()
map(LS_all, two_λs) do two_LS, (two_λa, two_λb)
    HRk = RecouplingLS(two_LS)
    Hij = ParityRecoupling(two_λa, two_λb, phase_factor)
    DecayChain(; k=1, two_j, Xlineshape, HRk, Hij)
end

