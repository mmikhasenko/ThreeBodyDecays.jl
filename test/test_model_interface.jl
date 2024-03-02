using Test
using ThreeBodyDecays
using Parameters



# II) possible_LS with spin and parity
ms = ThreeBodyMasses(1.1, 2.2, 3.3; m0=7.7)
(two_js, PC), (_, PV) =
    [ThreeBodySpinParities("1-", "1/2+", "0-"; jp0)
     for jp0 in ["1/2+", "1/2-"]]
tbs = ThreeBodySystem(masses, two_js)

chains = let
    i, j, k = ijk(3)
    two_ji, two_jj, two_jk, two_j0 = collect(two_js)[[i, j, k, 4]]
    mi, mj, mk, m0 = collect(ms)[[i, j, k, 4]]
    #
    jp = jp"3/2-"
    @unpack two_j = jp
    ls_all = possible_ls(jp, two_js, PC; k)
    LS_all = vcat(
        possible_LS(jp, two_js, PC; k),
        possible_LS(jp, two_js, PV; k))
    # 
    map(Iterators.product(LS_all, ls_all)) do (two_LS, two_ls)
        # vertex HRk
        (; two_LS, two_ls)
        HRk = RecouplingLS(two_LS)
        Hij = RecouplingLS(two_ls)
        L = div(two_LS[1], 2)
        l = div(two_ls[1], 2)
        # # 
        # Xlineshape =
        #     BW(m, Γ, ma, mb) *
        #     BlattWeisskopf{L}(dX)(σ -> breakup(m0, sqrt(σ), mk)) *
        #     BlattWeisskopf{l}(dR)(σ -> breakup(sqrt(σ), mi, mj))
        # #
        DecayChain(; k, two_j, Xlineshape=identity, HRk, Hij, tbs)
        # 
    end
end

# LS + parity basis
two_λs = possible_rehelicities()
map(LS_all, two_λs) do two_LS, (two_λa, two_λb)
    HRk = RecouplingLS(two_LS)
    Hij = ParityRecoupling(two_λa, two_λb, phase_factor)
    DecayChain(; k=1, two_j, Xlineshape, HRk, Hij)
end


# II) Reimplementation of Recoupling amplitude call
H = RecouplingLS(LS)
two_js = ThreeBodySpins(2, 1, 0; two_h0=1)
two_λs = two_js
@test "calling H_ls" begin
    @test amplitude(H, two_λs, two_js; k=1) = sqrt(two_j + 1)
end


let
