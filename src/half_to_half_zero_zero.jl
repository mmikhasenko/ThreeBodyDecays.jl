
export RkΛsτ
# couplings
# specific for 1/2 -> (1/2 0 0) final state
#
R1Λsτ(two_λp,two_Λ,two_τ,chains,W) =
    sum(v*HelicityRecoupling_doublearg([(two_j(chain),two_τ,1,two_λp)=>(two_J(chain),two_L(chain),two_S(chain))]) for (chain,v) in zip(chains,W))
R2Λsτ(two_λp,two_Λ,two_τ,chains,W) =
    sum(v*HelicityRecoupling_doublearg([(two_j(chain),two_τ,0,0)=>(two_J(chain),two_L(chain),two_S(chain)),
                                        (0,0,1,two_λp)=>(two_j(chain),two_l(chain),two_s(chain))]) for (chain,v) in zip(chains,W))
R3Λsτ(two_λp,two_Λ,two_τ,chains,W) =
    sum(v*HelicityRecoupling_doublearg([(two_j(chain),two_τ,0,0)=>(two_J(chain),two_L(chain),two_S(chain)),
                                        (1,two_λp,0,0)=>(two_j(chain),two_l(chain),two_s(chain))]) for (chain,v) in zip(chains,W))
function RkΛsτ(k,two_λp,two_Λ,two_τ,chains,W)
    (k==1) && return R1Λsτ(two_λp,two_Λ,two_τ,chains,W)
    (k==2) && return R2Λsτ(two_λp,two_Λ,two_τ,chains,W)
    (k==3) && return R3Λsτ(two_λp,two_Λ,two_τ,chains,W)
    return zero(W[1])
end

# angular functions
# specific for 1/2 -> (1/2 0 0) final state
# with reference to p <--o--> K*
#
function Z1Λsτ(two_λ,two_Λ,two_s,two_τ,σ3,σ1,tbs)
    σ2 = gσ2(σ3,σ1,tbs)
    return two_Λ!=(two_τ-two_λ) ? 0.0 :
        sqrt(2)*sqrt(two_s+1)*wignerd_doublearg(two_s,two_τ,0,cosθ23(σ2,σ3,tbs))
end
#
function Z2Λsτ(two_λ,two_Λ,two_s,two_τ,σ3,σ1,tbs)
    σ2 = gσ2(σ3,σ1,tbs)
    return sqrt(2) * wignerd_doublearg(1,two_Λ,two_τ,cosθhat31(σ3,σ1,tbs)) *
        ((two_τ+two_λ) % 4 == 2 ? -1 : 1)*
        sqrt(two_s+1)*wignerd_doublearg(two_s,two_τ,-two_λ,cosθ31(σ3,σ1,tbs))
end
#
function Z3Λsτ(two_λ,two_Λ,two_s,two_τ,σ3,σ1,tbs)
    σ2 = gσ2(σ3,σ1,tbs)
    return sqrt(2) * wignerd_doublearg(1,two_Λ,two_τ,cosθhat12(σ1,σ2,tbs)) *
        sqrt(two_s+1)*wignerd_doublearg(two_s,two_τ,-two_λ,cosθ12(σ1,σ2,tbs))
end
#
function ZkΛsτ(k,two_λ,two_Λ,two_s,two_τ,σ3,σ1,tbs)
    (k==1) && return Z1Λsτ(two_λ,two_Λ,two_s,two_τ,σ3,σ1,tbs)
    (k==2) && return Z2Λsτ(two_λ,two_Λ,two_s,two_τ,σ3,σ1,tbs)
    (k==3) && return Z3Λsτ(two_λ,two_Λ,two_s,two_τ,σ3,σ1,tbs)
    return 0.0
end

```
Full amplitude
 - specific for 1/2 -> (1/2 0 0) final state
 - with reference to p <--o--> K*
```
function amp_b2bzz(two_λ,two_Λ,σ3,σ1,CS,tbs,Cs)
    σ2 = gσ2(σ3,σ1,tbs)
    # Wigner rotations
    D1 = [two_λ==two_λp ? 1.0 : 0.0 for two_λp=-1:2:1]
    D2 = [((two_λp-two_λ) % 4 == 2 ? -1 : 1) *
        wignerd_doublearg(1,two_λp,two_λ,cosθtilde1to2_for1(σ1,σ2,tbs))
            for two_λp=-1:2:1]
    D3 = [wignerd_doublearg(1,two_λp,two_λ,cosθtilde3to1_for1(σ3,σ1,tbs))
            for two_λp=-1:2:1]
    #
    σs = (σ1,σ2,σ3)
    Ds = (D1,D2,D3)
    #
    val = zero(Cs[1][1,1]+0.0im)
    for (cs,W) in zip(CS, Cs) # coupling scheme and couplings
        (k,ξ,chains) = cs
        length(chains) == 0 && continue
        decay_amp = zero(Cs[1][1,1]+0.0im)
        for two_τ in -two_j(chains[1]):2:two_j(chains[1]), two_λp = -1:2:1
            Wc = [v1+v2*1im for (v1,v2) in zip(W[:,1],W[:,1])]
            decay_amp += RkΛsτ(k,two_λp,two_Λ,two_τ,chains, Wc) *
                    ZkΛsτ(k,two_λp,two_Λ,two_j(chains[1]),two_τ,σ3,σ1,tbs) *
                    Ds[k][div(two_λp+3,2)]
        end
        val += decay_amp * amp(σs[k],ξ)
    end
    return val;
end

```
Squared amplitude summed over polarization of initial and final state
 - specific for 1/2 -> (1/2 0 0) final state
 - with reference to p <--o--> K*
```
function rate_b2bzz(σ3,σ1,CS,tbs,Cs)
    return sum(
        abs2(amp_b2bzz(two_λ,two_Λ,σ3,σ1,CS,tbs,Cs))
            for two_λ in -1:2:1
                for two_Λ in -1:2:1)
end

```
A(C1)A*(C1) summed over polarization of initial and final state
 - specific for 1/2 -> (1/2 0 0) final state
 - with reference to p <--o--> K*
```
function rateCC_b2bzz(σ3,σ1,CS,tbs,Cs1,Cs2)
    return sum(
        amp_b2bzz(two_λ,two_Λ,σ3,σ1,CS,tbs,Cs1) *
            conj(amp_b2bzz(two_λ,two_Λ,σ3,σ1,CS,tbs,Cs2))
            for two_λ in -1:2:1
                for two_Λ in -1:2:1)
end

```
A(Λ1)A*(Λ2) summed over the final-state polarization
 - specific for 1/2 -> (1/2 0 0) final state
 - with reference to p <--o--> K*
```
function rateΛΛ_b2bzz(two_Λ1,two_Λ2,σ3,σ1,CS,tbs)
    return sum(
        amp_b2bzz(two_λ,two_Λ1,σ3,σ1,CS,tbs,Cs) *
            conj(amp_b2bzz(two_λ,two_Λ2,σ3,σ1,CS,tbs,Cs))
            for two_λ in -1:2:1)
end

function polarization_vector(M)
    P0 = real((M[1,1]+M[2,2]) / 2.0)
    Pz = real((M[1,1]-M[2,2]) / 2.0 / P0)
    Px =  real(M[1,2]) / P0
    Py = -imag(M[1,2]) / P0
    return (Px,Py,Pz,P0)
end
```
Gives a vector of the the polarization sensetivity
 - specific for 1/2 -> (1/2 0 0) final state
 - with reference to p <--o--> K*
```
function polSens_b2bzz(σ3,σ1,CS,tbs,Cs)
    M = [rateΛΛ_b2bzz(two_Λ1,two_Λ2,σ3,σ1,CS,tbs,Cs) for two_Λ1=-1:2:1, two_Λ2=-1:2:1]
    P = polarization_vector(M)
    return P
end

```
Gives a vector of the the polarization sensetivity
 - specific for 1/2 -> (1/2 0 0) final state
 - with reference to p <--o--> K*
```
function polSens_b2bzz(CS,tbs, Cs; gridN::Int = 100)
    σ1v = LinRange(tbs.mthsq[1],tbs.sthsq[1],gridN)
    σ3v = LinRange(tbs.mthsq[3],tbs.sthsq[3],gridN)
    #
    M11 = sum(Kibble31(σ3,σ1,tbs) > 0.0 ? 0.0 : rateΛΛ_b2bzz(-1,-1,σ3,σ1,CS,tbs,Cs) for σ3 in σ3v, σ1 in σ1v)
    M12 = sum(Kibble31(σ3,σ1,tbs) > 0.0 ? 0.0 : rateΛΛ_b2bzz( 1,-1,σ3,σ1,CS,tbs,Cs) for σ3 in σ3v, σ1 in σ1v)
    M22 = sum(Kibble31(σ3,σ1,tbs) > 0.0 ? 0.0 : rateΛΛ_b2bzz( 1, 1,σ3,σ1,CS,tbs,Cs) for σ3 in σ3v, σ1 in σ1v)
    #
    P = polarization_vector([M11 M12; conj(M12) M22])
    return P
end


# #############################################################################
#
# ```
# Get complex couplings of the decay chain
#
# Comment: for the b2bzz, the decay chain set up is transfered to the amplitude
# in a specific form
#     Vector{Tuple{Int,Lineshape,Vector{Pair{twochain,Number}}}}
# ```
# function get_couplings(Cls)
#     couplins = vcat(map(x->map(y->y[2],x[3]),Cls)...)
#     return couplins
# end
#
# ```
# Set the complex couplings to the decay chain
#
# Comment: for the b2bzz, the decay chain set up is transfered to the amplitude
# in a specific form
#     Vector{Tuple{Int,Lineshape,Vector{Pair{twochain,Number}}}}
# ```
# function set_couplings(Cls, p)
#     new_Cls = Vector{Tuple{Int,Lineshape,Vector{Pair{twochain,Number}}}}(undef,0);
#     ip=1;
#     for (k,f,couplings) in Cls
#         new_couplings = Vector{Pair{twochain,Number}}(undef,0)
#         for (chain,v) in couplings
#             push!(new_couplings,chain=>p[ip]); ip+=1
#         end
#         push!(new_Cls,(k,f,new_couplings))
#     end
#     return new_Cls
# end