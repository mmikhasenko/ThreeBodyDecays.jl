
#        _|                                                        _|                  _|
#    _|_|_|    _|_|      _|_|_|    _|_|_|  _|    _|        _|_|_|  _|_|_|      _|_|_|      _|_|_|
#  _|    _|  _|_|_|_|  _|        _|    _|  _|    _|      _|        _|    _|  _|    _|  _|  _|    _|
#  _|    _|  _|        _|        _|    _|  _|    _|      _|        _|    _|  _|    _|  _|  _|    _|
#    _|_|_|    _|_|_|    _|_|_|    _|_|_|    _|_|_|        _|_|_|  _|    _|    _|_|_|  _|  _|    _|
#                                                _|
#                                            _|_|

abstract type Recoupling end
@with_kw struct NoRecoupling <: Recoupling
    two_λa::Int
    two_λb::Int
end

amplitude(cs::NoRecoupling, (two_λa, two_λb), (two_j, two_ja, two_jb)) =
    (cs.two_λa == two_λa) *
    (cs.two_λb == two_λb)

@with_kw struct ParityRecoupling <: Recoupling
    two_λa::Int
    two_λb::Int
    ηηηphaseisplus::Bool
end

ParityRecoupling(two_λa::Int, two_λb::Int, ηηηphasesign::Char) = ParityRecoupling(two_λa, two_λb, ηηηphasesign == '+')
function ParityRecoupling(two_λa::Int, two_λb::Int, (jp, (jp1, jp2))::TwoBodyTopologySpinParity)
    ηηη = jp1.p ⊗ jp2.p ⊗ jp.p
    ηηηphase = (2 * (ηηη == '+') - 1) * x"-1"^(div(jp.two_j - jp1.two_j - jp2.two_j, 2))
    return ParityRecoupling(two_λa, two_λb, ηηηphase == 1)
end

function amplitude(cs::ParityRecoupling, (two_λa, two_λb), (two_j, two_ja, two_jb))
    (cs.two_λa == two_λa) * (cs.two_λb == two_λb) && return 1
    (cs.two_λa == -two_λa) * (cs.two_λb == -two_λb) && return 2 * cs.ηηηphaseisplus - 1
    return 0
end

@with_kw struct RecouplingLS <: Recoupling
    two_ls::Tuple{Int,Int}
end

amplitude(cs::RecouplingLS, (two_λa, two_λb), (two_j, two_ja, two_jb)) =
    jls_coupling(two_ja, two_λa, two_jb, two_λb, two_j, cs.two_ls[1], cs.two_ls[2])

abstract type AbstractDecayChain end

@with_kw struct DecayChain{X,T} <: AbstractDecayChain
    k::Int
    #
    two_j::Int # isobar spin
    #
    Xlineshape::X # lineshape
    #
    HRk::Recoupling
    Hij::Recoupling
    #
    tbs::T # the structure with masses and spins
end

iterate(x::AbstractDecayChain) = (x, nothing)
length(x::AbstractDecayChain) = 1
spins(d::DecayChain) = d.tbs.two_js
masses(d::DecayChain) = d.tbs.ms

"""
	DecayChainLS(;
		k, # chain is specified by the spectator index k
		Xlineshape, # lambda function for lineshape
		jp, # the spin-parity of the resonance, e.g. jp"1/2-"
		Ps, # need parities, e.g. Ps=ThreeBodyParities('+','+','+'; P0='+')
		tbs) # give three-body-system structure

	Returns the decay chain with the smallest LS, ls
"""
function DecayChainLS(;
    k, Xlineshape,
    tbs=error("give three-body-system structure"),
    jp=error("give spin-parity quantum numbers of the resonance"),
    Ps=error("need parities"))
    # 
    _jp = SpinParity(jp)
    two_lsLS = vcat(possible_lsLS(_jp, tbs.two_js, Ps; k)...)
    length(two_lsLS) == 0 && error("there are no possible LS couplings")
    # 
    two_lsLS_sorted = sort(two_lsLS, by=x -> x.two_LS[1])
    @unpack two_ls, two_LS = two_lsLS_sorted[1]
    #
    i, j = ij_from_k(k)
    return DecayChain(; k, Xlineshape, tbs, _jp.two_j,
        Hij=RecouplingLS(two_ls),
        HRk=RecouplingLS(two_LS))
end

"""
	DecayChainsLS(;
		k, # chain is specified by the spectator index k
		Xlineshape, # lambda function for lineshape
		jp, # the spin-parity of the resonance, e.g. jp"1/2-"
		Ps, # need parities, e.g. Ps=ThreeBodyParities('+','+','+'; P0='+')
		tbs) # give three-body-system structure

	Returns an array of the decay chains with all possible couplings
"""
function DecayChainsLS(;
    k,
    Xlineshape,
    jp=error("give spin-parity quantum numbers of the resonance"),
    Ps=error("need parities"),
    tbs=error("give three-body-system structure, tbs=..."))
    # 
    _jp = SpinParity(jp)
    LSlsv = possible_lsLS(_jp, tbs.two_js, Ps; k)
    return [DecayChain(;
        k, Xlineshape, tbs, _jp.two_j,
        Hij=RecouplingLS(two_ls),
        HRk=RecouplingLS(two_LS))
            for (two_ls, two_LS) in LSlsv]
end


wignerd_doublearg_sign(two_j, two_λ1, two_λ2, cosθ, ispositive) =
    (ispositive ? 1 : x"-1"^div(two_λ1 - two_λ2, 2)) *
    wignerd_doublearg(two_j, two_λ1, two_λ2, cosθ)

wignerd_doublearg_sign(two_j, cosθ, ispositive) =
    wignerd_doublearg_sign.(two_j, -two_j:2:two_j, transpose(-two_j:2:two_j), cosθ, ispositive)


function aligned_amplitude(dc::DecayChain, σs)
    @unpack k, tbs, two_j, HRk, Hij = dc
    i, j = ij_from_k(k)
    # 
    ms² = tbs.ms^2
    cosθ = cosθij(σs, ms²; k)
    d_θ = wignerd_doublearg_sign(two_j, cosθ, true)
    # 
    two_js = tbs.two_js
    two_js_Hij = (two_j, two_js[i], two_js[j])
    two_js_HRk = (two_js[4], two_j, two_js[k])
    # 
    VRk = [amplitude(HRk, (two_m1, two_m2), two_js_HRk) *
           phase(two_js[k] - two_m2) # particle-2 convention
           for two_m1 in -two_j:2:two_j, two_m2 in -two_js[k]:2:two_js[k]]
    # 
    Vij = [amplitude(Hij, (two_m1, two_m2), two_js_Hij) *
           phase(two_js[j] - two_m2) # particle-2 convention
           for two_m1 in -two_js[i]:2:two_js[i], two_m2 in -two_js[j]:2:two_js[j]]
    #
    Δ_zk = div(two_j - two_js[4] - two_js[k], 2) - 1
    Δ_ij = div(two_j - two_js[i] + two_js[j], 2) + 1
    #
    lineshape = dc.Xlineshape(σs[k])
    F0 = zeros(typeof(lineshape), Tuple(two_js) .+ 1)
    F = permutedims(F0, (i, j, k, 4))
    # 
    for itr in CartesianIndices(F)
        _zk = itr[4] + itr[3] + Δ_zk
        _ij = itr[1] - itr[2] + Δ_ij
        # 
        !(1 <= _zk <= two_j + 1) && continue
        !(1 <= _ij <= two_j + 1) && continue
        # 
        F[itr] = VRk[_zk, itr[3]] * d_θ[_zk, _ij] * Vij[itr[1], itr[2]]
    end
    # 
    one_T = one(typeof(two_js[1]))
    F .*= sqrt(two_j * one_T + 1) * # normalization
          lineshape # same for all amplutudes in the chain
    # 
    permutedims!(F0, F, invperm((i, j, k, 4)))
    return F0
end

function amplitude(dc::DecayChain, σs, two_λs; refζs=(1, 2, 3, 1))
    @unpack k, tbs, two_j = dc
    ms² = tbs.ms^2
    two_js = tbs.two_js

    F0 = aligned_amplitude(dc, σs)

    # alignment rotations
    d_ζs = map(enumerate(zip(two_js, refζs))) do (l, (_two_j, _refζ))
        _w = wr(k, _refζ, mod(l, 4))
        _cosζ = cosζ(_w, σs, ms²)
        wignerd_doublearg_sign(_two_j, _cosζ, ispositive(_w))
    end
    #
    ind = map(zip(two_js, two_λs)) do (_two_j, _two_λ)
        div(_two_j + _two_λ, 2) + 1
    end
    # 
    f = sum(
        d_ζs[4][ind[4], itr[4]] *
        F0[itr] *
        d_ζs[1][itr[1], ind[1]] *
        d_ζs[2][itr[2], ind[2]] *
        d_ζs[3][itr[3], ind[3]]
        # 
        for itr in CartesianIndices(F0)
    )
    return f
end

#
amplitude(dc::AbstractDecayChain, dpp; kw...) = amplitude(dc, dpp.σs, dpp.two_λs; kw...)
#
summed_over_polarization(fn, two_js) = σs -> sum(fn(σs, two_λs) for two_λs in itr(two_js))
#
itr(two_js) = Iterators.ProductIterator(Tuple([-two_j:2:two_j for two_j in two_js]))
