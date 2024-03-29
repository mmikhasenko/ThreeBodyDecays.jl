
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

function amplitude(dc::DecayChain, σs, two_λs; refζs=(1, 2, 3, 1))

    @unpack k, tbs, two_j, HRk, Hij = dc
    two_js = tbs.two_js
    # 
    i, j = ij_from_k(k)
    two_js_Hij = (two_j, two_js[i], two_js[j])
    two_js_HRk = (two_js[4], two_j, two_js[k])
    #
    ms² = tbs.ms^2
    # 
    w0 = wr(k, refζs[4], 0)
    wi = wr(k, refζs[i], i)
    wj = wr(k, refζs[j], j)
    wk = wr(k, refζs[k], k)
    # 
    cosζ0 = cosζ(w0, σs, ms²)
    cosζi = cosζ(wi, σs, ms²)
    cosζj = cosζ(wj, σs, ms²)
    cosζk = cosζ(wk, σs, ms²)
    #
    cosθ = cosθij(σs, ms²; k)
    #
    T = typeof(two_λs[1])
    T1 = one(T)
    # 
    itr_two_λs′ = itr(SVector{4,T}(two_js[1], two_js[2], two_js[3], two_js[4]))
    lineshape = dc.Xlineshape(σs[k])
    f = zero(lineshape)
    for two_λs′ in itr_two_λs′
        # two_λs′[4] = two_τ - two_λs′[k]
        # two_τ = two_λs′[4] + two_λs′[k]
        f +=
            wignerd_doublearg_sign(two_js[4], two_λs[4], two_λs′[4], cosζ0, ispositive(w0)) *
            # 
            amplitude(HRk, (two_λs′[4] + two_λs′[k], two_λs′[k]), two_js_HRk) * phase(two_js[k] - two_λs′[k]) * # particle-2 convention
            sqrt(two_j * T1 + 1) * wignerd_doublearg_sign(two_j, two_λs′[4] + two_λs′[k], two_λs′[i] - two_λs′[j], cosθ, true) *
            amplitude(Hij, (two_λs′[i], two_λs′[j]), two_js_Hij) * phase(two_js[j] - two_λs′[j]) * # particle-2 convention
            # 
            wignerd_doublearg_sign(two_js[i], two_λs′[i], two_λs[i], cosζi, ispositive(wi)) *
            wignerd_doublearg_sign(two_js[j], two_λs′[j], two_λs[j], cosζj, ispositive(wj)) *
            wignerd_doublearg_sign(two_js[k], two_λs′[k], two_λs[k], cosζk, ispositive(wk))
    end
    return f * lineshape
end
#
amplitude(dc::AbstractDecayChain, dpp) = amplitude(dc, dpp.σs, dpp.two_λs)
#
summed_over_polarization(fn, two_js) = σs -> sum(fn(σs, two_λs) for two_λs in itr(two_js))
#
itr(two_js) = Iterators.ProductIterator(Tuple([-two_j:2:two_j for two_j in two_js]))
