#                                  _|
#    _|_|_|  _|    _|    _|_|_|  _|_|_|_|    _|_|    _|_|_|  _|_|
#  _|_|      _|    _|  _|_|        _|      _|_|_|_|  _|    _|    _|
#      _|_|  _|    _|      _|_|    _|      _|        _|    _|    _|
#  _|_|_|      _|_|_|  _|_|_|        _|_|    _|_|_|  _|    _|    _|
#                  _|
#              _|_|

const MassTuple{T} = NamedTuple{(:m1, :m2, :m3, :m0),NTuple{4,T}}
function ThreeBodyMasses(m1, m2, m3; m0)
    tm0 = typeof(m0)
    tm0 <: Number && (m0 < m1 + m2 + m3) && error("m₀ should be bigger than m₁+m₂+m₃")
    MassTuple{tm0}((m1, m2, m3, m0))
end
# 
ThreeBodyMasses(; m1, m2, m3, m0) = ThreeBodyMasses(m1, m2, m3; m0)

function lims(ms::MassTuple; k::Int)
    i, j = ij_from_k(k)
    ((ms[i] + ms[j])^2, (ms[4] - ms[k])^2)
end
lims(k::Int, ms::MassTuple) = lims(ms; k)
lims1(ms::MassTuple) = lims(ms; k=1)
lims2(ms::MassTuple) = lims(ms; k=2)
lims3(ms::MassTuple) = lims(ms; k=3)
#
import Base: getindex, ^, length, iterate
^(ms::MassTuple, i::Int) = Tuple(ms) .^ i

# -----------------------------------------------------
const SpinTuple = NamedTuple{(:two_h1, :two_h2, :two_h3, :two_h0),NTuple{4,Int}}

function ThreeBodySpins(two_h1_or_h1, two_h2_or_h2, two_h3_or_h3;
    h0=-1, # unrealistic by default
    two_h0=-1) # unrealistic by default

    two_h0 == -1 && h0 == -1 && return error("Use either two_h0=..., or h0=... key word argument.")
    two_hs = (h0 == -1 && two_h0 != -1) ?
             SpinTuple(Tuple(Int[two_h1_or_h1, two_h2_or_h2, two_h3_or_h3, two_h0])) :
             SpinTuple(Tuple([two_h1_or_h1, two_h2_or_h2, two_h3_or_h3, h0] .|> x2))
    # 
    @unpack two_h1, two_h2, two_h3 = two_hs
    isodd(two_h1 + two_h2 + two_h3 + two_hs.two_h0) ?
    error("baryon number is not conserved") :
    return two_hs
end
# 
# 
@with_kw struct ThreeBodySystem{T,K}
    ms::T
    two_js::K = ThreeBodySpins(0, 0, 0; two_h0=0)
end
#
# convenient constructors
ThreeBodySystem(ms::MassTuple) = ThreeBodySystem(ms=ms)
ThreeBodySystem(m1, m2, m3; m0, two_js=ThreeBodySpins(0, 0, 0; two_h0=0)) =
    ThreeBodySystem(ThreeBodyMasses(m1, m2, m3; m0=m0), two_js)
#
masses(tbs::ThreeBodySystem) = tbs.ms
spins(tbs::ThreeBodySystem) = tbs.two_js

# -----------------------------------------------------

@with_kw struct SpinParity
    two_j::Int
    p::Char
end

SpinParity(s::String) = str2jp(s)

length(jp1::SpinParity) = 0

# dealing with spin 1/2
x2(v) = @. Int(2v)
x2(v::AbstractString) = Meta.parse(v) |> eval |> x2

_d2(two_s::Int) = iseven(two_s) ? "$(div(two_s,2))" : "$(two_s)/2"
d2(v) = _d2.(v)


⊗(p1::Char, p2::Char) = p1 == p2 ? '+' : '-'
⊗(jp1::SpinParity, jp2::SpinParity) =
    [SpinParity(two_j, ⊗(jp1.p, jp2.p))
     for two_j in abs(jp1.two_j - jp2.two_j):2:abs(jp1.two_j + jp2.two_j)]

function str2jp(pin::AbstractString)
    p = filter(x -> x != '^', pin)
    !(contains(p, '/')) && return SpinParity(x2(p[1:end-1]), p[end])
    p[end-2:end-1] != "/2" && error("the string should be `x/2±`, while it is $(p)")
    two_j = Meta.parse(p[1:end-3])
    !(typeof(two_j) <: Int) && error("the string should be `x/2±`, while it is $(p)")
    return SpinParity(two_j, p[end])
end

macro jp_str(p)
    return str2jp(p)
end

# -----------------------------------------------------

const ParityTuple = NamedTuple{(:P1, :P2, :P3, :P0),NTuple{4,Char}}
#
ThreeBodyParities(P1, P2, P3;
    P0=error("used the format ThreeBodyParities('+','-','+'; P0='±')")) =
    ParityTuple((P1, P2, P3, P0))

# -----------------------------------------------------

function ThreeBodySpinParities(jp1::SpinParity, jp2::SpinParity, jp3::SpinParity;
    jp0::SpinParity=error("Provide jp0 as a key argument"))
    two_js = ThreeBodySpins(jp1.two_j, jp2.two_j, jp3.two_j; two_h0=jp0.two_j)
    Ps = ThreeBodyParities(jp1.p, jp2.p, jp3.p; P0=jp0.p)
    return (two_js, Ps)
end

function ThreeBodySpinParities(jp1::AbstractString, jp2::AbstractString, jp3::AbstractString;
    jp0::AbstractString=error("Provide jp0 as a key argument"))
    ThreeBodySpinParities(str2jp(jp1), str2jp(jp2), str2jp(jp3); jp0=str2jp(jp0))
end

# Dynamic variables
const MandestamTuple{T} = NamedTuple{(:σ1, :σ2, :σ3),NTuple{3,T}}

"""
	Invariants(ms::MassTuple{T}; σ1, σ2)
	Invariants(ms::MassTuple{T}; σ1, σ3)
	Invariants(ms::MassTuple{T}; σ2, σ3)

Construct a tuple of (σ1, σ2, σ3) from just two invariants and the mass tuple.
"""
function Invariants(ms::MassTuple{T};
    σ1=-one(ms.m0), σ2=-one(ms.m0), σ3=-one(ms.m0)) where {T}
    # 
    !((σ1 == -one(ms.m0)) || (σ2 == -one(ms.m0)) || (σ3 == -one(ms.m0))) &&
        error("the method works with TWO invariants given: $((σ1,σ2,σ3))")
    # 
    σ3 == -one(ms.m0) && return MandestamTuple{T}((σ1, σ2, σ3=sum(ms^2) - σ1 - σ2))
    σ1 == -one(ms.m0) && return MandestamTuple{T}((sum(ms^2) - σ2 - σ3, σ2, σ3))
    return MandestamTuple{T}((σ1=σ1, σ2=sum(ms^2) - σ3 - σ1, σ3=σ3))
end
Invariants(; σ1, σ2, σ3) = MandestamTuple{typeof(σ1)}((σ1, σ2, σ3))
Invariants(σ1, σ2, σ3) = MandestamTuple{typeof(σ1)}((σ1, σ2, σ3))

# -----------------------------------------------------

circleorigin(k, t) = (t[mod(k, 3)+1], t[mod(k + 1, 3)+1], t[mod(k - 1, 3)+1])

fitin(y, (a, b)) = a + y * (b - a)
function x2σs(x, ms::MassTuple; k::Int)
    σk = fitin(x[1], lims(ms; k))
    σj = σjofk(2x[2] - 1, σk, ms^2; k)
    σi = sum(ms^2) - σk - σj
    σt = circleorigin(-k, (σi, σj, σk))
    return MandestamTuple{typeof(ms.m0)}(σt)
end

function y2σs(y, ms::MassTuple; k::Int=last(findmin(Tuple(ms))))
    i, j = ij_from_k(k)
    σi = fitin(y[1], lims(ms; k=i))
    σj = fitin(y[2], lims(ms; k=j))
    σk = sum(ms^2) - σi - σj
    σt = circleorigin(-k, (σi, σj, σk))
    return MandestamTuple{typeof(ms.m0)}(σt)
end

randomPoint(ms::MassTuple) = x2σs(rand(2), ms; k=3)
randomPoint(two_js::SpinTuple) = SpinTuple([rand(-j:2:j) for j in two_js])

@with_kw struct DalitzPlotPoint{I,S}
    σs::I
    two_λs::S
end

function randomPoint(tbs::ThreeBodySystem)
    DalitzPlotPoint(
        σs=randomPoint(tbs.ms),
        two_λs=randomPoint(tbs.two_js))
end

#                                                                _|
#    _|_|_|    _|_|    _|_|_|      _|_|    _|  _|_|    _|_|_|  _|_|_|_|    _|_|
#  _|    _|  _|_|_|_|  _|    _|  _|_|_|_|  _|_|      _|    _|    _|      _|_|_|_|
#  _|    _|  _|        _|    _|  _|        _|        _|    _|    _|      _|
#    _|_|_|    _|_|_|  _|    _|    _|_|_|  _|          _|_|_|      _|_|    _|_|_|
#        _|
#    _|_|

"""
	polardalitz2invariants(θ, expansion_point)

For given polar angle θ, it returns an (σ1,σ2,σ3) Tuple of polynomials of radius r(θ) around the expansion point.
The polynomial works as a function of the r coordinate.
"""
polardalitz2invariants(θ, expansion_point::Tuple) =
    (Polynomial([0, -cos(θ)]),
        Polynomial([0, cos(θ + π / 3)]),
        Polynomial([0, cos(θ - π / 3)])) .+ expansion_point

function border(ms::MassTuple{T}; Nx::Int=300) where {T}
    # 
    expansion_point = let
        f = 0.5
        z = 0.0
        σ1 = (ms[2] + ms[3])^2 + f * ((ms[4] - ms[1])^2 - (ms[2] + ms[3])^2)
        σ3 = σ3of1(z, σ1, ms^2)
        Invariants(ms; σ1, σ3)
    end
    # 
    σs(θ) = polardalitz2invariants(θ, expansion_point |> Tuple)
    ϕ0 = Kibble(expansion_point, ms^2)
    ϕ(σs) = Kibble(σs, ms^2)
    # 
    function rborder(θ)
        _roots = PolynomialRoots.roots(coeffs(ϕ(σs(θ))))
        filter(_roots) do r
            (abs(imag(r)) < 1e-10) && real(r) > 0.0
        end |> real |> minimum
    end
    function σsborder(θ) # evaluate the polynomials
        r = rborder(θ)
        return map(P -> P(r), σs(θ))
    end
    θs = range(-π / 9, 2π - π / 9, length=Nx)
    σs_tuple = σsborder.(θs)
    return MandestamTuple{T}.(σs_tuple)
end

# border13, border12, border21, border23, border32
for (i, j) in ((1, 2), (2, 1), (2, 3), (3, 2), (3, 1), (1, 3))
    eval(quote
        $(Symbol(:border, i, j))(ms; Nx::Int=300) =
            NamedTuple{$(Symbol(:σ, i), Symbol(:σ, j))}.(border(ms; Nx))
    end)
end

#
inrange(x, r) = r[1] < x < r[2]
inphrange(σs::MandestamTuple, ms::MassTuple) = Kibble(σs, ms^2) < 0 &&
                                               inrange(σs[1], lims1(ms)) && inrange(σs[2], lims2(ms)) && inrange(σs[3], lims3(ms))


