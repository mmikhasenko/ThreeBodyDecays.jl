#                                  _|
#    _|_|_|  _|    _|    _|_|_|  _|_|_|_|    _|_|    _|_|_|  _|_|
#  _|_|      _|    _|  _|_|        _|      _|_|_|_|  _|    _|    _|
#      _|_|  _|    _|      _|_|    _|      _|        _|    _|    _|
#  _|_|_|      _|_|_|  _|_|_|        _|_|    _|_|_|  _|    _|    _|
#                  _|
#              _|_|

"""
    MassTuple{T}

A named tuple representing the masses of a three-body system.
Contains masses m₁, m₂, m₃ of the decay products and m₀ of the parent particle.
"""
const MassTuple{T} = NamedTuple{(:m1, :m2, :m3, :m0), NTuple{4, T}}

"""
    ThreeBodyMasses(m1, m2, m3; m0)
    ThreeBodyMasses(; m1, m2, m3, m0)

Construct a MassTuple for a three-body system.

# Arguments
- `m1`, `m2`, `m3`: Masses of the decay products
- `m0`: Mass of the parent particle

# Returns
- `MassTuple{T}`: A named tuple containing the masses

# Throws
- `ErrorException` if m₀ is less than the sum of m₁, m₂, and m₃
"""
function ThreeBodyMasses(m1, m2, m3; m0)
    tm0 = typeof(m0)
    tm0 <: Number && (m0 < m1 + m2 + m3) && error("m₀ should be bigger than m₁+m₂+m₃")
    MassTuple{tm0}((m1, m2, m3, m0))
end

ThreeBodyMasses(; m1, m2, m3, m0) = ThreeBodyMasses(m1, m2, m3; m0)

"""
    lims(ms::MassTuple; k::Int)
    lims(k::Int, ms::MassTuple)

Calculate the kinematic limits for the k-th Mandelstam variable.

# Arguments
- `ms::MassTuple`: Masses of the system
- `k::Int`: Index of the variable (1, 2, or 3)

# Returns
- `Tuple{T,T}`: Lower and upper limits for the variable
"""
function lims(ms::MassTuple; k::Int)
    i, j = ij_from_k(k)
    ((ms[i] + ms[j])^2, (ms[4] - ms[k])^2)
end
lims(k::Int, ms::MassTuple) = lims(ms; k)
lims1(ms::MassTuple) = lims(ms; k = 1)
lims2(ms::MassTuple) = lims(ms; k = 2)
lims3(ms::MassTuple) = lims(ms; k = 3)

import Base: getindex, ^, length, iterate
^(ms::MassTuple, i::Int) = Tuple(ms) .^ i

# -----------------------------------------------------

"""
    SpinTuple

A named tuple representing the spins of a three-body system.
Contains twice the helicities of the particles (2h₁, 2h₂, 2h₃, 2h₀).
"""
const SpinTuple = NamedTuple{(:two_h1, :two_h2, :two_h3, :two_h0), NTuple{4, Int}}

"""
    ThreeBodySpins(two_h1_or_h1, two_h2_or_h2, two_h3_or_h3; h0=nothing, two_h0=nothing)

Construct a SpinTuple for a three-body system. Depending on the key arguments, `h0` or `two_h0`,
the position arguments can be either twice the helicity or helicity itself.

# Arguments
- `two_h1_or_h1`: Twice the helicity or helicity of first particle
- `two_h2_or_h2`: Twice the helicity or helicity of second particle
- `two_h3_or_h3`: Twice the helicity or helicity of third particle
- `h0`: Key argument for helicity of parent particle (optional)
- `two_h0`: Key argument for twice the helicity of parent particle (optional)

# Returns
- `SpinTuple`: A named tuple containing the spins

# Throws
- `ErrorException` if neither `h0` nor `two_h0` is provided
- `ErrorException` if baryon number is not conserved
"""
function ThreeBodySpins(
    two_h1_or_h1,
    two_h2_or_h2,
    two_h3_or_h3;
    h0 = nothing,  # default to nothing
    two_h0 = nothing, # default to nothing
)
    if isnothing(h0) && isnothing(two_h0)
        error("Use either `two_h0=...`, or `h0=...` keyword argument.")
    end

    two_hs =
        (isnothing(h0) && !isnothing(two_h0)) ?
        SpinTuple(Tuple(Int[two_h1_or_h1, two_h2_or_h2, two_h3_or_h3, two_h0])) :
        SpinTuple(Tuple([two_h1_or_h1, two_h2_or_h2, two_h3_or_h3, h0] .|> x2))

    @unpack two_h1, two_h2, two_h3 = two_hs
    isodd(two_h1 + two_h2 + two_h3 + two_hs.two_h0) ?
    error("baryon number is not conserved") : return two_hs
end

"""
    ThreeBodySystem{T,K}

A structure representing a three-body system with masses and spins.

# Fields
- `ms::T`: Masses of the system (MassTuple)
- `two_js::K`: Spins of the system (SpinTuple)
"""
@with_kw struct ThreeBodySystem{T, K}
    ms::T
    two_js::K = ThreeBodySpins(0, 0, 0; two_h0 = 0)
end

# convenient constructors
ThreeBodySystem(ms::MassTuple) = ThreeBodySystem(ms = ms)
ThreeBodySystem(m1, m2, m3; m0, two_js = ThreeBodySpins(0, 0, 0; two_h0 = 0)) =
    ThreeBodySystem(ThreeBodyMasses(m1, m2, m3; m0 = m0), two_js)

masses(tbs::ThreeBodySystem) = tbs.ms
spins(tbs::ThreeBodySystem) = tbs.two_js

# -----------------------------------------------------

"""
    SpinParity

A structure representing spin and parity of a particle.

# Fields
- `two_j::Int`: Twice the spin value
- `p::Char`: Parity ('+' or '-')
"""
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

"""
    ⊗(p1::Char, p2::Char)
    ⊗(jp1::SpinParity, jp2::SpinParity)

Compute the tensor product of two parities or spin-parity states.

# Arguments
- `p1`, `p2`: Parity characters ('+' or '-')
- `jp1`, `jp2`: SpinParity objects

# Returns
- For parities: The combined parity
- For spin-parity: Array of possible combined states
"""
⊗(p1::Char, p2::Char) = p1 == p2 ? '+' : '-'
⊗(jp1::SpinParity, jp2::SpinParity) = [
    SpinParity(two_j, ⊗(jp1.p, jp2.p)) for
    two_j ∈ abs(jp1.two_j - jp2.two_j):2:abs(jp1.two_j + jp2.two_j)
]

"""
    str2jp(pin::AbstractString)

Convert a string representation of spin-parity to a SpinParity object.

# Arguments
- `pin::AbstractString`: String in format "x±" or "x/2±"

# Returns
- `SpinParity`: The corresponding spin-parity object

# Throws
- `ErrorException` if the string format is invalid
"""
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

"""
    ParityTuple

A named tuple representing parities of a three-body system.
Contains parities P₁, P₂, P₃ of the decay products and P₀ of the parent particle.
"""
const ParityTuple = NamedTuple{(:P1, :P2, :P3, :P0), NTuple{4, Char}}

"""
    ThreeBodyParities(P1, P2, P3; P0)

Construct a ParityTuple for a three-body system.

# Arguments
- `P1`, `P2`, `P3`: Parities of the decay products
- `P0`: Parity of the parent particle

# Returns
- `ParityTuple`: A named tuple containing the parities
"""
ThreeBodyParities(
    P1,
    P2,
    P3;
    P0 = error("used the format ThreeBodyParities('+','-','+'; P0='±')"),
) = ParityTuple((P1, P2, P3, P0))

# -----------------------------------------------------

"""
    ThreeBodySpinParities(jp1, jp2, jp3; jp0)

Construct spin and parity information for a three-body system using shortcuts like "x/2±", or jp"x/2±".

# Arguments
- `jp1`, `jp2`, `jp3`: SpinParity objects for decay products
- `jp0`: SpinParity object for parent particle

# Returns
- `Tuple{SpinTuple,ParityTuple}`: Tuple of spin structure and parity structure
"""
function ThreeBodySpinParities(
    jp1::SpinParity,
    jp2::SpinParity,
    jp3::SpinParity;
    jp0::SpinParity = error("Provide jp0 as a key argument"),
)
    two_js = ThreeBodySpins(jp1.two_j, jp2.two_j, jp3.two_j; two_h0 = jp0.two_j)
    Ps = ThreeBodyParities(jp1.p, jp2.p, jp3.p; P0 = jp0.p)
    return (two_js, Ps)
end

function ThreeBodySpinParities(
    jp1::AbstractString,
    jp2::AbstractString,
    jp3::AbstractString;
    jp0::AbstractString = error("Provide jp0 as a key argument"),
)
    ThreeBodySpinParities(str2jp(jp1), str2jp(jp2), str2jp(jp3); jp0 = str2jp(jp0))
end

# Dynamic variables
"""
    MandelstamTuple{T}

A named tuple representing Mandelstam variables `(; σ1, σ2, σ3)` for a three-body system.
"""
const MandelstamTuple{T} = NamedTuple{(:σ1, :σ2, :σ3), NTuple{3, T}}

"""
    Invariants(ms::MassTuple{T}; σ1, σ2, σ3)

Construct a tuple of Mandelstam invariants from given invariants and the mass tuple.

# Arguments
- `ms::MassTuple{T}`: Masses of the system
- `σ1`, `σ2`, `σ3`: Optional invariants (two or three must be provided)

# Returns
- `MandelstamTuple{T}`: The complete set of invariants

# Throws
- `ErrorException` if not exactly two or three invariants are provided
- `ErrorException` if the invariants violate mass constraints
"""
function Invariants(
    ms::MassTuple{T};
    σ1 = -one(ms.m0),
    σ2 = -one(ms.m0),
    σ3 = -one(ms.m0),
) where {T}
    # Count how many invariants are provided
    n_provided = count(x -> x != -one(ms.m0), (σ1, σ2, σ3))

    if n_provided == 0 || n_provided == 1
        error("At least two invariants must be provided")
    elseif n_provided == 2
        # Calculate the third invariant
        if σ3 == -one(ms.m0)
            σ3 = sum(ms^2) - σ1 - σ2
        elseif σ1 == -one(ms.m0)
            σ1 = sum(ms^2) - σ2 - σ3
        else
            σ2 = sum(ms^2) - σ3 - σ1
        end
    end

    # Validate the invariants
    if !isphysical((σ1, σ2, σ3), ms)
        error("The provided invariants violate physical constraints")
    end

    return MandelstamTuple{T}((σ1, σ2, σ3))
end
Invariants(; σ1, σ2, σ3) = MandelstamTuple{typeof(σ1)}((σ1, σ2, σ3))
Invariants(σ1, σ2, σ3) = MandelstamTuple{typeof(σ1)}((σ1, σ2, σ3))

# -----------------------------------------------------

circleorigin(k, t) = (t[mod(k, 3)+1], t[mod(k + 1, 3)+1], t[mod(k - 1, 3)+1])

fitin(y, (a, b)) = a + y * (b - a)




"""
    x2σs(x, ms::MassTuple; k::Int = last(findmin(Tuple(ms))))

Maps a pair of variables `x` to a physical set of mandelstam invariants `σs` using a linear transformation.
x[1] is mapped to σ[k], and x[2] is mapped into cosθij, that is used to compute σ[i] and σ[j].
Uniform distribution of x does not lead to the phase space distribution of σs.

## Arguments

- `x` : a pair of numbers
- `ms` : masses of the system as a `MassTuple`, see `ThreeBodyMasses`.
- `k` : the index for which the variable is not generated. By default, the function picks the coordinates where the Dalitz plot has the closest shape to the squared fitting box.

## Returns
- an instance of `MandelstamTuple` with the squared masses.

## Example

The phase space sample with 100 points can be generated as follows:
````julia
σs0 = x2σs([0.3, 0.2], ms)
````

See also [`y2σs`](@ref).
"""
function x2σs(x, ms::MassTuple; k::Int)
    σk = fitin(x[1], lims(ms; k))
    σj = σjofk(2x[2] - 1, σk, ms^2; k)
    σi = sum(ms^2) - σk - σj
    σt = circleorigin(-k, (σi, σj, σk))
    return MandelstamTuple{typeof(ms.m0)}(σt)
end



"""
    y2σs(y, ms::MassTuple; k::Int = last(findmin(Tuple(ms))))

Maps a pair of variables to the plane of squared masses using a linear transformation,
by shifting and scaling `y[1]` to `σ[i]`, and `y[2]` to `σ[j]`.
The indices of variables are determined by the argument `k`.
The mapping does not guarantee that values are physical, however,
physical values of σs are phase-space distributed for uniform y.

## Arguments

- `y` : a pair of numbers
- `ms` : masses of the system as a `MassTuple`, see `ThreeBodyMasses`.
- `k` : the index for which the variable is not generated. By default, the function picks the coordinates
where the Dalitz plot has the closest shape to the squared fitting box.

## Returns

- an instance of `MandelstamTuple` with the squared masses.

## Example

The phase space sample with 100 points can be generated as follows:
````julia
data = let
    N = 100
    # map random variables to dalitz
    _data = mapslices(rand(N, 2); dims=2) do xy
        y2σs(xy, ms)
    end[:, 1]
    # select physical
    filter!(_data) do σs
        isphysical(σs, ms)
    end
    _data
end
````

See also [`x2σs`](@ref).
"""
function y2σs(y, ms::MassTuple; k::Int = last(findmin(Tuple(ms))))
    i, j = ij_from_k(k)
    σi = fitin(y[1], lims(ms; k = i))
    σj = fitin(y[2], lims(ms; k = j))
    σk = sum(ms^2) - σi - σj
    σt = circleorigin(-k, (σi, σj, σk))
    return MandelstamTuple{typeof(ms.m0)}(σt)
end

"""
    randomPoint(ms::MassTuple)
    randomPoint(two_js::SpinTuple)
    randomPoint(tbs::ThreeBodySystem)

Generate a random point in the phase space.

# Arguments
- `ms::MassTuple`: Masses of the system
- `two_js::SpinTuple`: Spins of the system
- `tbs::ThreeBodySystem`: Complete three-body system

# Returns
- For masses: Random Mandelstam variables
- For spins: Random spin configuration
- For system: Random DalitzPlotPoint
"""
randomPoint(ms::MassTuple) = x2σs(rand(2), ms; k = 3)
randomPoint(two_js::SpinTuple) = SpinTuple([rand(-j:2:j) for j in two_js])

@with_kw struct DalitzPlotPoint{I, S}
    σs::I
    two_λs::S
end

function randomPoint(tbs::ThreeBodySystem)
    DalitzPlotPoint(σs = randomPoint(tbs.ms), two_λs = randomPoint(tbs.two_js))
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

For given polar angle θ, returns an (σ1,σ2,σ3) Tuple of polynomials of radius r(θ) around the expansion point.
The polynomial works as a function of the r coordinate.

# Arguments
- `θ`: Polar angle
- `expansion_point`: Tuple of expansion point coordinates

# Returns
- `Tuple{Polynomial,Polynomial,Polynomial}`: Polynomials for each invariant
"""
polardalitz2invariants(θ, expansion_point::Tuple) =
    (
        Polynomial([0, -cos(θ)]),
        Polynomial([0, cos(θ + π / 3)]),
        Polynomial([0, cos(θ - π / 3)]),
    ) .+ expansion_point

"""
    border(ms::MassTuple{T}; Nx::Int=300) where T

Calculate the border of the Dalitz plot.

# Arguments
- `ms::MassTuple{T}`: Masses of the system
- `Nx::Int`: Number of points to generate

# Returns
- `Vector{MandelstamTuple{T}}`: Points on the border
"""
function border(ms::MassTuple{T}; Nx::Int = 300) where {T}
    # Calculate a physically valid expansion point
    expansion_point = let
        f = 0.5
        z = 0.0
        σ1 = (ms[2] + ms[3])^2 + f * ((ms[4] - ms[1])^2 - (ms[2] + ms[3])^2)
        σ3 = σ3of1(z, σ1, ms^2)
        Invariants(ms; σ1, σ3)
    end

    σs(θ) = polardalitz2invariants(θ, expansion_point |> Tuple)
    ϕ0 = Kibble(expansion_point, ms^2)
    ϕ(σs) = Kibble(σs, ms^2)

    function rborder(θ)
        _roots = PolynomialRoots.roots(coeffs(ϕ(σs(θ))))
        filter(_roots) do r
            (abs(imag(r)) < 1e-10) && real(r) > 0.0
        end |> real |> minimum
    end

    function σs_border(θ)
        r = rborder(θ)
        return map(P -> P(r), σs(θ))
    end

    θs = range(-π / 9, 2π - π / 9, length = Nx)
    σs_tuple = σs_border.(θs)
    return MandelstamTuple{T}.(σs_tuple)
end

# border13, border12, border21, border23, border32
for (i, j) in ((1, 2), (2, 1), (2, 3), (3, 2), (3, 1), (1, 3))
    eval(
        quote
            $(Symbol(:border, i, j))(ms; Nx::Int = DEFAULT_BORDER_POINTS) =
                NamedTuple{$(Symbol(:σ, i), Symbol(:σ, j))}.(border(ms; Nx))
        end,
    )
end

"""
    inrange(x, r)
    inphrange(σs::MandelstamTuple, ms::MassTuple)
    isphysical(σs::MandelstamTuple, ms::MassTuple)

Check if values are within physical ranges.

# Arguments
- `x`: Value to check
- `r`: Range to check against
- `σs`: Mandelstam variables
- `ms`: Masses of the system

# Returns
- `Bool`: Whether the values are physical
"""
inrange(x, r) = r[1] < x < r[2]
inphrange(σs, ms::MassTuple) = isphysical(σs, ms)
isphysical(σs, ms::MassTuple) =
    Kibble(σs, ms^2) < 0 &&
    inrange(σs[1], lims1(ms)) &&
    inrange(σs[2], lims2(ms)) &&
    inrange(σs[3], lims3(ms))
