"""
    SpinTuple

A named tuple representing the spins of a three-body system.
Contains twice the helicities of the particles (2h₁, 2h₂, 2h₃, 2h₀).
"""
const SpinTuple = NamedTuple{(:two_h1, :two_h2, :two_h3, :two_h0),NTuple{4,Int}}

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
