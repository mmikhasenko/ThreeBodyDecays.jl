"""
    ParityTuple

A named tuple representing parities of a three-body system.
Contains parities P₁, P₂, P₃ of the decay products and P₀ of the parent particle.
"""
const ParityTuple = NamedTuple{(:P1, :P2, :P3, :P0),NTuple{4,Char}}

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
