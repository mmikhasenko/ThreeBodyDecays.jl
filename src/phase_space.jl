"""
    Kallen(x, y, z)

Calculate the Källén function λ(x,y,z) = x² + y² + z² - 2xy - 2yz - 2zx.
This function appears frequently in relativistic kinematics calculations.

# Arguments
- `x`, `y`, `z`: Real numbers or complex values

# Returns
- The value of the Källén function

# Examples
```jldoctest; setup = :(using ThreeBodyDecays)
julia> Kallen(10.0, 1.0, 1.0)  # For a decay A -> B + C, calculate λ(mA², mB², mC²)
60.0

julia> Kallen(5.0, 1.0, 1.0)  # Another example with different masses
5.0
```

See also [`sqrtKallenFact`](@ref), [`Kibble`](@ref).
"""
Kallen(x, y, z) = x^2 + y^2 + z^2 - 2x * y - 2y * z - 2z * x

"""
    sqrtKallenFact(a, b, c)

Calculate the square root of the Källén function λ(a,b,c) in factorized form.
The Källén function is defined as λ(a,b,c) = a² + b² + c² - 2ab - 2bc - 2ca.
This function returns √λ(a,b,c) in a factorized form to improve numerical stability:
√(a-(b+c)) √(a-(b-c)) √(a+(b+c)) √(a+(b-c)).

# Arguments
- `a`, `b`, `c`: Input values for the Källén function

# Returns
- Square root of the Källén function in factorized form

# Examples
```jldoctest; setup = :(using ThreeBodyDecays)
julia> sqrtKallenFact(10.0, 1.0, 1.0)  # √λ(10,1,1) = √9600
97.97958971132714

julia> sqrtKallenFact(4.0, 1.0, 1.0)  # √λ(4,1,1) = √192
13.856406460551018
```

See also [`Kallen`](@ref), [`Kibble`](@ref).
"""
sqrtKallenFact(a, b, c) =
    sqrt(a - (b + c)) * sqrt(a - (b - c)) * sqrt(a + (b + c)) * sqrt(a + (b - c))

"""
    Kibble(σs, msq)

Calculate the Kibble function φ(σ₁,σ₂,σ₃) for three-body decays, which is defined in terms of Källén functions.
The Kibble function determines the physical boundaries of the three-body phase space.

The function is calculated as:
φ(σ₁,σ₂,σ₃) = λ(λ₁,λ₂,λ₃)
where λᵢ = λ(M², mᵢ², σᵢ) are Källén functions of the total mass squared M²,
the i-th particle mass squared mᵢ², and the corresponding Mandelstam variable σᵢ.

# Arguments
- `σs`: Tuple of Mandelstam variables (σ₁,σ₂,σ₃)
- `msq`: Tuple of squared masses (m₁²,m₂²,m₃²,M²)

# Returns
- Value of the Kibble function

# Examples
```jldoctest; setup = :(using ThreeBodyDecays)
julia> msq = (1.0, 1.0, 1.0, 16.0);  # squared masses: m₁², m₂², m₃², M²

julia> σs = (2.0, 2.0, 2.0);         # Mandelstam variables

julia> Kibble(σs, msq)               # returns the Kibble function value
-77763.0

julia> Kibble(σs, msq) < 0           # physical configuration if negative
true
```

See also [`Kallen`](@ref), [`isphysical`](@ref).
"""
Kibble(σs, msq) = Kallen(
    Kallen(msq[4], msq[1], σs[1]),
    Kallen(msq[4], msq[2], σs[2]),
    Kallen(msq[4], msq[3], σs[3]),
)


"""
    σjofk(z,σi,msq; k::Int)

Computes invariant σj = (p0-pj)² from
the scattering angle z=cosθij in the rest from of (i,j),
given the mass of the system m(i,j)² = σk

Explicit forms: `σ3of1`, `σ1of2`, `σ2of3`.

See also [`σiofk`](@ref)
"""
function σjofk(z, σk, msq; k::Int)
    i, j = ij_from_k(k)
    #
    s = msq[4]
    # σj = (p0-pj)² in the rest frame of (i,j)
    EE4σ = (σk + msq[j] - msq[i]) * (σk + s - msq[k])
    p²q²4σ = Kallen(σk, msq[i], msq[j]) * Kallen(s, σk, msq[k])
    p²q²4σ = (p²q²4σ < 0) ? 0.0 : p²q²4σ # for numerical errors
    σi = s + msq[j] - (EE4σ - sqrt(p²q²4σ) * z) / (2σk)
    return σi
end
σ3of1(z, σ1, msq) = σjofk(z, σ1, msq; k = 1)
σ1of2(z, σ2, msq) = σjofk(z, σ2, msq; k = 2)
σ2of3(z, σ3, msq) = σjofk(z, σ3, msq; k = 3)


"""
    σiofk(k,z,σj,msq)

Computes invariant σi = (p0 - pi)² from
the scattering angle z=cosθij in the rest from of (i,j),
given the mass of the system m(i,j)² = σk

Explicit forms: `σ3of2`, `σ1of3`, `σ2of1`.

See also [`σjofk`](@ref)
"""
function σiofk(z, σk, msq; k::Int)
    σj = σjofk(z, σk, msq; k)
    sum(msq) - σj - σk
end
σ3of2(z, σ2, msq) = σiofk(z, σ2, msq; k = 2)
σ1of3(z, σ3, msq) = σiofk(z, σ3, msq; k = 3)
σ2of1(z, σ1, msq) = σiofk(z, σ1, msq; k = 1)


# Scattering angle

"""
cosθij(k,σs,msq)

Isobar decay angle for the chain-k, i.e.
an angle of between vectors pi and -pk in the (ij) rest frame.

Explicit forms: `cosθ23`, `cosθ31`, `cosθ12`.
"""
function cosθij(σs, msq; k::Int)
    i, j = ij_from_k(k)
    #
    s = msq[4]
    EE4σ = (σs[k] + msq[i] - msq[j]) * (s - σs[k] - msq[k])
    pp4σ = sqrt(Kallen(σs[k], msq[i], msq[j]) * Kallen(s, σs[k], msq[k]))
    rest = σs[j] - msq[k] - msq[i]
    return (2σs[k] * rest - EE4σ) / pp4σ
end
cosθ23(σs, msq) = cosθij(σs, msq; k = 1)
cosθ31(σs, msq) = cosθij(σs, msq; k = 2)
cosθ12(σs, msq) = cosθij(σs, msq; k = 3)


# momentum
"""
    breakup(m, m1, m2)
    breakup_Rk(σk, ms; k)
    breakup_ij(σk, ms; k)

Calculate the breakup momentum for a two-body system.
- `breakup`: General form for any two-body decay
- `breakup_Rk`: For decay of the total system into a subsystem (ij), and a particle k
- `breakup_ij`: For decay of the subsystem (ij) into particles i and j

# Arguments
- `m`: Mass of the decaying system
- `m1`, `m2`: Masses of the decay products
- `σk`: Invariant mass squared of the isobar
- `ms`: MassTuple containing all masses
- `k`: Spectator index

# Returns
- Break-up momentum magnitude

# Examples
```jldoctest
julia> using ThreeBodyDecays  # hide

julia> ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0=4.0);

julia> p = breakup_Rk(2.5^2, ms; k=1)
0.897587913521567

julia> q = breakup_ij(2.5^2, ms; k=1)
0.75

julia> breakup(4.0, 2.5, 1.0)  # direct use of breakup function
0.897587913521567
```
"""
breakup(m, m1, m2) =
    sqrt((m - (m1 + m2)) * (m + (m1 + m2)) * (m - (m1 - m2)) * (m + (m1 - m2))) / 2m
breakup_Rk(σk, ms; k) = breakup(ms[4], sqrt(σk), ms[k])
function breakup_ij(σk, ms; k)
    (i, j) = ij_from_k(k)
    return breakup(sqrt(σk), ms[i], ms[j])
end
breakup_Rk(ms; k) = σ -> breakup(ms[4], sqrt(σ), ms[k])
function breakup_ij(ms; k)
    (i, j) = ij_from_k(k)
    return σ -> breakup(sqrt(σ), ms[i], ms[j])
end

"""
    aligned_four_vectors(σs,ms; k::Int)

Computes the four-momenta of the three particles in the center of momentum frame aligning the k-th particle with the -z-axis.

## Arguments
- `σs`: Tuple of mandelstam variables,
- `ms`: Tuple of masses in the order m1, m2, m3, m0.
- `k`: Index of the particle to be aligned with the -z-axis.

## Returns
A tuple of three four-momenta in the form of (px, py, pz, E).

## Examples
```jldoctest; setup = :(using ThreeBodyDecays)
julia> ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0=4.0);

julia> σs = x2σs([0.5, 0.5], ms; k=2);

julia> p1, p2, p3 = aligned_four_vectors(σs, ms; k=1);

julia> p1[3] ≈ -breakup_Rk(σs[1], ms; k=1)  # first particle aligned with -z axis
true

julia> all(p -> abs(p[2]) < 1e-10, (p1, p2, p3))  # all momenta in x-z plane
true
```
"""
function aligned_four_vectors(σs, ms; k::Int)
    #
    i, j = ij_from_k(k)
    m0 = ms[4]
    s = m0^2
    #
    mi, mj, mk = ms[i], ms[j], ms[k]
    σi, σj, σk = σs[i], σs[j], σs[k]
    misq, mjsq, mksq = (mi, mj, mk) .^ 2
    #
    Ei = (s + misq - σi) / (2m0)
    Ej = (s + mjsq - σj) / (2m0)
    Ek = (s + mksq - σk) / (2m0)
    #
    p3i = sqrt(Kallen(s, misq, σi)) / (2m0)
    p3j = sqrt(Kallen(s, mjsq, σj)) / (2m0)
    p3k = sqrt(Kallen(s, mksq, σk)) / (2m0)
    #
    cosθ = -cosζ(wr(k, i, 0), σs, ms^2) # pi-angle(1,2) in rest of 0
    sinθ = sqrt(1 - cosθ^2)
    #
    pk = (0, 0, -p3k, Ek)
    pi = (p3i * sinθ, 0, p3i * cosθ, Ei)
    pj = (-p3i * sinθ, 0, p3k - p3i * cosθ, Ej)
    #
    return circleorigin(-k, (pi, pj, pk))
end

inrange(x, r) = r[1] < x < r[2]
inphrange(σs, ms::MassTuple) = isphysical(σs, ms)

"""
    isphysical(σs::MandelstamTuple, ms::MassTuple) -> Bool
    inphrange(σs::MandelstamTuple, ms::MassTuple) -> Bool
    inrange(x, r) -> Bool

Check if a set of Mandelstam variables represents a physically valid configuration in the three-body phase space.
A configuration is physical if:
1. The Kibble function is negative (φ(σ₁,σ₂,σ₃) < 0)
2. All Mandelstam variables are within their kinematically allowed ranges

# Arguments
- `σs::MandelstamTuple`: Tuple of Mandelstam variables (σ₁,σ₂,σ₃)
- `ms::MassTuple`: Tuple of masses (m₁,m₂,m₃,M)
- `x`: Value to check (for `inrange`)
- `r`: Range tuple (min,max) to check against (for `inrange`)

# Returns
- `Bool`: `true` if the configuration is physical, `false` otherwise

# Examples
```jldoctest; setup = :(using ThreeBodyDecays)
julia> ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0=4.0);

julia> σs = (6.25, 6.25, 6.5);

julia> isphysical(σs, ms)  # checks if this configuration is physically possible
true

julia> σs_unphysical = (20.0, 20.0, 20.0);  # values outside allowed ranges

julia> isphysical(σs_unphysical, ms)  # this configuration is not physical
false
```

See also [`Kibble`](@ref), [`lims`](@ref).
"""
isphysical(σs, ms::MassTuple) =
    Kibble(σs, ms^2) < 0 &&
    inrange(σs[1], lims1(ms)) &&
    inrange(σs[2], lims2(ms)) &&
    inrange(σs[3], lims3(ms))
