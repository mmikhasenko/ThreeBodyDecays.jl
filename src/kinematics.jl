
Kallen(x, y, z) = x^2 + y^2 + z^2 - 2x * y - 2y * z - 2z * x
sqrtKallenFact(a, b, c) =
    sqrt(a - (b + c)) * sqrt(a - (b - c)) * sqrt(a + (b + c)) * sqrt(a + (b - c))
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

See also `σiofk(z,σj,msq; k)`
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

See also `σjofk(z,σk,msq; k)`
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

## Example
````julia
ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0=4.0)
σs = x2σs([0.5, 0.5], ms; k=2)
p1, p2, p3 = aligned_four_vectors(σs, ms; k=1)
````
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
