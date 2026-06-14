"""
    Recoupling

Abstract supertype for recoupling schemes used inside [`Vertex`](@ref).

Concrete subtypes implement `amplitude(::Recoupling, (two_λa, two_λb), (two_j, two_ja, two_jb))`,
which provides the spin/helicity-dependent factor for a given vertex.
"""
abstract type Recoupling end

"""
    NoRecoupling(two_λa, two_λb)

Trivial recoupling: select a single helicity configuration.

The corresponding `amplitude` is 1 only when the requested helicities match the stored
`two_λa`, `two_λb`, and 0 otherwise.
"""
@with_kw struct NoRecoupling <: Recoupling
    two_λa::Int
    two_λb::Int
end

amplitude(cs::NoRecoupling, (two_λa, two_λb), (two_j, two_ja, two_jb)) =
    (cs.two_λa == two_λa) * (cs.two_λb == two_λb)

"""
    ParityRecoupling(two_λa, two_λb, ηηηphasesign)
    ParityRecoupling(two_λa, two_λb, topology::Pair{SpinParity,Tuple{SpinParity,SpinParity}})

Parity-related recoupling that connects `(λa, λb)` with `(-λa, -λb)`
according to the intrinsic-parity phase of the two-body vertex.

This is commonly used to enforce parity constraints in helicity amplitudes.

The amplitude for input helicities ``λ₁, λ₂`` (in half-integer units, use `two_λa = 2λ₁`, etc.) is:

```math
A(λ₁, λ₂) = δ_{λ₁,λ_a} δ_{λ₂,λ_b} + η \\, δ_{λ₁,-λ_a} δ_{λ₂,-λ_b}
```

where ``η = ±1`` is the parity phase (``η = +1`` when `ηηηphaseisplus` is true).
Only the configurations ``(λ₁,λ₂) = (λ_a,λ_b)`` or ``(-λ_a,-λ_b)`` give a non-zero result.
"""
@with_kw struct ParityRecoupling <: Recoupling
    two_λa::Int
    two_λb::Int
    ηηηphaseisplus::Bool
end

ParityRecoupling(two_λa::Int, two_λb::Int, ηηηphasesign::Char) =
    ParityRecoupling(two_λa, two_λb, ηηηphasesign == '+')
function ParityRecoupling(
    two_λa::Int,
    two_λb::Int,
    (jp, (jp1, jp2))::TwoBodyTopologySpinParity,
)
    ηηη = jp1.p ⊗ jp2.p ⊗ jp.p
    ηηηphase = (2 * (ηηη == '+') - 1) * x"-1"^(div(jp.two_j - jp1.two_j - jp2.two_j, 2))
    return ParityRecoupling(two_λa, two_λb, ηηηphase == 1)
end

function amplitude(cs::ParityRecoupling, (two_λa, two_λb), (two_j, two_ja, two_jb))
    (cs.two_λa == two_λa) * (cs.two_λb == two_λb) && return 1
    (cs.two_λa == -two_λa) * (cs.two_λb == -two_λb) && return 2 * cs.ηηηphaseisplus - 1
    return 0
end

"""
    RecouplingLS(two_ls)

LS-recoupling for a two-body vertex.

Stores `(two_l, two_s)` (i.e. twice the orbital angular momentum and twice the total spin)
and evaluates the corresponding LS/helicity recoupling coefficient via [`jls_coupling`](@ref).

For input helicities ``λ₁, λ₂`` (particles with spins ``j_a, j_b`` coupling to total ``j``), the amplitude is:

```math
A(\\lambda_1, \\lambda_2) = \\sqrt{\\frac{2l+1}{2j+1}} \\;
\\langle j_a, \\lambda_1; j_b, {-}\\lambda_2 \\mid s, \\lambda_1{-}\\lambda_2 \\rangle
\\langle l, 0; s, \\lambda_1{-}\\lambda_2 \\mid j, \\lambda_1{-}\\lambda_2 \\rangle
```

where ``(l, s)`` are the orbital and total spin from `two_ls`, and the kets use the usual
Clebsch–Gordan convention (twice-angular-momentum arguments in the code).
"""
@with_kw struct RecouplingLS <: Recoupling
    two_ls::Tuple{Int,Int}
end

amplitude(cs::RecouplingLS, (two_λa, two_λb), (two_j, two_ja, two_jb)) =
    jls_coupling(two_ja, two_λa, two_jb, two_λb, two_j, cs.two_ls[1], cs.two_ls[2])


"""
    NoFormFactor()

Trivial form-factor functor used by default inside [`Vertex`](@ref).

Calling `NoFormFactor()(...)` always returns 1, i.e. no kinematic suppression.
"""
struct NoFormFactor end

"""
    Vertex{R<:Recoupling,F}

Local vertex payload for a two-body sub-decay: a [`Recoupling`](@ref) scheme and a form-factor functor.

Spin/helicity structure enters through `h` (via [`amplitude`](@ref) on the recoupling type).
Kinematic dependence enters through `ff`, called as `ff(m0², m1², m2²)` with masses of the
three particles involved in that vertex. The default [`NoFormFactor`](@ref) is constant.

# Constructors
```julia
Vertex(h::Recoupling)        # trivial form factor (always 1)
Vertex(h::Recoupling, ff::F) # custom form-factor functor
```

Vertices appear in [`DecayChain`](@ref) as `HRk` (production ``0 \\to Rk``) and `Hij` (decay ``R \\to ij``).

!!! note "Renamed from `VertexFunction`"
    `VertexFunction` is a deprecated alias for `Vertex` and will be removed in a future release.
"""
struct Vertex{R<:Recoupling,F}
    h::R
    ff::F
end
Vertex(h::Recoupling) = Vertex(h, NoFormFactor())

if VERSION >= v"1.11-"
    Base.@deprecate_binding VertexFunction Vertex
else
    const VertexFunction = Vertex
end
(ff::NoFormFactor)(m0², m1², m2²) = one(typeof(m0²))
