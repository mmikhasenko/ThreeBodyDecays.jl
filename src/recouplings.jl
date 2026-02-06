"""
    Recoupling

Abstract supertype for recoupling schemes used inside [`VertexFunction`](@ref).

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
"""
@with_kw struct RecouplingLS <: Recoupling
    two_ls::Tuple{Int,Int}
end

amplitude(cs::RecouplingLS, (two_λa, two_λb), (two_j, two_ja, two_jb)) =
    jls_coupling(two_ja, two_λa, two_jb, two_λb, two_j, cs.two_ls[1], cs.two_ls[2])


"""
    NoFormFactor()

Trivial form-factor functor used by default inside [`VertexFunction`](@ref).

Calling `NoFormFactor()(...)` always returns 1, i.e. no kinematic suppression.
"""
struct NoFormFactor end

"""
    VertexFunction{R<:Recoupling,F}

A struct that contains a recoupling and a form factor.
There are two constructors:
```julia
VertexFunction(h::Recoupling)        # translates to a trivial form factor
VertexFunction(h::Recoupling, ff::F) # with a form factor
# amplitude(V::VertexFunction{<:Recoupling, yourF}, args...) where {yourF} needs to be defined
```

"""
struct VertexFunction{R<:Recoupling,F}
    h::R
    ff::F
end
VertexFunction(h::Recoupling) = VertexFunction(h, NoFormFactor())
(ff::NoFormFactor)(m0², m1², m2²) = one(typeof(m0²))
