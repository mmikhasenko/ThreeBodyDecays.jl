abstract type Recoupling end
@with_kw struct NoRecoupling <: Recoupling
    two_λa::Int
    two_λb::Int
end

amplitude(cs::NoRecoupling, (two_λa, two_λb), (two_j, two_ja, two_jb)) =
    (cs.two_λa == two_λa) * (cs.two_λb == two_λb)

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

@with_kw struct RecouplingLS <: Recoupling
    two_ls::Tuple{Int,Int}
end

amplitude(cs::RecouplingLS, (two_λa, two_λb), (two_j, two_ja, two_jb)) =
    jls_coupling(two_ja, two_λa, two_jb, two_λb, two_j, cs.two_ls[1], cs.two_ls[2])


struct NoFormFactor end

"""
    EnergyDependentFormFactor{F}

A form factor that depends on energy (invariant mass) and particle masses.
The form factor function `f` should accept:
- `energy`: The invariant mass squared (energy squared)
- `masses`: A `MassTuple` containing the particle masses

# Example
```julia
# Blatt-Weisskopf form factor
ff = EnergyDependentFormFactor((energy, masses) -> 1.0 / (1.0 + energy / masses.m0^2))
```
"""
struct EnergyDependentFormFactor{F}
    f::F
end

"""
    MassDependentFormFactor{F}

A form factor that depends only on particle masses (not energy).
The form factor function `f` should accept:
- `masses`: A `MassTuple` containing the particle masses

# Example
```julia
# Mass-dependent coupling
ff = MassDependentFormFactor((masses) -> masses.m0 / (masses.m1 + masses.m2))
```
"""
struct MassDependentFormFactor{F}
    f::F
end

"""
    VertexFunction{R<:Recoupling,F}

A struct that contains a recoupling and a form factor.
Two constructors:
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
amplitude(V::VertexFunction{R,NoFormFactor}, args...) where {R} = amplitude(V.h, args...)

# Form factor amplitude methods
function amplitude(V::VertexFunction{R,EnergyDependentFormFactor}, 
                  (two_λa, two_λb), (two_j, two_ja, two_jb), energy, masses) where {R}
    base_amplitude = amplitude(V.h, (two_λa, two_λb), (two_j, two_ja, two_jb))
    form_factor = V.ff.f(energy, masses)
    return base_amplitude * form_factor
end

function amplitude(V::VertexFunction{R,MassDependentFormFactor}, 
                  (two_λa, two_λb), (two_j, two_ja, two_jb), masses) where {R}
    base_amplitude = amplitude(V.h, (two_λa, two_λb), (two_j, two_ja, two_jb))
    form_factor = V.ff.f(masses)
    return base_amplitude * form_factor
end
