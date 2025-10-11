# Energy dependence of the amplitude

Any helicity amplitude depends on a set of discrete indices (the spin projections) and energy kinematic variables that determine the energies, momenta, and angles of particles involved in the decay.

While dependence on helicity indices is entirely determined by the quantum numbers, energy dependence is subject to modeling.

The energy dependence for the decay amplitude is separated into three parts:

- **Production form-factor** -- a function of three arguments (m0, σk, mk²)
- **Propagator** -- a function of a single argument (σk)
- **Decay form-factor** -- a function of three arguments (σk, mi², mj²)

where σk = (pi + pk)² = mij² is the invariant mass squared of the intermediate resonance.

When building the model, these dependencies are specified as follows:

```julia
tbs = ThreeBodySystem(
    ThreeBodyMasses(1.0, 2.0, 3.0; m0=16.0),
    ThreeBodySpins(0, 0, 0; two_h0=0))

L = 1  # orbital angular momentum for production
l = 1  # orbital angular momentum for decay

chain1 = DecayChain(;
    k = 1,  # topology (23)1
    jp = jp"1-",
    HRk = VertexFunction(  # production form-factor
        RecouplingLS((2L, 0)),  # (two_l, two_s)
        BlattWeisskopf{L}(1.5)),  # {l}
    Xlineshape = BW(1.5195, 0.0156),  # propagator
    Hij = VertexFunction(  # decay form-factor
        RecouplingLS((2l, 0)),  # (two_l, two_s)
        BlattWeisskopf{l}(1.5)),  # {l}
    tbs)
```

Many common models for the propagators and form factors are available in `HadronicLineshapes.jl`.
Alternatively, you can use lambda functions for simple cases, define any custom objects with dispatch methods as follows:

```julia
# Custom recoupling
struct MyCoupling <: Recoupling
    # custom fields
end
amplitude(cs::MyCoupling, (two_λa, two_λb), (two_j, two_ja, two_jb))  # to defined

# Custom form factor
struct MyFormFactor
    # custom fields
end
(ff::MyFormFactor)(m0², m1², m2²)  # to defined

# Custom propagator
struct MyPropagator
    # custom fields
end
(prop::MyPropagator)(σk)  # to defined
```
