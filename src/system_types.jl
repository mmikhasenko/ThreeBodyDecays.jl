"""
    ThreeBodySystem{T,K}

A structure representing a three-body system with masses and spins.

# Fields
- `ms::T`: Masses of the system (MassTuple)
- `two_js::K`: Spins of the system (SpinTuple)
"""
@with_kw struct ThreeBodySystem{T,K}
    ms::T
    two_js::K = ThreeBodySpins(0, 0, 0; two_h0 = 0)
end

# convenient constructors
ThreeBodySystem(ms::MassTuple) = ThreeBodySystem(ms = ms)
ThreeBodySystem(m1, m2, m3; m0, two_js = ThreeBodySpins(0, 0, 0; two_h0 = 0)) =
    ThreeBodySystem(ThreeBodyMasses(m1, m2, m3; m0 = m0), two_js)

masses(tbs::ThreeBodySystem) = tbs.ms
spins(tbs::ThreeBodySystem) = tbs.two_js

"""
    DalitzPlotPoint{I, S}

A structure representing a point in the Dalitz plot.

# Fields
- `σs::I`: Mandelstam variables
- `two_λs::S`: Spin configuration
"""
@with_kw struct DalitzPlotPoint{I,S}
    σs::I
    two_λs::S
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

function randomPoint(tbs::ThreeBodySystem)
    DalitzPlotPoint(σs = randomPoint(tbs.ms), two_λs = randomPoint(tbs.two_js))
end
