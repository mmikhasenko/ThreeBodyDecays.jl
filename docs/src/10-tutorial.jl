#md # # Tutorial
#nb # # ThreeBodyDecays.jl Tutorial
#jl # # ThreeBodyDecays.jl Tutorial

# This tutorial guides the user through the basic example of
# building a model for a cascade decay of a particle into three particles.
# The chosen example is the decay of Λb ⟶ Jψ p K,
# where the pentaquarks candidates has been seen for the first time in 2015.
# The tutorial shows how to build a decay chain object,
# and how to calculate the amplitude, the matrix element of the transition.

using ThreeBodyDecays # import the module
using Plots
using QuadGK
#
theme(
    :wong,
    frame = :box,
    lab = "",
    minorticks = true,
    guidefontvalign = :top,
    guidefonthalign = :right,
    xlim = (:auto, :auto),
    ylim = (:auto, :auto),
    grid = false,
)

# decay Λb ⟶ Jψ p K
constants = Dict("mJψ" => 3.09, "mp" => 0.938, "mK" => 0.49367, "mLb" => 5.62) # masses of the particles

# `ThreeBodySystem` creates an immutable structure that describes the setup.
# Two work with particles with non-integer spin, the doubled quantum numbers are stored.

ms = ThreeBodyMasses(   # masses m1,m2,m3,m0
    constants["mJψ"],
    constants["mp"],
    constants["mK"];
    m0 = constants["mLb"],
)

# create two-body system
tbs = ThreeBodySystem(; ms, two_js = ThreeBodySpins(2, 1, 0; two_h0 = 1)) # twice spin

Conserving = ThreeBodyParities('-', '+', '-'; P0 = '+')
Violating = ThreeBodyParities('-', '+', '-'; P0 = '-')

# - the invariant variables, `σs = [σ₁,σ₂,σ₃]`,
# - helicities `two_λs = [λ₁,λ₂,λ₃,λ₀]`
# - and complex couplings `cs = [c₁,c₂,...]`

# ## Decay chains

# The following code creates six possible decay channels.

# The lineshape of the intermediate resonances is specified
# by the second argument.
# It is a simple Breit-Wigner function in this example.

struct BW
    m::Float64
    Γ::Float64
end
(bw::BW)(σ::Number) = 1 / (bw.m^2 - σ - 1im * bw.m * bw.Γ)


# chains-1, i.e. (2+3): Λs with the lowest ls, LS
Λ1520 = DecayChainLS(
    k = 1,
    Xlineshape = BW(1.5195, 0.0156),
    jp = jp"3/2+",
    Ps = Conserving,
    tbs = tbs,
)
Λ1690 = DecayChainLS(
    k = 1,
    Xlineshape = BW(1.685, 0.050),
    jp = jp"1/2+",
    Ps = Conserving,
    tbs = tbs,
)
Λ1810 = DecayChainLS(
    k = 1,
    Xlineshape = BW(1.80, 0.090),
    jp = jp"5/2+",
    Ps = Conserving,
    tbs = tbs,
)
Λs = (Λ1520, Λ1690, Λ1810);

# chains-3, i.e. (1+2): Pentaquarks with the lowest ls, LS
Pc4312 = DecayChainLS(
    k = 3,
    Xlineshape = BW(4.312, 0.015),
    jp = jp"1/2+",
    Ps = Conserving,
    tbs = tbs,
)
Pc4440 = DecayChainLS(
    k = 3,
    Xlineshape = BW(4.440, 0.010),
    jp = jp"1/2+",
    Ps = Conserving,
    tbs = tbs,
)
Pc4457 = DecayChainLS(
    k = 3,
    Xlineshape = BW(4.457, 0.020),
    jp = jp"3/2+",
    Ps = Conserving,
    tbs = tbs,
)
Pcs = (Pc4312, Pc4440, Pc4457);

# ## Unpolarized intensity

# The full model is the vector of decay chains, and
# couplings for each decay chain.

const model = ThreeBodyDecay(
    ["Λ1520", "Λ1690", "Λ1810", "Pc4312", "Pc4440", "Pc4457"] .=> zip(
        [2, 2.1, 1.4im, 0.4, 0.3im, -0.8im],
        [Λ1520, Λ1690, Λ1810, Pc4312, Pc4440, Pc4457],
    ),
)

# just a random point of the Dalitz Plot
σs = randomPoint(tbs.ms)
two_λs = randomPoint(tbs.two_js)

# Model is build, one can compute unpolarized intensity with it
@show amplitude(model, σs, two_λs)

# gives a real number - probability
@show unpolarized_intensity(model, σs)

# ## Plotting API
#
# A natural way to visualize the three-body decay with two degrees of freedom
# is a correlation plot of the subchannel invariant masses squared.
# Kinematic limits can visualized using the `border` function.
# Plot in the σ₁σ₃ variables is obtained by

plot(
    lab = "",
    grid = false,
    size = (600, 300),
    plot(border31(masses(model)), xlab = "σ₁ (GeV²)", ylab = "σ₃ (GeV²)"),
    plot(border12(masses(model)), xlab = "σ₂ (GeV²)", ylab = "σ₁ (GeV²)"),
)

# The matrix element as a function of the kinematic variables
# can be visualized by passing an amplitude function and the kinematic mass object.

plot(masses(model), σs -> abs2(amplitude(Pc4312, σs, (2, -1, 0, 1))))
plot(masses(model), Base.Fix1(unpolarized_intensity, model); iσx = 1, iσy = 3)

# The projections of the Dalitz Plot can be computed numerically using integration routine, e.g. `QuadGK`.

plot(4.2, 4.6) do e1
    I = Base.Fix1(unpolarized_intensity, model)
    e1 * quadgk(projection_integrand(I, masses(model), e1^2; k = 3), 0, 1)[1]
end
