#md # # Tutorial
#nb # # ThreeBodyDecays.jl Tutorial
#jl # # ThreeBodyDecays.jl Tutorial

# This tutorial guides the user through the basic example of
# building a model for a cascade decay of a picticle into three particles.
# The chosen example is the decay of Λb ⟶ Jψ p K,
# where the pentaquarks candidates has been seen for the first time in 2015.
# The tutorial shows how to build a decay chain object,
# and how to calculate the amplitude, the matrix element of the transition.

using ThreeBodyDecays # import the module
using Plots

# decay Λb ⟶ Jψ p K
ms = (Jψ=3.09, p=0.938, K=0.49367, Lb=5.62) # masses of the particles

# `ThreeBodySystem` creates an immutable structure that describes the setup.
# Two work with particles with non-integer spin, the doubled quantum numbers are stored.

# create two-body system
tbs = ThreeBodySystem(ms.Jψ, ms.p, ms.K; m0=ms.Lb,   # masses m1,m2,m3,m0
    two_js=ThreeBodySpins(2, 1, 0; two_h0=1)) # twice spin

Concerving = ThreeBodyParities('-', '+', '-'; P0='+')
Violating = ThreeBodyParities('-', '+', '-'; P0='-')

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
Λ1520 = DecayChainLS(1, BW(1.5195, 0.0156); two_j=3 / 2 |> x2, parity='+', Ps=Concerving, tbs=tbs)
Λ1690 = DecayChainLS(1, BW(1.685, 0.050); two_j=1 / 2 |> x2, parity='+', Ps=Concerving, tbs=tbs)
Λ1810 = DecayChainLS(1, BW(1.80, 0.090); two_j=5 / 2 |> x2, parity='+', Ps=Concerving, tbs=tbs)
Λs = (Λ1520, Λ1690, Λ1810)

# chains-3, i.e. (1+2): Pentaquarks with the lowest ls, LS
Pc4312 = DecayChainLS(3, BW(4.312, 0.015); two_j=1 / 2 |> x2, parity='+', Ps=Concerving, tbs=tbs)
Pc4440 = DecayChainLS(3, BW(4.440, 0.010); two_j=1 / 2 |> x2, parity='+', Ps=Concerving, tbs=tbs)
Pc4457 = DecayChainLS(3, BW(4.457, 0.020); two_j=3 / 2 |> x2, parity='+', Ps=Concerving, tbs=tbs)
Pcs = (Pc4312, Pc4440, Pc4457)

# ## Unpolarized intensity

# The full model is the vector of decay chains, and 
# couplings for each decay chain. 
struct Model
    chains::Vector{DecayChain}
    couplings::Vector{Complex}
end

model = Model(
    [Λ1520, Λ1690, Λ1810, Pc4312, Pc4440, Pc4457],
    [1, 1.1, 0.4im, 2.2, 2.1im, -0.3im])

# Amplitudes for the decay chains are added coherently with complex constants.

# The intensity (and probability) is a squared amplitude summed over
# the summed over helicities for the case the decay particle is unpolarized.

I(model::Model, σs) = sum(abs2,
    transpose(model.couplings) * amplitude.(model.chains, Ref(σs), Ref(two_λs))
    for two_λs in itr(tbs.two_js))

# just a random point of the Dalitz Plot
σs = randomPoint(tbs.ms) 

# gives a real number - probability
@show I(model, σs)

# ## Plotting API
# 
# A natural way to visualize the three-body decay with two degrees of freedom
# is a correlation plot of the subchannel invariant masses squared.
# Kinematic limits can visualized using the `border` function.
# Plot in the σ₁σ₃ variables is obtained by

plot(lab="", grid=false, size=(600, 300),
    plot(border31(tbs.ms), xlab="σ₁ (GeV²)", ylab="σ₃ (GeV²)"),
    plot(border12(tbs.ms), xlab="σ₂ (GeV²)", ylab="σ₁ (GeV²)"))

# The matrix element as a function of the kinematic variables
# can be visualized by passing an amplitude function and the kinematic mass object.

plot(tbs.ms, σs -> abs2(amplitude(Pc4312, σs, (2, -1, 0, 1))))
plot(tbs.ms, Base.Fix1(I, model))