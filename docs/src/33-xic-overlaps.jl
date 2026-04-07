#md # # Ξc⁺ → Ξ⁻π⁺π⁺ overlaps
#nb # # Ξc⁺ → Ξ⁻π⁺π⁺ overlaps
#jl # # Ξc⁺ → Ξ⁻π⁺π⁺ overlaps

# This tutorial builds a compact model for
# Ξc⁺ → Ξ⁻ π⁺ π⁺ with three prominent Ξ* resonances:
# Ξ(1530), Ξ(1620), and Ξ(1690).
#
# The goal is not a tuned physics model. Instead, we use the model as a
# transparent object: we inspect its LS content, build an integral
# representation of the amplitude overlaps, and then reduce the result to
# quantities that are often discussed in amplitude analyses.

using ThreeBodyDecays
using HadronicLineshapes: BreitWigner
using LinearAlgebra
using Plots
using Random
using Statistics

theme(:wong, frame = :box, grid = false);

# ## Kinematics and quantum numbers
#
# Particle labels are fixed by the mass tuple:
#
# - particle 1: Ξ⁻
# - particle 2: π⁺
# - particle 3: π⁺
#
# Therefore σ₂ is the invariant mass squared of particles (3,1), and σ₃
# is the invariant mass squared of particles (1,2). These are the two
# Ξ⁻π⁺ combinations. The pions are identical in charge, so a realistic
# model should include both pairings.

const mΞ = 1.32171
const mπ = 0.13957039
const mΞc = 2.46771

ms = ThreeBodyMasses(mΞ, mπ, mπ; m0 = mΞc);

two_js, parities_PC = ThreeBodySpinParities("1/2+", "0-", "0-"; jp0 = "1/2+");
_, parities_PV = ThreeBodySpinParities("1/2+", "0-", "0-"; jp0 = "1/2-");

tbs = ThreeBodySystem(ms, two_js);

# The `parities_PC` and `parities_PV` objects are a small bookkeeping trick.
# The weak decay does not force a single parity-conserving LS family, so we
# include both. This gives two production LS couplings per resonance and per
# pion pairing in this example.

const resonances = [
    (; name = "Ξ(1530)", jp = "3/2+", mass = 1.5318, width = 0.0091),
    (; name = "Ξ(1620)", jp = "1/2-", mass = 1.6200, width = 0.0320),
    (; name = "Ξ(1690)", jp = "1/2-", mass = 1.6900, width = 0.0300),
];

# ## Building the model
#
# `HadronicLineshapes.BreitWigner` can be passed directly as the
# energy-dependence of a `DecayChainsLS` object. The rest is a small map over
# resonances, the two Ξπ pairings, and the two parity choices.

function chain_records(resonance)
    map((2, 3)) do k
        map((("PC", parities_PC), ("PV", parities_PV))) do (parity_label, Ps)
            chains = DecayChainsLS(;
                k,
                Xlineshape = BreitWigner(resonance.mass, resonance.width),
                jp = resonance.jp,
                Ps,
                tbs,
            )
            qns = possible_l_s_L_S(SpinParity(resonance.jp), two_js, Ps; k)
            map(zip(chains, qns)) do (chain, qn)
                (;
                    name = resonance.name,
                    parity = parity_label,
                    k,
                    label = "$(resonance.name) k=$k $parity_label",
                    qn,
                    chain,
                )
            end
        end
    end |>
    Iterators.flatten |>
    Iterators.flatten |>
    collect
end;

records = vcat(chain_records.(resonances)...);

# A first validation step is simply to ask what the model contains. The LS
# labels use ordinary angular momenta, while the implementation stores doubled
# quantum numbers internally.

[(r.label, r.qn) for r in records]

# Now attach a reproducible set of complex couplings. In a fit these would be
# floating parameters. Bose symmetry constrains the two Ξπ pairings: the
# coefficient for a `k=2` chain and the matching `k=3` chain must be the same
# up to the exchange sign. For two π⁺ bosons we use the symmetric sign.

coupling_by_resonance_and_parity = Dict(
    ("Ξ(1530)", "PC") => 1.00 + 0.00im,
    ("Ξ(1530)", "PV") => 0.35 - 0.20im,
    ("Ξ(1620)", "PC") => 0.75 - 0.15im,
    ("Ξ(1620)", "PV") => -0.25 + 0.30im,
    ("Ξ(1690)", "PC") => 0.55 - 0.25im,
    ("Ξ(1690)", "PV") => -0.20 - 0.20im,
);

exchange_sign = +1;
couplings = [
    (r.k == 2 ? 1 : exchange_sign) * coupling_by_resonance_and_parity[(r.name, r.parity)] for r in records
];

[(r.label, couplings[i]) for (i, r) in pairs(records)]

model = ThreeBodyDecay(
    getproperty.(records, :label) .=> zip(couplings, getproperty.(records, :chain)),
);

# A pointwise sanity check is cheap and helpful: the unpolarized intensity is
# real and positive.

σs0 = randomPoint(ms);
unpolarized_intensity(model, σs0)

# ## Phase-space sample
#
# The hit-and-miss step is deliberately simple: draw random points in a square,
# map them to Mandelstam variables with `y2σs`, and keep only the points inside
# the physical phase space.

function random_mandelstam_candidates(ms, nsamples; seed = 12)
    rng = MersenneTwister(seed)
    map(1:nsamples) do _
        y2σs(rand(rng, 2), ms)
    end
end;

function physical_samples(ms, nsamples = 3_000; seed = 12)
    filter(σs -> isphysical(σs, ms), random_mandelstam_candidates(ms, nsamples; seed))
end;

dalitz_candidates = random_mandelstam_candidates(ms, 600; seed = 31);
dalitz_hits = filter(σs -> isphysical(σs, ms), dalitz_candidates);

scatter(
    getindex.(dalitz_candidates, 2),
    getindex.(dalitz_candidates, 3);
    lab = "rejected",
    c = :lightgray,
    alpha = 0.35,
    markerstrokewidth = 0,
    markersize = 2,
    aspect_ratio = 1,
    xlab = "σ₂ = m²(Ξ⁻π⁺) [GeV²]",
    ylab = "σ₃ = m²(Ξ⁻π⁺) [GeV²]",
    title = "Hit-and-miss phase-space sample",
)
scatter!(
    getindex.(dalitz_hits, 2),
    getindex.(dalitz_hits, 3);
    lab = "accepted",
    c = 2,
    alpha = 0.7,
    markerstrokewidth = 0,
    markersize = 2.5,
)
plot!(border23(ms); lab = "", c = :black, lw = 2)

# ## Integral representation of overlaps
#
# Write the model as a coherent sum
#
# ```math
# A(\sigma) = \sum_i c_i\, a_i(\sigma).
# ```
#
# The integrated intensity is then
#
# ```math
# I = \sum_{ij} c_i^*\, \rho_{ij}\, c_j = c^\dagger \rho c,
# \qquad
# \rho_{ij} = \int d\Phi_3\, a_i(\sigma)^* a_j(\sigma).
# ```
#
# The overlap API works with event samples. We generate a small reproducible
# phase-space sample, cache the stripped chain amplitudes, and use
# `event_overlap_contributions` to build the spin-summed overlap objects.

samples = physical_samples(ms);
cache = chain_amplitudes(model, samples);

ρ_overlap = mean(event_overlap_contributions(cache));
ρ = ρ_overlap.matrix;

I_from_matrix = real(model.couplings' * Matrix(ρ) * model.couplings)

# The same object can be made physical by applying the coefficients. The
# package then knows how to sum diagonal and interference pieces into the
# total sample intensity.

physical_ρ = physical_overlap(ρ_overlap, model.couplings);
total_intensity(physical_ρ)

# As a validation, this agrees with direct evaluation of the model over the
# same sample. The equality is exact up to floating-point summation order.

I_direct = mean(σs -> unpolarized_intensity(model, σs), samples)

isapprox(I_from_matrix, I_direct; rtol = 1e-12)

# ## Normalized overlap matrix
#
# The normalized overlap matrix removes the individual norms:
#
# ```math
# \tilde{\rho}_{ij} =
# \frac{\rho_{ij}}{\sqrt{\rho_{ii}\rho_{jj}}}.
# ```
#
# If one of its eigenvalues is close to zero, the basis contains a near-null
# direction and the coefficients can become ambiguous. If the eigenvalues are
# comfortably away from zero, this particular basis is numerically healthy.

ρdiag = real.(diag(ρ));
ρtilde = Hermitian(Matrix(ρ) ./ sqrt.(ρdiag * ρdiag'));

round.(real.(eigvals(ρtilde)); digits = 4)

# The diagonal elements are one by construction. Looking at the largest
# off-diagonal entries tells us which amplitudes overlap most strongly after
# normalization.

function largest_overlaps(ρtilde, labels; n = 6)
    pairs = [
        (absρ = abs(ρtilde[i, j]), ρij = ρtilde[i, j], i, j) for i in axes(ρtilde, 1)
        for j = 1:(i-1)
    ]
    sort!(pairs, by = x -> x.absρ, rev = true)
    return [
        (labels[p.i], labels[p.j], round(p.ρij; digits = 3)) for
        p in pairs[1:min(n, length(pairs))]
    ]
end;

largest_overlaps(ρtilde, model.names)

# ## Fit-fraction and interference matrix
#
# A real fit-fraction/interference matrix can be defined by absorbing the
# couplings into the overlap matrix:
#
# ```math
# N_{ij} = \operatorname{Re}\left[c_i^*\, \rho_{ij}\, c_j\right],
# \qquad
# I = \mathbf{1}^T N \mathbf{1}.
# ```
#
# The API does this with `physical_overlap`: its matrix is the physical
# numerator matrix in the chain basis. We only take `real` and divide by the
# total when we want a conventional fraction matrix.

N = real.(Matrix(physical_ρ.matrix));

total_intensity(physical_ρ), I_from_matrix

Nfrac = N ./ I_from_matrix;

[(model.names[i], round(Nfrac[i, i]; digits = 4)) for i in eachindex(model.names)]

fit_fractions(physical_ρ; percent = false, sort = false)

# The diagonal entries do not have to sum to one: coherent interference can
# add or remove intensity. That is a useful sense check rather than a bug.

sum(diag(Nfrac)), sum(Nfrac)

# ## Reduced matrix by resonance
#
# Often we want a resonance-level summary instead of an LS-amplitude-level
# summary. `group_overlap` performs the partial sum over all entries that
# share a resonance name. That includes both Ξπ pairings and both parity
# structures for each resonance.

resonance_ρ = group_overlap(physical_ρ, getproperty.(records, :name));
resonance_names = resonance_ρ.labels;
Nres = real.(Matrix(resonance_ρ.matrix)) ./ total_intensity(resonance_ρ);

round.(Nres; digits = 4)

Nres_plot = round.(Nres; digits = 3);
nres = length(resonance_names);
pNres = heatmap(
    1:nres,
    1:nres,
    Nres_plot;
    xticks = (1:nres, resonance_names),
    yticks = (1:nres, resonance_names),
    xrotation = 30,
    aspect_ratio = 1,
    c = :viridis,
    xlab = "resonance j",
    ylab = "resonance i",
    title = "Reduced fit-fraction/interference matrix",
    colorbar_title = "Nᵢⱼ",
)
for i in axes(Nres_plot, 1), j in axes(Nres_plot, 2)
    annotate!(pNres, j, i, text(string(Nres_plot[i, j]), 8, :white))
end
pNres

# The diagonal elements are the resonance fit fractions in this reduced basis.

[(resonance_names[i], round(Nres[i, i]; digits = 4)) for i in eachindex(resonance_names)]

fit_fractions(resonance_ρ; percent = false, sort = false)

interference_terms(resonance_ρ; percent = false, sort = false)

# This is also easy to validate: the resonance-level matrix still sums to the
# same total fraction, because reduction is just a partial sum over indices.

sum(Nres), sum(Nfrac)

# ## Quiz
#
# **Q1.** Why do we include both `k=2` and `k=3` chains for the same Ξ* resonance?
#
# ```@raw html
# <details>
# <summary>Answer</summary>
#
# The two π⁺ particles are identical in charge. Either pion can pair with the
# Ξ⁻ to form the Ξ* isobar, so the model includes both Ξ⁻π⁺ combinations.
#
# </details>
# ```
#
# **Q2.** What does a very small eigenvalue of `ρtilde` mean?
#
# ```@raw html
# <details>
# <summary>Answer</summary>
#
# It means that some linear combination of amplitudes has a tiny integrated
# norm. The corresponding couplings can be weakly determined or ambiguous.
#
# </details>
# ```
#
# **Q3.** Why are the `k=2` and `k=3` couplings tied together?
#
# ```@raw html
# <details>
# <summary>Answer</summary>
#
# The two final-state pions are identical bosons. Exchanging them maps one
# Ξ⁻π⁺ pairing into the other, so the paired amplitudes are not independent.
# In this symmetric example their coefficients are equal; in a convention with
# an exchange phase, they would differ only by that fixed sign.
#
# </details>
# ```
#
# **Q4.** Why can the sum of diagonal fit fractions be different from one?
#
# ```@raw html
# <details>
# <summary>Answer</summary>
#
# The diagonal terms are incoherent contributions. The total intensity also
# contains off-diagonal interference terms, which can be positive or negative.
#
# </details>
# ```
#
# **Q5.** Which API choice makes the overlap calculation transparent in this page?
#
# ```@raw html
# <details>
# <summary>Answer</summary>
#
# The workflow is split into visible stages: `chain_amplitudes` caches the
# stripped chain amplitudes, `event_overlap_contributions` or
# `chain_overlap_matrix` builds the chain-basis overlaps, `physical_overlap`
# applies couplings, and `group_overlap` reduces the basis. Each object still
# exposes labels and a matrix, so it is easy to inspect the calculation.
#
# </details>
# ```
