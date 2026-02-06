---
title: 'ThreeBodyDecays.jl: Dalitz-plot decomposition for three-body decays in Julia'
tags:
  - Julia
  - particle physics
  - amplitude analysis
  - Dalitz plot
  - helicity formalism
  - Wigner rotations
authors:
  - name: Mikhail Mikhasenko
    orcid: 0000-0002-6969-2063
    affiliation: 1
affiliations:
  - name: Ruhr University Bochum, Germany
    index: 1
    ror: 04tsk2644
date: 7 February 2026
bibliography: paper.bib
---

# Summary

`ThreeBodyDecays.jl` is a Julia package for constructing helicity amplitudes
of three-body particle decays using the Dalitz-plot decomposition (DPD) formalism
introduced in @Mikhasenko:2019xse.
The package factorizes the full decay amplitude into an overall rotation---encoding the
polarization of the parent particle---and a *Dalitz-plot function* that depends only on two
Lorentz-invariant variables (the subchannel invariant masses).
This separation is achieved by expressing all angular dependence through
Wigner $d$-functions with arguments given as explicit, real-valued functions of the
Mandelstam variables $\sigma_1$, $\sigma_2$, $\sigma_3$.
The Wigner rotation angles that match spin states across different isobar frames are
computed using the algorithm generalized in @Habermann:2024sxs.
`ThreeBodyDecays.jl` handles particles with arbitrary spin, supports both helicity and
$LS$-coupling parametrizations, and leaves the isobar lineshape as a user-supplied function,
making the framework lineshape-agnostic.
Written in Julia [@Bezanson:2017julia], the package leverages multiple dispatch,
static arrays, and composability with the broader Julia ecosystem for
automatic differentiation, optimization, and statistical inference.

# Statement of need

Amplitude analysis of multi-body hadronic decays is the primary tool for
discovering and characterizing exotic hadron states such as pentaquarks
[@Aaij:2015tga; @Aaij:2019vzc] and tetraquark candidates [@Karliner:2017qhf; @Guo:2017jvc],
as well as for precision measurements of CKM matrix elements and CP violation.
In the isobar model, a three-body decay $0 \to 1\,2\,3$ is described as a sum of
cascade processes $0 \to R\,k \to (ij)\,k$, each producing an intermediate resonance $R$ in
one of the three subchannels.

A persistent difficulty in implementing the helicity formalism for these decays is the
correct treatment of Wigner rotations: when two decay chains share a final-state particle
with spin, the helicity quantization axes defined in different isobar rest frames must be
matched via a Lorentz transformation that induces a rotation of the spin states.
These Wigner rotations are a known source of subtle implementation errors, as
discussed in detail by @Habermann:2024sxs, @Mikhasenko:2017rkh, and @Pilloni:2018kwm.
The DPD formalism eliminates this difficulty by working in an *aligned* configuration where
all Wigner rotations reduce to rotations about the $y$-axis, whose angles are given as closed-form
expressions in the Mandelstam variables.

No existing registered package in the Julia ecosystem provides a type-safe implementation of the
DPD formalism for arbitrary spins.
`ThreeBodyDecays.jl` fills this gap, targeting experimental and phenomenological physicists
at facilities such as LHCb, BESIII, Belle II, and GlueX, as well as theory groups in the
JPAC collaboration.
The package also serves as a building block for extensions to general multi-body cascade
decays; the companion package `DecayAngle` [@Habermann:2024decayangle] generalizes the
Wigner rotation computation to arbitrary decay topologies.

# State of the field

A comprehensive inventory of amplitude-analysis software is maintained at @PWAPages.
Established frameworks include:
**Laura++** [@Back:2017zqt], a C++ Dalitz-plot fitter used extensively by LHCb for
$B$-meson analyses;
**AmpGen** [@AmpGen], which generates C++ amplitude code symbolically for fast evaluation;
**ComPWA** [@ComPWA], a Python-based project (QRules + AmpForm + TensorWaves) offering
symbolic model building with backend-agnostic numerical evaluation;
**AmpTools** [@AmpTools], a C++ framework with GPU acceleration used by GlueX and BESIII;
and **TF-PWA** [@TFPWA], a TensorFlow-based tool focused on fitting performance.

`ThreeBodyDecays.jl` occupies a distinct niche among these tools.
First, it is the sole Julia entry in the PWA software inventory, making it natively composable with
Julia's ecosystem for automatic differentiation, Bayesian inference, and GPU computing.
Second, rather than constructing amplitudes from generic Lorentz-covariant or helicity building
blocks, it directly implements the DPD master equation (Eq.\ 8 of @Mikhasenko:2019xse),
in which the only model-dependent input is the isobar lineshape function.
This formalism-native design means the mathematical structure of the paper maps one-to-one onto the code,
reducing the risk of transcription errors and simplifying validation.
Third, the twice-integer spin convention and static-array storage enable both
exact angular-momentum arithmetic and stack-allocated, cache-friendly data layouts.

# Software design

The architecture of `ThreeBodyDecays.jl` is organized in layers that mirror the physics, as
shown in \autoref{fig:decay}.

![Diagram of the three-body decay $0 \to 1\,2\,3$ through an isobar $R$ in the $(ij)$ subchannel and a spectator $k$. This topology is the basic building block of `DecayChain`.\label{fig:decay}](figures/general_decay.pdf){ width=40% }

**Type system.**
Particle masses, spins (in the twice-integer convention $2j$), and parities are stored in named tuples
(`MassTuple`, `SpinTuple`, `ParityTuple`) collected into a `ThreeBodySystem`.
Kinematic configurations are represented by `MandelstamTuple` holding $(\sigma_1, \sigma_2, \sigma_3)$,
with the physical region enforced via the Kibble function.

**Coupling and recoupling.**
Allowed $(l,s)$ and $(L,S)$ quantum numbers for each vertex are enumerated by `possible_lsLS`,
which enforces angular-momentum and parity conservation.
Each vertex is wrapped in a `VertexFunction` that pairs a `Recoupling` (helicity or $LS$-coupling
coefficient) with an optional form factor.

**Aligned amplitude.**
The core computation lives in `aligned_amplitude`, which evaluates the product

$$F_0 = \mathcal{N}_J \; V_{0 \to Rk} \; d^s(\theta_{ij}) \; V_{R \to ij} \; \mathcal{L}(\sigma_k), \label{eq:aligned}$$

where $d^s$ is a Wigner small-$d$ matrix and $\mathcal{L}(\sigma_k) = X_s(\sigma_k)\,B_L'\,B_{l'}'$ is the
user-supplied lineshape times Blatt--Weisskopf form factors.
The full helicity amplitude is then

$$A = d^{j_0}(\zeta_0) \; F_0 \; d^{j_1}(\zeta_1)\, d^{j_2}(\zeta_2)\, d^{j_3}(\zeta_3), \label{eq:full}$$

with the Wigner rotation angles $\zeta_i$ computed from the Mandelstam variables by the
`cosζ` family of functions.
This structure---**d (parent) $\times$ (V$\cdot$d$\cdot$V) $\times$ d$\cdot$d$\cdot$d (daughters)**---directly
mirrors Eq.\ 8 of @Mikhasenko:2019xse.

![Three aligned center-of-momentum configurations corresponding to the three spectator choices $(ijk) \in \{(123),(231),(312)\}$. In the DPD, all three are related by known rotations, eliminating complex phases from Wigner rotations.\label{fig:aligned}](figures/decomposition_comparison.pdf){ width=90% }

**Decay model.**
A coherent sum of chains with complex couplings is represented by `ThreeBodyDecay`.
Computing the unpolarized intensity amounts to summing $|A|^2$ over all helicity configurations,
which the function `unpolarized_intensity` performs efficiently using `Tullio.jl` for
Einstein-summation contractions.

**Key design trade-offs.**
The lineshape is kept as an arbitrary callable ($\sigma \to \mathbb{C}$) rather than a fixed
parametrization, so users can supply Breit--Wigner, $K$-matrix, or dispersive representations without
modifying the framework.
`StaticArrays` are used for the helicity-index arrays, enabling stack allocation and
compile-time size inference.
Plotting is provided through `RecipesBase`, which defines Dalitz-plot and projection recipes
without introducing a hard dependency on any particular plotting backend.

# Research impact statement

The DPD formalism paper @Mikhasenko:2019xse has accumulated over 60 citations on INSPIRE-HEP.
The formalism has been adopted in LHCb amplitude analyses, including the
recent five-body cascade analysis of $B^+ \to \psi(2S)\,K^+\pi^+\pi^-$ [@LHCb:2024amr].
A dedicated follow-up paper on Wigner rotations for general cascade reactions
[@Habermann:2024sxs] directly cites `ThreeBodyDecays.jl` (via its Zenodo archive)
and extends the DPD approach to arbitrary multi-body topologies, demonstrating
the package's role as a reference implementation.
The companion package `DecayAngle` [@Habermann:2024decayangle] was developed in tandem
for computing Wigner rotation angles in general decay graphs.

`ThreeBodyDecays.jl` is listed in the community PWA software inventory [@PWAPages]
as the sole Julia entry among approximately twenty amplitude-analysis frameworks worldwide.
The JPAC collaboration maintains a complementary C++ repository
(`DalitzPlotDecomposition`) that shares the same formalism.
The package has been used in CERN Summer Student workshops on amplitude analysis
and in university courses on hadron spectroscopy at Ruhr University Bochum.

# AI usage disclosure

Generative AI tools (Cursor with Claude) were used to assist in drafting and editing
this manuscript.
All content was reviewed, verified, and revised by the author.
No AI tools were used in the development of the `ThreeBodyDecays.jl` software itself.

# Acknowledgements

The author thanks Kai Habermann for the collaboration on the Wigner rotation formalism
and the `DecayAngle` package, and the members of the JPAC collaboration---in particular
M. Albaladejo, A. Pilloni, V. Mathieu, and A. P. Szczepaniak---for
the development of the underlying theoretical framework.
This work is supported by the Deutsche Forschungsgemeinschaft (DFG)
through the Research Unit FOR 2926 "Next Generation Perturbative QCD
for Hadron Structure: Preparing for the Electron-Ion Collider."

# References
