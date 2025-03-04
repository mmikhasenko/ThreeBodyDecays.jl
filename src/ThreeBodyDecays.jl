module ThreeBodyDecays

using StaticArrays
using PartialWaveFunctions
using Parameters
using RecipesBase
using Polynomials
using PolynomialRoots
using Tullio

import Base: getindex, iterate, length, vcat
import PartialWaveFunctions: wignerD_doublearg

# Types
export MassTuple, ThreeBodyMasses, lims, lims1, lims2, lims3
include("mass_types.jl")

export SpinTuple, ThreeBodySpins, SpinParity, ⊗, str2jp, @jp_str, x2, d2
include("spin_types.jl")

export ParityTuple, ThreeBodyParities, ThreeBodySpinParities
include("parity_types.jl")

export ThreeBodySystem, DalitzPlotPoint, masses, spins, randomPoint
include("system_types.jl")

# Kinematics
export MandelstamTuple, Invariants, x2σs, y2σs, circleorigin
include("mandelstam.jl")

export Kallen, sqrtKallenFact, Kibble, inphrange, isphysical
export phase, ijk, ij_from_k
export cosθij, cosθ12, cosθ23, cosθ31
export σiofk, σjofk, σ1of2, σ2of3, σ3of1, σ1of3, σ2of1, σ3of2
export breakup, breakup_Rk, breakup_ij, aligned_four_vectors
include("phase_space.jl")

export border31, border12, border23, border13, border21, border32, border
include("dalitz.jl")

# Angular Momentum
export ispositive, TrivialWignerRotation, wr, cosζ
export cosζ12_for0, cosζ23_for0, cosζ31_for0
export cosζ21_for1, cosζ21_for2
export cosζ13_for1, cosζ13_for3
export cosζ32_for3, cosζ32_for2
export cosζ12_for3, cosζ23_for1, cosζ31_for2
export cosζk1_for1, cosζk2_for2, cosζk3_for3
include("wigner_rotations.jl")

# No exports from this file
include("wigner_d_matrix.jl")

# Utilities
export minusone, x_str, shift_by_half
include("utils.jl")

export letterL
export possible_ls, possible_ls_ij, possible_ls_Rk, possible_lsLS, possible_l_s_L_S
export possible_coupling_schemes, complete_l_s_L_S
include("coupling_schemes.jl")

# Decay Models
export Recoupling, NoRecoupling, ParityRecoupling, RecouplingLS
include("recouplings.jl")

export AbstractDecayChain, DecayChain, DecayChainLS, DecayChainsLS
export amplitude, itr, summed_over_polarization, system
export unpolarized_intensity
include("decay_channel.jl")

export ThreeBodyDecay
include("decay_model.jl")

export change_basis_3from1, change_basis_1from2, change_basis_2from3
include("cross_channel.jl")

# Integration
export phase_space_integrand, projection_integrand
include("integrand.jl")

export DalitzPlot # just a docstring

end  # module ThreeBodyDecays
