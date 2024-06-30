module ThreeBodyDecays

using StaticArrays
using PartialWaveFunctions
using Parameters
using RecipesBase
using Polynomials
using PolynomialRoots

import Base: getindex, iterate, length, vcat


export Kallen, sqrtKallenFact, Kibble,
    inphrange,
    lims, lims1, lims2, lims3
export phase
export ijk, ij_from_k
# 
export cosθij
export cosθ12, cosθ23, cosθ31
# 
export σiofk, σjofk
export σ1of2, σ2of3, σ3of1,
    σ1of3, σ2of1, σ3of2
export breakup, breakup_Rk, breakup_ij
include("kinematics.jl")

export x_str, minusone, shift_by_half
include("utils.jl")

export ispositive
export TriavialWignerRotation
export wr
export cosζ
# 
export cosζ12_for0, cosζ23_for0, cosζ31_for0
# 
export cosζ21_for1, cosζ21_for2
export cosζ13_for1, cosζ13_for3
export cosζ32_for3, cosζ32_for2
# 
export cosζ12_for3, cosζ23_for1, cosζ31_for2
# 
export cosθhatk1, cosθhatk2, cosθhatk3
export cosζk1_for1, cosζk2_for2, cosζk3_for3
# 
include("wignerrotations.jl")


export str2jp
export @jp_str
export ThreeBodyMasses, ThreeBodySpins, ThreeBodyParities
export SpinParity, ThreeBodySpinParities
export MandestamTuple, SpinTuple, ParityTuple
export ThreeBodySystem
export DalitzPlotPoint
export Invariants
export x2σs, y2σs
export randomPoint
export nt, x2, d2
export possible_helicities
export border31, border12, border23
export border13, border21, border32
export border
export circleorigin
export unpolarized_intensity
include("tbs_struct.jl")

export letterL
export ⊗
export possible_ls, possible_ls_ij, possible_ls_Rk, possible_lsLS
export possible_coupling_schemes
include("coupling_scheme.jl")

export change_basis_3from1,
    change_basis_1from2,
    change_basis_2from3
include("cross_channel_relations.jl")

export Recoupling
export NoRecoupling, ParityRecoupling, RecouplingLS
export AbstractDecayChain
export DecayChain, DecayChainLS, DecayChainsLS
export amplitude
export itr, summed_over_polarization
include("decay_channel.jl")

export masses, spins, system
export ThreeBodyDecay
include("model.jl")

export phase_space_integrand
export projection_integrand
include("integrand.jl")

include("dalitzplotsrecipe.jl")

end  # module ThreeBodyDecays
