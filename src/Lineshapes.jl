module Lineshapes
# 
using Parameters
# 

export AbstractFlexFunc
export Compose
export Scale
export Product
import Base: *
include(joinpath("lineshapes", "operations.jl"))

export BlattWeisskopf
include(joinpath("lineshapes", "formfactors.jl"))

export BreitWigner
breakup(m, m1, m2) = sqrt((m - (m1 + m2)) * (m + (m1 + m2)) * (m - (m1 - m2)) * (m + (m1 - m2))) / 2m
include(joinpath("lineshapes", "shapes.jl"))

end