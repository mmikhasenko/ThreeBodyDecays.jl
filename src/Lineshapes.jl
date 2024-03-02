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
include(joinpath("lineshapes", "shapes.jl"))

end