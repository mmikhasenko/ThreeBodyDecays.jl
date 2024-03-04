module Lineshapes
# 
using Parameters
# 

export AbstractFlexFunc
export UnionFlexFunc
export WrapFlexFunction
export ComposeFlexFunc
export ScaleFlexFunc
export ProductFlexFunc
import Base: *
include(joinpath("lineshapes", "operations.jl"))

export BlattWeisskopf
include(joinpath("lineshapes", "formfactors.jl"))

export BreitWigner
include(joinpath("lineshapes", "shapes.jl"))

end