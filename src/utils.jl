struct minusone end
import Base: ^
^(x::minusone, n::Number) = isodd(n) ? -1 : 1
macro x_str(s::String)
    minusone()
end

function letterL(l::Int)
    waves = ['S', 'P', 'D', 'F', 'G', 'H']
    return l < 6 ? waves[l+1] : string(l)[1]
end

letterL(l::String) = letterL(Meta.parse(l))