struct minus_one end
import Base: ^
^(x::minus_one, n::Number) = isodd(n) ? -1 : 1
macro x_str(s::String)
    minus_one()
end

function letterL(l::Int)
    waves = ['S', 'P', 'D', 'F', 'G', 'H']
    return l < 6 ? waves[l+1] : string(l)[1]
end

letterL(l::String) = letterL(Meta.parse(l))

shift_by_half(v) = v[1:end-1] .+ first(diff(v)) / 2
