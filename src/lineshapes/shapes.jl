# Lineshapes
# 

BW(σ, m, Γ) = 1 / (m^2 - σ - 1im * m * Γ)

@with_kw struct BreitWigner <: AbstractFlexFunc
    m::Float64
    channels::AbstractVector{<:NamedTuple{(:gsq, :ma, :mb, :l, :d)}}
end
function (bw::BreitWigner)(σ::Number)
    m0 = bw.m
    Γ = sum(bw.channels) do channel
        @unpack gsq, ma, mb, l, d = channel
        FF = BlattWeisskopf{l}(d)
        _p = breakup(sqrt(σ), ma, mb)
        gsq * 2_p / sqrt(σ) * FF(_p)^2
    end
    BW(σ, m0, Γ)
end

BreitWigner(m::Float64, Γ::Float64) = BreitWigner(m, [(gsq=Γ, ma=0.0, mb=0.0, l=0, d=1.1)])
function BreitWigner(m::Real, Γ::Real, ma::Number, mb::Number, l::Int, d::Real)
    _p0 = breakup(m, ma, mb)
    FF = BlattWeisskopf{l}(d)
    gsq = Γ / (2_p0) * m / FF(_p0)^2
    return BreitWigner(m, [(; gsq, ma, mb, l, d)])
end
