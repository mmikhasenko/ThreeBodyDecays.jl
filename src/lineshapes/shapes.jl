# Lineshapes
# 
@with_kw struct BreitWigner{T} <: AbstractFlexFunc
    m::Float64
    Γ::Float64
    ma::T = 0.0
    mb::T = 0.0
end
(bw::BreitWigner)(σ::Number) = 1 / (bw.m^2 - σ - 1im * bw.m * bw.Γ)
