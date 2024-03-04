# Operations
# 
abstract type AbstractFlexFunc end

const UnionFlexFunc = Union{Function,AbstractFlexFunc}
struct WrapFlexFunction{T<:UnionFlexFunc} <: AbstractFlexFunc
    F::T
end
(L::WrapFlexFunction)(σ::Number) = L.F(σ)
struct ComposeFlexFunc{T,X} <: AbstractFlexFunc
    A::T
    F::X
end
(f::ComposeFlexFunc)(σ::Number) = f.F(f.A(σ))
(f::AbstractFlexFunc)(a::Function) = ComposeFlexFunc(a, f)
(f::AbstractFlexFunc)(a::AbstractFlexFunc) = ComposeFlexFunc(a, f)
# 
struct ScaleFlexFunc{T<:UnionFlexFunc,X<:Number} <: AbstractFlexFunc
    F::T
    S::X
end
*(f::AbstractFlexFunc, x::Number) = ScaleFlexFunc(f, x)
*(x::Number, f::AbstractFlexFunc) = ScaleFlexFunc(f, x)
(p::ScaleFlexFunc)(σ::Number) = p.F(σ) * p.S
# 
struct ProductFlexFunc{T1<:UnionFlexFunc,T2<:UnionFlexFunc} <: AbstractFlexFunc
    F1::T1
    F2::T2
end
*(F1::AbstractFlexFunc, F2::Function) = ProductFlexFunc(F1, F2)
*(F1::Function, F2::AbstractFlexFunc) = ProductFlexFunc(F1, F2)
*(F1::AbstractFlexFunc, F2::AbstractFlexFunc) = ProductFlexFunc(F1, F2)
(p::ProductFlexFunc)(σ::Number) = p.F1(σ) * p.F2(σ)
# 