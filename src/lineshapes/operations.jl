# Operations
# 
abstract type AbstractFlexFunc end
struct Compose{T,X} <: AbstractFlexFunc
    A::T
    F::X
end
(f::Compose)(σ::Number) = f.F(f.A(σ))
(f::AbstractFlexFunc)(a::Function) = Compose(a, f)
(f::AbstractFlexFunc)(a::AbstractFlexFunc) = Compose(a, f)
# 
struct Scale{T<:AbstractFlexFunc,X<:Number} <: AbstractFlexFunc
    F::T
    S::X
end
*(f::AbstractFlexFunc, x::Number) = Scale(f, x)
*(x::Number, f::AbstractFlexFunc) = Scale(f, x)
(p::Scale)(σ::Number) = p.F(σ) * p.S
# 
struct Product{T1<:AbstractFlexFunc,T2<:AbstractFlexFunc} <: AbstractFlexFunc
    F1::T1
    F2::T2
end
*(F1::AbstractFlexFunc, F2::AbstractFlexFunc) = Product(F1, F2)
(p::Product)(σ::Number) = p.F1(σ) * p.F2(σ)
# 