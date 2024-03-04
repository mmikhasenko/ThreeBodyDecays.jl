# Form Factors
# 
struct BlattWeisskopf{L} <: AbstractFlexFunc
    d::Float64
end
breakup(m, m1, m2) = sqrt((m - (m1 + m2)) * (m + (m1 + m2)) * (m - (m1 - m2)) * (m + (m1 - m2))) / 2m
(ff::BlattWeisskopf{L})(p::Number) where {L} = error("BlattWeisskopf{L} is not defined for L>7")
(ff::BlattWeisskopf{0})(p::Number) = one(p)
(ff::BlattWeisskopf{1})(p::Number) = (z² = (ff.d * p)^2; sqrt(z² / chi_poly1(z²)))
(ff::BlattWeisskopf{2})(p::Number) = (z² = (ff.d * p)^2; sqrt(z²^2 / chi_poly2(z²)))
(ff::BlattWeisskopf{3})(p::Number) = (z² = (ff.d * p)^2; sqrt(z²^3 / chi_poly3(z²)))
(ff::BlattWeisskopf{4})(p::Number) = (z² = (ff.d * p)^2; sqrt(z²^4 / chi_poly4(z²)))
(ff::BlattWeisskopf{5})(p::Number) = (z² = (ff.d * p)^2; sqrt(z²^5 / chi_poly5(z²)))
(ff::BlattWeisskopf{6})(p::Number) = (z² = (ff.d * p)^2; sqrt(z²^6 / chi_poly6(z²)))
(ff::BlattWeisskopf{7})(p::Number) = (z² = (ff.d * p)^2; sqrt(z²^7 / chi_poly7(z²)))
# 
chi_poly0(z²) = one(z²)
chi_poly1(z²) = 1 + z²
chi_poly2(z²) = 9 + 3z² + z²^2
chi_poly3(z²) = 225 + 45z² + 6z²^2 + z²^3
chi_poly4(z²) = 11025 + 1575z² + 135z²^2 + 10z²^3 + z²^4
chi_poly5(z²) = 893025 + 99225z² + 6300z²^2 + 315z²^3 + 15z²^4 + z²^5
chi_poly6(z²) = 108056025 + 9823275z² + 496125z²^2 + 18900z²^3 + 630z²^4 + 21z²^5 + z²^6
chi_poly7(z²) = 18261468225 + 1404728325z² + 58939650z²^2 + 1819125z²^3 + 47250z²^4 + 1134z²^5 + 28z²^6 + z²^7
