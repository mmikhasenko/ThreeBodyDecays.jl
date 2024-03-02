# Form Factors
# 
struct BlattWeisskopf{L} <: AbstractFlexFunc
    d::Float64
end
breakup(m, m1, m2) = sqrt((m - (m1 + m2)) * (m + (m1 + m2)) * (m - (m1 - m2)) * (m + (m1 - m2))) / 2m
(bw::BlattWeisskopf{L})(p::Number) where {L} = error("BlattWeisskopf{L} is not defined for L>2")
(bw::BlattWeisskopf{0})(p::Number) = one(p)
(bw::BlattWeisskopf{1})(p::Number) = (z² = (bw.d * p)^2;
sqrt(z²) / sqrt(1 + z²))
(bw::BlattWeisskopf{2})(p::Number) = (z² = (bw.d * p)^2;
sqrt(z²^2) / sqrt(9 + 3 * z² + z²^2))
(bw::BlattWeisskopf{3})(p::Number) = (z² = (bw.d * p)^2;
sqrt(z²^3) / sqrt(225 + 45z² + 6z²^2 + z²^3))
# 