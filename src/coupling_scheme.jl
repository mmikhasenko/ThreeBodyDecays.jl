
struct SpinParity
    two_j::Int
    p::Char
end

SpinParity(s::String) = str2jp(s)

length(jp1::SpinParity) = 0

# dealing with spin 1/2
x2(v) = @. Int(2v)
x2(v::AbstractString) = Meta.parse(v) |> eval |> x2

_d2(two_s::Int) = iseven(two_s) ? "$(div(two_s,2))" : "$(two_s)/2"
d2(v) = _d2.(v)


⊗(p1::Char, p2::Char) = p1 == p2 ? '+' : '-'
⊗(jp1::SpinParity, jp2::SpinParity) = [jp(j, ⊗(jp1.p, jp2.p)) for j in abs(jp1.two_j - jp2.two_j):abs(jp1.two_j + jp2.two_j)]

function str2jp(pin::AbstractString)
    p = filter(x -> x != '^', pin)
    !(contains(p, '/')) && return SpinParity(x2(p[1:end-1]), p[end])
    p[end-2:end-1] != "/2" && error("the string should be `x/2±`, while it is $(p)")
    two_j = Meta.parse(p[1:end-3])
    !(typeof(two_j) <: Int) && error("the string should be `x/2±`, while it is $(p)")
    return SpinParity(two_j, p[end])
end

macro jp_str(p)
    return str2jp(p)
end

possible_ls((jp, (jp1, jp2))::Pair{SpinParity,Tuple{SpinParity,SpinParity}}) = possible_ls(jp1, jp2; jp)
possible_ls(jp1::AbstractString, jp2::AbstractString; jp::AbstractString) =
    possible_ls(str2jp(jp1), str2jp(jp2); jp=str2jp(jp))

function possible_ls(jp1::SpinParity, jp2::SpinParity; jp::SpinParity)
    two_ls = Vector{Tuple{Int,Int}}(undef, 0)
    for two_s in abs(jp1.two_j - jp2.two_j):abs(jp1.two_j + jp2.two_j)
        for two_l in abs(jp.two_j - s):abs(jp.two_j + s)
            if jp1.p ⊗ jp2.p ⊗ jp.p == (isodd(l) ? '-' : '+')
                push!(two_ls, (two_l, two_s))
            end
        end
    end
    return sort(ls, by=x -> x[1])
end

function possible_lsLS(k::Int, jpR::SpinParity, jps::AbstractVector{SpinParity})
    i, j = ij_from_k(k)
    lsv = possible_ls(jps[i], jps[j]; jp=jpR)
    LSv = possible_ls(jpR, jps[k]; jp=jps[4])
    return [(ls=ls, LS=LS) for (ls, LS) in Iterators.product(lsv, LSv)]
end
function possible_lsLS(k::Int, two_s, parity::Char, two_js, Ps)
    jpR = jp(two_s // 2, parity)
    jps = jp.(zip(Tuple(two_js) .// 2, Ps))
    return possible_lsLS(k, jpR, jps)
end

function jls_coupling(two_j1, two_λ1, two_j2, two_λ2, two_j, two_l, two_s)
    T1 = one(two_λ1) # type consistency
    return sqrt((two_l * T1 + 1) / (two_j * T1 + 1)) *
           CG_doublearg(two_j1, two_λ1, two_j2, -two_λ2, two_s, two_λ1 - two_λ2) *
           CG_doublearg(two_l, zero(two_λ1 - two_λ2), two_s, two_λ1 - two_λ2, two_j, two_λ1 - two_λ2)
end
