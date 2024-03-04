function possible_ls(jp1::SpinParity, jp2::SpinParity; jp::SpinParity)
    two_ls = Vector{Tuple{Int,Int}}(undef, 0)
    for two_s in abs(jp1.two_j - jp2.two_j):2:abs(jp1.two_j + jp2.two_j)
        for two_l in abs(jp.two_j - two_s):2:abs(jp.two_j + two_s)
            if jp1.p ⊗ jp2.p ⊗ jp.p == (isodd(div(two_l, 2)) ? '-' : '+')
                push!(two_ls, (two_l, two_s))
            end
        end
    end
    return sort(two_ls, by=x -> x[1])
end


const TwoBodyTopologySpinParity = Pair{SpinParity,Tuple{SpinParity,SpinParity}}
possible_ls((jp, (jp1, jp2))::TwoBodyTopologySpinParity) = possible_ls(jp1, jp2; jp)
possible_ls(jp1::AbstractString, jp2::AbstractString; jp::AbstractString) =
    possible_ls(str2jp(jp1), str2jp(jp2); jp=str2jp(jp))
function possible_ls_ij(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k::Int)
    i, j = ij_from_k(k)
    return possible_ls(SpinParity(two_js[i], Ps[i]), SpinParity(two_js[j], Ps[j]); jp)
end
possible_ls_Rk(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k::Int) =
    possible_ls(jp, SpinParity(two_js[k], Ps[k]); jp=SpinParity(two_js[4], Ps[4]))

function possible_lsLS(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k::Int)
    lsv = possible_ls_ij(jp, two_js, Ps; k)
    LSv = possible_ls_Rk(jp, two_js, Ps; k)
    return [(; two_ls, two_LS) for (two_ls, two_LS) in Iterators.product(lsv, LSv)]
end

function jls_coupling(two_j1, two_λ1, two_j2, two_λ2, two_j, two_l, two_s)
    T1 = one(two_λ1) # type consistency
    return sqrt((two_l * T1 + 1) / (two_j * T1 + 1)) *
           CG_doublearg(two_j1, two_λ1, two_j2, -two_λ2, two_s, two_λ1 - two_λ2) *
           CG_doublearg(two_l, zero(two_λ1 - two_λ2), two_s, two_λ1 - two_λ2, two_j, two_λ1 - two_λ2)
end
