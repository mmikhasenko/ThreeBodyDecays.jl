function possible_ls(jp1::SpinParity, jp2::SpinParity; jp::SpinParity)
    two_ls = Vector{Tuple{Int, Int}}(undef, 0)
    for two_s ∈ abs(jp1.two_j - jp2.two_j):2:abs(jp1.two_j + jp2.two_j)
        for two_l ∈ abs(jp.two_j - two_s):2:abs(jp.two_j + two_s)
            if jp1.p ⊗ jp2.p ⊗ jp.p == (isodd(div(two_l, 2)) ? '-' : '+')
                push!(two_ls, (two_l, two_s))
            end
        end
    end
    return sort(two_ls, by = x -> x[1])
end


const TwoBodyTopologySpinParity = Pair{SpinParity, Tuple{SpinParity, SpinParity}}
possible_ls((jp, (jp1, jp2))::TwoBodyTopologySpinParity) = possible_ls(jp1, jp2; jp)
possible_ls(jp1::AbstractString, jp2::AbstractString; jp::AbstractString) =
    possible_ls(str2jp(jp1), str2jp(jp2); jp = str2jp(jp))
function possible_ls_ij(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k::Int)
    i, j = ij_from_k(k)
    return possible_ls(SpinParity(two_js[i], Ps[i]), SpinParity(two_js[j], Ps[j]); jp)
end
possible_ls_Rk(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k::Int) =
    possible_ls(jp, SpinParity(two_js[k], Ps[k]); jp = SpinParity(two_js[4], Ps[4]))

function possible_lsLS(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k::Int)
    lsv = possible_ls_ij(jp, two_js, Ps; k)
    LSv = possible_ls_Rk(jp, two_js, Ps; k)
    return [(; two_ls, two_LS) for (two_ls, two_LS) in Iterators.product(lsv, LSv)]
end

function possible_l_s_L_S(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k)
    _lsLS = possible_lsLS(jp, two_js, Ps; k) |> vec
    map(_lsLS) do (two_ls, two_LS)
        l, s = two_ls .|> d2
        L, S = two_LS .|> d2
        (; l, s, L, S)
    end
end

function jls_coupling(two_j1, two_λ1, two_j2, two_λ2, two_j, two_l, two_s)
    T1 = one(two_λ1) # type consistency
    return sqrt((two_l * T1 + 1) / (two_j * T1 + 1)) *
           CG_doublearg(two_j1, two_λ1, two_j2, -two_λ2, two_s, two_λ1 - two_λ2) *
           CG_doublearg(
               two_l,
               zero(two_λ1 - two_λ2),
               two_s,
               two_λ1 - two_λ2,
               two_j,
               two_λ1 - two_λ2,
           )
end


function check_and_attach(data, df, df_all)
    if size(df, 1) == 0
        error("""
        Decay description
        $(data)
        is inconsistent. Pick one of all possible,
        $(df_all)
        """)
    end
    if size(df, 1) > 1
        for q in [:L, :l, :S, :s]
            possible_q = unique(getproperty.(df, q))
            !(q in keys(data)) && error("""
            Decay description
            $(data)
            is not complete. Add $q = $(possible_q[1]), or other from list $(Tuple(possible_q))
            """)
        end
    end
    @unpack l, s, L, S = df[1]
    return (; data..., l, s, L, S)
end

filter_qn!(possible_qn, constraints::NamedTuple) =
    for var_name in [:L, :S, :l, :s]
        if var_name in keys(constraints)
            value = getproperty(constraints, var_name)
            value_str = value |> x2 |> d2
            filter!(x -> getproperty(x, var_name) == value_str, possible_qn)
        end
    end
#
function complete_l_s_L_S(jp::SpinParity,
    two_js::SpinTuple,
    Ps::ParityTuple,
    constraints; k)
    #
    qn = possible_l_s_L_S(jp, two_js, Ps; k)
    all_qn = copy(qn)
    filter_qn!(qn, constraints)
    complete_data = check_and_attach(constraints, qn, all_qn)
    return complete_data
end

function complete_l_s_L_S(
    jp::SpinParity,
    two_js::SpinTuple,
    several_Ps::AbstractArray{ParityTuple},
    constraints; k)
    #
    qn = vcat([possible_l_s_L_S(jp, two_js, Ps; k) for Ps in several_Ps]...)
    all_qn = copy(qn)
    filter_qn!(qn, constraints)
    complete_data = check_and_attach(constraints, qn, all_qn)
    return complete_data
end
