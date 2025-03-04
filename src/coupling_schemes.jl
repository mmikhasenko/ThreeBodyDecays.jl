"""
    possible_ls(jp1::SpinParity, jp2::SpinParity; jp::SpinParity)
    possible_ls_ij(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k::Int)
    possible_ls_Rk(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k::Int)

Calculate possible orbital angular momentum (l) and spin (s) combinations for a two-body system.
- `possible_ls`: For any two particles with given spin-parities
- `possible_ls_ij`: For the (i,j) pair in a three-body decay
- `possible_ls_Rk`: For the resonance-spectator system

# Arguments
- `jp1`, `jp2`: Spin-parity quantum numbers of the particles
- `jp`: Total spin-parity of the system
- `two_js`: Spins of all particles in the system
- `Ps`: Parities of all particles
- `k`: Spectator index

# Returns
- Vector of tuples (l,s) containing allowed combinations

# Example
```julia
# For a decay 1⁻ → 1/2⁺ + 1/2⁺
ls = possible_ls(jp"1/2+", jp"1/2+"; jp=jp"1-")
# For a three-body system
ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0=4.0)
two_js = ThreeBodySpins(1, 1, 0; h0=2)
Ps = ThreeBodyParities('+','+','+'; P0='-')
ls_ij = possible_ls_ij(jp"1-", two_js, Ps; k=1)
```
"""
function possible_ls(jp1::SpinParity, jp2::SpinParity; jp::SpinParity)
    two_ls = Vector{Tuple{Int,Int}}(undef, 0)
    for two_s ∈ abs(jp1.two_j - jp2.two_j):2:abs(jp1.two_j + jp2.two_j)
        for two_l ∈ abs(jp.two_j - two_s):2:abs(jp.two_j + two_s)
            if jp1.p ⊗ jp2.p ⊗ jp.p == (isodd(div(two_l, 2)) ? '-' : '+')
                push!(two_ls, (two_l, two_s))
            end
        end
    end
    return sort(two_ls, by = x -> x[1])
end

const TwoBodyTopologySpinParity = Pair{SpinParity,Tuple{SpinParity,SpinParity}}
possible_ls((jp, (jp1, jp2))::TwoBodyTopologySpinParity) = possible_ls(jp1, jp2; jp)
possible_ls(jp1::AbstractString, jp2::AbstractString; jp::AbstractString) =
    possible_ls(str2jp(jp1), str2jp(jp2); jp = str2jp(jp))

function possible_ls_ij(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k::Int)
    i, j = ij_from_k(k)
    return possible_ls(SpinParity(two_js[i], Ps[i]), SpinParity(two_js[j], Ps[j]); jp)
end

possible_ls_Rk(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k::Int) =
    possible_ls(jp, SpinParity(two_js[k], Ps[k]); jp = SpinParity(two_js[4], Ps[4]))

"""
    possible_lsLS(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k::Int)

Calculate all possible combinations of orbital angular momenta and spins for both
the isobar decay (l,s) and the total system decay (L,S).

# Arguments
- `jp`: Spin-parity of the resonance
- `two_js`: Spins of all particles
- `Ps`: Parities of all particles
- `k`: Spectator index

# Returns
- Vector of named tuples containing with all possible (two_ls, two_LS) combinations

# Example
```julia
ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0=4.0)
two_js = ThreeBodySpins(1, 1, 0; h0=2)
Ps = ThreeBodyParities('+','+','+'; P0='-')
lsLS = possible_lsLS(jp"1-", two_js, Ps; k=1)
```
"""
function possible_lsLS(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k::Int)
    lsv = possible_ls_ij(jp, two_js, Ps; k)
    LSv = possible_ls_Rk(jp, two_js, Ps; k)
    return [(; two_ls, two_LS) for (two_ls, two_LS) in Iterators.product(lsv, LSv)]
end

"""
    possible_l_s_L_S(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k)

Convert the output of `possible_lsLS` to a more readable format with half-integer values.

# Arguments
Same as `possible_lsLS`

# Returns
- Vector of named tuples containing (l,s,L,S) as strings representing half-integer values

# Example
```julia
ms = ThreeBodyMasses(1.0, 1.0, 1.0; m0=4.0)
two_js = ThreeBodySpins(1, 1, 0; h0=2)
Ps = ThreeBodyParities('+','+','+'; P0='-')
qn = possible_l_s_L_S(jp"1-", two_js, Ps; k=1)
# Each element has fields like qn[1].l == "0", qn[1].s == "1/2", etc.
```
"""
function possible_l_s_L_S(jp::SpinParity, two_js::SpinTuple, Ps::ParityTuple; k)
    _lsLS = possible_lsLS(jp, two_js, Ps; k) |> vec
    map(_lsLS) do (two_ls, two_LS)
        l, s = two_ls .|> d2
        L, S = two_LS .|> d2
        (; l, s, L, S)
    end
end

"""
    jls_coupling(two_j1, two_λ1, two_j2, two_λ2, two_j, two_l, two_s)

Calculate the LS-coupling coefficient for two particles coupling to total angular momentum j.
Uses Clebsch-Gordan coefficients to compute the coupling of orbital angular momentum l
and total spin s to total j, with specified helicities.

# Arguments
- `two_j1`, `two_j2`: Twice the spins of the particles
- `two_λ1`, `two_λ2`: Twice the helicities
- `two_j`: Twice the total angular momentum
- `two_l`: Twice the orbital angular momentum
- `two_s`: Twice the total spin

# Returns
- The coupling coefficient value

# Example
```julia
# For j₁=1/2, λ₁=1/2, j₂=1/2, λ₂=-1/2, j=1, l=1, s=1
julia> ThreeBodyDecays.jls_coupling(1, 1, 1, -1, 2, 2, 2)
-0.7071067811865477
```
"""
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
function complete_l_s_L_S(
    jp::SpinParity,
    two_js::SpinTuple,
    Ps::ParityTuple,
    constraints;
    k,
)
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
    constraints;
    k,
)
    #
    qn = vcat([possible_l_s_L_S(jp, two_js, Ps; k) for Ps in several_Ps]...)
    all_qn = copy(qn)
    filter_qn!(qn, constraints)
    complete_data = check_and_attach(constraints, qn, all_qn)
    return complete_data
end
