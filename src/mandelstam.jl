"""
    MandelstamTuple{T}

A named tuple representing Mandelstam variables `(; σ1, σ2, σ3)` for a three-body system.

Primarily created via [`Invariants`](@ref).
"""
const MandelstamTuple{T} = NamedTuple{(:σ1, :σ2, :σ3),NTuple{3,T}}

"""
    Invariants(ms::MassTuple{T}; σ1, σ2, σ3)

Construct a tuple of Mandelstam invariants from given invariants and the mass tuple.

# Arguments
- `ms::MassTuple{T}`: Masses of the system
- `σ1`, `σ2`, `σ3`: Optional invariants (two or three must be provided)

# Returns
- `MandelstamTuple{T}`: The complete set of invariants

# Throws
- `ErrorException` if not exactly two or three invariants are provided
- `ErrorException` if the invariants violate mass constraints
"""
function Invariants(
    ms::MassTuple{T};
    σ1 = -one(ms.m0),
    σ2 = -one(ms.m0),
    σ3 = -one(ms.m0),
) where {T}
    # Count how many invariants are provided
    n_provided = count(x -> x != -one(ms.m0), (σ1, σ2, σ3))

    if n_provided == 0 || n_provided == 1
        error("At least two invariants must be provided")
    elseif n_provided == 2
        # Calculate the third invariant
        if σ3 == -one(ms.m0)
            σ3 = sum(ms^2) - σ1 - σ2
        elseif σ1 == -one(ms.m0)
            σ1 = sum(ms^2) - σ2 - σ3
        else
            σ2 = sum(ms^2) - σ3 - σ1
        end
    elseif n_provided == 3
        # Check if the invariants are physical
        if !isphysical((σ1, σ2, σ3), ms)
            error("The invariants violate mass constraints")
        end
    end

    return MandelstamTuple{T}((σ1, σ2, σ3))
end
Invariants(; σ1, σ2, σ3) = MandelstamTuple{typeof(σ1)}((σ1, σ2, σ3))
Invariants(σ1, σ2, σ3) = MandelstamTuple{typeof(σ1)}((σ1, σ2, σ3))

"""
    x2σs(x, ms::MassTuple; k::Int = last(findmin(Tuple(ms))))

Maps a pair of variables `x` to a physical set of mandelstam invariants `σs` using a linear transformation.
x[1] is mapped to σ[k], and x[2] is mapped into cosθij, that is used to compute σ[i] and σ[j].
Uniform distribution of x does not lead to the phase space distribution of σs.

## Arguments

- `x` : a pair of numbers
- `ms` : masses of the system as a `MassTuple`, see `ThreeBodyMasses`.
- `k` : the index for which the variable is not generated. By default, the function picks the coordinates where the Dalitz plot has the closest shape to the squared fitting box.

## Returns
- an instance of `MandelstamTuple` with the squared masses.

## Example

The phase space sample with 100 points can be generated as follows:
````julia
σs0 = x2σs([0.3, 0.2], ms)
````

See also [`y2σs`](@ref).
"""
function x2σs(x, ms::MassTuple; k::Int)
    σk = fitin(x[1], lims(ms; k))
    σj = σjofk(2x[2] - 1, σk, ms^2; k)
    σi = sum(ms^2) - σk - σj
    σt = circleorigin(-k, (σi, σj, σk))
    return MandelstamTuple{typeof(ms.m0)}(σt)
end

"""
    y2σs(y, ms::MassTuple; k::Int = last(findmin(Tuple(ms))))

Maps a pair of variables to the plane of squared masses using a linear transformation,
by shifting and scaling `y[1]` to `σ[i]`, and `y[2]` to `σ[j]`.
The indices of variables are determined by the argument `k`.
The mapping does not guarantee that values are physical, however,
physical values of σs are phase-space distributed for uniform y.

## Arguments

- `y` : a pair of numbers
- `ms` : masses of the system as a `MassTuple`, see `ThreeBodyMasses`.
- `k` : the index for which the variable is not generated. By default, the function picks the coordinates
where the Dalitz plot has the closest shape to the squared fitting box.

## Returns

- an instance of `MandelstamTuple` with the squared masses.

## Example

The phase space sample with 100 points can be generated as follows:
````julia
data = let
    N = 100
    # map random variables to dalitz
    _data = mapslices(rand(N, 2); dims=2) do xy
        y2σs(xy, ms)
    end[:, 1]
    # select physical
    filter!(_data) do σs
        isphysical(σs, ms)
    end
    _data
end
````

See also [`x2σs`](@ref).
"""
function y2σs(y, ms::MassTuple; k::Int = last(findmin(Tuple(ms))))
    i, j, _ = ij_from_k(k)
    σi = fitin(y[1], lims(ms; k = i))
    σj = fitin(y[2], lims(ms; k = j))
    σk = sum(ms^2) - σi - σj
    σt = circleorigin(-k, (σi, σj, σk))
    return MandelstamTuple{typeof(ms.m0)}(σt)
end

"""
    circleorigin(k, t)

Apply a cyclic permutation to a 3-tuple `t` with an origin/handedness convention used across
the package.

This is a small utility used to re-order invariant tuples `(σ1, σ2, σ3)`
depending on which index is treated as “special” (often a spectator index `k`).

# Arguments
- `k`: An integer shift parameter (typically `±1`, `±2`, `±3` in this package).
- `t`: A 3-tuple.

# Returns
- A 3-tuple containing the permuted elements of `t`.
"""
circleorigin(k, t) = (t[mod(k, 3)+1], t[mod(k+1, 3)+1], t[mod(k-1, 3)+1])

fitin(y, (a, b)) = a + y * (b - a)
