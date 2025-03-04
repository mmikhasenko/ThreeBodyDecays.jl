"""
    ijk(k::Int)
    ij_from_k(k::Int)

Return a tuple of indices (i,j,k) for a three-body system, where i,j are ordered cyclically from k.
Used extensively in three-body decay calculations to determine the spectator particle (k) and the pair (i,j).

# Arguments
- `k::Int`: Index of the spectator particle (1, 2, or 3)

# Returns
- `Tuple{Int,Int,Int}`: A tuple of indices (i,j,k) where i,j are ordered cyclically from k

# Example
```julia
ijk(1) # returns (2, 3, 1)
ijk(2) # returns (3, 1, 2)
ijk(3) # returns (1, 2, 3)
```
"""
ijk(k::Int) = (k + 1, k + 2, k) |> x -> mod.(x, Ref(Base.OneTo(3)))
ij_from_k(k::Int) = ijk(k)

"""
    struct minusone
    @x_str(s::String)

A utility type and macro for handling -1 powers in calculations.
The macro creates a `minusone` instance that returns -1 when raised to odd powers and 1 when raised to even powers.

# Example
```julia
x"-1"^3  # returns -1
x"-1"^2  # returns 1
```
"""
struct minusone end
import Base: ^
^(x::minusone, n::Number) = isodd(n) ? -1 : 1
macro x_str(s::String)
    minusone()
end

"""
    letterL(l::Int)
    letterL(l::String)

Convert an angular momentum quantum number to its spectroscopic notation.
For l = 0,1,2,3,4,5 returns S,P,D,F,G,H respectively. For l ≥ 6 returns the first character of the string representation.

# Arguments
- `l`: Angular momentum quantum number (as Int or String)

# Returns
- `Char`: The spectroscopic notation for the given angular momentum

# Example
```julia
letterL(0)  # returns 'S'
letterL(1)  # returns 'P'
letterL("2")  # returns 'D'
```
"""
function letterL(l::Int)
    waves = ['S', 'P', 'D', 'F', 'G', 'H']
    return l < 6 ? waves[l+1] : string(l)[1]
end

letterL(l::String) = letterL(Meta.parse(l))

"""
    shift_by_half(v)

Shift a vector of values by half the step size between elements.
Used in plotting and binning operations to center values between grid points.

# Arguments
- `v`: Vector of values

# Returns
- Vector of values shifted by half the step size

# Example
```julia
v = [1, 2, 3, 4]
shift_by_half(v)  # returns [1.5, 2.5, 3.5]
```
"""
shift_by_half(v) = v[1:end-1] .+ first(diff(v)) / 2
