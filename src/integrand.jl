"""
    phase_space_integrand(function_σs, ms; k::Int)

Calculate the phase space integrand for a given function `function_σs`,
and mass tuple `ms`. The key argument `k` specifies the mapping: `σk->[0,1]`, zk->[0,1].
It returns an integrand function of x, x ∈ [0,1]x[0,1] domain to pass to a numerical integrator.

# Arguments
- `function_σs`: A function that takes a MandelstamTuple and returns a scalar.
- `ms`: A scalar representing the mass.
- `k`: An integer representing the mapping index.

# Usage
```julia
integrand = phase_space_integrand(function_σs, ms; k)
```

# See also
- `x2σs`
- `projection_integrand`

"""
function phase_space_integrand(function_σs, ms; k::Int)
    mssq = ms^2
    i, j = ij_from_k(k)
    misq, mjsq, mksq, m0sq = mssq[i], mssq[j], mssq[k], mssq[4]
    σk_min, σk_max = lims(ms; k)
    #
    function integrand(x)
        σs = x2σs(x, ms; k)
        σk = σs[k]
        jac_z = sqrt(Kallen(m0sq, σk, mksq) * Kallen(σk, misq, mjsq)) / σk
        jac_σk = (σk_max - σk_min)
        value = function_σs(σs) * jac_σk * jac_z / m0sq
        return value
    end
    return integrand
end

"""
projection_integrand(function_σs, ms, σk; k)

Calculate the projection integrand for a given function `function_σs`,
mass tuple `ms`, and Mandelstam variable `σk`, with `k` specified by a keyword argument.
It returns an integrand function of x, x ∈ [0,1] to pass to a numerical integrator.

# Arguments
- `function_σs`: A function that takes a MandelstamTuple and returns a scalar.
- `ms`: A scalar representing the mass.
- `σk`: A scalar representing the Mandelstam variable.
- `k`: A scalar representing the momentum transfer (optional).

# Usage
```julia
plot(4.2, 4.6) do e1
I = Base.Fix1(unpolarized_intensity, model)
integrand = projection_integrand(I, masses(model), e1^2; k = 3)
e1 * quadgk(integrand, 0, 1)[1]
end
```
"""
function projection_integrand(function_σs, ms, σk; k)
    l, h = lims(ms; k)
    !(l < σk < h) && return x -> 0.0
    σj_lims = σjofk.([-1, 1], Ref(σk), Ref(ms^2); k)
    function integrand(x)
        σj = σj_lims[1] + x[1] * (σj_lims[2] - σj_lims[1])
        σi = sum(ms^2) - σk - σj
        σt = circleorigin(-k, (σi, σj, σk))
        σs = MandelstamTuple{typeof(ms.m0)}(σt)
        return function_σs(σs) * (σj_lims[2] - σj_lims[1])
    end
    return integrand
end

"""
    project_cosθij_intergand(fs, ms, z; k)

Calculate the integrand for projecting onto cos(θij) for a given function `fs`,
mass tuple `ms`, and cosine value `z`, with `k` specifying the spectator.
Returns an integrand function of x ∈ [0,1] to pass to a numerical integrator.

# Arguments
- `fs`: A function that takes a MandelstamTuple and returns a scalar.
- `ms`: A MassTuple containing the masses.
- `z`: The cosine value to project onto.
- `k`: The spectator index.

# Returns
- An integrand function for numerical integration.

# Example
```julia
let Nb = 100
    [sum(range(-1, 1, Nb)) do z
        integrand = project_cosθij_intergand(I, ms, z; k)
        quadgk(integrand, 0, 1)[1] * 2/Nb
    end for k in 1:3]
end
```
"""
function project_cosθij_intergand(fs, ms, z; k)
    i, j = ij_from_k(k)
    mi, mj, mk, m0 = ms[i], ms[j], ms[k], ms[4]

    function integrand(xσk)
        σs = x2σs([xσk, (z + 1) / 2], ms; k)
        σk = σs[k]
        jac = sqrt(Kallen(σk, mk^2, m0^2) * Kallen(σk, mi^2, mj^2)) / σk
        return fs(σs) / 2 * jac * ((m0 - mk)^2 - (mi + mj)^2)
    end
    return integrand
end
