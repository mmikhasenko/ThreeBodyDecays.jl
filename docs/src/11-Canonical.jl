#md # # Canonical formalism
#nb # # Canonical formalism
#jl # # Canonical formalism

# This notebook demonstrates the equivalence of the canonical formulation of the decay to the helicity formalism with LS couplings,
# for a single chain

using Markdown
using InteractiveUtils
using Random
Random.seed!(1212);

# ## Implementation

using ThreeBodyDecays
using ThreeBodyDecays.Parameters
using ThreeBodyDecays.PartialWaveFunctions

# The `CC` function calculates the Clebsch-Gordan coefficients for the given quantum numbers.
function CC(two_m12t, two_j12t, two_ls)
    two_j1, two_j2, two_j = two_j12t
    two_m1, two_m2, two_mj = two_m12t
    two_l, two_s = two_ls

    two_ms = two_m1 + two_m2
    two_ml = two_mj - two_ms

    clebschgordan_doublearg(two_j1, two_m1, two_j2, two_m2, two_s, two_ms) *
    clebschgordan_doublearg(two_l, two_ml, two_s, two_ms, two_j, two_mj)
end;

# The `Y_lm_doublearg` function calculates the spherical harmonics for the given quantum numbers and angles.
Y_lm_doublearg(two_l, two_m, ϕ, cosθ) =
    sqrt((two_l + 1) / (4π)) * conj(wignerD_doublearg(two_l, two_m, 0, ϕ, cosθ, 0));

# The `A_node_canonical` function calculates the amplitude for a node in the canonical formalism.
function A_node_canonical(angles, two_m12t, two_j12t, two_ls)
    two_m1, two_m2, two_mj = two_m12t
    two_ms = two_m1 + two_m2
    two_ml = two_mj - two_ms
    two_l, _ = two_ls
    @unpack ϕ, cosθ = angles

    return CC(two_m12t, two_j12t, two_ls) *
           Y_lm_doublearg(two_l, two_ml, ϕ, cosθ) *
           sqrt(4π)
end;

# The `A_canonical` function calculates the amplitude for the entire decay chain in the canonical formalism.
function A_canonical(chain, angles, two_ms)
    @unpack Hij, HRk = chain
    two_LS = HRk.two_ls
    two_ls = Hij.two_ls
    @unpack angles_Hij, angles_HRk = angles

    @unpack k = chain
    i, j = ij_from_k(k)

    @unpack two_js = chain.tbs
    two_ji, two_jj, two_jk, two_j0 = two_js[i], two_js[j], two_js[k], two_js[4]
    two_mi, two_mj, two_mk, two_m0 = two_ms[i], two_ms[j], two_ms[k], two_ms[4]

    @unpack two_j = chain
    value = sum(-two_j:2:two_j) do two_m
        A_node_canonical(angles_HRk, (two_m, two_mk, two_m0), (two_j, two_jk, two_j0), two_LS) * A_node_canonical(
            angles_Hij,
            (two_mi, two_mj, two_m),
            (two_ji, two_jj, two_j),
            two_ls,
        )
    end

    matching = 1 / sqrt(two_j0 + 1) # ThreeBodyDecays implementation does not have sqrt(two_j0+1) factor
    return value * matching
end;

# The `unpolarized_intensity_canonical` function calculates the unpolarized intensity for the canonical formalism.
unpolarized_intensity_canonical(chain, angles) =
    sum(itr(chain.tbs.two_js)) do two_ms
        A = A_canonical(chain, angles, two_ms)
        abs2(A)
    end;

#md # ## Example

# In this example, we use a decay process, `0` → `1` `2` `3` with spins `3/2` → `1/2` `1` `0`, respectively.
# A decay chain `(23)1` with `jp=2-` permit 8 different helicity couplings.
# The equivalent between formalism is checked separately for each of them.

# Create a three-body system with specified masses and spins.
tbs = ThreeBodySystem(
    ThreeBodyMasses(1, 2, 3; m0 = 6.32397),
    ThreeBodySpins(1, 2, 0; two_h0 = 3),
);

# Define a decay chain with LS couplings.
decay_chains = DecayChainsLS(;
    k = 1,
    Xlineshape = x -> 1.0,
    jp = jp"2-",
    Ps = ThreeBodyParities('+', '+', '+'; P0 = '+'), # need parities
    tbs,
); # give three-body-system structure

# Generate a random point in the phase space.
dpp = randomPoint(tbs);
(masses = dpp.σs, spin_projections = dpp.two_λs)

# Define the angles for the decay.
# Here, there is something cool 🔥: for canonical formulation, `RBR⁻¹` is applied when going to a subchannel rest frame.
# Hence, the θ_ij is that is passed to Y_lm is actually the sum of both polar helicity angle, θ_Rk and θ_ij.
# Going over π might be subtle.

angles = let
    cosθ_Rh = 0.7
    cosθ_ij = cos(acos(cosθ_Rh) + acos(cosθ23(dpp.σs, tbs.ms^2)))

    angles_HRk = (ϕ = 0.3, cosθ = cosθ_Rh) # ϕ is arbitrary
    angles_Hij = (ϕ = 0.5, cosθ = cosθ_ij) # ϕ is arbitrary
    (; angles_Hij, angles_HRk)
end

# Calculate and compare the unpolarized intensities.

let
    Ic = [unpolarized_intensity_canonical(dc, angles) for dc in decay_chains]
    Ih = [unpolarized_intensity(dc, dpp.σs) for dc in decay_chains]
    Ic ./ Ih
end
