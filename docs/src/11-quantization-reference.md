# [Quantization reference frames](@id quantization_reference)

Particles with spin quantized using the helicity basis require an additional specification: the **reference frame** in which the spin projection (helicity) is defined. This subtlety is unique to multi-channel amplitude analyses where several decay topologies are added coherently to form the full amplitude.

This section explains the concept of the reference kinematics and its implementation in `ThreeBodyDecays.jl`.
The formalism is presented in detail in the paper
*"Wigner rotations for cascade reactions"* by K. Habermann and M. Mikhasenko, ([INSPIRE](https://inspirehep.net/literature/2827198)).

## The problem: helicity is frame-dependent

In the helicity formalism, the spin state of a particle is described by its projection onto the direction of motion. However, for a given three-body decay ``0 \to 1\,2\,3``, the "direction of motion" of a final-state particle depends on *which rest frame* one considers.

Each decay topology ``k = 1, 2, 3``, referring to ``(23)1``, ``(31)2``, ``(12)3`` respectively, defines a specific sequence of Lorentz boosts from the parent rest frame to the rest frames of the final-state particles. These sequences are different for different topologies, so the same particle can have its helicity defined with respect to different quantization axes depending on which topology is used.

When amplitudes from different decay chains are added together, it is essential that the helicity states of the final-state particles refer to the **same quantization frame**. If they do not, the Wigner rotation matrix ``d^j_{\lambda'\lambda}(\zeta)`` must be applied to rotate from one frame to another.

## Natural reference frames for each topology

Each decay topology ``k`` naturally defines a reference frame for every particle in the decay. These frames are reached by a specific chain of Lorentz boosts from the parent rest frame, following the Jacob-Wick particle-2 convention.

### Topology ``k = 1``: decay ``0 \to (23)\, 1``

For the reference topology ``k = 1``, the natural quantization frames for each particle are:

| Particle | Path from parent rest frame |
|:--------:|:----------------------------|
| **2**    | ``(1,2,3) \to (23) \to 2`` |
| **3**    | ``(1,2,3) \to (23) \to R_y(\pi) \to 3`` |
| **1**    | ``(1,2,3) \to R_y(\pi) \to 1`` |
| **0**    | Parent rest frame (trivial) |

Here, the arrow ``\to`` stands for the helicity transformation ``[Rz Ry B]^{-1}`` and ``R_y(\pi)`` denotes the rotation by ``\pi`` about the ``y``-axis, which implements the approach to the particle state using the **Jacob-Wick particle-2 convention**: instead of going directly to the particle's rest frame, one first applies a ``\pi``-rotation to align with the convention where particle 2 in a two-body decay defines the positive ``z``-direction.

The same logic applies to topologies ``k = 2`` and ``k = 3`` with cyclic permutations of the particle indices.

## Aligned amplitude

Using the helicities defined in the natural frame of a given topology is called the **aligned amplitude** in this package.
The function [`aligned_amplitude`](@ref) computes the decay amplitude for a single chain ``k`` with all spin projections quantized along the natural axes of that chain:

```julia
F = aligned_amplitude(dc, σs)
```

This returns a multidimensional array ``F_{\lambda_1' \lambda_2' \lambda_3' \lambda_0'}`` where the primed helicities are defined in the natural frame of the decay chain `dc`.

Since each chain has its *own* natural frame, the aligned amplitudes from different chains cannot be directly added. To combine them coherently, one must rotate all helicities to a common reference frame.

## The full amplitude with Wigner rotations

The full helicity amplitude [`amplitude`](@ref) performs this rotation automatically.
Given a decay chain with topology ``k`` and a choice of reference frames specified by `refζs`, the amplitude is:

```math
A_{\lambda_1 \lambda_2 \lambda_3 \lambda_0} =
\sum_{\lambda_0' \lambda_1' \lambda_2' \lambda_3'}
d^{j_0}_{\lambda_0 \lambda_0'}(\zeta_0) \;
F_{\lambda_1' \lambda_2' \lambda_3' \lambda_0'} \;
d^{j_1}_{\lambda_1' \lambda_1}(\zeta_1) \;
d^{j_2}_{\lambda_2' \lambda_2}(\zeta_2) \;
d^{j_3}_{\lambda_3' \lambda_3}(\zeta_3)
```

The Wigner ``d``-matrices rotate from the aligned (primed) helicities of chain ``k`` to the helicities defined in the reference frames specified by `refζs`. The rotation angles ``\zeta`` depend on the kinematics and are computed automatically.

## The `refζs` parameter

The keyword argument `refζs` in the [`amplitude`](@ref) function specifies which topology defines the quantization frame for each particle:

```julia
amplitude(dc, σs, two_λs; refζs = (1, 1, 1, 1))
```

The tuple `refζs = (r₁, r₂, r₃, r₀)` means:

- Particle 1 has its helicity quantized in the natural frame of topology ``r_1``
- Particle 2 has its helicity quantized in the natural frame of topology ``r_2``
- Particle 3 has its helicity quantized in the natural frame of topology ``r_3``
- Particle 0 (parent) has its helicity quantized in the natural frame of topology ``r_0``

### Common choices

| `refζs` | Description |
|:--------|:------------|
| `(1, 1, 1, 1)` | All particles quantized using topology ``k=1``. All helicities refer to the frames defined by the ``0 \to (23)\,1`` chain. |
| `(1, 2, 3, 1)` | Each final-state particle uses its own "natural" topology (particle ``n`` uses topology ``n``). The parent uses topology 1. This is the **default** in `ThreeBodyDecays.jl`. |

The choice `refζs = (1, 1, 1, 1)` is often the most intuitive: it defines a single, common reference topology for all particles. This is equivalent to saying that the helicity states ``|\lambda_1, \lambda_2, \lambda_3, \lambda_0\rangle`` are defined as if the decay always proceeds through the ``k=1`` chain.

!!! note "Consistency across chains"
    The critical requirement is that `refζs` must be the **same** for all chains that are added coherently in a model. The specific choice of `refζs` does not affect physical observables (such as the unpolarized intensity), but it changes the meaning of the individual helicity couplings. Different choices of `refζs` correspond to different helicity bases.

## Example

```julia
using ThreeBodyDecays

ms = ThreeBodyMasses(1.0, 2.0, 3.0; m0 = 7.0)
tbs = ThreeBodySystem(; ms, two_js = ThreeBodySpins(1, 2, 0; two_h0 = 1))
Ps = ThreeBodyParities('+', '+', '+'; P0 = '+')
Xlineshape = σ -> 1.0

# Two chains with different topologies
dc1 = DecayChainLS(; k = 1, Xlineshape, jp = jp"1/2-", Ps, tbs)
dc3 = DecayChainLS(; k = 3, Xlineshape, jp = jp"1/2-", Ps, tbs)

σs = randomPoint(ms)
two_λs = (1, 0, 0, -1)

# All helicities quantized in topology-1 frames
A1 = amplitude(dc1, σs, two_λs; refζs = (1, 1, 1, 1))
A3 = amplitude(dc3, σs, two_λs; refζs = (1, 1, 1, 1))

# The sum is meaningful because both use the same reference
A_total = A1 + A3
```

## Further reading

The concept of Wigner rotations for cascade decays, including alternative conventions (minus-phi and canonical), is discussed thoroughly in:

- K. Habermann, M. Mikhasenko, *"Wigner rotations for cascade reactions"*,
  [INSPIRE:2827198](https://inspirehep.net/literature/2827198)

The original Dalitz-plot decomposition formalism is described in:

- M. Mikhasenko et al. (JPAC), *"Dalitz-plot decomposition for three-body decays"*,
  [INSPIRE:1758460](https://inspirehep.net/record/1758460)
