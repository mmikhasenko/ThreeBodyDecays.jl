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

| Particle | Path from parent rest frame                        |
|:--------:|:---------------------------------------------------|
| **2**    | ``(1,2,3) \to (23) \to 2``                         |
| **3**    | ``(1,2,3) \to (23) \to R_y(\pi) \to 3``            |
| **1**    | ``(1,2,3) \to R_y(\pi) \to 1``                     |
| **0**    | Parent rest frame (trivial)                        |

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
Every chain is rotated to the single reference topology stored in the model.
The cyclic choices are `1 => ((2,3),1)`, `2 => ((3,1),2)`, and `3 => ((1,2),3)`.
All four helicity indices use that topology. Wigner rotations depend on the
invariants and are computed automatically.

## Constructing and evaluating a model

```julia
model = ThreeBodyDecay(chains, couplings, names; reference_k = 1)
reference_k(model) # 1
A = amplitude(model, σs)
I = unpolarized_intensity(model, σs)
```

The keyword is required: fitted coefficients must have an explicit reference.
Evaluation cannot override the stored reference. Slicing preserves it, and
concatenating models with different references throws an `ArgumentError`.

For a single chain, `amplitude(dc, σs; reference_k = 1)` uses the same common
reference convention. The chain-only default is 1. Choosing `reference_k = dc.k`
reproduces `aligned_amplitude(dc, σs)`.

For an oriented event, use `amplitude(model, angles, σs)`, where `angles` contains
`α`, `cosβ`, and `γ`. The canonical reference event has the spectator along
negative z, the isobar along positive z, and the cyclic first daughter in the
positive-x half-plane. Its orientation is `Rz(α) Ry(β) Rz(γ)`; the parent spin
index transforms with the conjugate Wigner D matrix. The corresponding
`unpolarized_intensity(model, angles, σs)` sums all helicity components.

## Migrating to 0.16

The evaluation keyword `refζs` has been removed. An old common reference tuple
`(k,k,k,k)` becomes `reference_k = k` **when constructing the model**, rather than
when evaluating it. Direct chain calls accept the integer keyword instead.
Custom `AbstractDecayChain` implementations used inside models must also accept
and honor `reference_k` in their `amplitude` methods.

The old default `(1,2,3,1)` mixed different final-state reference frames. It has
no literal single-topology replacement. In particular, choosing topology 1 does
not preserve the old individual polarized amplitudes. Audit the convention of
stored fit coefficients and polarization information before migrating; the
package does not transform fitted coefficients automatically. Record the
reference alongside serialized models and use the same topology when converting
four-vectors into angles.

Every model carries one reference and its couplings are tied to it. Changing a
fit's reference convention may require SU(2) phase changes in its couplings; no
convenience operation copies a model into a new reference.

!!! note "Parent rotation between topologies 2 and 3"
    The three parent rotations between topologies share one axis and their angles
    add up to ``2π``. For a half-integer parent spin they cannot all be taken as
    the short rotation and still compose consistently. This release keeps the
    previous rotations unchanged, so for a half-integer parent the same couplings
    give different spin-summed intensities for references 2 and 3; equivalently,
    one topology's couplings differ in sign. How to fix this convention is
    discussed in #127.

## Further reading

The concept of Wigner rotations for cascade decays, including alternative conventions (minus-phi and canonical), is discussed thoroughly in:

- K. Habermann, M. Mikhasenko, *"Wigner rotations for cascade reactions"*,
  [INSPIRE:2827198](https://inspirehep.net/literature/2827198)

The original Dalitz-plot decomposition formalism is described in:

- M. Mikhasenko et al. (JPAC), *"Dalitz-plot decomposition for three-body decays"*,
  [INSPIRE:1758460](https://inspirehep.net/record/1758460)
