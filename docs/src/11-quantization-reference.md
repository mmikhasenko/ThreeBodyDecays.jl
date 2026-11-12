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
ThreeBodyDecays.reference_k(model) # 1
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

## Reference conventions and couplings

The reference topology is part of the model's amplitude convention, together
with its coupling coefficients. Changing the reference can require corresponding
phase changes in the couplings, in addition to rotating the external spin basis.
Keeping the same numerical couplings while changing the reference therefore does
not generally describe the same physical model.

The package does not convert couplings between reference conventions automatically.
Record the reference alongside serialized models and use the same topology when
converting four-vectors into orientation angles. Custom `AbstractDecayChain`
implementations used inside a model must honor its `reference_k`.

<!-- cspell:ignore pmatrix mathrm nondegenerate spinor -->

## CD–TBD signs and reference conversion

For the ordered cyclic chains `1 => ((2,3),1)`, `2 => ((3,1),2)`, and
`3 => ((1,2),3)`, let ``j_a`` be the spin of final particle ``a`` and define
``a_a=(-1)^{2j_a}``. The sign relating CascadeDecays (CD) and the present
ThreeBodyDecays (TBD) convention is

```math
S(k,r)=\begin{cases}
1,&k=r,\\
(-1)^{2j_\ell},&\ell\notin\{k,r\},
\end{cases}
\qquad
S=\begin{pmatrix}
1&a_3&a_2\\
a_3&1&a_1\\
a_2&a_1&1
\end{pmatrix}.
```

Rows label the chain ``k``; columns label the common reference ``r``.
For matching local vertices, cyclic daughter order, particle-2 convention,
normalization, and external helicity labels, the full complex tensors obey

```math
A_{\mathrm{CD}}^{k;r}
=S(k,r)\,e^{i\lambda_r\gamma_r}A_{\mathrm{TBD}}^{k;r}.
```

Here ``\lambda_r`` is the **physical helicity**, not its doubled integer, of
reference spectator ``r``. The angle ``\gamma_r`` is the inner azimuth in the
reference topology's helicity construction (the third Euler angle in the
corresponding factorization). This phase is common to every chain at fixed
reference and helicities. It cancels from unpolarized intensities but must be
retained when comparing polarized tensors. See the
[Cascade comparison note](https://github.com/RUB-EP1/CascadeDecays.jl/blob/e72eb7c44e378be7a549f8e2921c60511eef6d79/docs/cascade-vs-dpd.qmd)
for the frame construction. In that revision, the explanatory inner conjugated-D
exponent should have a plus sign, ``+i\lambda_R\gamma_r``.

This relation uses the existing TBD principal-angle convention, without the
alternative parent full-turn prescription proposed in
[PR #127](https://github.com/mmikhasenko/ThreeBodyDecays.jl/pull/127).
It applies to massive, nondegenerate kinematics with consistent helicity-frame
lifts. At exact Euler poles or degenerate boundaries, use a coordinated limiting
frame; independent Euler-angle choices need not preserve the spinor sign.

For a coherent model, the coefficients therefore convert as

```math
c_{\mathrm{TBD},n}^{(r)}=S(k_n,r)c_{\mathrm{CD},n},\qquad
c_{\mathrm{TBD},n}^{(s)}=S(k_n,s)S(k_n,r)c_{\mathrm{TBD},n}^{(r)}.
```

The index ``n`` labels a chain, so the same sign applies to every resonance and
LS coupling with topology ``k_n``. The latter formula converts an existing TBD
model from reference ``r`` to ``s``. It preserves unpolarized intensities;
polarized tensors additionally require transforming their external spin bases.
They are not expected to agree component by component in different frames.

For example, the following explicit conversion constructs a new model without
changing the source model. It applies to the standard `DecayChain` convention;
custom chain implementations must establish their own convention first.

```julia
function convert_reference(model::ThreeBodyDecay, s::Integer)
    s in 1:3 && !(s isa Bool) || throw(ArgumentError("reference must be 1, 2, or 3"))
    r = ThreeBodyDecays.reference_k(model)
    two_js = spins(model) # doubled spins
    sign_for(k, ref) = k == ref ? 1 : (-1)^two_js[6-k-ref]
    couplings = map(model.chains, model.couplings) do chain, c
        sign_for(chain.k, s) * sign_for(chain.k, r) * c
    end
    return ThreeBodyDecay(model.chains, couplings, model.names; reference_k = s)
end

converted = convert_reference(model, 3)
# unpolarized_intensity(converted, σs) ≈ unpolarized_intensity(model, σs)
```

For example, with only particle 1 carrying half-integer spin, the signs for
references 1, 2, and 3 are respectively ``(+,+,+)``, ``(+,+,-)``, and ``(+,-,+)``
as vectors over chain topologies. Converting from reference 1 to reference 2
negates topology-3 couplings; converting to reference 3 negates topology-2
couplings. Relabeling the reference without these changes generally changes
interference.

### Why a full turn is not an additional default correction

The three signed principal parent rotations around the spectator directions
accumulate a full ``2\pi`` turn. Its spin representation is ``(-1)^{2j_0}``.
Consequently, the direct short rotation and the composition through another
topology can differ by that central sign: they reach the same spatial orientation
along different rotation paths.

PR #127 chose topology 1 as an anchor and added a parent full-turn sign on the
``2\leftrightarrow3`` pair. This gives consistently composed spin-frame maps,
but changes the coupling convention. For a half-integer parent it negates chain
3 at reference 2 and chain 2 at reference 3; reference 1 is unchanged.
Its CD conversion signs become ``g_k g_r``, where
``g=(1,(-1)^{2j_3},(-1)^{2j_2})``. Both conventions can describe the same physical
model when their couplings are converted.

The package retains the symmetric ``S(k,r)`` rule above and binds the reference
to the model. There is no extra full-turn factor for users to apply during
model evaluation. The old expectation that one could freely change the reference
while keeping all numerical couplings was too strong. Tests instead compare
models after the explicit coupling conversion.

## Further reading

The concept of Wigner rotations for cascade decays, including alternative conventions (minus-phi and canonical), is discussed thoroughly in:

- K. Habermann, M. Mikhasenko, *"Wigner rotations for cascade reactions"*,
  [INSPIRE:2827198](https://inspirehep.net/literature/2827198)

The original Dalitz-plot decomposition formalism is described in:

- M. Mikhasenko et al. (JPAC), *"Dalitz-plot decomposition for three-body decays"*,
  [INSPIRE:1758460](https://inspirehep.net/record/1758460)
