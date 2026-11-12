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

<!-- cspell:ignore pmatrix mathrm nondegenerate spinor spinors zeta perp lambda boldsymbol ispositive noncollinear mathtt leftarrow leftrightarrow notin -->

## Planar rotations and the parent-angle sum rule

Work in the parent rest frame, first in the canonical plane and then oriented
relative to external axes. The overall orientation below is a rotation, not a
boost of the parent into a laboratory frame.
Momentum conservation, ``\boldsymbol p_1+\boldsymbol p_2+\boldsymbol p_3=0``,
puts the three momenta in one plane. TBD takes this to be the ``xz`` plane,
with its normal along ``y``. For each cyclic topology ``k=((i,j),k)``, the
canonical configuration has the isobar along ``+z``, spectator ``k`` along
``-z``, and daughter ``i`` in the positive-``x`` half-plane.

The internal transformations preserve this common plane: rotations are about
``y`` and boosts have their velocities in the ``xz`` plane. The resulting
Wigner rotations also have their axis along ``y``. Thus every external-leg
alignment is described by a **signed angle** and a real small-``d`` matrix,

```math
R_y\!\left(\zeta^a_{k(r)}\right),\qquad
 d^{j_a}_{\lambda'\lambda}\!\left(\zeta^a_{k(r)}\right),
\qquad a=0,1,2,3.
```

Here ``a=0`` denotes the parent, ``k`` the chain, and ``r`` the reference.
In the implementation, `cosζ(wr(k,r,a), σs, ms^2)` supplies the cosine and
`ispositive(wr(k,r,a))` supplies the direction. Each signed angle is chosen in
``[-\pi,\pi]``. There is no independently accumulated winding number in this
pairwise prescription.

For the parent, the positive cyclic angles are the angles between spectator
momenta:

```math
\cos\zeta^0_{i(j)}
=\frac{\boldsymbol p_i\cdot\boldsymbol p_j}
 {|\boldsymbol p_i|\,|\boldsymbol p_j|},
\qquad (i,j)=(1,2),(2,3),(3,1).
```

For a noncollinear event they lie strictly between zero and ``\pi``. The three
momentum directions surround the origin, so their successive angular gaps obey

```math
\boxed{\zeta^0_{1(2)}+\zeta^0_{2(3)}+\zeta^0_{3(1)}=2\pi},
\qquad
\zeta^0_{j(i)}=-\zeta^0_{i(j)}.
```

No individual angle crosses ``2\pi``; the **cyclic sum is identically** ``2\pi``.
For example,

```math
\zeta^0_{1(2)}+\zeta^0_{2(3)}
=\zeta^0_{1(3)}+2\pi.
```

The direct principal rotation and the two-step rotation therefore reach the
same spatial axes, but differ by a full turn along their rotation paths. In a
spin representation this means

```math
d^{j_0}\!\left(\zeta^0_{1(2)}\right)
d^{j_0}\!\left(\zeta^0_{2(3)}\right)
=(-1)^{2j_0}\,
d^{j_0}\!\left(\zeta^0_{1(3)}\right).
```

This is the spin representation of the stated angle prescription, not a loss
of planarity. In particular, one must not replace a direct principal rotation
by a composed path while discarding its full-turn sign.

### Negative angles and Euler representations

Use the convention
``R(\phi,\theta,\psi)=R_z(\phi)R_y(\theta)R_z(\psi)``.
Writing ``\theta=|\zeta|``, the planar rotation has the Euler representations

```math
R_y(\zeta)=
\begin{cases}
R(0,\theta,0),&\zeta\geq0,\\
R(\pi,\theta,-\pi),&\zeta<0.
\end{cases}
```

The second equality is exact even for spinors:
``R_z(\pi)R_y(\theta)R_z(-\pi)=R_y(-\theta)``.
It introduces **no** extra full-turn sign. The apparent ``z`` rotations are
only an Euler representation of a pure ``y`` rotation; the net transformation
preserves the chosen plane normal.

By contrast, changing the last angle from ``-\pi`` to ``+\pi`` gives

```math
D^j(\pi,\theta,\pi)
=(-1)^{2j}D^j(\pi,\theta,-\pi)
=(-1)^{2j}d^j(-\theta).
```

These rotations agree as ordinary three-dimensional rotation matrices, but
have opposite signs for half-integer spin. Reducing Euler angles modulo
``2\pi`` without retaining this sign loses information needed by amplitudes.

## Overall orientation and the reference-spectator phase

Let ``r=\mathtt{reference\_k}`` and let
``\Omega_r=(\alpha_r,\beta_r,\gamma_r)`` orient its canonical planar event:

```math
R(\Omega_r)=R_z(\alpha_r)R_y(\beta_r)R_z(\gamma_r).
```

The angles are those of the **reference topology**, not those of each term in
the coherent sum. If ``\mathcal A^{k;r}`` is the planar chain tensor, TBD forms

```math
A_{\mathrm{TBD},\,m_0,\boldsymbol\lambda}^{k;r}(\Omega_r,\sigma)
=\sum_{\nu=-j_0}^{j_0}
D^{j_0*}_{m_0\nu}(\alpha_r,\beta_r,\gamma_r)
\mathcal A_{\boldsymbol\lambda,\nu}^{k;r}(\sigma).
```

Here ``m_0`` is the parent projection in the external frame, ``\nu`` its
projection in the reference planar frame, and
``\boldsymbol\lambda=(\lambda_1,\lambda_2,\lambda_3)`` the final helicities.
The API is `amplitude(model, (α=αr, cosβ=cos(βr), γ=γr), σs)`.
The overall rotation acts on the parent index; it is distinct from the
internal small-``d`` alignments.

To compare with the usual two-body helicity construction
``R(\phi_A,\theta_A,0)`` for ``R\to A B``, fix the Wigner convention

```math
D^j_{mn}(\alpha,\beta,\gamma)
=e^{-im\alpha}d^j_{mn}(\beta)e^{-in\gamma}.
```

Consider the reference chain ``0\to R+r``, ``R\to i+j``. Its ordinary
helicity construction contains the angular factors

```math
D^{j_0*}_{m_0,\lambda_R-\lambda_r}(\alpha_r,\beta_r,0)
D^{j_R*}_{\lambda_R,\lambda_i-\lambda_j}(\gamma_r,\theta_{ij},0).
```

The inner azimuth ``\gamma_r`` specifies the daughter plane after the first
rotation. With ``\mu=\lambda_R-\lambda_r``,

```math
D^{j_R*}_{\lambda_R,\lambda_i-\lambda_j}(\gamma_r,\theta_{ij},0)
=e^{+i\lambda_R\gamma_r}
 d^{j_R}_{\lambda_R,\lambda_i-\lambda_j}(\theta_{ij}),
```

whereas the third angle of the overall rotation contributes
``e^{+i\mu\gamma_r}``. Consequently,

```math
\begin{aligned}
&D^{j_0*}_{m_0,\lambda_R-\lambda_r}(\alpha_r,\beta_r,0)
D^{j_R*}_{\lambda_R,\lambda_i-\lambda_j}(\gamma_r,\theta_{ij},0)\\
&\quad=
D^{j_0*}_{m_0,\lambda_R-\lambda_r}(\alpha_r,\beta_r,\gamma_r)
e^{+i\lambda_r\gamma_r}
 d^{j_R}_{\lambda_R,\lambda_i-\lambda_j}(\theta_{ij}).
\end{aligned}
```

Thus the TBD overall ``D^*`` factor absorbs the inner azimuth except for
``e^{+i\lambda_r\gamma_r}``. This remaining phase involves the **reference
spectator helicity**, not the parent projection ``\nu`` or ``m_0``. All
helicities in these equations are physical spin projections, not the doubled
integers used in parts of the Julia API. TBD does not include this spectator
phase in its oriented-amplitude API; include it when converting to the CD
helicity convention.

## How cascade frame paths produce the discrete signs

TBD factors out the reference orientation once, then uses only signed planar
``R_y`` alignments. A cascade construction instead applies
``R(\phi_A,\theta_A,0)`` and a boost at each vertex, and aligns the resulting
final-state helicity frames with those of the reference chain. Its relative
Wigner matrices generally have three Euler angles.

The geometric decay plane is the same in both constructions. What can differ
is the orientation of the local axes attached to it: an azimuthal ``\pi``
rotation reverses both ``x`` and ``y``. Thus a cascade path can reverse a local
plane normal and later restore it. Products of such rotations can leave an
extra ``R_z(2\pi)`` (equivalently the central full turn about any axis) relative
to the planar path. It is invisible to spatial axes but acts as
``(-1)^{2j}`` on spin ``j``: minus for half-integer spin, plus for integer spin.
This statement concerns the chosen frame paths, not a physical reflection of
the event.

Before extracting these signs, the continuous spectator phases must also be
accounted for. An explicit matrix comparison makes the separation precise.
Let ``C_a(k\leftarrow r)`` be the CD conjugated Wigner matrix on final leg
``a``, with rows in the chain frame and columns in the reference frame, and set

```math
P_a(q)_{mm'}=\delta_{mm'}e^{i m\gamma_q\delta_{aq}}.
```

With ``d_a(k,r)`` denoting the TBD signed small-``d`` alignment in the same
index order, the corresponding path comparisons have the form

```math
\begin{aligned}
D^{j_0*}(\Omega_k)
 &=\epsilon_0^{\,2j_0}D^{j_0*}(\Omega_r)d_0(k,r),\\
P_a(k)C_a(k\leftarrow r)
 &=\epsilon_a^{\,2j_a}d_a(k,r)P_a(r),\qquad a=1,2,3,
\end{aligned}
\qquad \epsilon_a\in\{+1,-1\}.
```

The ``\epsilon_a`` are the residual signs of the spin-``1/2`` frame maps.
For a direct numerical check, compare the two sides at spin ``1/2`` before
inserting ``\epsilon_a``: their relative matrix must be ``+I`` or ``-I``.
A noncentral residual is not a sign ambiguity; it indicates unmatched frame
conventions, kinematics, or a numerical problem. Individual leg signs can
change with Euler branches; their physical product below does not.

### Example: chain ``(12)3`` with reference ``(23)1``

Here ``k=3`` and ``r=1``. Both calculations express the final helicities in
the ``(23)1`` reference frames. In CD the chain starts from its own two-vertex
angles, ``\Omega_3`` and ``\theta_{12}``, and its final legs are aligned by
``C_a(3\leftarrow1)``. In TBD the parent orientation is ``\Omega_1`` and the
chain is aligned by the planar ``d_a(3,1)`` matrices. The path identities above
relate these two ways of assembling the same event.

For a continuous Euler branch one can obtain
``(\epsilon_0,\epsilon_1,\epsilon_2,\epsilon_3)=(+,+,-,+)``:
the remaining full-turn sign sits on particle 2. A different branch can
redistribute the leg signs, so this assignment is not itself an invariant.
The complete chain relation is

```math
\boxed{A_{\mathrm{CD}}^{3;1}
=(-1)^{2j_2}e^{i\lambda_1\gamma_1}A_{\mathrm{TBD}}^{3;1}}.
```

The discrete sign is therefore negative when particle 2 has half-integer spin,
regardless of whether other particles also have half-integer spin. It is not
in general a sign determined by the parent alone.

## General CD–TBD conversion recipe

Match the cyclic daughter order, local vertex/LS conventions, particle-2
phases, line shapes, and normalization first. The product of the four external
spin signs is

```math
S(k,r)=\epsilon_0^{2j_0}\prod_{a=1}^3\epsilon_a^{2j_a}.
```

Angular-momentum coupling requires
``2j_0=2j_1+2j_2+2j_3\pmod 2``. Define the reduced leg signs
``q_a=\epsilon_0\epsilon_a``; then

```math
S(k,r)=\prod_{a=1}^3 q_a^{2j_a}.
```

Comparing the cyclic frame paths gives the following reduced signs. Unlike an
individual ``\epsilon_a``, these combinations are unchanged when a regular
Euler branch change flips all four spinor signs together.

| Chain/reference pair | ``q_1`` | ``q_2`` | ``q_3`` |
| :--------------------- | :-------: | :-------: | :-------: |
| ``1\leftrightarrow2`` | ``+1`` | ``+1`` | ``-1`` |
| ``1\leftrightarrow3`` | ``+1`` | ``-1`` | ``+1`` |
| ``2\leftrightarrow3`` | ``-1`` | ``+1`` | ``+1`` |
| ``k=r`` | ``+1`` | ``+1`` | ``+1`` |

Only the particle absent from the chain/reference pair has a negative reduced
sign. This gives the symmetric rule

```math
\boxed{S(k,r)=\begin{cases}
1,&k=r,\\
(-1)^{2j_\ell},&\ell\notin\{k,r\},
\end{cases}}
\qquad
S=\begin{pmatrix}
1&a_3&a_2\\
a_3&1&a_1\\
a_2&a_1&1
\end{pmatrix},\quad a_i=(-1)^{2j_i}.
```

Rows label chain ``k`` and columns reference ``r``. This covers multiple
half-integer spins as well as integer spins; it does not depend on the chosen
spin magnitudes or allowed LS couplings. For every chain at a fixed reference,

```math
\boxed{A_{\mathrm{CD}}^{k;r}
=S(k,r)e^{i\lambda_r\gamma_r}A_{\mathrm{TBD}}^{k;r}}.
```

To convert a coherent CD model, multiply each coefficient by the sign for its
chain topology. To change the reference of a TBD model, apply the signs for
both the old and new references:

```math
c_{\mathrm{TBD},n}^{(r)}=S(k_n,r)c_{\mathrm{CD},n},\qquad
\boxed{c_{\mathrm{TBD},n}^{(s)}
=S(k_n,s)S(k_n,r)c_{\mathrm{TBD},n}^{(r)}}.
```

The index ``n`` labels a chain, including its resonance and LS coupling. The
spectator phase is common to all chains at fixed external helicities, so it
cancels from unpolarized intensities. Keep it for polarized comparisons; when
changing reference, also transform the external spin bases. Tensors in
different reference frames need not agree component by component.

The following explicit conversion constructs a new model without changing the
source. It applies to the standard `DecayChain` convention; custom chain
implementations must establish their own convention first.

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

This prescription retains the principal planar rotations. No additional
full-turn factor is inserted during TBD evaluation. The coupling conversion
is essential: changing only `reference_k` generally changes interference.
The formulas assume massive, nondegenerate events and consistently tracked
spin-frame paths. At exact Euler poles or degenerate boundaries, a common
limiting-frame prescription is required.

The [Cascade comparison note](https://github.com/RUB-EP1/CascadeDecays.jl/blob/e72eb7c44e378be7a549f8e2921c60511eef6d79/docs/cascade-vs-dpd.qmd)
contains explicit path and angle examples. In that revision the displayed
inner conjugated-``D`` exponent has a sign typo: the derivation above uses
``+i\lambda_R\gamma_r``, consistently with the stated ``D`` convention.

## Further reading

The concept of Wigner rotations for cascade decays, including alternative conventions (minus-phi and canonical), is discussed thoroughly in:

- K. Habermann, M. Mikhasenko, *"Wigner rotations for cascade reactions"*,
  [INSPIRE:2827198](https://inspirehep.net/literature/2827198)

The original Dalitz-plot decomposition formalism is described in:

- M. Mikhasenko et al. (JPAC), *"Dalitz-plot decomposition for three-body decays"*,
  [INSPIRE:1758460](https://inspirehep.net/record/1758460)
