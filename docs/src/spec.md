# Specification: Dalitz to Unit Circle Conformal Mapping

## Kinematic Domain D

### Definition

The physical Dalitz region \(D\) is defined as the set of Mandelstam invariants \((\sigma_1, \sigma_2, \sigma_3)\) satisfying:

1. **Kibble condition**: \(\phi(\sigma_1, \sigma_2, \sigma_3) < 0\) where
   \[
   \phi(\sigma_1, \sigma_2, \sigma_3) = \lambda(\lambda_1, \lambda_2, \lambda_3)
   \]
   with \(\lambda_i = \lambda(M^2, m_i^2, \sigma_i)\) being Källén functions.

2. **Physical ranges**: Each \(\sigma_i\) lies in its kinematic limits:
   \[
   (m_j + m_k)^2 \leq \sigma_i \leq (M - m_i)^2
   \]
   where \((i,j,k)\) is a cyclic permutation of \((1,2,3)\).

3. **Sum rule**: \(\sigma_1 + \sigma_2 + \sigma_3 = \Sigma = m_1^2 + m_2^2 + m_3^2 + M^2\)

### Boundary Parameterization

The boundary \(\partial D\) is parameterized by angle \(\theta \in [0, 2\pi)\):

\[
\sigma_i(\theta) = \sigma_i^{(0)} + r(\theta) \cdot p_i(\theta)
\]

where:

- \(\sigma_i^{(0)}\) is the expansion point (computed from masses)
- \(p_i(\theta)\) are direction polynomials:
    - \(p_1(\theta) = -\cos(\theta)\)
    - \(p_2(\theta) = \cos(\theta + \pi/3)\)
    - \(p_3(\theta) = \cos(\theta - \pi/3)\)
- \(r(\theta)\) is the boundary radius, found by solving \(\phi(\sigma_1(\theta, r), \sigma_2(\theta, r), \sigma_3(\theta, r)) = 0\) for the smallest positive real root

## Democratic Coordinate System

### Forward Transformation

The democratic coordinate \(z \in \mathbb{C}\) is defined as:

\[
z = \frac{1}{\Lambda} \sum_{i=1}^{3} (\sigma_i - \bar{\sigma}) \omega^{i-1}
\]

where:

- \(\bar{\sigma} = \frac{\Sigma}{3}\) is the mean invariant
- \(\omega = e^{2\pi i/3} = -\frac{1}{2} + \frac{\sqrt{3}}{2}i\) is the cube root of unity
- \(\Lambda\) is a normalization parameter chosen such that the boundary in \(z\) has appropriate scale

### Inverse Transformation

Given \(z\) and \(\Sigma\), we recover the invariants:

From \(\Lambda z = (\sigma_1 - \bar{\sigma}) + (\sigma_2 - \bar{\sigma})\omega + (\sigma_3 - \bar{\sigma})\omega^2\), we extract:

\[
\begin{align}
\sigma_1 &= \frac{2\text{Re}(\Lambda z) + \Sigma}{3} \\
\sigma_2 &= \frac{\Sigma - \sigma_1 + 2\text{Im}(\Lambda z)/\sqrt{3}}{2} \\
\sigma_3 &= \frac{\Sigma - \sigma_1 - 2\text{Im}(\Lambda z)/\sqrt{3}}{2}
\end{align}
\]

This follows from:

- Real part: \(\text{Re}(\Lambda z) = (3\sigma_1 - \Sigma)/2\)
- Imaginary part: \(\text{Im}(\Lambda z) = \sqrt{3}(\sigma_2 - \sigma_3)/2\)
- Sum constraint: \(\sigma_1 + \sigma_2 + \sigma_3 = \Sigma\)

### Normalization Parameter \(\Lambda\)

The parameter \(\Lambda\) is chosen to ensure the democratic coordinate boundary has reasonable scale. One option is to normalize by the maximum distance from center:

\[
\Lambda = \max_{\theta \in [0, 2\pi)} \left| \sum_{i=1}^{3} (\sigma_i(\theta) - \bar{\sigma}) \omega^{i-1} \right|
\]

## Conformal Mapping: z → w

### Zipper Algorithm

The zipper algorithm maps the boundary polygon (approximation of curved boundary) to the unit circle:

\[
w = f_{\text{zipper}}(z)
\]

where \(f_{\text{zipper}}\) is constructed via iterative geodesic arcs.

**Implementation note**: The boundary is discretized into \(n\) points for the zipper algorithm, but the resulting mapping function is stored and can be evaluated at any point.

### Alternative: Boundary Integral Method

If zipper is not suitable, use the boundary integral method with Szegő kernel:

\[
w = \frac{1}{2\pi i} \oint_{\partial D} \frac{S(z, \zeta)}{z - \zeta} d\zeta
\]

where \(S(z, \zeta)\) is the Szegő kernel.

## Channel-Democratic Normalization

### Landmark Selection

Select three boundary landmarks \(p_1, p_2, p_3\) by maximizing each invariant:

\[
\theta_i = \arg\max_{\theta \in [0, 2\pi)} \sigma_i(\theta), \quad i = 1, 2, 3
\]

with deterministic tie-breaking if multiple maxima exist.

### Möbius Transformation

The Möbius transformation \(A: \mathbb{D} \to \mathbb{D}\) mapping landmarks to \(1, \omega, \omega^2\) is:

\[
A(w) = \frac{(w - w_1)(w_2 - w_3)}{(w - w_3)(w_2 - w_1)} \cdot \frac{1 - \omega}{1 - \omega^2}
\]

where:

- \(w_i = f_{\text{zipper}}(\text{democratic}(p_i))\) are the mapped landmarks
- The transformation ensures: \(A(w_1) = 1\), \(A(w_2) = \omega\), \(A(w_3) = \omega^2\)

### Final Mapping

The complete mapping from Dalitz region to normalized unit disk is:

\[
w = A(f_{\text{zipper}}(\text{democratic}(\sigma_1, \sigma_2, \sigma_3)))
\]

## Numerical Tolerances

### Boundary Closure

\[
|boundary(0) - boundary(2\pi)| < 10^{-10}
\]

### Conformality (Cauchy-Riemann)

\[
\left| \frac{\partial u}{\partial x} - \frac{\partial v}{\partial y} \right| < 10^{-8}, \quad
\left| \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right| < 10^{-8}
\]
where \(w = u + iv\) and \(z = x + iy\).

### Normalization Exactness

\[
|A(w_i) - \text{target}_i| < 10^{-6}
\]
where \(\text{target}_i \in \{1, \omega, \omega^2\}\).

### Convergence

The zipper algorithm error should decrease monotonically with discretization:

- \(n = 200\): baseline
- \(n = 400\): error reduced by factor \(\geq 1.5\)
- \(n = 800\): error reduced by factor \(\geq 2.0\)
- \(n = 1600\): error reduced by factor \(\geq 2.5\)

## Coordinate Choices Summary

1. **Mandelstam invariants**: \((\sigma_1, \sigma_2, \sigma_3)\) with \(\sigma_1 + \sigma_2 + \sigma_3 = \Sigma\)
2. **Democratic coordinates**: \(z \in \mathbb{C}\) via cube-root-weighted linear combination
3. **Zipper mapping**: \(w_{\text{raw}} \in \mathbb{D}\) via conformal mapping
4. **Normalized coordinates**: \(w \in \mathbb{D}\) with landmarks at \(1, \omega, \omega^2\)

## API Specification

### Type-Stable Structs

```julia
struct BoundaryFunction{T}
    ms::MassTuple{T}
    expansion_point::MandelstamTuple{T}
    msq::NTuple{4,T}
end

struct DemocraticMap{T}
    ms::MassTuple{T}
    Λ::T
    Σ::T
end

struct ZipperMap{T}
    # Structure TBD after research
    # Will store zipper algorithm state
end

struct MobiusTransform{T}
    a::Complex{T}
    b::Complex{T}
    c::Complex{T}
    d::Complex{T}
end

struct DalitzMap{T}
    ms::MassTuple{T}
    boundary::BoundaryFunction{T}
    democratic::DemocraticMap{T}
    zipper::ZipperMap{T}
    normalization::MobiusTransform{T}
end
```

### Function Signatures

```julia
# Boundary
(bf::BoundaryFunction)(θ::Real) -> MandelstamTuple{T}

# Democratic coordinates
(dm::DemocraticMap)(σs::MandelstamTuple) -> Complex{T}
(dm::DemocraticMap)(z::Complex{T}) -> MandelstamTuple{T}  # inverse

# Zipper
(zm::ZipperMap)(z::Complex{T}) -> Complex{T}

# Möbius
(mt::MobiusTransform)(w::Complex{T}) -> Complex{T}

# Top-level
dalitzmap(ms::MassTuple; method=:zipper) -> DalitzMap{T}
(map::DalitzMap)(σs::MandelstamTuple) -> Complex{T}
(map::DalitzMap)(s12, s23) -> Complex{T}  # convenience method
```

## References

- Kibble, T. W. B. (1960). "Lorentz-invariant phase-space factors and Dalitz plots"
- Marshall, D. E., & Rohde, S. (2007). "The Zipper Algorithm for Conformal Mapping"
- Dalitz, R. H. (1953). "On the analysis of τ-meson data and the nature of the τ-meson"
