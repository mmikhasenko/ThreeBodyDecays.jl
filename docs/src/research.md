# Research Documentation: Dalitz to Unit Circle Mapping

## Democratic Coordinates

Based on the requirement for "cube-root-weighted linear combination", the democratic coordinate transformation is:

\[
z = \frac{1}{\Lambda} \sum_{i=1}^{3} (\sigma_i - \bar{\sigma}) \omega^{i-1}
\]

where:

- \(\sigma_i\) are the Mandelstam invariants: \(\sigma_1 = (p_2 + p_3)^2\), \(\sigma_2 = (p_3 + p_1)^2\), \(\sigma_3 = (p_1 + p_2)^2\)
- \(\bar{\sigma} = \frac{\Sigma}{3}\) where \(\Sigma = \sigma_1 + \sigma_2 + \sigma_3 = m_1^2 + m_2^2 + m_3^2 + m_0^2\)
- \(\omega = e^{2\pi i/3}\) is the cube root of unity
- \(\Lambda\) is a normalization parameter (to be determined)

### Properties

1. **Permutation covariance**: Cyclic permutation \((\sigma_1, \sigma_2, \sigma_3) \to (\sigma_2, \sigma_3, \sigma_1)\) rotates \(z\) by \(\omega\)
2. **D₃ symmetry**: For equal masses, the boundary in \(z\) has dihedral group D₃ symmetry
3. **Inverse**: The inverse transformation is linear and can be solved analytically given \(\Sigma\)

### Inverse Transformation

Given \(z\) and \(\Sigma\), we can recover \((\sigma_1, \sigma_2)\):

- \(\sigma_3 = \Sigma - \sigma_1 - \sigma_2\)
- The transformation is linear, so the inverse is straightforward

### Potential Issues

- **Singular cases**: Need to verify for degenerate mass configurations
- **Normalization**: \(\Lambda\) must be chosen to ensure the boundary maps appropriately

## Zipper Algorithm

### Overview

The zipper algorithm (also known as geodesic algorithm) is a method for conformal mapping of polygons to the unit disk. It was developed by Don Marshall and others.

### Key References

- Marshall, D. E., & Rohde, S. (2007). "The Zipper Algorithm for Conformal Mapping"
- Need to find Julia implementations or reference code

### Algorithm Requirements

1. **Input**: Polygon vertices (ordered, counter-clockwise)
2. **Output**: Conformal map from polygon interior to unit disk
3. **Method**: Iterative process using geodesic arcs

### Potential Deadends

1. **Curved boundaries**: Dalitz boundary is curved (not a polygon)
   - **Solution**: Discretize boundary into polygon approximation
   - **Alternative**: Use boundary integral method (Schwarz-Christoffel or Szegő kernel)

2. **Parametric vs discrete**: Zipper typically requires discrete vertices
   - **Solution**: Discretize parametric boundary on-demand during zipper construction
   - Store the resulting mapping function, not the discretization

3. **Julia packages**: Need to check for existing implementations
   - Search for: `ConformalMaps.jl`, `Zipper.jl`, `SchwarzChristoffel.jl`

### Alternative Approaches

If zipper algorithm is not suitable:

1. **Schwarz-Christoffel mapping**: For polygons, but may require adaptation
2. **Boundary integral method**: Szegő kernel or Kerzman-Stein kernel
3. **Numerical conformal mapping**: Using boundary element methods

## Boundary Parameterization

### Existing Implementation

The codebase already has a parametric boundary representation in `src/dalitz.jl`:

- `polardalitz2invariants(θ, expansion_point)`: Returns polynomials in \(r\) for each \(\sigma_i\)
- `rborder(θ)`: Finds the boundary radius by solving \(\phi(\sigma_s(\theta, r)) = 0\) where \(\phi\) is the Kibble function
- `σs_border(θ)`: Evaluates the boundary point at angle \(\theta\)

### Stability

- Uses `PolynomialRoots.roots()` for root finding
- Filters for real, positive roots
- Takes minimum root (closest to expansion point)

### Alternative Construction

For cross-validation, can construct boundary via:

- Cosine scan: Use \(\cos\theta_{ij} = \pm 1\) in sequential frame
- This gives analytic \(s_2^{\pm}(s_1)\) branches

## Möbius Transformation for Normalization

To map three points \(w_1, w_2, w_3\) to \(1, \omega, \omega^2\):

\[
A(w) = \frac{(w - w_1)(w_2 - w_3)}{(w - w_3)(w_2 - w_1)} \cdot \frac{(\omega - \omega^2)}{(1 - \omega^2)} \cdot \frac{(1 - \omega)}{(\omega - \omega^2)}
\]

Simplified form:
\[
A(w) = \frac{(w - w_1)(w_2 - w_3)}{(w - w_3)(w_2 - w_1)} \cdot \frac{1 - \omega}{1 - \omega^2}
\]

### Well-Conditioned Check

- Need to verify when \(w_2 \approx w_1\) or \(w_3 \approx w_1\) (degenerate cases)
- For equal masses, landmarks should be well-separated due to symmetry

## Landmark Selection

### Method

Find three boundary points that maximize each \(\sigma_i\):

- \(\theta_1 = \arg\max_\theta \sigma_1(\theta)\)
- \(\theta_2 = \arg\max_\theta \sigma_2(\theta)\)
- \(\theta_3 = \arg\max_\theta \sigma_3(\theta)\)

### Uniqueness

- For generic masses: Should be unique
- For equal masses: May have symmetry, need deterministic tie-breaking
- Use optimization (e.g., Brent's method) along parametric boundary

## Numerical Tolerances

### Boundary Closure

- Tolerance: \(10^{-10}\) for \(|boundary(0) - boundary(2\pi)|\)

### Conformality

- Cauchy-Riemann residuals: \(< 10^{-8}\) on interior grid

### Normalization

- Landmark mapping: \(|A(w_i) - target_i| < 10^{-6}\) where target ∈ {1, ω, ω²}

### Convergence

- Zipper discretization: Error should decrease with n=200, 400, 800, 1600

## Deadend Analysis

### Critical Deadends

1. **Zipper with curved boundary**:
   - **Risk**: High - zipper is designed for polygons
   - **Mitigation**: Discretize boundary, verify convergence
   - **Fallback**: Boundary integral method

2. **Democratic coordinate inverse instability**:
   - **Risk**: Low - transformation is linear
   - **Mitigation**: Verify for all physical configurations

3. **Landmark selection ambiguity**:
   - **Risk**: Medium - for symmetric cases
   - **Mitigation**: Deterministic tie-breaking, optimization with constraints

4. **Near-threshold numerical issues**:
   - **Risk**: Medium - polynomial root finding may be unstable
   - **Mitigation**: Use existing `border()` which handles this

### Action Items

- [ ] Find Julia package for zipper or conformal mapping
- [ ] Verify democratic coordinate formula with literature
- [ ] Test boundary parameterization stability for edge cases
- [ ] Implement fallback plan if zipper fails
