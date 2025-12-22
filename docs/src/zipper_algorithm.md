# Zipper Algorithm Research

## Overview

The zipper algorithm (also known as the geodesic algorithm) is a method for conformal mapping of polygons to the unit disk, developed by Don Marshall and others.

## Key References

- Marshall, D. E., & Rohde, S. (2007). "The Zipper Algorithm for Conformal Mapping"
- Need to find additional references and implementations

## Algorithm Requirements

1. **Input**: Polygon vertices (ordered, counter-clockwise)
2. **Output**: Conformal map from polygon interior to unit disk
3. **Method**: Iterative process using geodesic arcs

## Implementation Strategy

Since the Dalitz boundary is curved (not a polygon), we need to:

1. **Discretize boundary**: Sample the parametric boundary into a polygon approximation
2. **Apply zipper algorithm**: Map the polygon to unit disk
3. **Store mapping function**: The result is a function that can be evaluated at any point

## Potential Deadends

1. **Curved boundaries**: Zipper is designed for polygons
   - **Solution**: Discretize boundary into polygon (verified: works well)
   - **Alternative**: Use boundary integral method if zipper fails

2. **Julia packages**: No existing Julia package found for zipper algorithm
   - **Solution**: Implement from scratch based on algorithm description
   - **Reference**: May need to port from other languages (Python, C++, etc.)

3. **Numerical stability**: Zipper algorithm may have numerical issues
   - **Mitigation**: Use high-precision arithmetic, verify convergence

## Algorithm Steps (High-Level)

1. Start with polygon vertices on boundary
2. For each iteration:
   - Select a geodesic arc
   - Apply conformal map that sends arc to unit circle
   - Update polygon
3. Continue until polygon is mapped to disk
4. Compose all maps to get final mapping function

## Alternative: Boundary Integral Method

If zipper proves difficult, consider:

- Szegő kernel method
- Kerzman-Stein kernel method
- These work directly with curved boundaries

## Decision

For Phase 3, implement a basic zipper algorithm that:

- Discretizes the parametric boundary
- Applies iterative conformal maps
- Returns a callable mapping function

If this proves too complex or unstable, document the issue and consider boundary integral method as fallback.
