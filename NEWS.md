# Release notes

## 0.16 (unreleased)

### Reference topology

Construct models with an explicit reference:
`ThreeBodyDecay(chains, couplings, names; reference_k = 1)`.
The model retains this choice for amplitude evaluation, slicing, and overlap
calculations. Evaluation cannot override it, and concatenation rejects models
with different references.

The former `refζs` evaluation keyword is removed. A common reference tuple
`(k,k,k,k)` becomes `reference_k = k` at model construction. The old default
`(1,2,3,1)` used mixed reference frames and has no literal single-topology
replacement. Review the reference convention of fitted couplings before reuse;
changing it may require coupling-phase changes, which are not applied automatically.
See the [reference-frame guide](docs/src/11-quantization-reference.md).
