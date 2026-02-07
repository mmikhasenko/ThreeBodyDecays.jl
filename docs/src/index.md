```@meta
CurrentModule = ThreeBodyDecays
```

# ThreeBodyDecays

Documentation for [ThreeBodyDecays](https://github.com/mmikhasenko/ThreeBodyDecays.jl).

`ThreeBodyDecays.jl` is a Julia package for building hadronic decay models using a cascade reaction. Decays with three particles is the main application of the approach, however, it is also found useful in multibody decays where transitions can be factorized to a product of sequential decays with ≤3 products. Particles can have arbitrary spin.

The implementation is based on a research paper, "Dalitz-plot decomposition for three-body decays" by JPAC Collaboration (M. Mikhasenko at al.) [(inspire reference)](http://inspirehep.net/record/1758460).

## Documentation

- **[Energy dependence of the amplitude](@ref energy_dependence)** - Understanding how helicity amplitudes depend on energy kinematic variables
- **[Quantization reference frames](@ref quantization_reference)** - Reference kinematics and Wigner rotations for particles with spin
- **[Tutorial](30-tutorial.md)** - Getting started with ThreeBodyDecays.jl
- **[Canonical examples](31-canonical.md)** - Common usage patterns and examples
- **[Contributing guidelines](90-contributing.md)** - How to contribute to the project
- **[Developer documentation](91-developer.md)** - Technical details for developers
- **[Reference](95-reference.md)** - Complete API reference

## Contributors

```@raw html
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
```
