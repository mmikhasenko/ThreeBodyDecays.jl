# ThreeBodyDecays.jl

![Build Status](https://github.com/mmikhasenko/ThreeBodyDecays.jl/actions/workflows/ci.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/mmikhasenko/ThreeBodyDecays.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mmikhasenko/ThreeBodyDecays.jl)
[![arXiv article](https://img.shields.io/badge/article-PRD%20101%2C%20034033-yellowgreen)](https://arxiv.org/abs/1910.04566)

## Description

ThreeBodyDecays.jl is a Julia package for building three-body decay models using a cascade reaction.
It is designed for researchers and scientists in the field of particle physics.

The main focus of the project is the three-particle amplitudes.
Decays with three particles is the main application of the approach, however, it is also found useful in multibody decays where
transitions can be factorized to a product of sequential decays with $≤3$ products. Particles can have arbitrary spin.

The implementation is based on a research paper,
"Dalitz-plot decomposition for three-body decays" by JPAC Collaboration (M Mikhasenko at al.) [(arxiv)](http://inspirehep.net/record/1758460).
The code mostly inherits notations of the paper:

- Particles are numbered 1,2,3, and 0 for the decay products and the mother particle, respectively.
- `m0` is a mass of the decay particle, and `m1`, `m2`, `m3` are masses the final-state particles.
- `σ` is a two-particle invariant mass squared, `σk = (pi+pj)²`,
- `θij` is a scattering angle, an angle between `vec pi` and `- vec pk`.
- `ζ⁰ₖ₍ⱼ₎` is the Wigner angle of the 0-particle, an angle of `vec pⱼ+pⱼ` with respect the the chain `j`.
- `ζᵏᵢ₍ⱼ₎` is the Wigner angle for the particle `k` (the angle in the rest frame of particle `k`) that is mismatched for the chain `i` with respect to the chain `j`.

See [example](docs/src/demo.jl) for a demomostration case.

## Installation
```julia
using Pkg
Pkg.add("ThreeBodyDecays")
```

## Contributing
We welcome contributions! Please submit an issue to start a discussion.

## License
This project is licensed under the [MIT License](LICENSE.md).

## Contact
For more information, please contact [mikhail.mikhasenko@cern.ch](mailto:mikhail.mikhasenko@cern.ch).
