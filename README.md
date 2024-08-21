# ThreeBodyDecays

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://mmikhasenko.github.io/ThreeBodyDecays.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://mmikhasenko.github.io/ThreeBodyDecays.jl/dev)
[![Build Status](https://github.com/mmikhasenko/ThreeBodyDecays.jl/workflows/Test/badge.svg)](https://github.com/mmikhasenko/ThreeBodyDecays.jl/actions)
[![Test workflow status](https://github.com/mmikhasenko/ThreeBodyDecays.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/ThreeBodyDecays.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/mmikhasenko/ThreeBodyDecays.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/ThreeBodyDecays.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/mmikhasenko/ThreeBodyDecays.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/ThreeBodyDecays.jl/actions/workflows/Docs.yml?query=branch%3Amain)

[![Coverage](https://codecov.io/gh/mmikhasenko/ThreeBodyDecays.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mmikhasenko/ThreeBodyDecays.jl)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![All Contributors](https://img.shields.io/github/all-contributors/mmikhasenko/ThreeBodyDecays.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors)

## Description

`ThreeBodyDecays.jl` is a Julia package for building hadronic decay models using a cascade reaction.
Decays with three particles is the main application of the approach, however, it is also found useful in multibody decays where transitions can be factorized to a product of sequential decays with $≤3$ products. Particles can have arbitrary spin.

The implementation is based on a research paper, "Dalitz-plot decomposition for three-body decays" by JPAC Collaboration (M. Mikhasenko at al.) [(inspire reference)](http://inspirehep.net/record/1758460).

The code mostly inherits notations of the paper:

- Particles are numbered 1,2,3, and 0 for the decay products and the mother particle, respectively.
- `m0` is a mass of the decay particle, and `m1`, `m2`, `m3` are masses the final-state particles.
- `σ` is a two-particle invariant mass squared, `σk = (pi+pj)²`,
- `θij` is a scattering angle, an angle between `vec pi` and `- vec pk`.
- `ζ_kj_for_0` is the Wigner angle of the 0-particle, an angle of `vec pⱼ+pⱼ` with respect the the chain `j`.
- `ζ_ij_for_k` is the Wigner angle for the particle `k` (the angle in the rest frame of particle `k`) that is mismatched for the chain `i` with respect to the chain `j`.

See [example](docs/src/demo.jl) for a demonstration case.

## Installation

```julia
using Pkg
Pkg.add("ThreeBodyDecays")
```

## How to Cite

If you use ThreeBodyDecays.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/mmikhasenko/ThreeBodyDecays.jl/blob/main/CITATION.cff).
The [PRD101 (2020) 3, 034033 paper](http://inspirehep.net/record/1758460) is available for academic citations.

## Contributing

If you want to make contributions of any kind, please first that a look into our [contributing guide directly on GitHub](docs/src/90-contributing.md) or the [contributing page on the website](https://mmikhasenko.github.io/ThreeBodyDecays.jl/dev/90-contributing/).

---

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
