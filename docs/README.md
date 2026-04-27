# Documentation sources (contributor notes)

This directory feeds the [Documenter](https://github.com/JuliaDocs/Documenter.jl) site via `make.jl`. The notes here are for people editing sources locally; they are **not** themselves published as a documentation page.

## Building

From the `docs` folder, with the `docs/Project.toml` environment:

```bash
julia --project=. make.jl
```

`make.jl` runs [Literate.jl](https://github.com/fredrikekre/Literate.jl) on every `*.jl` file under `src/`, writing matching `*.md` files, then calls `makedocs`.

Those generated `3*.md` files are listed in the root `.gitignore` (they are recreated in CI by `docs/make.jl`). Edit the `*.jl` sources only.

## Page map

| Doc | Source | Role |
| --- | --- | --- |
| ✓ | `src/index.md` | Site home |
| ✓ | `src/10-energy-dependence.md` | Energy dependence |
| ✓ | `src/11-quantization-reference.md` | Quantization / conventions |
| ✓ | `src/12-performance.md` | Performance |
| ✓ | `src/30-tutorial.jl` | Main tutorial (Literate) → `30-tutorial.md` in the build |
| ✓ | `src/31-canonical.jl` | Canonical vs helicity (Literate) → `31-canonical.md` |
| ✓ | `src/32-visualization.jl` | Visualization (Literate) → `32-visualization.md` |
| ✓ | `src/33-xic-overlaps.jl` | Ξc overlaps (Literate) → `33-xic-overlaps.md` |
| ✓ | `src/34-angular-distributions.jl` | Angular distributions (Literate) → `34-angular-distributions.md` |
| ✓ | `src/90-contributing.md` | Contributing |
| ✓ | `src/91-developer.md` | Developer docs |
| ✓ | `src/95-reference.md` | API reference |

The **Doc** column marks sources that become part of the published manual: `make.jl` collects `index.md` and every other `src/*.md` after Literate runs (so each `*.jl` row is included **via** its generated `*.md`). This `README.md` file is not part of that set.

## Overview notebook (not part of the official manual)

A Jupyter notebook can serve as one continuous tour of the same material as the numbered pages (interactive cells, optional extra plots). **It is not included in the Documenter build:** nothing points to it from `index.md` or `make.jl`, and it will not show up in the published site.

If you keep such a notebook in the repo, a conventional path is:

`docs/notebooks/overview.ipynb`
