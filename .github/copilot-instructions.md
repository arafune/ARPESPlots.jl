# Copilot Instructions for ARPESPlots.jl

## Build, Test, and Lint

```julia
# Run full test suite
julia --project=. -e "using Pkg; Pkg.test()"

# Run a single test file
julia --project=. test/crosshair_heatmap.jl

# Build documentation
julia --project=docs docs/make.jl

# Format code (JuliaFormatter with "blue" style)
julia -e 'using JuliaFormatter; format(".")'
```

> **Note:** The `ARPES` package is an unregistered dependency. In fresh environments, `Pkg.develop` it first (see CI.yml for the pattern).

## Architecture

This is a Makie-based visualization library for ARPES (Angle-Resolved Photoemission Spectroscopy) data. It is a companion to [ARPES.jl](https://github.com/arafune/ARPES.jl) and deliberately kept as a separate package (polyrepo approach).

**Two plot types are implemented:**

- `waterfall_dispersion` / `waterfall_dispersion!` (`src/waterfall.jl`): Stacked offset line/fill plots sliced along a chosen dimension of a `DimArray`. Follows the Makie `func` / `func!` convention — the `!` variant takes an existing `Axis`, the plain variant creates a `Figure`.
- `crosshair_heatmap` (`src/crosshair_heatmap.jl`): Interactive 2D heatmap with linked crosshair, side-profile line plots, and integration-thickness sliders. Mouse position is tracked via `events(scene).mouseposition`; reactive updates use Makie `Observable`s and `lift`.

All public functions accept `AbstractDimArray` (from DimensionalData.jl) as their primary data argument.

## Key Conventions

- **Makie `func` / `func!` pattern:** Mutating variants (`!`) accept an existing `AbstractAxis` and return `(plots, axis_right)`; non-mutating variants create and return a `Makie.FigureAxisPlot` (waterfall) or `Makie.Figure` (heatmap).
- **Keyword argument forwarding via NamedTuples:** Plot customisation is passed as `figure=(;)`, `axis=(;)`, `axis_right=(;)`, `heatmap_setting=(;)`, etc., and merged with defaults using `merge`.
- **DimensionalData dimensions drive slicing:** Use `eachslice` with a named dimension (not integer indices) to iterate over stacking slices.
- **Code style:** JuliaFormatter "blue" style (configured in `.JuliaFormatter.toml`).
- **Julia ≥ 1.12** is required (see `[compat]` in `Project.toml`).
- **Tests use CairoMakie** (headless backend) for rendering assertions; interactive backends are not used in CI.
