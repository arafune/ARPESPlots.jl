# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.2] - 2026-XX-XX

### Changed

- Add x- and y-labels to `crosshair_heatmap`

## [0.0.1] - 2026-04-07

### Added

- `crosshair_heatmap`: Interactive 2D heatmap with linked crosshair, side-profile line plots, and integration-thickness sliders. Supports both 2D and 3D `DimArray` input.
- `waterfall_dispersion` / `waterfall_dispersion!`: Stacked offset line/fill plots sliced along a chosen dimension of a `DimArray`.
- Makie `func` / `func!` pattern for both plot types.
- `color` keyword accepted in addition to `colormap` for the `cmap` parameter in `crosshair_heatmap`.
- `isfinite` filtering to handle NaN values in axis limit determination.
- Fix for blank/all-NaN data when computing `ylims!` in `waterfall_dispersion`.
- Documentation (docstrings) for both plot functions.
- Test suite using CairoMakie (headless) backend.
- CI workflow with code coverage via Codecov.
