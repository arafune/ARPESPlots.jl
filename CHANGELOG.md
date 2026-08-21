# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.5] - Unreleased

### Added

- `tarpes_evolution_mp4(A, outpath; ...)`: new helper to create animations by iterating the `stack_dim` (default `:delay`) frames and writing a video file. When `transparent=false` it uses `Makie.record` for in-process encoding; when `transparent=true` it writes PNG frames (preserving alpha) and invokes `ffmpeg` to encode to an alpha-capable container.
- `transparent` option and `encoder` choice (`:webm` (VP9) or `:mov` (ProRes 4444)) to preserve transparency. Note: MP4/H.264 does not support alpha — use `.webm` or `.mov` for transparent output.
- Export `tarpes_evolution_mp4` from `src/tarpes.jl`.
- Documentation: mention `ffmpeg` is required on PATH when `transparent=true`.

### Changed

- Add transparent keyward argument to tarpes_evolution_heatmaps to switch the background transparency.

## [0.0.4] - 2026-08-14

### Added

- Add `tarpes.jl`: `tarpes_evolution_heatmaps` for time-resolved (tr-)ARPES snapshots and temporal-evolution heatmaps.
- Tests: add `test/tarpes.jl` covering the new API and headless rendering of the heatmap figure.
- Documentation: add `tarps.md` (documentation and docstrings updated for the new functions).

## [0.0.3] - 2026-04-18

### Fixed

- Fix bug in `crosshair_heatmap` where the crosshair and line plots would not update correctly when the heatmap data contained NaN values. Now properly handles NaNs by ignoring them in the crosshair position and line plot updates.

## [0.0.2] - 2026-04-13

### Added

- `stitch_ui``: Graphical UI for`stitch_along`function. Only for 2D`AbstractDimArrays`.

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
