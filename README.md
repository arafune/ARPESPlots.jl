# ARPESPlots.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://arafune.github.io/ARPESPlots.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://arafune.github.io/ARPESPlots.jl/dev/)
[![Build Status](https://github.com/arafune/ARPESPlots.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/arafune/ARPESPlots.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/arafune/ARPESPlots.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/arafune/ARPESPlots.jl)

A supporting library for [ARPES.jl](https://arafune.github.io/ARPES.jl).
Currently, there are no plans to include ARPESPlots.jl in ARPES.jl itself. (ARPES.jl is being developed with a polyrepo approach.)

## Overview

**ARPESPlots.jl** is a library for plotting ARPES (Angle-Resolved Photoemission Spectroscopy) data, primarily built on top of Makie.
Plotting ARPES data often requires tuning plot styles and arguments depending on the situation.
This package provides reusable functions for commonly used, but not-so-trivial, plots—functions
that can otherwise be tedious to write from scratch.

The package is designed to allow as much fine-tuning of the plots as possible,
while still streamlining routine visualization tasks in ARPES data analysis.

## Features

- Convenient functions for common—but complex—ARPES plot types
- Highly flexible design for customization and fine adjustments
- Built using [Makie.jl](https://makie.juliaplots.org/)

## Installation

```julia
import Pkg
Pkg.add(url="https://github.com/arafune/ARPESPlots.jl")
```

## Usage

See examples and API documentation (coming soon) for usage instructions.

---

This library is intended to complement ARPES.jl and streamline common visualization tasks in ARPES research.
