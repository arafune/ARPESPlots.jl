using DimensionalData
using Makie
using ARPES: stitch_along

export stitch_ui

"""
    stitch_ui(A, B, dim, figure=(;); heatmap_kwargs...)
    stitch_ui(A, B; dim=:phi, figure=(;), kwargs...)
    stitch_ui(A::ARPESData, B::ARPESData, figure; heatmap_kwargs...)

Interactively stitch two 2-D `AbstractDimArray`s along `dim` and display the result
as a heatmap with live-updating sliders.

# Arguments
- `A`, `B`: Two 2-D `AbstractDimArray`s (or `ARPESData`) to stitch together.
- `dim`: Dimension along which to stitch. Accepts a `DimensionalData.Dimension` or a
  `Symbol` (e.g. `:phi`). Defaults to `:phi` in the keyword-argument overloads.
- `figure`: `NamedTuple` of keyword arguments forwarded to `Makie.Figure`.
  Defaults to `(size = (650, 450))`.

# Keyword Arguments
- `heatmap_kwargs...`: Additional keyword arguments forwarded to `Makie.heatmap!`.
  `colormap` defaults to `:turbo`.

# Sliders
- **Seam Ratio** (`0`–`1`, default `0.5`): Controls where the seam between `A` and `B`
  is placed. Forwarded to `ARPES.stitch_along` as `seam_ratio`.
- **Gain A** (`0`–`3`, default `1`): Intensity gain applied to `A` before stitching.
  Forwarded to `ARPES.stitch_along` as `gain_a`.

# Returns
A `Makie.Figure` containing the heatmap and the slider grid.

# Examples
```julia
fig = stitch_ui(data_a, data_b, :phi)
display(fig)
```
"""
function stitch_ui(
    A::AbstractDimArray{T,2},
    B::AbstractDimArray{T,2},
    dim::Union{DimensionalData.Dimensions.Dimension,Symbol},
    figure::NamedTuple = (;);
    heatmap_kwargs...,
) where {T}
    default_figure_setting = (size = (650, 450),)
    fig_kwargs = merge(default_figure_setting, figure)
    default_heatmap_setting = (colormap = :turbo,)
    heatmap_setting = merge(default_heatmap_setting, NamedTuple(heatmap_kwargs))

    fig = Figure(; fig_kwargs...)
    ax = Axis(
        fig[1, 1],
        xlabel = DimensionalData.Dimensions.label(dims(A, 1)),
        ylabel = DimensionalData.Dimensions.label(dims(A, 2)),
    )
    slider_layout = GridLayout(fig[2, 1], tellwidth = false)
    sg = SliderGrid(
        slider_layout[1, 1],
        (label = "Seam Ratio", range = 0:0.001:1, startvalue = 0.5, width = 200),
        (label = "Gain A(1st Arg.)", range = 0:0.002:3, startvalue = 1, width = 200),
        width = 500,
        tellheight = true,
    )
    seam_ratio_obs = sg.sliders[1].value
    gain_a_obs = sg.sliders[2].value
    C = lift(seam_ratio_obs, gain_a_obs) do r, g
        stitch_along(A, B, dim; seam_ratio = r, gain_a = g)
    end
    x_obs = lift(s -> lookup(s, 1), C)
    y_obs = lift(s -> lookup(s, 2), C)
    updated_C = lift(parent, C)
    heatmap!(ax, x_obs, y_obs, updated_C; heatmap_setting...)
    fig
end

stitch_ui(
    A::AbstractDimArray{T,2},
    B::AbstractDimArray{T,2};
    dim::Union{DimensionalData.Dimensions.Dimension,Symbol} = :phi,
    figure::NamedTuple = (;),
    kwargs...,
) where {T} = stitch_ui(A, B, dim, figure; kwargs...)

function stitch_ui(
    A::ARPESData{T,2},
    B::ARPESData{T,2},
    figure::NamedTuple;
    heatmap_kwargs...,
) where {T}
    stitch_ui(A, B, :phi, figure; heatmap_kwargs...)
end
