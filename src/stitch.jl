using DimensionalData
using Makie
using ARPES: stitch_along

export stitch_ui

function stitch_ui(
    A::AbstractDimArray{T,2},
    B::AbstractDimArray{T,2},
    dim::Union{DimensionalData.Dimensions.Dimension,Symbol},
    figure::NamedTuple = (;),
    heatmap_kwargs...,
) where {T}
    default_figure_setting = (size = (650, 450),)
    fig_kwargs = merge(default_figure_setting, figure)
    default_heatmap_setting = (colormap = :turbo,)
    heatmap_setting = merge(default_heatmap_setting, NamedTuple(heatmap_kwargs))

    fig = Figure(; fig_kwargs...)
    ax = Axis(fig[1, 1])
    slider_layout = GridLayout(fig[2, 1], tellwidth = false)
    sg = SliderGrid(
        slider_layout[1, 1],
        (label = "Seam Ratio", range = 0:0.001:1, startvalue = 0.5, width = 200),
        (label = "Gain A(1st Arg.)", range = 0:0.001:1, startvalue = 1, width = 200),
        width = 500,
        tellheight = true,
    )
    seam_ratio_obs = sg.sliders[1].value
    gain_a_obs = sg.sliders[2].value
    C = lift(seam_ratio_obs, gain_a_obs) do r, g
        stitch_along(A, B, dim; seam_ratio = r, gain = g)
    end
    x_obs = lift(s -> lookup(s, 1), C)
    y_obs = lift(s -> lookup(s, 2), C)
    updated_C = lift(parent, C)
    heatmap!(ax, x_obs, y_obs, updated_C, heatmap_setting...)
    fig
end
