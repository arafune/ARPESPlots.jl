using DimensionalData
using Statistics
using Makie
using Makie: AbstractAxis, AbstractPlot, Axis, Figure, lines!, colormap, ylims!
using ARPES

export crosshair_heatmap

"""
    crosshair_heatmap(A::AbstractDimArray{T,2} where T; kwargs...) -> Figure

Create an interactive 2D heatmap for ARPES data with linked crosshair cursors and
side-panel line profiles.

# Arguments
- `A::AbstractDimArray{T,2}`: The 2D dimensional array to plot. Must have valid
`DimensionalData` dimensions (e.g., `eV`, `phi`).

# Keyword Arguments
- `figure::NamedTuple = (;)`: Arguments passed to `Makie.Figure` (e.g., `size = (800, 600)`).
- `axis_top::NamedTuple = (;)`: Custom attributes for the top slice axis (X-profile).
- `axis_right::NamedTuple = (;)`: Custom attributes for the right slice axis (Y-profile).

# Features
- **Dynamic Slicing**: Moving the mouse over the main heatmap updates the top and right
  line plots in real-time.
- **Linked Axes**: Zooming or panning the main heatmap automatically stays in sync with the
  slice axes via `linkxaxes!` and `linkyaxes!`.
- **Live Labels**: Displays the current coordinate and intensity value in the top-right
  corner.

# Example

```julia
using DimensionalData, ARPESPlots, GLMakie
data = rand(100, 100)
spec = ARPESData(data, (phi(-10:10), eV(0:0.1:9.9)))
fig = crosshair_heatmap(spec)
display(fig)
```

"""
function crosshair_heatmap(
    A::AbstractDimArray{T,2} where {T};
    figure::NamedTuple = (;),
    heatmap_setting::NamedTuple = (;),
    axis_top::NamedTuple = (;),
    axis_right::NamedTuple = (;),
)
    # Get the dimensions and their corresponding axes
    default_figure_setting = (size = (900, 700),)
    default_top_axis_setting = (xticklabelsvisible = false, ylabel = "Intensity")
    default_right_axis_setting = (yticklabelsvisible = false, xlabel = "Intensity")
    default_heatmap_setting = (colormap = :turbo,)

    fig_kwargs = merge(default_figure_setting, figure)
    fig = Figure(; fig_kwargs...)
    nx, ny = size(A)
    # Create a figure and axis
    ax_main = Axis(fig[2, 1])
    ax_top = Axis(fig[1, 1]; merge(default_top_axis_setting, axis_top)...)
    ax_right = Axis(fig[2, 2]; merge(default_right_axis_setting, axis_right)...)
    #
    slider_layout = GridLayout(fig[3, 1:2], tellwidth = false)
    sg = SliderGrid(
        slider_layout[1, 1],
        (label = "Horizontal radius(pts)", range = 0:(nx÷2), startvalue = 0),
        (label = "Vertical radius(pts)", range = 0:(ny÷2), startvalue = 0),
        width = 500,
        tellheight = true,
    )
    crosshair_thick_x = sg.sliders[1].value
    crosshair_thick_y = sg.sliders[2].value
    #
    rowsize!(fig.layout, 1, Relative(0.25))
    colsize!(fig.layout, 2, Relative(0.25))
    rowgap!(fig.layout, 10)
    colgap!(fig.layout, 10)

    # Plot the heatmap
    heatmap!(ax_main, A; merge(default_heatmap_setting, heatmap_setting)...)
    #get the mouse event 
    pos = Observable(Point2f(median(lookup(A, 1)), median(lookup(A, 2))))
    on(events(ax_main.scene).mouseposition) do _
        if is_mouseinside(ax_main.scene)
            pos[] = mouseposition(ax_main.scene)
        end
    end

    v_line_width = lift(w -> w == 0 ? 1.5 : 0.0, crosshair_thick_x)
    h_line_width = lift(w -> w == 0 ? 1.5 : 0.0, crosshair_thick_y)
    v_line_style = :dash
    h_line_style = :dash
    v_span_color = lift(w -> w == 0 ? (:red, 0.0) : (:red, 0.2), crosshair_thick_x)
    h_span_color = lift(w -> w == 0 ? (:red, 0.0) : (:red, 0.2), crosshair_thick_y)

    vlines!(
        ax_main,
        lift(p -> p[1], pos),
        color = (:red, 0.3),
        linewidth = v_line_width,
        linestyle = v_line_style,
    )

    function get_range_indices(val, lookup_vals, n, thickness)
        idx = clamp(searchsortedfirst(lookup_vals, val), 1, n)
        return max(1, idx-thickness):min(n, idx+thickness)
    end

    vspan_range = lift(pos, crosshair_thick_x) do p, tx
        # get index
        r = get_range_indices(p[1], lookup(A, 1), nx, tx)
        return lookup(A, 1)[first(r)], lookup(A, 1)[last(r)]
    end

    vspan!(
        ax_main,
        lift(r -> r[1], vspan_range),
        lift(r -> r[2], vspan_range),
        color = v_span_color,
    )

    hspan_range = lift(pos, crosshair_thick_y) do p, ty
        r = get_range_indices(p[2], lookup(A, 2), ny, ty)
        return lookup(A, 2)[first(r)], lookup(A, 2)[last(r)]
    end

    hlines!(
        ax_main,
        lift(p -> p[2], pos),
        color = (:red, 0.3),
        linewidth = h_line_width,
        linestyle = h_line_style,
    )

    hspan!(
        ax_main,
        lift(r -> r[1], hspan_range),
        lift(r -> r[2], hspan_range),
        color = h_span_color,
    )

    function get_range(val, lookup_vals, n, thickness)
        idx = clamp(searchsortedfirst(lookup_vals, val), 1, n)
        return max(1, idx-thickness):min(n, idx+thickness)
    end


    # get index of the closest point in the heatmap to the mouse position
    function get_indices(p)
        x, y = p
        x_idx = clamp(searchsortedfirst(lookup(A, 1), x), 1, nx)
        y_idx = clamp(searchsortedfirst(lookup(A, 2), y), 1, ny)
        return x_idx, y_idx
    end

    line_top_data = lift(pos, crosshair_thick_y) do p, ty
        ry = get_range_indices(p[2], lookup(A, 2), ny, ty)
        sliced_data = vec(mean(parent(A[:, ry]), dims = 2))
        return Point2f.(lookup(A, 1), sliced_data)
    end
    lines!(ax_top, line_top_data, color = :blue)

    line_right_data = lift(pos, crosshair_thick_x) do p, tx
        rx = get_range_indices(p[1], lookup(A, 1), nx, tx)
        sliced_data = vec(mean(parent(A[rx, :]), dims = 1))
        return Point2f.(sliced_data, lookup(A, 2))
    end
    lines!(ax_right, line_right_data, color = :blue)

    label_text = lift(pos, crosshair_thick_x, crosshair_thick_y) do p, tx, ty
        ix = clamp(searchsortedfirst(lookup(A, 1), p[1]), 1, nx)
        iy = clamp(searchsortedfirst(lookup(A, 2), p[2]), 1, ny)
        # Calculate mean in the current integration box
        rx = max(1, ix-tx):min(nx, ix+tx)
        ry = max(1, iy-ty):min(ny, iy+ty)
        z = mean(parent(A[rx, ry]))
        return """
            $(name(dims(A, 1))): $(round(lookup(A, 1)[ix], digits=4)) (±$tx pts)
            $(name(dims(A, 2))): $(round(lookup(A, 2)[iy], digits=4)) (±$ty pts)
            Avg Intensity: $(round(z, digits=4))
            """
    end

    Label(
        fig[1, 2],
        label_text,
        tellwidth = false,
        halign = :left,
        font = :bold,
        fontsize = 16,
        padding = (10, 10, 10, 10),
        justification = :left,
    )

    linkxaxes!(ax_main, ax_top)
    linkyaxes!(ax_main, ax_right)
    zmin, zmax = minimum(A), maximum(A)
    ylims!(ax_top, zmin, zmax)
    xlims!(ax_right, zmin, zmax)

    return fig
end

