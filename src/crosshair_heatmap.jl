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

    @debug "size(A)" size(A)
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

function crosshair_heatmap(
    A::AbstractDimArray{T,3} where {T},
    stack_dim::Union{DimensionalData.Dimension,Symbol};
    figure::NamedTuple = (;),
    heatmap_setting::NamedTuple = (;),
    axis_top::NamedTuple = (;),
    axis_right::NamedTuple = (;),
)
    # Default settings for figure and axes
    default_figure_setting = (size = (900, 800),)
    default_top_axis_setting = (xticklabelsvisible = false, ylabel = "Intensity")
    default_right_axis_setting = (yticklabelsvisible = false, xlabel = "Intensity")
    default_heatmap_setting = (colormap = :turbo,)

    fig_kwargs = merge(default_figure_setting, figure)
    fig = Figure(; fig_kwargs...)

    # Resolve stack dimension (Symbol → Dimension)
    stack_dim = stack_dim isa Symbol ? dims(A, stack_dim) : stack_dim
    stack_lookup = lookup(A, stack_dim)
    n_stack = length(stack_lookup)
    stack_dim_idx = dimnum(A, stack_dim)
    # Observable for current stack index
    stack_idx = Observable(1)

    # Observable 2D slice extracted from 3D array
    A2 = lift(stack_idx) do i
        selectdim(A, stack_dim_idx, i)
    end

    xvals = lift(A2) do a
        lookup(a, 1)
    end

    yvals = lift(A2) do a
        lookup(a, 2)
    end

    # Get size of a representative 2D slice
    nx = size(selectdim(A, stack_dim_idx, 1), 1)
    ny = size(selectdim(A, stack_dim_idx, 1), 2)

    # Create axes
    ax_main = Axis(fig[2, 1])
    ax_top = Axis(fig[1, 1]; merge(default_top_axis_setting, axis_top)...)
    ax_right = Axis(fig[2, 2]; merge(default_right_axis_setting, axis_right)...)

    # Slider UI (horizontal radius, vertical radius, stack index)
    slider_layout = GridLayout(fig[3, 1:2], tellwidth = false)
    sg = SliderGrid(
        slider_layout[1, 1],
        (label = "Horizontal radius(pts)", range = 0:(nx÷2), startvalue = 0),
        (label = "Vertical radius(pts)", range = 0:(ny÷2), startvalue = 0),
        (label = "Stack index", range = 1:n_stack, startvalue = 1),
        width = 500,
        tellheight = true,
    )

    # Observables for crosshair thickness
    crosshair_thick_x = sg.sliders[1].value
    crosshair_thick_y = sg.sliders[2].value

    # Update stack index from slider
    on(sg.sliders[3].value) do v
        stack_idx[] = Int(round(v))
    end

    # Layout adjustments
    rowsize!(fig.layout, 1, Relative(0.25))
    colsize!(fig.layout, 2, Relative(0.25))
    rowgap!(fig.layout, 10)
    colgap!(fig.layout, 10)

    # Plot heatmap using observable 2D slice
    heatmap!(
        ax_main,
        xvals,
        yvals,
        lift(a -> parent(a), A2);
        merge(default_heatmap_setting, heatmap_setting)...,
    )

    # Mouse position observable (initialized at center of slice)
    pos = Observable(Point2f(median(xvals[]), median(yvals[])))

    on(A2) do a
        pos[] = Point2f(median(lookup(a, 1)), median(lookup(a, 2)))
    end

    # Update mouse position when cursor moves inside axis
    on(events(ax_main.scene).mouseposition) do _
        if is_mouseinside(ax_main.scene)
            pos[] = mouseposition(ax_main.scene)
        end
    end

    # Crosshair line styles
    v_line_width = lift(w -> w == 0 ? 1.5 : 0.0, crosshair_thick_x)
    h_line_width = lift(w -> w == 0 ? 1.5 : 0.0, crosshair_thick_y)
    v_span_color = lift(w -> w == 0 ? (:red, 0.0) : (:red, 0.2), crosshair_thick_x)
    h_span_color = lift(w -> w == 0 ? (:red, 0.0) : (:red, 0.2), crosshair_thick_y)

    # Draw vertical and horizontal crosshair lines
    vlines!(
        ax_main,
        lift(p -> p[1], pos),
        color = (:red, 0.3),
        linewidth = v_line_width,
        linestyle = :dash,
    )

    hlines!(
        ax_main,
        lift(p -> p[2], pos),
        color = (:red, 0.3),
        linewidth = h_line_width,
        linestyle = :dash,
    )

    # Helper: get index range around a value with given thickness
    function get_range_indices(val, lookup_vals, n, thickness)
        idx = clamp(searchsortedfirst(lookup_vals, val), 1, n)
        return max(1, idx-thickness):min(n, idx+thickness)
    end

    # Vertical integration span (x-direction)
    vspan_range = lift(A2, pos, crosshair_thick_x) do a, p, tx
        lx = lookup(a, 1)
        r = get_range_indices(p[1], lx, length(lx), tx)
        return lx[first(r)], lx[last(r)]
    end

    vspan!(
        ax_main,
        lift(r -> r[1], vspan_range),
        lift(r -> r[2], vspan_range),
        color = v_span_color,
    )

    # Horizontal integration span (y-direction)
    hspan_range = lift(A2, pos, crosshair_thick_y) do a, p, ty
        ly = lookup(a, 2)
        r = get_range_indices(p[2], ly, length(ly), ty)
        return ly[first(r)], ly[last(r)]
    end

    hspan!(
        ax_main,
        lift(r -> r[1], hspan_range),
        lift(r -> r[2], hspan_range),
        color = h_span_color,
    )

    # Top plot: average along y-direction (horizontal slice)
    line_top_data = lift(A2, pos, crosshair_thick_y) do a, p, ty
        ly = lookup(a, 2)
        ry = get_range_indices(p[2], ly, length(ly), ty)
        sliced = vec(mean(parent(a[:, ry]), dims = 2))
        #    @debug "length(sliced), length(lookup(a, 1))" length(sliced), length(lookup(a, 1))
        zmin, zmax = minimum(a), maximum(a)
        ylims!(ax_top, zmin, zmax)
        return Point2f.(lookup(a, 1), sliced)
    end
    lines!(ax_top, line_top_data, color = :blue)

    # Right plot: average along x-direction (vertical slice)
    line_right_data = lift(A2, pos, crosshair_thick_x) do a, p, tx
        lx = lookup(a, 1)
        rx = get_range_indices(p[1], lx, length(lx), tx)
        sliced = vec(mean(parent(a[rx, :]), dims = 1))
        return Point2f.(sliced, lookup(a, 2))
    end
    lines!(ax_right, line_right_data, color = :blue)

    # Label showing current position and averaged intensity
    label_text =
        lift(A2, pos, crosshair_thick_x, crosshair_thick_y, stack_idx) do a, p, tx, ty, si
            lx = lookup(a, 1)
            ly = lookup(a, 2)
            nx, ny = length(lx), length(ly)

            ix = clamp(searchsortedfirst(lx, p[1]), 1, nx)
            iy = clamp(searchsortedfirst(ly, p[2]), 1, ny)

            rx = max(1, ix-tx):min(nx, ix+tx)
            ry = max(1, iy-ty):min(ny, iy+ty)

            z = mean(parent(a[rx, ry]))

            return """
                $(name(dims(a, 1))): $(round(lx[ix], digits=4)) (±$tx pts)
                $(name(dims(a, 2))): $(round(ly[iy], digits=4)) (±$ty pts)
                $(name(stack_dim)): $(round(stack_lookup[si], digits=4))
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

    # Link axes between main and projections
    linkxaxes!(ax_main, ax_top)
    linkyaxes!(ax_main, ax_right)

    # Update projection axis limits based on current slice
    on(A2) do a
        zmin, zmax = minimum(a), maximum(a)
        ylims!(ax_top, zmin, zmax)
        xlims!(ax_right, zmin, zmax)
    end

    return fig
end


# =============================
# Internal helper functions (prefixed with `_`)
# =============================

# Compute the index range around a value with a given thickness
# Returns a UnitRange of indices
function _get_range_indices(
    val::Real,
    lookup_vals::AbstractVector{<:Real},
    n::Int,
    thickness::Int,
)::UnitRange{Int}
    idx = clamp(searchsortedfirst(lookup_vals, val), 1, n)
    return max(1, idx-thickness):min(n, idx+thickness)
end

# Get the closest indices in both x and y directions from a mouse position
function _get_indices(p::Point2f, A::AbstractDimArray{T,2}) where {T}
    nx, ny = size(A)
    x_idx = clamp(searchsortedfirst(lookup(A, 1), p[1]), 1, nx)
    y_idx = clamp(searchsortedfirst(lookup(A, 2), p[2]), 1, ny)
    return x_idx, y_idx
end

# Compute the horizontal slice (top line) as a vector of Point2f
# Averaging over a vertical range defined by thickness
function _compute_slice_line_top(
    A::AbstractDimArray{T,2},
    pos::Point2f,
    thickness::Int,
)::Vector{Point2f} where {T}
    ry = _get_range_indices(pos[2], lookup(A, 2), size(A, 2), thickness)
    sliced_data = vec(mean(parent(A[:, ry]), dims = 2))
    return Point2f.(lookup(A, 1), sliced_data)
end

# Compute the vertical slice (right line) as a vector of Point2f
# Averaging over a horizontal range defined by thickness
function _compute_slice_line_right(
    A::AbstractDimArray{T,2},
    pos::Point2f,
    thickness::Int,
)::Vector{Point2f} where {T}
    rx = _get_range_indices(pos[1], lookup(A, 1), size(A, 1), thickness)
    sliced_data = vec(mean(parent(A[rx, :]), dims = 1))
    return Point2f.(sliced_data, lookup(A, 2))
end

# Compute the vertical span range (x-direction) for the crosshair
# Returns a tuple (xmin, xmax) for the shaded region
function _compute_vspan_range(
    A::AbstractDimArray{T,2},
    pos::Point2f,
    thickness::Int,
)::Tuple{Float64,Float64} where {T}
    lx = lookup(A, 1)
    r = _get_range_indices(pos[1], lx, length(lx), thickness)
    return lx[first(r)], lx[last(r)]
end

# Compute the horizontal span range (y-direction) for the crosshair
# Returns a tuple (ymin, ymax) for the shaded region
function _compute_hspan_range(
    A::AbstractDimArray{T,2},
    pos::Point2f,
    thickness::Int,
)::Tuple{Float64,Float64} where {T}
    ly = lookup(A, 2)
    r = _get_range_indices(pos[2], ly, length(ly), thickness)
    return ly[first(r)], ly[last(r)]
end

# Generate the label string showing coordinates and averaged intensity
# Optionally includes stack dimension info for 3D data
function _make_label(
    A::AbstractDimArray{T,2},
    pos::Point2f,
    crosshair_thick_x::Int,
    crosshair_thick_y::Int;
    stack_val::Union{Nothing,Real} = nothing,
    stack_dim_name::Union{Nothing,String} = nothing,
    stack_lookup::Union{Nothing,AbstractVector{<:Real}} = nothing,
)::String where {T}

    nx, ny = size(A)
    ix, iy = _get_indices(pos, A)

    # Compute the integration range around the cursor
    rx = max(1, ix-crosshair_thick_x):min(nx, ix+crosshair_thick_x)
    ry = max(1, iy-crosshair_thick_y):min(ny, iy+crosshair_thick_y)
    z = mean(parent(A[rx, ry]))

    # Build the label string
    s = """
        $(name(dims(A, 1))): $(round(lookup(A, 1)[ix], digits=4)) (±$crosshair_thick_x pts)
        $(name(dims(A, 2))): $(round(lookup(A, 2)[iy], digits=4)) (±$crosshair_thick_y pts)
        Avg Intensity: $(round(z, digits=4))
        """
    # Add stack dimension info if provided (for 3D slices)

    if stack_val === nothing && stack_lookup !== nothing && stack_dim_name !== nothing
        stack_val = stack_lookup[stack_val]
    end
    return s
end
