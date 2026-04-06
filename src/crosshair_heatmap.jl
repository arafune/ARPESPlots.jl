using DimensionalData
using Statistics
using Makie
using Makie: AbstractAxis, AbstractPlot, Axis, Figure, lines!, colormap, ylims!
using ARPES

export crosshair_heatmap

"""
    crosshair_heatmap(A::AbstractDimArray{T,2} where T; kwargs...) -> Figure
    crosshair_heatmap(A::AbstractDimArray{T,3} where T,
                      stack_dim::Union{DimensionalData.Dimension, Symbol};
                      kwargs...) -> Figure

Create an interactive heatmap with a linked crosshair cursor and side-panel line
profiles for exploring ARPES data.

The 2D method plots `A` directly.  The 3D method treats `A` as a stack of 2D
heatmaps along `stack_dim` and adds a slider to step through the slices.

# Arguments
- `A`: 2D or 3D data array with `DimensionalData` dimensions.
- `stack_dim` *(3D only)*: Dimension along which to slice (e.g., `Dim{:hv}` or `:hv`).

# Keyword Arguments
- `figure::NamedTuple = (;)`: Passed to `Makie.Figure` (e.g., `size = (900, 600)`).
- `axis_top::NamedTuple = (;)`: Attributes for the top (X-profile) axis.
- `axis_right::NamedTuple = (;)`: Attributes for the right (Y-profile) axis.
- `heatmap_kwargs...`: Additional keyword arguments forwarded to `heatmap!`
  (e.g., `colormap = :inferno`, `colorrange = (0, 1)`).

# Layout
```
┌─────────────┬──────────┐
│  ax_top     │  label   │  row 1 (25 %)
├─────────────┼──────────┤
│  ax_main    │ ax_right │  row 2
├─────────────┴──────────┤
│  SliderGrid            │  row 3
└────────────────────────┘
```
- **ax_main**: central heatmap; crosshair updates as the mouse moves.
- **ax_top**: X-direction profile, averaged over the vertical integration band.
- **ax_right**: Y-direction profile, averaged over the horizontal integration band.
- **SliderGrid**: sliders control the half-width (in index points) of the integration
  bands.  The 3D method adds a third **Stack index** slider.

# Features
- Mouse movement over `ax_main` updates `ax_top` and `ax_right` in real-time.
- `ax_main` is linked to `ax_top` (x-axis) and `ax_right` (y-axis).
- A label in the top-right corner shows the current coordinates, integration range,
  and averaged intensity.  The 3D method also shows the current `stack_dim` value.
- When a slider value is non-zero, a shaded span is drawn on `ax_main` to visualise
  the integration region.

# Examples

**2D**
```julia
using DimensionalData, ARPESPlots, GLMakie

A = DimArray(rand(100, 80), (Dim{:phi}(range(-10, 10, 100)), Dim{:eV}(range(0, 8, 80))))
fig = crosshair_heatmap(A; colormap = :inferno)
display(fig)
```

**3D**
```julia
using DimensionalData, ARPESPlots, GLMakie

A = DimArray(
    rand(100, 80, 5),
    (Dim{:phi}(range(-10, 10, 100)), Dim{:eV}(range(0, 8, 80)), Dim{:hv}([100, 105, 110, 115, 120])),
)
fig = crosshair_heatmap(A, :hv)
display(fig)
```
"""
function crosshair_heatmap(
    A::AbstractDimArray{T,2} where {T};
    figure::NamedTuple = (;),
    axis_top::NamedTuple = (;),
    axis_right::NamedTuple = (;),
    heatmap_kwargs...,
)
    # Get the dimensions and their corresponding axes
    default_figure_setting = (size = (650, 450),)
    default_top_axis_setting = (xticklabelsvisible = false, ylabel = "Intensity")
    default_right_axis_setting = (yticklabelsvisible = false, xlabel = "Intensity")
    default_heatmap_setting = (colormap = :turbo,)
    heatmap_setting = merge(default_heatmap_setting, heatmap_kwargs)

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
    heatmap!(ax_main, A; heatmap_setting...)
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

    vspan_range = lift(pos, crosshair_thick_x) do p, tx
        _compute_vspan_range(A, p, tx)
    end

    vspan!(
        ax_main,
        lift(r -> r[1], vspan_range),
        lift(r -> r[2], vspan_range),
        color = v_span_color,
    )

    hspan_range = lift(pos, crosshair_thick_y) do p, ty
        _compute_hspan_range(A, p, ty)
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

    line_top_data = lift(pos, crosshair_thick_y) do p, ty
        _compute_slice_line_top(A, p, ty)
    end
    lines!(ax_top, line_top_data, color = :blue)

    line_right_data = lift(pos, crosshair_thick_x) do p, tx
        _compute_slice_line_right(A, p, tx)
    end
    lines!(ax_right, line_right_data, color = :blue)

    label_text = lift(pos, crosshair_thick_x, crosshair_thick_y) do p, tx, ty
        _make_label(A, p, tx, ty)
    end

    Label(
        fig[1, 2],
        label_text,
        tellwidth = false,
        halign = :left,
        font = :bold,
        fontsize = 16,
        padding = (5, 10, 10, 10),
        justification = :left,
    )

    linkxaxes!(ax_main, ax_top)
    linkyaxes!(ax_main, ax_right)
    zmin, zmax = extrema(parent(A)[isfinite.(parent(A))])
    ylims!(ax_top, zmin, zmax)
    xlims!(ax_right, zmin, zmax)

    return fig
end

"""
    crosshair_heatmap(A::AbstractDimArray{T,3} where T,
                      stack_dim::Union{DimensionalData.Dimension, Symbol};
                      kwargs...) -> Figure

Create an interactive 2D heatmap with a crosshair cursor, side-panel line profiles,
and a slider to step through slices of a 3D data array.

The 3D array is treated as a stack of 2D heatmaps along `stack_dim`.  All
crosshair and profile features from the 2D method are available; the extra
slider selects which 2D slice is displayed.

# Arguments
- `A::AbstractDimArray{T,3}`: 3D data array with `DimensionalData` dimensions.
- `stack_dim`: The dimension along which to slice (e.g., `Dim{:hv}` or `:hv`).

# Keyword Arguments
- `figure::NamedTuple = (;)`: Passed to `Makie.Figure`.
- `axis_top::NamedTuple = (;)`: Attributes for the top (X-profile) axis.
- `axis_right::NamedTuple = (;)`: Attributes for the right (Y-profile) axis.
- `heatmap_kwargs...`: Additional keyword arguments forwarded to `heatmap!`.

# Layout
Same 2×2 grid as the 2D method, plus a **SliderGrid** with three sliders:
1. **Horizontal radius (pts)** — half-width of the x-direction integration band.
2. **Vertical radius (pts)** — half-width of the y-direction integration band.
3. **Stack index** — selects the active 2D slice along `stack_dim`.

The label in the top-right corner also shows the current value of `stack_dim`.

# Example

```julia
using DimensionalData, ARPESPlots, GLMakie

A = DimArray(
    rand(100, 80, 5),
    (Dim{:phi}(range(-10, 10, 100)), Dim{:eV}(range(0, 8, 80)), Dim{:hv}([100, 105, 110, 115, 120])),
)
fig = crosshair_heatmap(A, :hv)
display(fig)
```

See also: [`crosshair_heatmap`](@ref) (2D method).
"""
function crosshair_heatmap(
    A::AbstractDimArray{T,3} where {T},
    stack_dim::Union{DimensionalData.Dimension,Symbol};
    figure::NamedTuple = (;),
    axis_top::NamedTuple = (;),
    axis_right::NamedTuple = (;),
    heatmap_kwargs...,
)
    # Default settings for figure and axes
    default_figure_setting = (size = (650, 500),)
    default_top_axis_setting = (xticklabelsvisible = false, ylabel = "Intensity")
    default_right_axis_setting = (yticklabelsvisible = false, xlabel = "Intensity")
    default_heatmap_setting = (colormap = :turbo,)
    heatmap_setting = merge(default_heatmap_setting, heatmap_kwargs)

    fig_kwargs = merge(default_figure_setting, figure)
    fig = Figure(; fig_kwargs...)

    # Resolve stack dimension (SHymbol → Dimension)
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
    heatmap!(ax_main, xvals, yvals, lift(a -> parent(a), A2); heatmap_setting...)

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

    # Vertical integration span (x-direction)
    vspan_range = lift(A2, pos, crosshair_thick_x) do a, p, tx
        _compute_vspan_range(a, p, tx)
    end

    vspan!(
        ax_main,
        lift(r -> r[1], vspan_range),
        lift(r -> r[2], vspan_range),
        color = v_span_color,
    )

    # Horizontal integration span (y-direction)
    hspan_range = lift(A2, pos, crosshair_thick_y) do a, p, ty
        _compute_hspan_range(a, p, ty)
    end

    hspan!(
        ax_main,
        lift(r -> r[1], hspan_range),
        lift(r -> r[2], hspan_range),
        color = h_span_color,
    )

    # Top plot: average along y-direction (horizontal slice)
    line_top_data = lift(A2, pos, crosshair_thick_y) do a, p, ty
        _compute_slice_line_top(a, p, ty)
    end
    lines!(ax_top, line_top_data, color = :blue)

    # Right plot: average along x-direction (vertical slice)
    line_right_data = lift(A2, pos, crosshair_thick_x) do a, p, tx
        _compute_slice_line_right(a, p, tx)
    end
    lines!(ax_right, line_right_data, color = :blue)

    # Label showing current position and averaged intensity
    label_text =
        lift(A2, pos, crosshair_thick_x, crosshair_thick_y, stack_idx) do a, p, tx, ty, si
            _make_label(
                a,
                p,
                tx,
                ty;
                stack_dim_name = string(name(stack_dim)),
                stack_val = stack_lookup[si],
            )
        end

    Label(
        fig[1, 2],
        label_text,
        tellwidth = false,
        halign = :left,
        font = :bold,
        fontsize = 16,
        padding = (5, 10, 10, 10),
        justification = :left,
    )

    # Link axes between main and projections
    linkxaxes!(ax_main, ax_top)
    linkyaxes!(ax_main, ax_right)

    # Update projection axis limits based on current slice.
    # update=true fires immediately on registration to set the initial limits.
    on(A2; update = true) do a
        valid_data = parent(a)[isfinite.(parent(a))]
        if isempty(valid_data)
            zmin, zmax = 0.0, 1.0
        else
            zmin, zmax = extrema(valid_data)
        end
        if zmin == zmax
            zmin -= 0.1
            zmax += 0.1
        end
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

# Generate the label string showing coordinates and averaged intensity.
# Pass stack_dim_name and stack_val together to include stack dimension info (3D data).
function _make_label(
    A::AbstractDimArray{T,2},
    pos::Point2f,
    crosshair_thick_x::Int,
    crosshair_thick_y::Int;
    stack_dim_name::Union{Nothing,AbstractString} = nothing,
    stack_val::Union{Nothing,Real} = nothing,
)::String where {T}
    nx, ny = size(A)
    ix, iy = _get_indices(pos, A)

    rx = max(1, ix-crosshair_thick_x):min(nx, ix+crosshair_thick_x)
    ry = max(1, iy-crosshair_thick_y):min(ny, iy+crosshair_thick_y)
    z = mean(parent(A[rx, ry]))

    lines = [
        "$(name(dims(A, 1))): $(round(lookup(A, 1)[ix], digits=4)) (±$crosshair_thick_x pts)",
        "$(name(dims(A, 2))): $(round(lookup(A, 2)[iy], digits=4)) (±$crosshair_thick_y pts)",
    ]
    if stack_dim_name !== nothing && stack_val !== nothing
        push!(lines, "$stack_dim_name: $(round(stack_val, digits=4))")
    end
    push!(lines, "Avg Intensity: $(round(z, digits=4))")
    return join(lines, "\n")
end
