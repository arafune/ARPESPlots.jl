using DimensionalData
using Makie
using Makie: AbstractAxis, AbstractPlot, Axis, Figure, lines!, colormap, ylims!
using ARPES

export waterfall_dispersion, waterfall_dispersion!

"""
    waterfall_dispersion!(ax, A, stack_dim, scale_factor, cmap, mode, alpha)
    waterfall_dispersion!(ax, A, stack_dim; scale_factor=1.0, cmap=:turbo, mode=:line,
alpha=0.5)

Create a waterfall plot in an existing Axis, showing multiple slices of a multidimensional
array with vertical offsets.

# Arguments
- `ax::AbstractAxis`: The Axis to plot into
- `A::AbstractDimArray`: The input dimensional array to plot
- `stack_dim::Union{DimensionalData.Dimension,Symbol}`: The dimension along which to
  slice the array
- `scale_factor::Real`: Scaling factor for the vertical offset between slices
- `cmap::Symbol`: Color map to use for coloring different slices (e.g., `:turbo`, `:viridis`)
- `mode::Symbol`: Plot mode (:line, :fill, :hide)
- `alpha::Real`: Transparency value for the `:fill` and `:hide` plot elements

# Returns
- `Vector{AbstractPlot}`: A vector containing all the plotted line objects

# Description
This function creates a waterfall-style dispersion plot where each slice along the specified
dimension is plotted with a vertical offset proportional to its position along that dimension.
Each slice is colored using the specified colormap, creating a visually appealing 3D-like
effect that helps visualize how the data changes across the stacking dimension.

The vertical offset for each slice is calculated as
  `scale_factor * abs(stack_axis[i] - bottom)`,
where `bottom` is the last value in the stacking dimension and `stack_axis[i]` is the
current position along that dimension.

# Examples
```julia
# Create test data
data = DimArray(rand(100, 20), (Ti(1:100), Freq(1:20)))

# Plot into existing axis
fig = Figure()
ax = Axis(fig[1, 1], title="Waterfall Plot")
plots = waterfall_dispersion!(ax, data, :Freq, scale_factor=0.1)

# Using keyword arguments
plots = waterfall_dispersion!(ax, data, :Freq,
                             scale_factor=0.2,
                             cmap=:viridis)
```

# Notes

- All plotted line objects are returned in a vector for further manipulation if needed
- The mode and alpha parameters are defined but not currently used in the implementation
- Each line object in the returned vector corresponds to one slice of the input array
"""
function waterfall_dispersion!(
    ax::AbstractAxis,
    A::AbstractDimArray,
    stack_dim::Union{DimensionalData.Dimension,Symbol},
    scale_factor::Real,
    cmap::Symbol,
    mode::Symbol,
    alpha::Union{Real,Nothing},
    axis_right::NamedTuple,
)

    stack_axis = collect(lookup(A, stack_dim))
    n = length(stack_axis)
    bottom = last(stack_axis)

    colors = [cgrad(cmap)[i] for i in range(0, stop = 1, length = n)]
    plotted_objects = AbstractPlot[]

    map_to_right_axis(y_left) =
        (first(stack_axis)-last(stack_axis)) /
        (scale_factor * abs(stack_axis[1] - bottom)) * y_left + last(stack_axis)

    for (i, a_dimarray) in enumerate(eachslice(A, dims = stack_dim))
        offset = scale_factor * abs(stack_axis[i] - bottom)
        offset_array = rebuild(a_dimarray, zeros(size(a_dimarray))) .+ offset
        if mode == :fill
            alpha = isnothing(alpha) ? 0.5 : alpha
            band!(
                ax,
                lookup(a_dimarray, 1),
                parent(offset_array),
                parent(a_dimarray) .+ offset,
                color = colors[i],
                alpha = alpha,
            )
        elseif mode == :hide
            alpha = isnothing(alpha) ? 1.0 : alpha
            band!(
                ax,
                lookup(a_dimarray, 1),
                parent(offset_array),
                parent(a_dimarray) .+ offset,
                color = "white",
                alpha = alpha,
            )
        end
        line_obj = lines!(ax, a_dimarray .+ offset, color = colors[i])
        push!(plotted_objects, line_obj)
    end
    ax_right = Axis(
        ax.parent[1, 1],
        yaxisposition = :right,
        xgridvisible = false,
        ygridvisible = false,
        axis_right...,
    )
    hidespines!(ax_right)

    on(ax.finallimits) do lims
        (xmin, ymin), (xmax, ymax) = extrema(lims)
        ymin_trans, ymax_trans = map_to_right_axis.((ymin, ymax))
        limits!(ax_right, BBox(xmin, xmax, ymin_trans, ymax_trans))
    end

    tick_func = LinearTicks(5)
    r_tick_values =
        Makie.get_tickvalues(tick_func, minimum(stack_axis), maximum(stack_axis))
    ax_right.yticks = r_tick_values

    return plotted_objects, ax_right
end

waterfall_dispersion!(
    ax::AbstractAxis,
    A::AbstractDimArray,
    stack_dim::Union{DimensionalData.Dimension,Symbol};
    scale_factor::Real = 1.0,
    cmap::Symbol = :turbo,
    mode::Symbol = :line,
    alpha::Union{Real,Nothing} = nothing,
    axis_right::NamedTuple = (;),
) = waterfall_dispersion!(ax, A, stack_dim, scale_factor, cmap, mode, alpha, axis_right)

"""
    waterfall_dispersion(A, stack_dim; scale_factor=1.0, cmap=:turbo, mode=:line, alpha=0.5,
figure=(;), axis=(;), axis_right=(;)

Create a waterfall plot showing multiple slices of a multidimensional array with vertical
  offsets.

# Arguments
- `A::AbstractDimArray`: The input dimensional array to plot.
- `stack_dim::Union{DimensionalData.Dimension,Symbol}`: The dimension along which to slice
  the array.
- `scale_factor::Real = 1.0`: Scaling factor for the vertical offset between slices.
- `cmap::Symbol = :turbo`: Color map to use for coloring different slices.
- `mode::Symbol = :line` : Plot mode (:line, :fill, :hide)
- `alpha::Union{Real, Nothing} = nothing`: Transparency value for the plot elements
  (0.0 to 1.0).
- `figure::NamedTuple = (;)`: Configuration parameters for the Makie `Figure`
  (e.g., `resolution=(800, 600)`).
- `axis::NamedTuple = (;)`: Configuration parameters for the Makie `Axis`
  (e.g., `xlabel="Time"`, `title="Waterfall"`).
- `axis_right::NamedTuple = (;)`: Configuration parameters for the Makie `Axis` for the
  right side

# Returns
- `Makie.FigureAxisPlot`: A Makie object containing the Figure, Axis, and the plot object.

# Description
This function creates a waterfall-style dispersion plot where each slice along the specified
dimension is plotted with a vertical offset proportional to its position along that
dimension. The slices are colored using the specified colormap, creating a visually
appealing 3D-like effect.

The vertical offset for each slice is calculated based on its position along the stacking
dimension, multiplied by the `scale_factor`.

# Examples
```julia
data = DimArray(rand(100, 50), (Ti(1:100), Freq(1:50)))

# Basic usage
fap = waterfall_dispersion(data, :Freq)

# Customizing figure and axis via NamedTuples
fap = waterfall_dispersion(
    data, :Freq; 
    scale_factor=0.2, 
    figure=(size=(1000, 400), backgroundcolor=:gray90),
    axis=(title="My Waterfall Plot", xlabel="Time (s)")
)
```
# Notes

- The alpha and mode parameters are defined in the signature but their usage depends on
  the underlying waterfall_dispersion! implementation.
- Returns a FigureAxisPlot which can be unpacked or displayed directly in a
  Makie-compatible environment.

"""
function waterfall_dispersion(
    A::AbstractDimArray,
    stack_dim::Union{DimensionalData.Dimension,Symbol},
    scale_factor::Real,
    cmap::Symbol,
    mode::Symbol,
    alpha::Union{Real,Nothing},
    figure::NamedTuple,
    axis::NamedTuple,
    axis_right::NamedTuple,
)
    fig = Figure(; figure...)
    ax = Axis(fig[1, 1]; axis...)
    plot_obj, _ =
        waterfall_dispersion!(ax, A, stack_dim, scale_factor, cmap, mode, alpha, axis_right)

    return Makie.FigureAxisPlot(fig, ax, first(plot_obj))
end

waterfall_dispersion(
    A::AbstractDimArray,
    stack_dim::Union{DimensionalData.Dimension,Symbol};
    scale_factor::Real = 1.0,
    cmap::Symbol = :turbo,
    mode::Symbol = :line,
    alpha::Union{Real,Nothing} = nothing,
    figure::NamedTuple = (;),
    axis::NamedTuple = (;),
    axis_right::NamedTuple = (;),
) = waterfall_dispersion(
    A,
    stack_dim,
    scale_factor,
    cmap,
    mode,
    alpha,
    figure,
    axis,
    axis_right,
)


