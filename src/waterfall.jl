using DimensionalData
using Makie
using Makie: AbstractAxis, AbstractPlot, Axis, Figure, lines!, colormap, ylims!
using ARPES

export waterfall_dispersion, waterfall_dispersion!

"""
    waterfall_dispersion!(ax::AbstractAxis, A::ARPESData; stack_dim=:phi, scale_factor=1.0, cmap=:turbo, mode=:line, alpha=nothing, axis_right=(;))
    waterfall_dispersion!(ax::AbstractAxis, A::AbstractDimArray, stack_dim; scale_factor=1.0, cmap=:turbo, mode=:line, alpha=nothing, axis_right=(;))
    waterfall_dispersion!(ax::Axis3, A::AbstractDimArray, stack_dim; cmap=:turbo, mode=:line, alpha=nothing)

Plot a waterfall dispersion into an existing axis, showing multiple slices of a 2D
`AbstractDimArray` stacked with vertical offsets (2D axis variants) or as 3D ribbon lines
(`Axis3` variant).

# Arguments
- `ax`: Target axis — an `AbstractAxis` (2D) or `Axis3` (3D).
- `A`: 2D dimensional array to plot.
- `stack_dim::Union{DimensionalData.Dimension,Symbol}`: Dimension along which to slice.
  Defaults to `:phi` when `A` is an `ARPESData`.
- `scale_factor::Real = 1.0`: Vertical offset multiplier between slices (2D only).
- `cmap = :turbo`: Colormap / Color applied across slices.
- `mode::Symbol = :line`: Rendering mode — `:line` (lines only), `:fill` (lines + colored
  bands), or `:hide` (lines + white bands to mask lower slices).
- `alpha::Union{Real,Nothing} = nothing`: Opacity for band elements; defaults to `0.5`
  (`:fill`) or `1.0` (`:hide`) when `nothing`.
- `axis_right::NamedTuple = (;)`: Keyword arguments forwarded to the secondary right-hand
  `Axis` (2D only; ignored for `Axis3`).
- `kwargs...`: Additional keyword arguments forwarded to `lines!` (e.g. `linewidth`,
  `linestyle`).

# Returns
- **2D variants**: `(Vector{AbstractPlot}, Axis)` — plotted line objects and the linked
  right-hand axis whose ticks reflect the original stacking-dimension values.
- **`Axis3` variant**: `Vector{AbstractPlot}` — plotted line objects.

# Description
Each slice along `stack_dim` is drawn as a line (and optionally a band) offset vertically by

    offset = scale_factor * abs(stack_axis[i] - stack_axis[end])

so that the first slice sits highest and the last slice sits at the baseline. Slices are
colored by sampling `cmap` uniformly across the number of slices.

For the 2D variants a secondary right-hand `Axis` is created and linked to the left axis,
with tick values drawn from the actual stacking-dimension coordinates.

# Examples
```julia
using DimensionalData, GLMakie

E  = range(-2.0, 0.0, length=200)
kx = range(-1.0, 1.0, length=40)
A  = DimArray(rand(200, 40), (X(E), Y(kx)))

fig = Figure()
ax  = Axis(fig[1, 1], xlabel="Energy (eV)", ylabel="Intensity (arb. u.)")
plots, ax_right = waterfall_dispersion!(ax, A, :Y; scale_factor=0.05, cmap=:viridis)
ax_right.ylabel = "kx (Å⁻¹)"
display(fig)
```
"""
function waterfall_dispersion!(
    ax::AbstractAxis,
    A::ARPESData;
    stack_dim::Union{DimensionalData.Dimension,Symbol} = :phi,
    scale_factor::Real = 1.0,
    cmap = :turbo,
    mode::Symbol = :line,
    alpha::Union{Real,Nothing} = nothing,
    axis_right::NamedTuple = (;),
    kwargs...,
)
    return waterfall_dispersion!(
        ax,
        A,
        stack_dim,
        scale_factor,
        cmap,
        mode,
        alpha,
        axis_right,
        kwargs...,
    )
end

function waterfall_dispersion!(
    ax::Axis3,
    A::AbstractDimArray{T,2} where {T},
    stack_dim::Union{DimensionalData.Dimension,Symbol};
    cmap = :turbo,
    mode::Symbol = :line,
    alpha::Union{Real,Nothing} = nothing,
    kwargs...,
)
    stack_dim = dims(A, stack_dim)
    xi = otherdims(A, stack_dim)[1]
    plotted_objects = AbstractPlot[]

    colors = _make_colors(cmap, length(stack_dim))

    for (i, a_dimarry) in enumerate(eachslice(A, dims = stack_dim))
        yi = fill(stack_dim[i], length(lookup(A, xi)))
        if mode == :fill
            alpha_val = isnothing(alpha) ? 0.5 : alpha
            band_obj = band!(
                ax,
                Point3.(lookup(a_dimarry, xi), yi, 0.0),
                Point3.(lookup(a_dimarry, xi), yi, parent(a_dimarry)),
                color = colors[i],
                alpha = alpha_val,
            )
            push!(plotted_objects, band_obj)
        elseif mode == :hide
            alpha_val = isnothing(alpha) ? 1.0 : alpha
            band_obj = band!(
                ax,
                Point3.(lookup(a_dimarry, xi), yi, 0.0),
                Point3.(lookup(a_dimarry, xi), yi, parent(a_dimarry)),
                color = :white,
                alpha = alpha_val,
            )
            push!(plotted_objects, band_obj)
        end
        line_obj = lines!(
            ax,
            lookup(a_dimarry, xi),
            yi,
            parent(a_dimarry);
            color = colors[i],
            kwargs...,
        )
        push!(plotted_objects, line_obj)
    end
    return plotted_objects
end

function waterfall_dispersion!(
    ax::AbstractAxis,
    A::AbstractDimArray{T,2} where {T},
    stack_dim::Union{DimensionalData.Dimension,Symbol},
    scale_factor::Real,
    cmap,
    mode::Symbol,
    alpha::Union{Real,Nothing},
    axis_right::NamedTuple,
    kwargs...,
)
    stack_axis = collect(lookup(A, stack_dim))
    bottom = last(stack_axis)

    colors = _make_colors(cmap, length(stack_axis))
    plotted_objects = AbstractPlot[]

    map_to_right_axis(y_left) =
        (first(stack_axis)-last(stack_axis)) /
        (scale_factor * abs(stack_axis[1] - bottom)) * y_left + last(stack_axis)

    for (i, a_dimarray) in enumerate(eachslice(A, dims = stack_dim))
        offset = scale_factor * abs(stack_axis[i] - bottom)
        if mode == :fill
            alpha = isnothing(alpha) ? 0.5 : alpha
            band!(
                ax,
                lookup(a_dimarray, 1),
                fill(offset, size(a_dimarray, 1)),
                parent(a_dimarray) .+ offset,
                color = colors[i],
                alpha = alpha,
            )
        elseif mode == :hide
            alpha = isnothing(alpha) ? 1.0 : alpha
            band!(
                ax,
                lookup(a_dimarray, 1),
                fill(offset, size(a_dimarray, 1)),
                parent(a_dimarray) .+ offset,
                color = :white,
                alpha = alpha,
            )
        end
        line_obj = lines!(ax, a_dimarray .+ offset; color = colors[i], kwargs...)
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
    A::AbstractDimArray{T,2} where {T},
    stack_dim::Union{DimensionalData.Dimension,Symbol};
    scale_factor::Real = 1.0,
    cmap = :turbo,
    mode::Symbol = :line,
    alpha::Union{Real,Nothing} = nothing,
    axis_right::NamedTuple = (;),
    kwargs...,
) = waterfall_dispersion!(
    ax,
    A,
    stack_dim,
    scale_factor,
    cmap,
    mode,
    alpha,
    axis_right,
    kwargs...,
)

"""
    waterfall_dispersion(A, stack_dim; scale_factor=1.0, cmap=:turbo, mode=:line, alpha=nothing, figure=(;), axis=(;), axis_right=(;))

Create a waterfall dispersion plot, returning a `Makie.FigureAxisPlot`.

# Arguments
- `A::AbstractDimArray`: 2D dimensional array to plot.
- `stack_dim::Union{DimensionalData.Dimension,Symbol}`: Dimension along which to slice.
- `scale_factor::Real = 1.0`: Vertical offset multiplier between slices.
- `cmap = :turbo`: Colormap /color applied across slices.
- `mode::Symbol = :line`: Rendering mode — `:line`, `:fill`, or `:hide`.
- `alpha::Union{Real,Nothing} = nothing`: Band opacity; defaults to `0.5` (`:fill`) or
  `1.0` (`:hide`) when `nothing`.
- `figure::NamedTuple = (;)`: Keyword arguments forwarded to `Figure` (e.g. `size=(800,600)`).
- `axis::NamedTuple = (;)`: Keyword arguments forwarded to the left `Axis`
  (e.g. `xlabel="Energy (eV)"`).
- `axis_right::NamedTuple = (;)`: Keyword arguments forwarded to the linked right-hand `Axis`.
- `kwargs...`: Additional keyword arguments forwarded to `lines!` (e.g. `linewidth`,
  `linestyle`).

# Returns
- `Makie.FigureAxisPlot`: Contains the `Figure`, left `Axis`, and the first plot object.
  The linked right-hand axis (with stacking-dimension tick labels) is accessible via
  `waterfall_dispersion!` if needed.

# Examples
```julia
using DimensionalData, GLMakie

E  = range(-2.0, 0.0, length=200)
kx = range(-1.0, 1.0, length=40)
A  = DimArray(rand(200, 40), (X(E), Y(kx)))

fap = waterfall_dispersion(A, :Y;
    scale_factor = 0.05,
    cmap         = :plasma,
    figure       = (size = (900, 500),),
    axis         = (xlabel = "Energy (eV)", ylabel = "Intensity (arb. u.)"),
    axis_right   = (ylabel = "kx (Å⁻¹)",),
)
display(fap)
```
"""
function waterfall_dispersion(
    A::AbstractDimArray,
    stack_dim::Union{DimensionalData.Dimension,Symbol},
    scale_factor::Real,
    cmap,
    mode::Symbol,
    alpha::Union{Real,Nothing},
    figure::NamedTuple,
    axis::NamedTuple,
    axis_right::NamedTuple,
    kwargs...,
)
    fig = Figure(; figure...)
    ax = Axis(fig[1, 1]; axis...)
    plot_obj, _ = waterfall_dispersion!(
        ax,
        A,
        stack_dim,
        scale_factor,
        cmap,
        mode,
        alpha,
        axis_right,
        kwargs...,
    )

    return Makie.FigureAxisPlot(fig, ax, first(plot_obj))
end

waterfall_dispersion(
    A::AbstractDimArray,
    stack_dim::Union{DimensionalData.Dimension,Symbol};
    scale_factor::Real = 1.0,
    cmap = :turbo,
    mode::Symbol = :line,
    alpha::Union{Real,Nothing} = nothing,
    figure::NamedTuple = (;),
    axis::NamedTuple = (;),
    axis_right::NamedTuple = (;),
    kwargs...,
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
    kwargs...,
)

function _make_colors(cmap, n::Integer)
    try
        return [cgrad(cmap)[i] for i in range(0, stop = 1, length = n)]
    catch
        try
            color = Makie.to_color(cmap)
            return fill(color, n)
        catch
            throw(ArgumentError("Invalid colormap or color: $cmap"))
        end
    end
end

