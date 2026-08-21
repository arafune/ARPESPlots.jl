using DimensionalData
using Makie
using Printf
using ARPES

export tarpes_evolution_heatmaps

"""
    tarpes_evolution_heatmaps(A, delay_time=0.0, evolution_at=0.0; kwargs...)
    tarpes_evolution_heatmaps(A, delay_index::Integer, evolution_at=0.0; ...)

Create a Makie Figure showing an ARPES snapshot and the corresponding temporal evolution.

Layout
- Left: heatmap of the ARPES snapshot at `delay_time`.
- Center: heatmap of the temporal evolution at `evolution_at`.
- Right: vertical colorbar.

Keyword arguments
- stack_dim::Symbol, vertical_dim::Symbol: dimension names (defaults shown above).
- transparent::Bool: If true, and Makie backend is CairoMakie, the background is set transparent.
- figure::NamedTuple: forwarded to `Figure`.
- arpes_kwargs::NamedTuple: forwarded to the left Axis.
- evolution_kwargs::NamedTuple: forwarded to the right Axis.
- heatmap_kwargs::NamedTuple: forwarded to `heatmap!` (e.g., colorrange, colormap). Both heatmaps share these settings. For log scale, use `colorscale=log10` with an explicit `colorrange`.
- colorbar_kwargs::NamedTuple: forwarded to `Colorbar`.
- full_temporal::Bool: If true, the evolution heatmap uses the full dataset time range; otherwise it is limited to times <= `delay_time`.

Returns
- fig::Makie.Figure: A figure ready for display or further modification.

Notes
- By default, the heatmap colorrange uses the extrema of finite values in A.
- When the current Makie backend supports transparency (checked via `Makie.current_backend()`),
  the figure and axes backgrounds are set to transparent by default.
"""
function tarpes_evolution_heatmaps(
    A::AbstractDimArray{T,3} where {T<:AbstractFloat},
    delay_time::Real = 0.0,
    evolution_at::Union{Real,Tuple{<:Real,<:Real}} = 0.0;
    stack_dim::Symbol = :delay,
    vertical_dim::Symbol = :eV,
    transparent::Bool = false,
    figure::NamedTuple = (;),
    arpes_kwargs::NamedTuple = (;),
    evolution_kwargs::NamedTuple = (;),
    heatmap_kwargs::NamedTuple = (;),  # ex.) colorrange=(vmin, vmax), colorscale=log10, colormap = :turbo
    colorbar_kwargs::NamedTuple = (;),
    full_temporal::Bool = false,
)

    non_dispersion_axis = otherdims(A, (vertical_dim, stack_dim)) |> first |> name
    heatmap_kwargs_default = (colorrange = extrema(filter(!isnan, A)),)

    # default settings for figure and axes
    default_figure_setting = (size = (650, 350),)
    default_left_axis_setting = (
        xlabel = string(non_dispersion_axis),
        ylabel = string(vertical_dim),
        xticks = WilkinsonTicks(3),
        xgridvisible = false,
        ygridvisible = false,
    )
    default_right_axis_setting = (
        xlabel = string(stack_dim),
        ylabel = string(vertical_dim),
        yticklabelsvisible = false,
        ylabelvisible = false,
        xgridvisible = false,
        ygridvisible = false,
    )
    default_colorbar_setting =
        (label = "Intensity", vertical = true, ticklabelsvisible = true)

    if nameof(Makie.current_backend()) == :CairoMakie && transparent
        default_figure_setting =
            merge(default_figure_setting, (backgroundcolor = :transparent,))
        default_left_axis_setting =
            merge(default_left_axis_setting, (backgroundcolor = :transparent,))
        default_right_axis_setting =
            merge(default_right_axis_setting, (backgroundcolor = :transparent,))
    end

    # keyword argument merging: user-specified settings override defaults
    fig_kwargs = merge(default_figure_setting, figure)
    left_axis_kwargs = merge(default_left_axis_setting, arpes_kwargs)
    right_axis_kwargs = merge(default_right_axis_setting, evolution_kwargs)
    heatmap_kwargs = merge(heatmap_kwargs_default, heatmap_kwargs)
    colorbar_kwargs = merge(default_colorbar_setting, colorbar_kwargs)

    fig = Figure(; fig_kwargs...)
    left_axis = Axis(fig[1, 1]; left_axis_kwargs...)
    right_axis = Axis(fig[1, 2]; right_axis_kwargs...)
    colsize!(fig.layout, 1, Relative(0.25))

    arpes_data_at_delay, temporal_evolution_data = tarpes_evolution(
        A,
        delay_time,
        evolution_at,
        stack_dim = stack_dim,
        vertical_dim = vertical_dim,
        full_temporal = full_temporal,
    )

    z_left = Observable(arpes_data_at_delay)
    z_right = Observable(temporal_evolution_data)

    xlims!(right_axis, dims(A, stack_dim) |> extrema)
    hm_left = heatmap!(left_axis, z_left; heatmap_kwargs...)
    hm_right = heatmap!(right_axis, z_right; heatmap_kwargs...)
    Colorbar(fig[1, 3], hm_left; colorbar_kwargs...)

    return fig
end

function tarpes_evolution_heatmaps(
    A::AbstractDimArray{T,3} where {T<:AbstractFloat},
    delay_index::Integer,
    evolution_at::Union{Real,Tuple{<:Real,<:Real}} = 0.0;
    stack_dim::Symbol = :delay,
    vertical_dim::Symbol = :eV,
    transparent::Bool = false,
    figure::NamedTuple = (;),
    arpes_kwargs::NamedTuple = (;),
    evolution_kwargs::NamedTuple = (;),
    heatmap_kwargs::NamedTuple = (;),
    colorbar_kwargs::NamedTuple = (;),
    full_temporal::Bool = false,
)
    delay_time = dims(A, stack_dim)[delay_index] |> float
    return tarpes_evolution_heatmaps(
        A,
        delay_time,
        evolution_at;
        stack_dim = stack_dim,
        vertical_dim = vertical_dim,
        transparent = transparent,
        figure = figure,
        arpes_kwargs = arpes_kwargs,
        evolution_kwargs = evolution_kwargs,
        heatmap_kwargs = heatmap_kwargs,
        colorbar_kwargs = colorbar_kwargs,
        full_temporal = full_temporal,
    )
end


"""
    tarpes_evolution_mp4(A, evolution_at=0.0; outpath="output.mp4", fps=25, ...)

Create an MP4 animation iterating the `stack_dim` (default :delay) frames.
`outpath` should end with .mp4. Uses Makie.record and updates two Observables
for the left (ARPES snapshot) and right (temporal evolution) heatmaps.

Keyword arguments
- outpath::AbstractString: output path for the MP4 file.
- fps::Integer: frames per second for the MP4 file.
- stack_dim::Symbol, vertical_dim::Symbol: dimension names (defaults shown above).
- figure::NamedTuple: forwarded to `Figure`.
- arpes_kwargs::NamedTuple: forwarded to the left Axis.
- evolution_kwargs::NamedTuple: forwarded to the center Axis.
- heatmap_kwargs::NamedTuple: forwarded to `heatmap!` (e.g., colorrange, colormap). Both heatmaps share these settings. For log scale, use `colorscale=log10` with an explicit `colorrange`.
- colorbar_kwargs::NamedTuple: forwarded to `Colorbar`.
- full_temporal::Bool: If true, the evolution heatmap uses the full dataset time range; otherwise it is limited to times <= `delay_time`.
- transparent::Bool: If true, produce alpha-preserving frames and (optionally) encode with ffmpeg.
- encoder::Symbol: :webm or :mov when transparent=true. Determines the ffmpeg encoding settings for the output video.


"""
function tarpes_evolution_mp4(
    A::AbstractDimArray{T,3} where {T},
    evolution_at::Union{Real,Tuple{<:Real,<:Real}} = 0.0;
    outpath::AbstractString = "output.mp4",
    fps::Integer = 25,
    stack_dim::Symbol = :delay,
    vertical_dim::Symbol = :eV,
    figure::NamedTuple = (;),
    arpes_kwargs::NamedTuple = (;),
    evolution_kwargs::NamedTuple = (;),
    heatmap_kwargs::NamedTuple = (;),
    colorbar_kwargs::NamedTuple = (;),
    full_temporal::Bool = false,
    transparent::Bool = false,    # if true, produce alpha-preserving frames and (optionally) encode with ffmpeg
    encoder::Symbol = :webm,     # :webm or :mov when transparent=true
)
    # frame indices across the stack dimension
    indices = 1:(dims(A, stack_dimlay)|>length)
    nframes = length(dims(A, stack_dim))

    # prepare initial data
    delay_time = dims(A, stack_dim)[1] |> float
    arpes_data_at_delay, temporal_evolution_data = tarpes_evolution(
        A,
        delay_time,
        evolution_at,
        stack_dim = stack_dim,
        vertical_dim = vertical_dim,
        full_temporal = full_temporal,
    )

    # copy defaults from tarpes_evolution_heatmaps
    non_dispersion_axis = otherdims(A, (vertical_dim, stack_dim)) |> first |> name
    heatmap_kwargs_default = (colorrange = extrema(filter(!isnan, A)),)

    default_figure_setting = (size = (650, 350),)
    default_left_axis_setting = (
        xlabel = string(non_dispersion_axis),
        ylabel = string(vertical_dim),
        xticks = WilkinsonTicks(3),
        xgridvisible = false,
        ygridvisible = false,
    )
    default_right_axis_setting = (
        xlabel = string(stack_dim),
        ylabel = string(vertical_dim),
        yticklabelsvisible = false,
        ylabelvisible = false,
        xgridvisible = false,
        ygridvisible = false,
    )
    default_colorbar_setting =
        (label = "Intensity", vertical = true, ticklabelsvisible = true)

    if nameof(Makie.current_backend()) == :CairoMakie && transparent
        default_figure_setting =
            merge(default_figure_setting, (backgroundcolor = :transparent,))
        default_left_axis_setting =
            merge(default_left_axis_setting, (backgroundcolor = :transparent,))
        default_right_axis_setting =
            merge(default_right_axis_setting, (backgroundcolor = :transparent,))
    end

    fig_kwargs = merge(default_figure_setting, figure)
    left_axis_kwargs = merge(default_left_axis_setting, arpes_kwargs)
    right_axis_kwargs = merge(default_right_axis_setting, evolution_kwargs)
    heatmap_kwargs = merge(heatmap_kwargs_default, heatmap_kwargs)
    colorbar_kwargs = merge(default_colorbar_setting, colorbar_kwargs)

    fig = Figure(; fig_kwargs...)
    left_axis = Axis(fig[1, 1]; left_axis_kwargs...)
    right_axis = Axis(fig[1, 2]; right_axis_kwargs...)
    colsize!(fig.layout, 1, Relative(0.25))

    # Observables so frames can be updated in-place
    z_left = Observable(arpes_data_at_delay)
    z_right = Observable(temporal_evolution_data)

    hm_left = heatmap!(left_axis, z_left; heatmap_kwargs...)
    hm_right = heatmap!(right_axis, z_right; heatmap_kwargs...)
    Colorbar(fig[1, 3], hm_left; colorbar_kwargs...)

    # Record frames and optionally encode with alpha-preserving encoder
    if !transparent
        # Record into MP4 (or extension-chosen) by updating observables per frame
        Makie.record(fig, outpath, 1:nframes; framerate = fps) do i
            i
            delay_time = dims(A, stack_dim)[i] |> float
            left_z, right_z = tarpes_evolution(
                A,
                delay_time,
                evolution_at,
                stack_dim = stack_dim,
                vertical_dim = vertical_dim,
                full_temporal = full_temporal,
            )
            z_left[] = left_z
            z_right[] = right_z
        end

    else
        # Save PNG frames (preserves alpha) then encode with ffmpeg to desired container
        mktempdir() do tmp
            for i = 1:nframes
                delay_time = dims(A, stack_dim)[i] |> float
                left_z, right_z = tarpes_evolution(
                    A,
                    delay_time,
                    evolution_at,
                    stack_dim = stack_dim,
                    vertical_dim = vertical_dim,
                    full_temporal = full_temporal,
                )
                z_left[] = left_z
                z_right[] = right_z
                framepath = joinpath(tmp, @sprintf("frame_%04d.png", i))
                save(framepath, fig)
            end

            # Choose ffmpeg command based on encoder
            if encoder == :webm
                cmd = `ffmpeg -y -framerate $fps -i $(joinpath(tmp, "frame_%04d.png")) -c:v libvpx-vp9 -pix_fmt yuva420p -auto-alt-ref 0 -b:v 0 $outpath`
            elseif encoder == :mov
                cmd = `ffmpeg -y -framerate $fps -i $(joinpath(tmp, "frame_%04d.png")) -c:v prores_ks -profile:v 4444 -pix_fmt yuva444p10le $outpath`
            else
                error("Unsupported encoder: $encoder. Use :webm or :mov")
            end

            try
                run(cmd)
            catch e
                error("ffmpeg encoding failed: ", e)
            end
        end
    end

    return outpath
end
