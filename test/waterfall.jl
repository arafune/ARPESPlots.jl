using Test
using CairoMakie
using DimensionalData
using DimensionalData: @dim
using Statistics
using Makie
using ARPES

# Force headless backend
CairoMakie.activate!()

@testset "waterfall_dispersion!" begin
    # -------------------------------
    # Test data
    # -------------------------------
    t = 1:50
    f = 1:10
    @dim Freq
    data = DimArray(rand(length(t), length(f)), (Ti(t), Freq(f)))

    fig = Figure()
    ax = Axis(fig[1, 1])

    # -------------------------------
    # Basic execution (no error)
    # -------------------------------
    plots, ax_right =
        waterfall_dispersion!(ax, data, :Freq; scale_factor = 0.1, cmap = :viridis)

    @test length(plots) == length(f)
    @test ax_right isa Axis

    # -------------------------------
    # Each plot is a Makie plot object
    # -------------------------------
    @test all(p -> p isa AbstractPlot, plots)

    # -------------------------------
    # Offset is actually applied
    # -------------------------------
    stack_axis = collect(lookup(data, :Freq))
    bottom = last(stack_axis)

    for (i, p) in enumerate(plots)
        # Extract y data from plotted object
        pts = p[1][]
        y = getindex.(pts, 2)  # 
        expected_offset = 0.1 * abs(stack_axis[i] - bottom)

        # Check mean shift roughly matches
        original = parent(data[:, i])
        @test isapprox(mean(y .- original), expected_offset; atol = 1e-6)
    end

    # -------------------------------
    # Right axis ticks exist
    # -------------------------------
    @test !isnothing(ax_right.yticks)

    # -------------------------------
    # Mode = :fill
    # -------------------------------
    fig2 = Figure()
    ax2 = Axis(fig2[1, 1])

    plots2, _ = waterfall_dispersion!(ax2, data, :Freq; mode = :fill, alpha = 0.3)

    @test length(plots2) == length(f)

    # -------------------------------
    # Mode = :hide
    # -------------------------------
    fig3 = Figure()
    ax3 = Axis(fig3[1, 1])

    plots3, _ = waterfall_dispersion!(ax3, data, :Freq; mode = :hide)

    @test length(plots3) == length(f)
end


@testset "waterfall_dispersion! (3D)" begin
    t = 1:50
    f = 1:10
    @dim Freq
    data = DimArray(rand(length(t), length(f)), (Ti(t), Freq(f)))
    fig_3d_1 = Figure()
    ax3d_1 = Axis3(fig_3d_1[1,1])
    plot3d_1 = waterfall_dispersion!(ax3d_1, data, :Freq)
    @test length(plot3d_1) == length(f)

    fig_3d_2 = Figure()
    ax3d_2 = Axis(fig_3d_2[1, 1])
    plots3d_2 = waterfall_dispersion!(ax3d_2, data, :Freq; mode = :hide)
    @test length(plots3d_2) == 2

    fig_3d_3 = Figure()
    ax3d_3 = Axis(fig_3d_3[1, 1])
    plots3d_3 = waterfall_dispersion!(ax3d_3, data, :Freq; mode = :fill)
    @test length(plots3d_3) == 2

end


@testset "waterfall_dispersion (wrapper)" begin
    t = 1:30
    f = 1:5
    data = DimArray(rand(length(t), length(f)), (Ti(t), Freq(f)))

    fap = waterfall_dispersion(data, :Freq)

    @test fap isa Makie.FigureAxisPlot

    fig, ax, plt = fap

    @test fig isa Figure
    @test ax isa Axis
    @test plt isa AbstractPlot
end
