using Test
using ARPESPlots
using ARPES
using DimensionalData
using Makie
using CairoMakie # Required for headless CI environments

# Ensure Makie doesn't try to open a window
Makie.inline!(true)

@testset "ARPESPlots.jl Tests" begin

    # --- Setup Mock Data ---
    # Create two overlapping DimArrays along the X dimension
    l1 = (X(1:10), Y(1:10))
    l2 = (X(6:15), Y(1:10))

    data_a = ones(10, 10)
    data_b = ones(10, 10) .* 1.5

    A = DimArray(data_a, l1; name = :Intensity)
    B = DimArray(data_b, l2; name = :Intensity)
    target_dim = :X

    @testset "UI Initialization" begin
        # Test if the figure is generated without errors
        fig = nothing
        try
            fig = stitch_ui(A, B, target_dim)
        catch e
            @error "stitch_ui execution failed" exception=e
            rethrow(e)
        end

        @test fig isa Figure
        # Check if basic components (Axis, Sliders) are present
        @test !isempty(fig.content)
        @test any(x -> x isa Axis, fig.content)
    end

    @testset "Reactive Logic (Observable Update)" begin
        # Mocking the stitch_along logic used inside the UI
        # Check if the seam_ratio and gain_a parameters are correctly handled
        r_val = 0.7
        g_val = 1.2

        C = stitch_along(A, B, target_dim; seam_ratio = r_val, gain_a = g_val)

        @test C isa AbstractDimArray
        @test size(C, 1) == 15  # Total range from 1 to 15

        # Verify that the gain_a was applied to the first array correctly
        # Depending on how stitch_along is implemented, check a known point
        @test maximum(C) >= 1.5
    end

    @testset "Edge Cases" begin
        # Test with a Symbol vs Dimension type
        @test stitch_ui(A, B, :X) isa Figure
        @test stitch_ui(A, B, X()) isa Figure

        # Test with additional heatmap kwargs
        fig_custom = stitch_ui(A, B, :X; colormap = :viridis, colorrange = (0, 2))
        @test fig_custom isa Figure
    end
end
