using Test
using DimensionalData
using DimensionalData: @dim
using CairoMakie
using ARPES
using ARPES: ARPESData, kx, ky, kz, phi, psi, eV
using ARPESPlots

@testset "ARPESPlots Interactive Heatmap Tests" begin

    # Create Dummy Data
    data = rand(Float64, 40, 60)
    A = ARPESData(data, (phi(range(-10, 10, 40)), eV(range(0, 5, 60))))

    @testset "Figure Generation" begin
        fig = nothing
        # Ensure no errors occur during figure creation
        @test_nowarn fig = crosshair_heatmap(A)
        @test fig isa Makie.Figure
    end

    @testset "Argument Customization" begin
        # Test if custom named tuples are merged correctly
        @test_nowarn crosshair_heatmap(
            A;
            figure = (size = (500, 500),),
            axis_top = (title = "Top View",),
        )
    end

    @testset "Headless Rendering" begin
        fig = crosshair_heatmap(A)
        # Verify the figure can actually be saved to a file (renders the pipeline)
        mktempdir() do tmp
            path = joinpath(tmp, "test_render.png")
            save(path, fig)
            @test isfile(path)
            @test filesize(path) > 0
        end
    end
end
