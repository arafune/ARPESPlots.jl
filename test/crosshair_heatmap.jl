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
    data[3, 3] = NaN
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

    @testset "Log-scaled projections" begin
        data_log = copy(data)
        data_log[1, 1] = 0.0
        A_log = ARPESData(data_log, (phi(range(-10, 10, 40)), eV(range(0, 5, 60))))
        @test_nowarn crosshair_heatmap(
            A_log; axis_top = (yscale = log10,), axis_right = (xscale = log10,)
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


@testset "crosshair_heatmap 3D Tests" begin
    # Create Dummy Data
    data = rand(Float64, 40, 60, 30)
    data[3, 3, 3] = NaN
    A = ARPESData(
        data,
        (phi(range(-10, 10, 40)), eV(range(0, 5, 60)), psi(range(-5, 5, 30))),
    )
    @testset "Figure Generation" begin
        fig = nothing
        # Ensure no errors occur during figure creation
        @test_nowarn fig = crosshair_heatmap(A, :psi)
        @test fig isa Makie.Figure
    end

    @testset "Argument Customization" begin
        # Test if custom named tuples are merged correctly
        @test_nowarn crosshair_heatmap(
            A,
            :psi;
            figure = (size = (500, 500),),
            axis_top = (title = "Top View",),
        )
    end

    @testset "Headless Rendering" begin
        fig = crosshair_heatmap(A, :psi)
        # Verify the figure can actually be saved to a file (renders the pipeline)
        mktempdir() do tmp
            path = joinpath(tmp, "test_render.png")
            save(path, fig)
            @test isfile(path)
            @test filesize(path) > 0
        end
    end

end
