using ARPESPlots
using Test

@testset "ARPESPlots.jl" begin
    @testset "test for crosshair_heatmap" begin    # Write your tests here.
        include("./crosshair_heatmap.jl")
        include("./waterfall.jl")
    end
end
