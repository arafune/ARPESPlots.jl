using ARPESPlots
using Documenter

DocMeta.setdocmeta!(ARPESPlots, :DocTestSetup, :(using ARPESPlots); recursive=true)

makedocs(;
    modules=[ARPESPlots],
    authors="Ryuichi Arafune",
    sitename="ARPESPlots.jl",
    format=Documenter.HTML(;
        canonical="https://arafune.github.io/ARPESPlots.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/arafune/ARPESPlots.jl",
    devbranch="main",
)
