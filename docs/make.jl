using FractionalDiffEq
using Documenter

DocMeta.setdocmeta!(FractionalDiffEq, :DocTestSetup, :(using FractionalDiffEq); recursive=true)

makedocs(;
    modules=[FractionalDiffEq],
    authors="Qingyu Qu",
    repo="https://github.com/SciFracX/FractionalDiffEq.jl/blob/{commit}{path}#{line}",
    sitename="FractionalDiffEq.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://SciFracX.github.io/FractionalDiffEq.jl",
        assets=["assets/favicon.ico"],
    ),
    pages=[
        "FractionalDiffEq.jl" => "index.md",
        "Get Started" => "get_start.md",
        "Detailed Models" => [
            "Relaxation Oscillation" => "relaxation-oscillation.md",
            "Bagley Torvik" => "bagley-torvik.md",
            "Diffusion" => "diffusion.md"
        ],
        "Multi-term FDE" => "multiterm.md",
        "Examples" => "example/example.md",
        "FractionalDiffEq APIs" => "APIs.md"
    ],
)

deploydocs(;
    repo="github.com/SciFracX/FractionalDiffEq.jl",
)
