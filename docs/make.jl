using FractionalDiffEq, Documenter

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
        "Get Started" => "get_started.md",
        "Multi-term Fractional Differential Equations" => "multiterm.md",
        "System of Fractional Differential Equations" => "system_of_FDE.md",
        "Fractional Partial Differential Equations" => "fpde.md",
        "Fractional Delay Differential Equations" => "fdde.md",
        "Distributed Order Differential Equations" => "dode.md",
        "Detailed Models" => "models.md",
        "Examples" => "example.md",
        "Algorithms" => "algorithms.md",
        "Auxiliary Functions" => "auxiliary.md",
        "Comparison" => "comparison.md",
        "FractionalDiffEq.jl APIs" => "APIs.md"
    ],
)

deploydocs(;
    repo="github.com/SciFracX/FractionalDiffEq.jl",
)
