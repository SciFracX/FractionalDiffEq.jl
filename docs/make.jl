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
        "Multi-term FDE" => "multiterm.md",
        "System of FDE" => "system_of_FDE.md",
        "Detailed Models" => "models.md",
        "Examples" => "example.md",
        "FractionalDiffEq APIs" => "APIs.md"
    ],
)

deploydocs(;
    repo="github.com/SciFracX/FractionalDiffEq.jl",
)
