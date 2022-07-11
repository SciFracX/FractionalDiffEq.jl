using FractionalDiffEq, Documenter, DemoCards, JSON

ChaosGallery, postprocess_cb, ChaosGallery_assets = makedemos("ChaosGallery")

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
    pages = [
        "FractionalDiffEq.jl" => "index.md",
        "Get Started" => "get_started.md",
        "Multi-terms Fractional Differential Equations" => "multiterm.md",
        "System of Fractional Differential Equations" => "fdesystem.md",
        "Fractional Partial Differential Equations" => "fpde.md",
        "Fractional Delay Differential Equations" => "fdde.md",
        "Distributed Order Differential Equations" => "dode.md",
        "Fractional Differences Equations" => "differences.md",
        "Fractal-Fractional Differential Equations" => "ffode.md",
        "Fractional Integral Equations" => "fie.md",
        "Ordinary Differential Equations" => "ode.md",
        ChaosGallery,
        "Problems" => "problems.md",
        "Algorithms" => "algorithms.md",
        "Mittag Leffler Function" => "mittagleffler.md",
        "Comparison" => "comparison.md"
    ],
)

postprocess_cb()

deploydocs(;
    repo="github.com/SciFracX/FractionalDiffEq.jl",
)
