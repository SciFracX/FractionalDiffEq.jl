# Set plot globals
ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "nul"

using FractionalDiffEq, Documenter
gr()
default(show=false, size=(800,450))

dir = joinpath(dirname(pathof(FractionalDiffEq)), "..")
cd(dir)

DocMeta.setdocmeta!(FractionalDiffEq, :DocTestSetup, :(using FractionalDiffEq); recursive=true)


# Copied from Documenter/src/Document.jl, modified to remove # hide lines
Markdown = Documenter.Documents.Markdown
function Documenter.Documents.doctest_replace!(block::Markdown.Code)
    startswith(block.language, "jldoctest") || return false
    # suppress output for `#output`-style doctests with `output=false` kwarg
    if occursin(r"^# output$"m, block.code) && occursin(r";.*output\h*=\h*false", block.language)
        input = first(split(block.code, "# output\n", limit = 2))
        block.code = rstrip(input)
    end
    # Remove # hide lines
    block.code = Documenter.Expanders.droplines(block.code)
    # correct the language field
    block.language = occursin(r"^julia> "m, block.code) ? "julia-repl" : "julia"
    return false
end

println("Making docs")


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
    strict=[
        :doctest, 
        :linkcheck, 
        :parse_error, 
        # Other available options are
        # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block, :footnote, :meta_block, :missing_docs, :setup_block
    ],
    pages=[
        "FractionalDiffEq.jl" => "index.md",
        "Get Started" => "get_start.md",
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
