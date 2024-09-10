using Documenter
using SPRFittingPaper2023

makedocs(
    sitename = "SPRFittingPaper2023.jl",
    authors = "Samuel Isaacson",
    format = Documenter.HTML(prettyurls = (get(ENV, "CI", nothing) == "true")),
    modules = [SPRFittingPaper2023],
    doctest = false,
    clean = true,
    pages = Any[
        "Home" => "index.md",
        "Forward Model" => "forward_simulation.md",
        "Surrogate Construction" => "surrogate_construction.md",
        "Fitting Data" => "fitting.md",
        "API" => "SPRFitting_api.md"
    ],
    warnonly = [:missing_docs]
)

deploydocs(
   repo = "github.com/isaacsas/SPRFittingPaper2023.jl.git";
   push_preview = true
)
