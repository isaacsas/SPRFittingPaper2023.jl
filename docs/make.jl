using Documenter
using SPRFitting

makedocs(
    sitename = "SPRFitting.jl",
    authors = "Samuel Isaacson",
    format = Documenter.HTML(mathengine=Documenter.Writers.HTMLWriter.MathJax(), prettyurls = (get(ENV, "CI", nothing) == "true")),
    modules = [SPRFitting],
    doctest = false,
    clean = true,
    pages = Any[
        "Home" => "index.md",
        "API" => "SPRFitting_api.md"
    ]
)

# deploydocs(
#    repo = "github.com/isaacsas/SPRFitting.jl.git";
#    push_preview = true
# )