using Documenter
using ODEInterfaceDiffEq

makedocs(;
    modules = [ODEInterfaceDiffEq],
    sitename = "ODEInterfaceDiffEq.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://docs.sciml.ai/ODEInterfaceDiffEq/stable/",
        edit_link = "master",
    ),
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo = "github.com/SciML/ODEInterfaceDiffEq.jl",
    devbranch = "master",
)
