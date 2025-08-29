using Documenter
using PioneerEntrapment

DocMeta.setdocmeta!(PioneerEntrapment, :DocTestSetup, :(using PioneerEntrapment); recursive=true)

makedocs(
    modules = [PioneerEntrapment],
    sitename = "PioneerEntrapment.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/nathanwamsley/PioneerEntrapment.jl.git",
    push_preview = true,
)
