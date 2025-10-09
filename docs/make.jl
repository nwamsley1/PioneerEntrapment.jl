using Documenter
using PioneerEntrapment

DocMeta.setdocmeta!(PioneerEntrapment, :DocTestSetup, :(using PioneerEntrapment); recursive=true)

makedocs(
    modules = [PioneerEntrapment],
    sitename = "PioneerEntrapment.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    pages = [
        "Home" => "index.md",
        "Entrapment EFDR" => "entrapment_computation.md",
    ],
)

deploydocs(
    repo = "github.com/nwamsley1/PioneerEntrapment.jl.git",
    push_preview = true,
)
