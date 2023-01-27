using SafeSelfTriggered
using Documenter

DocMeta.setdocmeta!(SafeSelfTriggered, :DocTestSetup, :(using SafeSelfTriggered); recursive=true)

makedocs(;
    modules=[SafeSelfTriggered],
    authors="Arvind Adimoolam <asarvind.adimoolam@gmail.com>",
    repo="https://github.com/asarvind/SafeSelfTriggered.jl/blob/{commit}{path}#{line}",
    sitename="SafeSelfTriggered.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://asarvind.github.io/SafeSelfTriggered.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/asarvind/SafeSelfTriggered.jl",
    devbranch="main",
)
