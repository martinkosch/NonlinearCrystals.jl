using NonlinearCrystals
using Documenter

DocMeta.setdocmeta!(NonlinearCrystals, :DocTestSetup, :(using NonlinearCrystals); recursive=true)

makedocs(;
    modules=[NonlinearCrystals],
    authors="Martin Kosch <martin.kosch@gmail.com> and contributors",
    sitename="NonlinearCrystals.jl",
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical="https://martinkosch.github.io/NonlinearCrystals.jl",
        edit_link="main",
        assets=String[],
        repolink = "https://github.com/martinkosch/NonlinearCrystals.jl"
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo = "github.com/martinkosch/NonlinearCrystals.jl",
    devbranch="main",
)
