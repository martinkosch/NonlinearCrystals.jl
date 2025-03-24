using NonlinearCrystals
using Documenter

DocMeta.setdocmeta!(NonlinearCrystals, :DocTestSetup, :(using NonlinearCrystals); recursive=true)

makedocs(;
    modules=[NonlinearCrystals],
    authors="Martin Kosch <martin.kosch@gmail.com> and contributors",
    sitename="NonlinearCrystals.jl",
    format=Documenter.HTML(;
        canonical="https://martinkosch.github.io/NonlinearCrystals.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo = "github.com/martinkosch/NonlinearCrystals.jl",
    devbranch="main",
)
