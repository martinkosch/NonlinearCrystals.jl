using NonlinearCrystals
using Documenter

DocMeta.setdocmeta!(NonlinearCrystals, :DocTestSetup, :(using NonlinearCrystals); recursive=true)

makedocs(;
    modules=[NonlinearCrystals],
    authors="Martin Kosch <martin.kosch@gmail.com> and contributors",
    sitename="NonlinearCrystals.jl",
    checkdocs=:none,
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://martinkosch.github.io/NonlinearCrystals.jl",
        edit_link="main",
        assets=String[],
        repolink="https://github.com/martinkosch/NonlinearCrystals.jl",
        size_threshold = 1000 * 2^10,
    ),
    pages=[
        "Introduction" => "index.md",
        "API" => [
            "Refractive Index" => "refractive_index.md",
            "Phasematching" => "phasematching.md",
            "Utils" => "utils.md",
            ]
    ],
)

deploydocs(;
    repo="github.com/martinkosch/NonlinearCrystals.jl",
    devbranch="main",
)
