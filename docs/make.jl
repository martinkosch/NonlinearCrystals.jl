using NonlinearCrystals
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(NonlinearCrystals, :DocTestSetup, :(using NonlinearCrystals); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "src", "bibliography.bibtex"))

makedocs(;
    modules=[NonlinearCrystals],
    authors="Martin Kosch <martin.kosch@gmail.com> and contributors",
    sitename="NonlinearCrystals.jl",
    checkdocs=:none,
    plugins=[bib],
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://martinkosch.github.io/NonlinearCrystals.jl",
        edit_link="main",
        assets=String["assets/citations.css"],
        repolink="https://github.com/martinkosch/NonlinearCrystals.jl",
        size_threshold=1000 * 2^10,
    ),
    pages=[
        "Introduction" => "index.md",
        "Examples" => [
            "Visualizing refractive indices" => "examples/vis_refr_indices.md",
            "Adding new nonlinear crystals" => "examples/adding_new_crystal.md",
        ],
        "API" => [
            "Refractive Index" => "refractive_index.md",
            "Phasematching" => "phasematching.md",
            "Utils" => "utils.md",
        ],
        "Bibliography" => "bibliography.md",
    ],
)

deploydocs(;
    repo="github.com/martinkosch/NonlinearCrystals.jl",
    devbranch="main",
)
