# push!(LOAD_PATH, "../src/")

using Documenter, MoleculeSpectrum

DocMeta.setdocmeta!(MoleculeSpectrum, :DocTestSetup, :(using MoleculeSpectrum); recursive=true)

makedocs(
    sitename = "MoleculeSpectrum documentation",
    modules = [MoleculeSpectrum],
    doctest = true,
    clean = true,
    format = Documenter.HTML(prettyurls = false),
    pages = Any[
        "Introduction" => "index.md",
        "Test" => "man/basics.md",
        "API" => Any[
            "Indexing molecular states" => "api/state.md",
        ]
    ],
    # strict = true
)
