push!(LOAD_PATH,"../src/")

using Documenter, MoleculeSpectrum

makedocs(
    sitename="MoleculeSpectrum Documentation",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "MoleculeSpectrum" => "index.md",
        "API Documentation" => "api.md"
    ]
)