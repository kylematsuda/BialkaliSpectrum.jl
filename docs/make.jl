# push!(LOAD_PATH, "../src/")

using Documenter, BialkaliSpectrum

DocMeta.setdocmeta!(BialkaliSpectrum,
    :DocTestSetup,
    :(using BialkaliSpectrum);
    recursive=true
)

makedocs(
    sitename = "BialkaliSpectrum documentation",
    modules = [BialkaliSpectrum],
    doctest = true,
    clean = true,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = Any[
        "Introduction" => "index.md",
        "Test" => "man/basics.md",
        "API" => Any[
            "Setting up a calculation" => "api/setting_up.md",
            "Analyzing the results" => "api/analyzing_spectrum.md"
        ]
    ],
    # strict = true
)

deploydocs(
    repo = "github.com/kylematsuda/BialkaliSpectrum.jl.git",
    devbranch = "main",
    push_preview = true,
)


