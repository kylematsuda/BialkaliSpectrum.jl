```@meta
CurrentModule = MoleculeSpectrum
```

# DataFrame

```@index
Pages = ["dataframe.md"]
```

Utilities for working with the output of [`calculate_spectrum`](@ref) and
[`calculate_spectra_vs_fields`](@ref), both of which return a `DataFrame`.

These are mainly simple filters or transforms, defined for convenience.
Anything more complicated should use the methods in the
[`DataFrames.jl`](https://dataframes.juliadata.org/stable/) library.

## Methods
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["dataframe.jl"]
Private = false
Order = [:function]
```