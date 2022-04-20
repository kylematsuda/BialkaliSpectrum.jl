```@meta
CurrentModule = MoleculeSpectrum
```

# Plotting
```@index
Pages = ["plotting.md"]
```

Methods for plotting data from [`calculate_spectrum`](@ref) and
[`calculate_spectra_vs_fields`](@ref) using
[`CairoMakie`](https://makie.juliaplots.org/stable/documentation/backends/cairomakie/).

In a REPL session, entering `using ElectronDisplay` will plot the figures in a
pop-up window.

## Methods
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["plotting.jl"]
Private = false
Order = [:function]
```