```@meta
CurrentModule = MoleculeSpectrum
```

# Analyzing the spectrum

## Contents

```@contents
Pages = ["analyzing_spectrum.md"]
```

## Index

```@index
Pages = ["analyzing_spectrum.md"]
```

## Spectrum

These are the main functions for calculating the molecular spectrum.

In most cases, you should just call [`get_spectra`](@ref) to calculate
the spectrum at each point in a `field_scan::Vector{ExternalFields}`.

A common pattern is to put the resulting `spectra` into
[`transform_spectra`](@ref) to calculate quantities of interest at each
point in the `field_scan`. [`transform_spectra`](@ref) groups `spectra`
(the default is by `:fields`), then applies a transformation `f` to each
subgroup. See the source for [`transitions`](@ref) for an example of this.

[`get_spectrum`](@ref) can be called to calculate the spectrum at a single
value of `external_fields`.

TODO: get rid of `get_spectrum`?? When is it ever helpful? Then the utility functions
might need to be extended to work on grouped dataframes as well?

### Methods
```@docs
get_spectrum
find_closest_eigenstate
get_energy
get_energy_difference
get_spectra
transform_spectra
```

## DataFrame helpers

Utilities for working with the output of [`get_spectrum`](@ref) and
[`get_spectra`](@ref), both of which return a `DataFrame`.

These are mainly simple filters or transforms, defined for convenience.
Anything more complicated should use the methods in the
[`DataFrames.jl`](https://dataframes.juliadata.org/stable/) library.

### Methods
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["dataframe.jl"]
Private = false
Order = [:function]
```

## Analysis

Functions for analyzing the output of [`get_spectrum`](@ref) or
[`get_spectra`](@ref).

### Methods
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["analysis.jl"]
Private = false
Order = [:function]
```

## Plotting

Methods for plotting data from [`get_spectrum`](@ref) and
[`get_spectra`](@ref) using
[`CairoMakie`](https://makie.juliaplots.org/stable/documentation/backends/cairomakie/).

In a REPL session, entering `using ElectronDisplay` will plot the figures in a
pop-up window.

### Methods
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["plotting.jl"]
Private = false
Order = [:function]
```