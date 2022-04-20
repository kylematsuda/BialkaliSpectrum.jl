```@meta
CurrentModule = MoleculeSpectrum
```

# Spectrum
```@index
Pages = ["spectrum.md"]
```

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

## Methods
```@docs
get_spectrum
find_closest_eigenstate
get_energy
get_energy_difference
get_spectra
transform_spectra
```