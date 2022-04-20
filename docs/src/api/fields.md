```@meta
CurrentModule = MoleculeSpectrum
```

# Fields
```@index
Pages = ["fields.md"]
```

Utilities for defining the external magnetic, electric, and optical fields
experienced by the molecules.

The [`ExternalFields`](@ref) holds the external fields, which are each
defined by a [`SphericalVector`](@ref) with a `magnitude`, polar angle
``\theta``, and azimuthal angle ``\phi``. Note that the molecular basis
is defined with respect to the ``z``-axis (``\theta = 0``).

A `Vector{ExternalFields}` is expected by [`get_spectra`](@ref),
and can be generated manually or with [`generate_fields_scan`](@ref).

## Types
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["fields.jl"]
Private = false
Order = [:type]
```

## Methods
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["fields.jl"]
Private = false
Order = [:function]
```