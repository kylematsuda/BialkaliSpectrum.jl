```@meta
CurrentModule = MoleculeSpectrum
```

# Setting up the calculation

## Contents

```@contents
Pages = ["setting_up.md"]
```

## Index

```@index
Pages = ["setting_up.md"]
```

## Molecular parameters

Structs for organizing the coupling constants, dipole moment, and nuclear spins
of the molecule.

### Types
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["molecular_parameters.jl"]
Private = false
Order = [:type]
```

### Constants
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["molecular_parameters.jl"]
Private = false
Order = [:constant]
```

## Hamiltonian

Functions for building the molecular Hamiltonian.

The [`HamiltonianParts`](@ref) struct is used as an input to several analysis
methods, e.g. computing dipole matrix elements, since it holds the individual
terms in the Hamiltonian separately.

### Types
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["hamiltonian.jl"]
Private = false
Order = [:type]
```

### Methods
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["hamiltonian.jl"]
Private = false
Order = [:function]
```

## Fields

Utilities for defining the external magnetic, electric, and optical fields
experienced by the molecules.

The [`ExternalFields`](@ref) holds the external fields, which are each
defined by a [`SphericalVector`](@ref) with a `magnitude`, polar angle
``\theta``, and azimuthal angle ``\phi``. Note that the molecular basis
is defined with respect to the ``z``-axis (``\theta = 0``).

A `Vector{ExternalFields}` is expected by [`get_spectra`](@ref),
and can be generated manually or with [`generate_fields_scan`](@ref).

### Types
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["fields.jl"]
Private = false
Order = [:type]
```

### Methods
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["fields.jl"]
Private = false
Order = [:function]
```

## Basis states

Utilities for working with molecular states in the uncoupled basis
``|N, m_N, m_{i1}, m_{i2} \rangle``.

The basis states are indexed with the quantum numbers changing fastest on the right.
For example, with ``I_1 = 1/2`` and ``I_2 = 1/2``, the first few states are:

1. ``|0, 0, -1/2, -1/2⟩``
2. ``|0, 0, -1/2, +1/2⟩``
3. ``|0, 0, +1/2, -1/2⟩``
4. ``|0, 0, +1/2, +1/2⟩``
5. ``|1, -1, -1/2, -1/2⟩``

etc.

The methods [`basis_state`](@ref) and [`basis_index`](@ref) are used to switch between an `index`
and the corresponding basis [`State`](@ref).

### Types
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["state.jl"]
Private = false
Order = [:type]
```

### Standard library interfaces
```@docs
convert
show
```

### Methods
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["state.jl"]
Private = false
Order = [:function]
```

## Molecule-specific definitions

Constants and method defs for working with specific molecule species.

### ``{}^{40}\text{K}^{87}\text{Rb}`` (`K40Rb87`)
```@autodocs
Modules = [MoleculeSpectrum.K40Rb87]
```
