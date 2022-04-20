```@meta
CurrentModule = MoleculeSpectrum
```

# Hamiltonian

```@index
Pages = ["hamiltonian.md"]
```

Functions for building the molecular Hamiltonian.

The [`HamiltonianParts`](@ref) struct is used as an input to several analysis
methods, e.g. computing dipole matrix elements, since it holds the individual
terms in the Hamiltonian separately.

## Types
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["hamiltonian.jl"]
Private = false
Order = [:type]
```

## Methods
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["hamiltonian.jl"]
Private = false
Order = [:function]
```