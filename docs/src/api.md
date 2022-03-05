# MoleculeSpectrum.jl API Documentation

```@meta
    CurrentModule = MoleculeSpectrum
    DocTestSetup  = quote
        using MoleculeSpectrum
    end
```

## Contents
```@contents
    Pages = ["api.md"]
    Depth = 4
```

## Molecular Parameters
### Molecular Parameters
```@docs
MolecularParameters
KRb_Parameters_Neyenhuis
KRb_Parameters_Ospelkaus
DEFAULT_MOLECULAR_PARAMETERS
TOY_MOLECULE_PARAMETERS
```
### Polarizability Parameters
```@docs
Polarizability
KRb_Polarizability
```
### Zeeman Parameters
```@docs
ZeemanParameters
KRb_Zeeman
```
### Nuclear Parameters
```@docs
NuclearParameters
KRb_Nuclear_Neyenhuis
KRb_Nuclear_Ospelkaus
```

## Vectors and Tensors
### Spherical vectors
```@docs
SphericalVector
VectorX
VectorY
VectorZ
```
### Unit vectors
```@docs
SphericalUnitVector
UnitVectorX
UnitVectorY
UnitVectorZ
```
### Tensors
```@docs
T⁽¹⁾
T⁽²⁾
get_tensor_component
tensor_dot
```

## States
### Types
```@docs
State
KRbState
```
### Functions
```@docs
index_to_state
state_to_index
order_by_overlap_with
max_overlap_with
find_closest_basis_state
decompose_to_basis_states
```

## Fields
```@docs
ExternalFields
DEFAULT_FIELDS
TEST_FIELDS
```

## Hamiltonians
```@docs
hamiltonian
HamiltonianParts
make_hamiltonian_parts
make_krb_hamiltonian_parts
```

## Spectra and Analysis
```@docs
Spectrum
calculate_spectrum
get_energy
get_energy_difference
find_transition_strengths
plot_transition_strengths
calculate_dipolar_interaction
```

## Index
```@index
```