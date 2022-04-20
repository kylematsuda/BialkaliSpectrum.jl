```@meta
CurrentModule = MoleculeSpectrum
```

# State
```@index
Pages = ["state.md"]
```

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

## Types
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["state.jl"]
Private = false
Order = [:type]
```

## Standard library interfaces
```@docs
convert
show
```

## Methods
```@autodocs
Modules = [MoleculeSpectrum]
Pages = ["state.jl"]
Private = false
Order = [:function]
```