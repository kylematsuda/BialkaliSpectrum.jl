# BialkaliSpectrum.jl

Calculate the energy levels of diatomic
$${}^{1} \Sigma^+$$
molecules in magnetic, electric, and optical fields.

Example calculation: ${}^{40}\text{K}^{87}\text{Rb}$ transitions from a particular hyperfine state as a function of electric field. 
![](krb_transitions_vs_E_bypol.png)

## Manual
```@contents
Pages = ["man/basics.md",
    "man/worked_example.md",
]
Depth = 2
```

## Public API
```@contents
Pages = [
    "api/setting_up.md",
    "api/analyzing_spectrum.md",
]
Depth = 3
```

## Internals
```@contents
Pages = ["internals/constants.md",
    "internals/state.md",
    "internals/matrix_elements.md",
    "internals/hamiltonian.md",
    "internals/fields.md",
    "internals/analysis.md",
]
Depth = 2
```

## Index