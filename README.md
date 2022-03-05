# MoleculeSpectrum.jl

[![CI Actions Status](https://github.com/cal-miller-harvard/MoleculeSpectrum/actions/workflows/ci.yml/badge.svg)](https://github.com/cal-miller-harvard/MoleculeSpectrum/actions)

Compute the energy levels of diatomic ${}^{1} \Sigma$ molecules in magnetic, electric, and optical fields.

Documentation is available [here](https://cal-miller-harvard.github.io/MoleculeSpectrum).

# Todo
- [x] Add transition strengths
- [ ] Dump to file, and/or dataframes?
- [ ] Graphing
- [x] Improve code organization
- [x] Documentation
    - [x] Make it build!
    - [ ] Finish documenting all functions
    - [ ] Document units of all numeric quantities
- [ ] Remove unnecessary types in function arguments (?)
- [ ] API ergonomics
    - [x] ExternalFields
    - [ ] Get dipole matrix elements
    - [ ] More generally... inspect any matrix element?
    - [ ] Make more `calculate` functions polymorphic over `State`, eigenstate index, etc
    - [ ] ...
- [x] Switch to `get_eigenstates` function
- [x] Explicitly type fields of all structs
- [ ] Many more tests!!!