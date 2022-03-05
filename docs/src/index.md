# MoleculeSpectrum.jl

Computes the energy levels of diatomic ${}^{1} \Sigma$ molecules in magnetic, electric, and optical fields.

# Todo
- [x] Add transition strengths
- [ ] Dump to file, and/or dataframes?
- [ ] Graphing
- [x] Improve code organization
- [ ] Documentation
    - [ ] Make it build!
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