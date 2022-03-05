"""
```
struct ZeemanParameters
    "Rotational g factor"
    gᵣ::Float64
    "Nuclear g factor"
    gᵢ::SVector{2,Float64}
    "Nuclear shielding factor"
    σᵢ::SVector{2,Float64}
end
```

Contains the g-factors and nuclear shielding factors
for computing Zeeman shifts.
"""
struct ZeemanParameters
    "Rotational g factor\n"
    gᵣ::Float64
    "Nuclear g factor\n"
    gᵢ::SVector{2,Float64}
    "Nuclear shielding factor\n"
    σᵢ::SVector{2,Float64}
end

"""
```
struct NuclearParameters
    "Nuclear electric quadrupole (MHz)"
    eqQᵢ::SVector{2,Float64}
    "Nuclear spin-rotation interaction (MHz)"
    cᵢ::SVector{2,Float64}
    "Nuclear spin-spin scalar interaction (MHz)"
    c₄::Float64
end
```

Contains the nuclear electric quadrupole moments,
nuclear spin-rotation couplings, and nuclear
spin-spin coupling for calculating the hyperfine
Hamiltonian.
"""
struct NuclearParameters
    "Nuclear electric quadrupole (MHz)\n"
    eqQᵢ::SVector{2,Float64}
    "Nuclear spin-rotation interaction (MHz)\n"
    cᵢ::SVector{2,Float64}
    "Nuclear spin-spin scalar interaction (MHz)\n"
    c₄::Float64
end

"""
```
struct Polarizability
    "Parallel ac polarizability (MHz / (W / cm^2))"
    α_par::Float64
    "Perpendicular ac polarizability (MHz / (W / cm^2))"
    α_perp::Float64
end
```

Contains the parallel and perpendicular ac
polarizabilities at a particular optical wavelength λ.

So far, it is assumed that λ is far-detuned from
any electronic transitions, such that the polarizability
does not depend on rotational state ``N``.
"""
struct Polarizability
    "Parallel ac polarizability (MHz / (W / cm^2))\n"
    α_par::Float64
    "Perpendicular ac polarizability (MHz / (W / cm^2))\n"
    α_perp::Float64
end

"""
```
struct MolecularParameters
    "Permanent dipole moment (Debye)"
    dₚ::Float64
    "Rotational constant (MHz)"
    Bᵣ::Float64
    "Nuclear angular momenta"
    I::SVector{2,HalfInt}
    "Zeeman parameters"
    zeeman::ZeemanParameters
    "Nuclear Parameters"
    nuclear::NuclearParameters
    "Molecular polarizability at the trapping wavelength"
    α::Polarizability
end
```

Contains the coupling constants for the molecular
Hamiltonian and the nuclear angular momenta (needed to 
construct the basis states).
"""
struct MolecularParameters
    "Permanent dipole moment (Debye)\n"
    dₚ::Float64
    "Rotational constant (MHz)\n"
    Bᵣ::Float64
    "Nuclear angular momenta\n"
    I::SVector{2,HalfInt}
    "Zeeman parameters\n"
    zeeman::ZeemanParameters
    "Nuclear Parameters\n"
    nuclear::NuclearParameters
    "Molecular polarizability at the trapping wavelength\n"
    α::Polarizability
end

"
    KRb_Zeeman = ZeemanParameters(0.014, [-0.324, 1.834], [1321e-6, 3469e-6])

[`ZeemanParameters`](@ref) with theoretical values from [Aldegunde et al., PRA 78, 033434 (2008)](https://doi.org/10.1103/PhysRevA.78.033434)\n"
const KRb_Zeeman = ZeemanParameters(0.014, [-0.324, 1.834], [1321e-6, 3469e-6])

"
    KRb_Polarizability = Polarizability(10.0e-5, 3.3e-5)

[`Polarizability`](@ref) with experimental values from [Neyenhuis et al., PRL 109, 230403 (2012)](https://doi.org/10.1103/PhysRevLett.109.230403)\n"
const KRb_Polarizability = Polarizability(10.0e-5, 3.3e-5)

"
    KRb_Nuclear_Neyenhuis = NuclearParameters([0.45, -1.308], [-24.1e-6, 420.1e-6], -2030.4e-6)

[`NuclearParameters`](@ref) with experimental values from [Neyenhuis et al., PRL 109, 230403 (2012)](https://doi.org/10.1103/PhysRevLett.109.230403)\n"
const KRb_Nuclear_Neyenhuis =
    NuclearParameters([0.45, -1.308], [-24.1e-6, 420.1e-6], -2030.4e-6)

"
    KRb_Nuclear_Ospelkaus = NuclearParameters([0.45, -1.41], [-24.1e-6, 420.1e-6], -2030.4e-6)

[`NuclearParameters`](@ref) with experimental values from [Ospelkaus et al., PRL 104, 030402 (2010)](https://doi.org/10.1103/PhysRevLett.104.030402)\n"
const KRb_Nuclear_Ospelkaus =
    NuclearParameters([0.45, -1.41], [-24.1e-6, 420.1e-6], -2030.4e-6)

"
    KRb_Parameters_Neyenhuis = MolecularParameters(
        0.574,
        1113.9514,
        [HalfInt(4), HalfInt(3 / 2)],
        KRb_Zeeman,
        KRb_Nuclear_Neyenhuis,
        KRb_Polarizability,
    )

[`MolecularParameters`](@ref) with experimental values from [Neyenhuis et al., PRL 109, 230403 (2012)](https://doi.org/10.1103/PhysRevLett.109.230403)\n"
const KRb_Parameters_Neyenhuis = MolecularParameters(
    0.574,
    1113.9514,
    [HalfInt(4), HalfInt(3 / 2)],
    KRb_Zeeman,
    KRb_Nuclear_Neyenhuis,
    KRb_Polarizability,
)

"
    KRb_Parameters_Ospelkaus = MolecularParameters(
        0.574,
        1113.950,
        [HalfInt(4), HalfInt(3 / 2)],
        KRb_Zeeman,
        KRb_Nuclear_Ospelkaus,
        KRb_Polarizability,
    )

[`MolecularParameters`](@ref) with experimental values from [Ospelkaus et al., PRL 104, 030402 (2010)](https://doi.org/10.1103/PhysRevLett.104.030402)\n"
const KRb_Parameters_Ospelkaus = MolecularParameters(
    0.574,
    1113.950,
    [HalfInt(4), HalfInt(3 / 2)],
    KRb_Zeeman,
    KRb_Nuclear_Ospelkaus,
    KRb_Polarizability,
)

"""
    DEFAULT_MOLECULAR_PARAMETERS

Default molecular parameters; alias for [`KRb_Parameters_Neyenhuis`](@ref)
"""
const DEFAULT_MOLECULAR_PARAMETERS = KRb_Parameters_Neyenhuis

"""
    TOY_MOLECULE_PARAMETERS

Toy model values with dipole = 1 D, no hyperfine structure, etc.
Intended for testing.

Note that the formulas break down for `I = 0`, which is why we use
`I = 1` here.
"""
const TOY_MOLECULE_PARAMETERS = MolecularParameters(
    1.0,
    1000.0,
    [HalfInt(1), HalfInt(1)], # Some of the matrix elements don't make sense for I = 0
    ZeemanParameters(0.0, [0.0, 0.0], [0.0, 0.0]),
    NuclearParameters([0.0, 0.0], [0.0, 0.0], 0.0),
    Polarizability(0.0, 0.0),
)
