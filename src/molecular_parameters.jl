using StaticArrays
import HalfIntegers: HalfInt

"""
    ZeemanParameters

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
    NuclearParameters

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
    Polarizability

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
    MolecularParameters

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

"Theoretical values from Aldegunde et al. PRA (2008)\n" # TODO: Put exact citation
const KRb_Zeeman = ZeemanParameters(0.014, [-0.324, 1.834], [1321e-6, 3469e-6])

"Experimental values from Neyenhuis et al., PRL 109, 230403 (2012)\n"
const KRb_Polarizability = Polarizability(10.0e-5, 3.3e-5)

"Experimental values from Neyenhuis et al., PRL 109, 230403 (2012)\n"
const KRb_Nuclear_Neyenhuis =
    NuclearParameters([0.45, -1.308], [-24.1e-6, 420.1e-6], -2030.4e-6)

"Experimental values from Ospelkaus et al., PRL 104, 030402 (2010)\n"
const KRb_Nuclear_Ospelkaus =
    NuclearParameters([0.45, -1.41], [-24.1e-6, 420.1e-6], -2030.4e-6)

"Experimental values from Neyenhuis et al., PRL 109, 230403 (2012)\n"
const KRb_Parameters_Neyenhuis = MolecularParameters(
    0.574,
    1113.9514,
    [HalfInt(4), HalfInt(3 / 2)],
    KRb_Zeeman,
    KRb_Nuclear_Neyenhuis,
    KRb_Polarizability,
)

"Experimental values from Ospelkaus et al., PRL 104, 030402 (2010)\n"
const KRb_Parameters_Ospelkaus = MolecularParameters(
    0.574,
    1113.950,
    [HalfInt(4), HalfInt(3 / 2)],
    KRb_Zeeman,
    KRb_Nuclear_Ospelkaus,
    KRb_Polarizability,
)

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