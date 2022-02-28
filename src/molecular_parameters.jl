struct ZeemanParameters
    "Rotational g factor"
    gᵣ::Float64
    "Nuclear g factor"
    gᵢ::SVector{2,Float64}
    "Nuclear shielding factor"
    σᵢ::SVector{2,Float64}
end

struct NuclearParameters
    "Nuclear electric quadrupole (MHz)"
    eqQᵢ::SVector{2,Float64}
    "Nuclear spin-rotation interaction (MHz)"
    cᵢ::SVector{2,Float64}
    "Nuclear spin-spin scalar interaction (MHz)"
    c₄::Float64
end

struct Polarizability
    "Parallel ac polarizability (MHz / (W / cm^2))"
    α_par::Float64
    "Perpendicular ac polarizability (MHz / (W / cm^2))"
    α_perp::Float64
end

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

const KRb_Zeeman = ZeemanParameters(0.014, [-0.324, 1.834], [1321e-6, 3469e-6])
const KRb_Polarizability = Polarizability(10.0e-5, 3.3e-5)

const KRb_Nuclear_Neyenhuis =
    NuclearParameters([0.45, -1.308], [-24.1e-6, 420.1e-6], -2030.4e-6)
const KRb_Nuclear_Ospelkaus =
    NuclearParameters([0.45, -1.41], [-24.1e-6, 420.1e-6], -2030.4e-6)

const KRb_Parameters_Neyenhuis = MolecularParameters(
    0.574,
    1113.9514,
    [HalfInt(4), HalfInt(3 / 2)],
    KRb_Zeeman,
    KRb_Nuclear_Neyenhuis,
    KRb_Polarizability,
)
const KRb_Parameters_Ospelkaus = MolecularParameters(
    0.574,
    1113.950,
    [HalfInt(4), HalfInt(3 / 2)],
    KRb_Zeeman,
    KRb_Nuclear_Ospelkaus,
    KRb_Polarizability,
)

const DEFAULT_MOLECULAR_PARAMETERS = KRb_Parameters_Neyenhuis

const TOY_MOLECULE_PARAMETERS = MolecularParameters(
    1.0,
    1000.0,
    [HalfInt(1), HalfInt(1)], # Some of the matrix elements don't make sense for I = 0
    ZeemanParameters(0.0, [0.0, 0.0], [0.0, 0.0]),
    NuclearParameters([0.0, 0.0], [0.0, 0.0], 0.0),
    Polarizability(0.0, 0.0),
)
