module K40Rb87

import HalfIntegers: HalfInt
import ..BialkaliSpectrum: ZeemanParameters, Polarizability, NuclearParameters,
    MolecularParameters, make_hamiltonian_parts, State

export KRb_Zeeman, KRb_Nuclear_Neyenhuis, KRb_Nuclear_Ospelkaus, KRb_Polarizability
export KRb_Parameters_Neyenhuis, KRb_Parameters_Ospelkaus
export KRbState, make_krb_hamiltonian_parts

"""
    KRb_Zeeman

Theoretical values from Aldegunde et al. PRA (2008)
"""
const KRb_Zeeman = ZeemanParameters(0.014, [-0.324, 1.834], [1321e-6, 3469e-6])

"""
    KRb_Polarizability

Experimental values from Neyenhuis et al., PRL 109, 230403 (2012)
"""
const KRb_Polarizability = Polarizability(10.0e-5, 3.3e-5)

"""
    KRb_Nuclear_Neyenhuis

Experimental values from Neyenhuis et al., PRL 109, 230403 (2012)
"""
const KRb_Nuclear_Neyenhuis =
    NuclearParameters([0.45, -1.308], [-24.1e-6, 420.1e-6], -2030.4e-6)

"""
    KRb_Nuclear_Ospelkaus

Experimental values from Ospelkaus et al., PRL 104, 030402 (2010)
"""
const KRb_Nuclear_Ospelkaus =
    NuclearParameters([0.45, -1.41], [-24.1e-6, 420.1e-6], -2030.4e-6)

"""
    KRb_Parameters_Neyenhuis

Experimental values from Neyenhuis et al., PRL 109, 230403 (2012)
"""
const KRb_Parameters_Neyenhuis = MolecularParameters(
    0.574,
    1113.9514,
    [HalfInt(4), HalfInt(3 / 2)],
    KRb_Zeeman,
    KRb_Nuclear_Neyenhuis,
    KRb_Polarizability,
)

"""
    KRb_Parameters_Ospelkaus

Experimental values from Ospelkaus et al., PRL 104, 030402 (2010)
"""
const KRb_Parameters_Ospelkaus = MolecularParameters(
    0.574,
    1113.950,
    [HalfInt(4), HalfInt(3 / 2)],
    KRb_Zeeman,
    KRb_Nuclear_Ospelkaus,
    KRb_Polarizability,
)

"""
    make_krb_hamiltonian_parts(N_max::Int)

Construct all parts of the ``{}^{40}\\text{K}^{87}\\text{Rb}`` Hamiltonian
that do not depend on external fields.

The rotational states `0:N_max` are included. This is a shortcut method that
replaces [`make_hamiltonian_parts`](@ref) for KRb.

See also [`make_hamiltonian_parts`](@ref).
"""
function make_krb_hamiltonian_parts(N_max::Int)
    return make_hamiltonian_parts(KRb_Parameters_Neyenhuis, N_max)
end

"""
    KRbState(N, mₙ, mK, mRb)

Creates a basis state ``|N, m_n, m_{\\text{K}}, m_{\\text{Rb}}⟩`` for ``{}^{40}\\text{K}^{87}\\text{Rb}``.

This is a wrapper around [`State`](@ref) to avoid having to specify the nuclear spins ``I_k`` each time.

See also [`State`](@ref).
"""
KRbState(N, mₙ, mK, mRb) =
    State(N, mₙ, KRb_Parameters_Neyenhuis.I, [HalfInt(mK) HalfInt(mRb)])

end # module
using .K40Rb87

module Toy

import HalfIntegers: HalfInt
import ..BialkaliSpectrum: ZeemanParameters, Polarizability, NuclearParameters,
    MolecularParameters, State,
    generate_basis, SparseHamiltonian, HamiltonianParts,
    h_rotation, h_dipole,
    h_diagonal, h_rank_1, h_rank_2,
    make_hamiltonian_parts

export make_toy_hamiltonian_parts, State

"""
    TOY_PARAMETERS

Toy model values with dipole = 1 D and no hyperfine structure.
Intended for testing.
"""
const TOY_PARAMETERS = MolecularParameters(
    1.0, # d_p
    1.0, # B_r
    [HalfInt(0), HalfInt(0)],
    ZeemanParameters(0.0, [0.0, 0.0], [0.0, 0.0]),
    NuclearParameters([0.0, 0.0], [0.0, 0.0], 0.0),
    Polarizability(0.0, 0.0),
)

"""
    make_toy_hamiltonian_parts(N_max::Int)

Construct all parts of the toy molecule Hamiltonian
that do not depend on external fields.

The rotational states `0:N_max` are included. This is a shortcut method that
replaces [`make_hamiltonian_parts`](@ref) for the toy molecule.

See also [`make_hamiltonian_parts`](@ref).
"""
function make_toy_hamiltonian_parts(N_max::Int)::HamiltonianParts
    basis = generate_basis(TOY_PARAMETERS, N_max)

    rotation = TOY_PARAMETERS.Bᵣ * h_rotation(basis)
    dipole_relative = h_dipole(basis)
    dipole = (-1) * TOY_PARAMETERS.dₚ * h_dipole(basis)

    return HamiltonianParts(
        basis,
        rotation,
        dipole,
        dipole_relative,
        h_diagonal(basis, (_, _) -> 0),
        h_rank_1(basis, (_, _, _) -> 0),
        h_diagonal(basis, (_, _) -> 0),
        h_rank_2(basis, (_, _, _) -> 0),
    )
end

"""
    State(N, mₙ)

Creates a basis state ``|N, m_n⟩`` for the toy molecule.

This is a wrapper around [`State`](@ref) to avoid having to specify the nuclear spins ``I_k`` each time.

See also [`State`](@ref).
"""
State(N, m_n) = State(N, m_n, [0, 0], [0, 0])

end # module