module MoleculeSpectrum

import WignerSymbols: wigner3j
import HalfIntegers: HalfInt
using LinearAlgebra, SparseArrays, StaticArrays, Test

export ZeemanParameters, NuclearParameters, Polarizability, MolecularParameters
export KRb_Zeeman, KRb_Nuclear_Neyenhuis, KRb_Nuclear_Ospelkaus, KRb_Polarizability
export KRb_Parameters_Neyenhuis, KRb_Parameters_Ospelkaus, DEFAULT_MOLECULAR_PARAMETERS

export SphericalVector, VectorX, VectorY, VectorZ
export SphericalUnitVector, UnitVectorX, UnitVectorY, UnitVectorZ
export T⁽¹⁾, T⁽²⁾, get_tensor_component, tensor_dot
export ExternalFields, DEFAULT_FIELDS, TEST_FIELDS

export State, index_to_state, state_to_index
export order_by_overlap_with, max_overlap_with, get_energy, get_energy_difference

export HamiltonianParts, make_hamiltonian_parts, hamiltonian, make_krb_hamiltonian_parts

export calculate_spectrum, transition_strengths

module Constants
    const μN = 7.622593285e-4 # MHz/G
    const DToSI = 3.33564e-30
    const h = 6.62607004e-34
    const DVcm⁻¹ToMHz = (DToSI / h) * 1e-4
end # module

include("molecular_parameters.jl")
include("fields.jl")

struct State
    N::Int
    mₙ::Int
    I::SVector{2, HalfInt} # [K, Rb]
    mᵢ::SVector{2, HalfInt}
end

State(N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂) = State(N, mₙ, SVector(I₁, I₂), SVector(mᵢ₁, mᵢ₂))
State(N, mₙ, mᵢ₁::Number, mᵢ₂::Number) = State(N, mₙ, DEFAULT_MOLECULAR_PARAMETERS.I, [HalfInt(mᵢ₁) HalfInt(mᵢ₂)])

include("utility.jl")
include("matrix_elements.jl")
include("hamiltonian.jl")

struct Spectrum
    hamiltonian_parts
    energies
    eigenstates
end

function make_krb_hamiltonian_parts(N_max::Int)
    return make_hamiltonian_parts(KRb_Parameters_Neyenhuis, N_max)
end

function calculate_spectrum(hamiltonian_parts::HamiltonianParts, external_fields::ExternalFields)::Spectrum
    h = hamiltonian(hamiltonian_parts, external_fields)
    energies = eigvals(h)
    eigenstates = eigvecs(h)
    return Spectrum(hamiltonian_parts, energies, eigenstates)
end

function transition_strengths(spectrum::Spectrum, g::State, f_low, f_high; polarization::SphericalUnitVector = Unpolarized())
    eigenstates = spectrum.eigenstates
    energies = spectrum.energies

    (overlap, index_g) = max_overlap_with(g, eigenstates)
    if overlap < 0.5
        @warn "The best overlap with your requested ground state is < 0.5."
    end
    E_g = energies[index_g]
    state_range = (searchsortedfirst(energies, f_low + E_g), searchsortedlast(energies, f_high + E_g))

    T1ϵ = T⁽¹⁾(polarization)
    h_dipole = tensor_dot(T1ϵ, spectrum.hamiltonian_parts.dipole_relative) / sqrt(3)

    transitions = [(energies[k] - E_g, abs(eigenstates[:, index_g]' * h_dipole * eigenstates[:, k])) for k in state_range[1]:state_range[2]]
    return sort!(transitions, by=t->t[2], rev=true)
end

end # module
