module MoleculeSpectrum

import WignerSymbols: wigner3j
import HalfIntegers: HalfInt
import Gadfly
using LinearAlgebra, SparseArrays, StaticArrays, Test

export ZeemanParameters, NuclearParameters, Polarizability, MolecularParameters
export KRb_Zeeman, KRb_Nuclear_Neyenhuis, KRb_Nuclear_Ospelkaus, KRb_Polarizability
export KRb_Parameters_Neyenhuis, KRb_Parameters_Ospelkaus, DEFAULT_MOLECULAR_PARAMETERS, TOY_MOLECULE_PARAMETERS

export SphericalVector, VectorX, VectorY, VectorZ
export SphericalUnitVector, UnitVectorX, UnitVectorY, UnitVectorZ
export T⁽¹⁾, T⁽²⁾, get_tensor_component, tensor_dot
export ExternalFields, DEFAULT_FIELDS, TEST_FIELDS

export State, KRbState, index_to_state, state_to_index
export order_by_overlap_with, max_overlap_with, get_energy, get_energy_difference

export HamiltonianParts, make_hamiltonian_parts, hamiltonian, make_krb_hamiltonian_parts

export Spectrum, calculate_spectrum, transition_strengths, plot_transition_strengths

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

KRbState(N, mₙ, mᵢ₁::Number, mᵢ₂::Number) = State(N, mₙ, KRb_Parameters_Neyenhuis.I, [HalfInt(mᵢ₁) HalfInt(mᵢ₂)])

struct Spectrum
    hamiltonian_parts
    energies
    eigenstates
end

include("utility.jl")
include("matrix_elements.jl")
include("hamiltonian.jl")

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

    (overlap, index_g) = max_overlap_with(spectrum, g)
    if overlap < 0.5
        @warn "The best overlap with your requested ground state is < 0.5."
    end
    E_g = energies[index_g]
    g_state = eigenstates[:, index_g]

    state_range = searchsortedfirst(energies, f_low + E_g):searchsortedlast(energies, f_high + E_g)
    states = [eigenstates[:, k] for k in state_range]
    frequencies = [energies[k] - E_g for k in state_range]

    T1ϵ = T⁽¹⁾(polarization)
    # Normalize to d/sqrt(3), which is the largest transition dipole (between |0,0> and |1,0>)
    h_dipole = tensor_dot(T1ϵ, spectrum.hamiltonian_parts.dipole_relative) / (1/sqrt(3))

    strengths = [abs(g_state' * h_dipole * e) for e in states]
    closest_basis_states = map(e -> find_closest_basis_state(spectrum, e), state_range)
    
    out = [x for x in zip(frequencies, strengths, closest_basis_states)]
    return sort!(out, by=t->t[2], rev=true)
end

function plot_transition_strengths(spectrum::Spectrum, g::State, f_low, f_high; polarization::SphericalUnitVector = Unpolarized())
    transitions = transition_strengths(spectrum, g, f_low, f_high; polarization=polarization)

    freqs = [t[1] for t in transitions]
    strengths = [t[2] for t in transitions]

    Gadfly.plot(
        x=freqs,
        y=strengths,
        Gadfly.Geom.hair,
        Gadfly.Geom.point,
        Gadfly.Guide.xlabel("Frequency (MHz)"),
        Gadfly.Guide.ylabel("Transition strength")
    )
end

end # module
