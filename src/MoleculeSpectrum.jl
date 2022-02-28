module MoleculeSpectrum

import WignerSymbols: wigner3j
import HalfIntegers: HalfInt
import Gadfly
using LinearAlgebra, SparseArrays, StaticArrays, Test

export ZeemanParameters, NuclearParameters, Polarizability, MolecularParameters
export KRb_Zeeman, KRb_Nuclear_Neyenhuis, KRb_Nuclear_Ospelkaus, KRb_Polarizability
export KRb_Parameters_Neyenhuis,
    KRb_Parameters_Ospelkaus, DEFAULT_MOLECULAR_PARAMETERS, TOY_MOLECULE_PARAMETERS

export SphericalVector, VectorX, VectorY, VectorZ
export SphericalUnitVector, UnitVectorX, UnitVectorY, UnitVectorZ, Unpolarized
export T⁽¹⁾, T⁽²⁾, get_tensor_component, tensor_dot
export ExternalFields, DEFAULT_FIELDS, TEST_FIELDS

export State, KRbState, index_to_state, state_to_index
export order_by_overlap_with, max_overlap_with, find_closest_basis_state
export get_energy, get_energy_difference

export HamiltonianParts, make_hamiltonian_parts, hamiltonian, make_krb_hamiltonian_parts

export Spectrum, calculate_spectrum, transition_strengths, plot_transition_strengths

module Constants
"Nuclear magneton in MHz/G\n"
const μN = 7.622593285e-4
"Factor to convert from Debye to C m\n"
const DToSI = 3.33564e-30
"Planck's constant (SI)\n"
const h = 6.62607004e-34
"Convert from D*(V/cm) to MHz\n"
const DVcm⁻¹ToMHz = (DToSI / h) * 1e-4
end # module

include("molecular_parameters.jl")
include("fields.jl")

"""
    State

Represents a molecular state in the uncoupled basis.
"""
struct State
    N::Int
    mₙ::Int
    I::SVector{2,HalfInt} # [K, Rb]
    mᵢ::SVector{2,HalfInt}
end

"""
    State(N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂)

Creates a basis state ``|N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂⟩``.

"""
State(N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂) = State(N, mₙ, SVector(HalfInt(I₁), HalfInt(I₂)), SVector(HalfInt(mᵢ₁), HalfInt(mᵢ₂)))

"""
    KRbState(N, mₙ, mK, mRb)

Creates a basis state ``|N, m_n, m_{\\text{K}}, m_{\\text{Rb}}⟩`` for ``{}^{40}\\text{K}^{87}\\text{Rb}``.

This is a wrapper around [`State`](@ref) to avoid having to specify the nuclear spins ``I_k`` each time.

See also [`State`](@ref).
"""
KRbState(N, mₙ, mK, mRb) =
    State(N, mₙ, KRb_Parameters_Neyenhuis.I, [HalfInt(mK) HalfInt(mRb)])

struct Spectrum
    hamiltonian_parts::Any
    energies::Any
    eigenstates::Any
end

include("utility.jl")
include("matrix_elements.jl")
include("hamiltonian.jl")

"""
    calculate_spectrum(hamiltonian_parts, external_fields)

Compute the energies and eigenstates under the external fields.

To avoid reconstructing the Hamiltonian each time, `hamiltonian_parts` can be reused over calls
to `calculate_spectrum`. The output [`Spectrum`](@ref) object is used as an input for
further analysis, for example in [`transition_strengths`](@ref).

See also [`make_hamiltonian_parts`](@ref), [`make_krb_hamiltonian_parts`](@ref),
[`Spectrum`](@ref).
"""
function calculate_spectrum(
    hamiltonian_parts::HamiltonianParts,
    external_fields::ExternalFields,
)::Spectrum
    h = hamiltonian(hamiltonian_parts, external_fields)
    energies = eigvals(h)
    eigenstates = eigvecs(h)
    return Spectrum(hamiltonian_parts, energies, eigenstates)
end

"""
    transition_strengths(spectrum, g, frequency_range; polarization::SphericalUnitVector=Unpolarized())

Compute electric dipole transitions out of `g` with energy between `frequency_range[1]` and `frequency_range[2]`.

The output is a `Vector` of tuples `(frequency, strength, closest_basis_state)`, produced in
decreasing order of transition strength. The `strength` is the absolute value of the dipole matrix element,
normalized by ``D/\\sqrt{3}`` (the maximum transition dipole between ``N = 0`` and ``N = 1``). We use
the absolute value of the matrix element, rather than its square, so the results are proportional to
Rabi frequency ``Ω``.

The `polarization` keyword argument can be used to choose a specific microwave polarization. This defaults to
[`Unpolarized()`](@ref), and expects a [`SphericalUnitVector`](@ref).

There is a convenience method [`plot_transition_strengths`](@ref) that calls `transition_strengths` internally
and immediately produces a plot.

See also [`plot_transition_strengths`](@ref), [`make_hamiltonian_parts`](@ref), [`make_krb_hamiltonian_parts`](@ref),
[`Spectrum`](@ref).
"""
function transition_strengths(
    spectrum::Spectrum,
    g::State,
    frequency_range;
    polarization::SphericalUnitVector = Unpolarized(),
)
    eigenstates = spectrum.eigenstates
    energies = spectrum.energies

    (overlap, index_g) = max_overlap_with(spectrum, g)
    if overlap < 0.5
        @warn "The best overlap with your requested ground state is < 0.5."
    end
    E_g = energies[index_g]
    g_state = eigenstates[:, index_g]

    state_range =
        searchsortedfirst(energies, frequency_range[1] + E_g):searchsortedlast(energies, frequency_range[2] + E_g)
    states = [eigenstates[:, k] for k in state_range]
    frequencies = [energies[k] - E_g for k in state_range]

    T1ϵ = T⁽¹⁾(polarization)
    # Normalize to d/sqrt(3), which is the largest transition dipole (between |0,0> and |1,0>)
    h_dipole = tensor_dot(T1ϵ, spectrum.hamiltonian_parts.dipole_relative) / (1 / sqrt(3))

    strengths = [abs(g_state' * h_dipole * e) for e in states]
    closest_basis_states = map(e -> find_closest_basis_state(spectrum, e), state_range)

    out = [x for x in zip(frequencies, strengths, closest_basis_states)]
    return sort!(out, by = t -> t[2], rev = true)
end

"""
    plot_transition_strengths(spectrum, g, frequency_range; polarization::SphericalUnitVector=Unpolarized())

Plot the frequencies and strengths of electric dipole transitions out of `g`,
with energy between `frequency_range[1]` and `frequency_range[2]`.

The `polarization` keyword argument can be used to choose a specific microwave polarization. This defaults to
[`Unpolarized()`](@ref), and expects a [`SphericalUnitVector`](@ref).

This method calls [`transition_strengths`](@ref) internally.

See also [`transition_strengths`](@ref), [`make_hamiltonian_parts`](@ref), [`make_krb_hamiltonian_parts`](@ref),
[`Spectrum`](@ref).
"""
function plot_transition_strengths(
    spectrum::Spectrum,
    g::State,
    frequency_range,
    polarization::SphericalUnitVector = Unpolarized(),
)
    transitions =
        transition_strengths(spectrum, g, frequency_range; polarization = polarization)

    freqs = [t[1] for t in transitions]
    strengths = [t[2] for t in transitions]

    Gadfly.plot(
        x = freqs,
        y = strengths,
        Gadfly.Geom.hair,
        Gadfly.Geom.point,
        Gadfly.Guide.xlabel("Frequency (MHz)"),
        Gadfly.Guide.ylabel("Transition strength"),
    )
end

end # module
