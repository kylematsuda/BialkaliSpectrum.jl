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
export order_by_overlap_with, max_overlap_with, find_closest_basis_state, decompose_to_basis_states
export get_energy, get_energy_difference

export HamiltonianParts, make_hamiltonian_parts, hamiltonian, make_krb_hamiltonian_parts

export Spectrum, calculate_spectrum
export find_transition_strengths, plot_transition_strengths
export calculate_dipolar_interaction

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

include("matrix_elements.jl")
include("hamiltonian.jl")

struct Spectrum
    hamiltonian_parts::HamiltonianParts
    energies::Vector{Float64}
    eigenstates::Any
end

get_eigenstate(spectrum, k) = spectrum.eigenstates[:, k]

include("utility.jl")

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
    find_transition_strengths(spectrum::Spectrum, g::State, frequency_range; polarization::Union{Int, SphericalUnitVector, Nothing}=nothing)

Compute electric dipole transitions out of `g` with energy between `frequency_range[1]` and `frequency_range[2]`.

The output is a `Vector` of tuples `(frequency, strength, closest_basis_state, eigenstate_index)`, produced in
decreasing order of transition strength. The `strength` is the absolute value of the dipole matrix element,
normalized by ``D/\\sqrt{3}`` (the maximum transition dipole between ``N = 0`` and ``N = 1``). We use
the absolute value of the matrix element, rather than its square, so the results are proportional to
Rabi frequency ``Ω``.

The `polarization` keyword argument can be used to choose a specific microwave polarization. This defaults to `nothing`,
which returns the incoherent sum over ``σ-``, ``π``, and ``σ+``. If `polarization` is an `Int`, then it is interpreted as
the spherical component `p = -1:1` of the dipole operator (`p == -1` corresponds to ``σ-`` polarization). If
`polarization` is a `SphericalUnitVector`, then the polarization is interpreted as linear along that axis.

There is a convenience method [`plot_transition_strengths`](@ref) that immediately produces a plot from the result.

See also [`plot_transition_strengths`](@ref), [`make_hamiltonian_parts`](@ref), [`make_krb_hamiltonian_parts`](@ref),
[`Spectrum`](@ref).
"""
function find_transition_strengths(
    spectrum::Spectrum,
    g::State,
    frequency_range;
    polarization::Union{Int, SphericalUnitVector, Nothing} = nothing
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

    if polarization === nothing
        strengths = [calculate_transition_strength_incoherent(spectrum, g_state, e) for e in states]
    else
        strengths = [abs(calculate_transition_strength_coherent(spectrum, g_state, e, polarization)) for e in states]
    end
    closest_basis_states = map(e -> find_closest_basis_state(spectrum, e), state_range)

    out = [x for x in zip(frequencies, strengths, closest_basis_states, state_range)]
    return sort!(out, by = t -> t[2], rev = true)
end


"""
    calculate_transition_strength_incoherent(spectrum::Spectrum, g::Vector{ComplexF64}, e::Vector{ComplexF64})

Compute the transition strength from `g` to `e` for an even incoherent mixture of polarizations.

The incoherent sum is formed by calculating the squared matrix elements of ``|⟨g|H_p|e⟩|^2`` for each
spherical component `p = -1:1` of the dipole Hamiltonian, summing the three values, and then taking the square root.
"""
function calculate_transition_strength_incoherent(spectrum::Spectrum, g::Vector{ComplexF64}, e::Vector{ComplexF64})
    h_dipole = spectrum.hamiltonian_parts.dipole_relative
    intensities = [abs2(e' * h_dipole[p] * g) for p in 1:3]

    # Normalize to d/sqrt(3), which is the largest transition dipole (between |0,0> and |1,0>)
    strength = sqrt(reduce(+, intensities)) / (1 / sqrt(3))
    return strength
end

"""
    calculate_transition_strength_coherent(spectrum::Spectrum, g::Vector{ComplexF64}, e::Vector{ComplexF64}, polarization::Int)
    calculate_transition_strength_coherent(spectrum::Spectrum, g::Vector{ComplexF64}, e::Vector{ComplexF64}, polarization::SphericalUnitVector)

Compute the transition strength from `g` to `e`, driven by a field with the given `polarization`.

If `polarization` is an `Int`, then it is interpreted as the spherical component `p = -1:1` of the dipole operator (`p == -1` corresponds to 
``σ-`` polarization). If `polarization` is a `SphericalUnitVector`, then the polarization is interpreted as linear along that axis.

"""
function calculate_transition_strength_coherent(spectrum::Spectrum, g::Vector{ComplexF64}, e::Vector{ComplexF64}, polarization::Int)::ComplexF64
    @assert polarization <= 1 && polarization >= -1

    index = polarization + 2 # components are p = -1, 0, 1
    h_dipole = spectrum.hamiltonian_parts.dipole_relative[index]
    return e' * h_dipole * g / (1 / sqrt(3))
end

function calculate_transition_strength_coherent(spectrum::Spectrum, g::Vector{ComplexF64}, e::Vector{ComplexF64}, polarization::SphericalUnitVector)::ComplexF64
    h_dipole = tensor_dot(T⁽¹⁾(polarization), spectrum.hamiltonian_parts.dipole_relative)
    return e' * h_dipole * g / (1 / sqrt(3))
end

"""
    plot_transition_strengths(spectrum::Spectrum, g::State, frequency_range; polarization::Union{Int, SphericalUnitVector, Nothing}=nothing)

Plot the frequencies and strengths of electric dipole transitions out of `g`,
with energy between `frequency_range[1]` and `frequency_range[2]`.

The `polarization` keyword argument can be used to choose a specific microwave polarization. This defaults to `nothing`,
which returns the incoherent sum over ``σ-``, ``π``, and ``σ+``. If `polarization` is an `Int`, then it is interpreted as
the spherical component `p = -1:1` of the dipole operator (`p == -1` corresponds to ``σ-`` polarization). If
`polarization` is a `SphericalUnitVector`, then the polarization is interpreted as linear along that axis.

This method calls [`find_transition_strengths`](@ref) internally.

See also [`find_transition_strengths`](@ref), [`make_hamiltonian_parts`](@ref), [`make_krb_hamiltonian_parts`](@ref),
[`Spectrum`](@ref).
"""
function plot_transition_strengths(
    spectrum::Spectrum,
    g::State,
    frequency_range;
    polarization::Union{Int, SphericalUnitVector, Nothing}=nothing
)
    transitions =
        find_transition_strengths(spectrum, g, frequency_range; polarization = polarization)

    freqs = [t[1] for t in transitions]
    strengths = [abs(t[2]) for t in transitions]

    Gadfly.plot(
        x = freqs,
        y = strengths,
        Gadfly.Geom.hair,
        Gadfly.Geom.point,
        Gadfly.Guide.xlabel("Frequency (MHz)"),
        Gadfly.Guide.ylabel("Transition strength"),
    )
end

function calculate_dipolar_interaction(
    spectrum::Spectrum,
    g::State,
    e::State;
    p::Int = 0
)
    (overlap, index_g) = max_overlap_with(spectrum, g)
    if overlap < 0.5
        @warn "The best overlap with your requested ground state is < 0.5."
    end

    (overlap, index_e) = max_overlap_with(spectrum, e)
    if overlap < 0.5
        @warn "The best overlap with your requested excited state is < 0.5."
    end

    return calculate_dipolar_interaction(
        spectrum,
        get_eigenstate(spectrum, index_g),
        get_eigenstate(spectrum, index_e);
        p=p
    )
end

function calculate_dipolar_interaction(
    spectrum::Spectrum,
    g::Vector{ComplexF64},
    e::Vector{ComplexF64};
    p::Int = 0
)
    d_1 = [Complex(calculate_transition_strength_coherent(spectrum, g, e, pol)) for pol=-1:1]
    d_2 = [Complex(calculate_transition_strength_coherent(spectrum, e, g, pol)) for pol=-1:1]
    return get_tensor_component(p, T⁽²⁾(d_1, d_2)) * sqrt(6) / 2
end

end # module
