module MoleculeSpectrum

using LinearAlgebra, SparseArrays, StaticArrays, Test
import WignerSymbols: wigner3j
import HalfIntegers: HalfInt
import DataFrames
import ProgressMeter
using CairoMakie

export ZeemanParameters, NuclearParameters, Polarizability, MolecularParameters
export KRb_Zeeman, KRb_Nuclear_Neyenhuis, KRb_Nuclear_Ospelkaus, KRb_Polarizability
export KRb_Parameters_Neyenhuis,
    KRb_Parameters_Ospelkaus, DEFAULT_MOLECULAR_PARAMETERS, TOY_MOLECULE_PARAMETERS

export SphericalVector, VectorX, VectorY, VectorZ
export SphericalUnitVector, UnitVectorX, UnitVectorY, UnitVectorZ, Unpolarized
export T⁽¹⁾, T⁽²⁾, get_tensor_component, tensor_dot
export ExternalFields, DEFAULT_FIELDS, TEST_FIELDS, generate_fields_scan

export State, KRbState, index_to_state, state_to_index
export order_by_overlap_with,
    max_overlap_with, find_closest_basis_state, decompose_to_basis_states
export get_energy, get_energy_difference, get_row_by_state

export HamiltonianParts, make_hamiltonian_parts, hamiltonian, make_krb_hamiltonian_parts

export calculate_spectrum, calculate_spectra_vs_fields
export calculate_transition_strengths, plot_transition_strengths
export calculate_transitions_vs_E, plot_transitions_vs_E

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
include("state.jl")

include("utility.jl")
include("matrix_elements.jl")
include("hamiltonian.jl")

"""
    calculate_spectrum(hamiltonian_parts, external_fields)

Compute the energies and eigenstates under the external fields.

To avoid reconstructing the Hamiltonian each time, `hamiltonian_parts` can be reused over calls
to `calculate_spectrum`. The output is a [`DataFrames.DataFrame`](@ref), with the following fields:

`fields`:       value of `external_fields`

`index`:        index of the eigenenergy (from lowest to highest energy)

`energy`:       energy of the state (MHz)

`eigenstate`:   vector of state amplitudes

`basis_index`:  the index of nearest basis state (in wavefunction overlap) from `hamiltonian_parts.basis`

`N`:            rotational number `N` of the nearest basis state

`m_n`:          rotational projection `m_n` of the nearest basis state

`I_1`:          nuclear angular momentum `I_1` of the nearest basis state

`m_i1`:         nuclear projection `m_i1` of the nearest basis state

`I_2`:          nuclear angular momentum `I_2` of the nearest basis state

`m_i2`:         nuclear projection `m_i2` of the nearest basis state

See also [`make_hamiltonian_parts`](@ref), [`make_krb_hamiltonian_parts`](@ref), [`ExternalFields`](@ref).
"""
function calculate_spectrum(
    hamiltonian_parts::HamiltonianParts,
    external_fields::ExternalFields,
)
    h = hamiltonian(hamiltonian_parts, external_fields)
    es = eigen(h)
    states = [c for c in eachcol(es.vectors)]
    df = DataFrames.DataFrame(
        fields = external_fields,
        index = 1:length(es.values),
        energy = es.values,
        eigenstate = states,
    )

    closest_basis_states = map(
        s -> (basis_index = state_to_index(s), state_to_named_tuple(s)...),
        map(s -> find_closest_basis_state(hamiltonian_parts, s).state, states),
    )
    basis_states = DataFrames.DataFrame(closest_basis_states)

    return DataFrames.hcat(df, basis_states)
end

"""
    calculate_spectra_vs_fields(hamiltonian_parts, fields_scan, df_transform)

Compute the energies and eigenstates at each point in `fields_scan`, applying
`df_transform` to each intermediate spectrum.

The `fields_scan` can be conveniently generated with [`generate_fields_scan`](@ref).
The closure `df_transform`, which must have the signature `DataFrame -> DataFrame`,
can be used to filter away unneeded rows (typically large `N` states), or to do further
analysis.

Internally, this method calls [`calculate_spectrum`](@ref) for each point in `fields_scan`,
calls `df_transform` on each point (if provided), and vertically concatenates the results.
The output is a [`DataFrames.DataFrame`](@ref), see [`calculate_spectrum`](@ref) for details
on the dataframe columns.

See also [`make_hamiltonian_parts`](@ref), [`make_krb_hamiltonian_parts`](@ref),
[`generate_fields_scan`](@ref), [`calculate_spectrum`](@ref).
"""
function calculate_spectra_vs_fields(
    hamiltonian_parts::HamiltonianParts,
    fields_scan::Vector{ExternalFields},
    df_transform::Union{Function,Nothing} = nothing,
)
    out = DataFrames.DataFrame()
    ProgressMeter.@showprogress for field in fields_scan
        df = calculate_spectrum(hamiltonian_parts, field)

        if df_transform !== nothing
            df = df_transform(df)
        end
        DataFrames.append!(out, df)
    end

    return out
end

include("analysis.jl")
include("plotting.jl")

end # module
