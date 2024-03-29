module BialkaliSpectrum

using LinearAlgebra, SparseArrays, StaticArrays, Test
import WignerSymbols: wigner3j
import HalfIntegers: HalfInt
import DataFrames
import ProgressMeter
using CairoMakie

export ZeemanParameters, NuclearParameters, Polarizability, MolecularParameters

export SphericalVector, VectorX, VectorY, VectorZ
export ExternalFields, DEFAULT_FIELDS, TEST_FIELDS, generate_fields_scan

export State, basis_state, basis_index, closest_basis_state

export HamiltonianParts, make_hamiltonian_parts, hamiltonian

export get_spectrum, get_spectra
export find_closest_eigenstate, get_energy, get_energy_difference

export filter_rotational, filter_rotational!, filter_hyperfine, filter_hyperfine!,
    filter_basis_state, filter_basis_state!, expand_fields!
export wide_format

export induced_dipole_moments, transitions, adiabatic
export plot_transition_strengths, plot_induced_dipole, plot_transition_dipole
export plot_states_adiabatic, plot_states_adiabatic_weighted, plot_transitions_adiabatic

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

include("matrix_elements.jl")
include("hamiltonian.jl")

include("dataframe.jl")

"""
    get_spectrum(
        hamiltonian_parts::HamiltonianParts,
        external_fields::ExternalFields,
    )

Compute the energies and eigenstates under the external fields.

To avoid reconstructing the Hamiltonian each time, `hamiltonian_parts` can be reused over calls
to `get_spectrum`. The output is a `DataFrame`, with the following fields:

| Field name | Description |
| ---------- | :----------- |
| `fields`      | value of `external_fields` |
| `B`           | magnitude of `external_fields.B` |
| `E`           | magnitude of `external_fields.E` |
| `index`       | index of the eigenenergy (from lowest to highest energy) |
| `energy`      | energy of the state (MHz) |
| `eigenstate`  | vector of state amplitudes |
| `basis_index` | the index of nearest basis state (in wavefunction overlap) from `hamiltonian_parts.basis` |
| `N`           | rotational number `N` of the nearest basis state | 
| `m_n`         | rotational projection `m_n` of the nearest basis state |
| `I_1`         | nuclear angular momentum `I_1` of the nearest basis state |
| `m_i1`        | nuclear projection `m_i1` of the nearest basis state |
| `I_2`         | nuclear angular momentum `I_2` of the nearest basis state |
| `m_i2`        | nuclear projection `m_i2` of the nearest basis state |

See also [`get_spectra`](@ref), [`make_hamiltonian_parts`](@ref),
[`make_krb_hamiltonian_parts`](@ref), [`ExternalFields`](@ref).
"""
function get_spectrum(
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
        s -> (basis_index = basis_index(s), convert(NamedTuple, s)...),
        map(s -> closest_basis_state(hamiltonian_parts, s).state, states),
    )
    basis_states = DataFrames.DataFrame(closest_basis_states)

    out = DataFrames.hcat(df, basis_states)
    expand_fields!(out)
    return out
end

"""
    find_closest_eigenstate(spectrum, basis_state::State; tol=0.5)

Find the row in `spectrum` whose `:eigenstate` has the highest overlap with `basis_state`.

Example!!!
"""
function find_closest_eigenstate(spectrum, basis_state::State; tol=0.5)
    state_index = basis_index(basis_state)
    get_weight(e) = abs2(e[state_index])

    states =
        DataFrames.transform(spectrum, :eigenstate => (e -> map(get_weight, e)) => :weight)
    DataFrames.sort!(states, DataFrames.order(:weight, rev=true))
    out = DataFrames.first(states)

    if out.weight < tol
        @warn "The best overlap with your requested state is lower than $tol."
    end

    return out
end

"""
    get_energy(spectrum, basis_state::State; tol = 0.5)

Find the row in `spectrum` whose `:eigenstate` has the highest overlap with `basis_state`.

Example!!!
"""
function get_energy(spectrum, basis_state::State; tol = 0.5)
    closest = find_closest_eigenstate(spectrum, basis_state; tol = tol)
    return closest.energy
end

"""
    get_energy_difference(spectrum, basis_g::State, basis_e::State; tol = 0.5)

Find the row in `spectrum` whose `:eigenstate` has the highest overlap with `basis_state`.

Example!!!
"""
function get_energy_difference(spectrum, basis_g::State, basis_e::State; tol = 0.5)
    return find_closest_eigenstate(spectrum, basis_e; tol = tol).energy -
           find_closest_eigenstate(spectrum, basis_g; tol = tol).energy
end


"""
    get_spectra(hamiltonian_parts, fields_scan, df_transform)

Compute the energies and eigenstates at each point in `fields_scan`, applying
`df_transform` to each intermediate spectrum.

The `fields_scan` can be conveniently generated with [`generate_fields_scan`](@ref).
The closure `df_transform`, which must have the signature `DataFrame -> DataFrame`,
can be used to filter away unneeded rows (typically large `N` states), or to do further
analysis.

Internally, this method calls [`get_spectrum`](@ref) for each point in `fields_scan`,
calls `df_transform` on each point (if provided), and vertically concatenates the results.
The output is a `DataFrame`, see [`get_spectrum`](@ref) for details
on the dataframe columns.

See also [`make_hamiltonian_parts`](@ref), [`make_krb_hamiltonian_parts`](@ref),
[`generate_fields_scan`](@ref), [`get_spectrum`](@ref).
"""
function get_spectra(
    hamiltonian_parts::HamiltonianParts,
    fields_scan::Vector{ExternalFields},
    df_transform::Union{Function,Nothing} = nothing,
)
    out = DataFrames.DataFrame()
    ProgressMeter.@showprogress for field in fields_scan
        df = get_spectrum(hamiltonian_parts, field)

        if df_transform !== nothing
            df = df_transform(df)
        end
        DataFrames.append!(out, df)
    end

    return out
end

"""
    transform_spectra(spectra, f; groupby=:fields)

A generic function for transforming the output of [`get_spectra`](@ref).

Returns the result of grouping `spectra` by `groupby` and applying `f` to each group,
then combining the results into a new `DataFrame`. The signature of `f` must be
`DataFrame -> DataFrame`.

Example??
"""
function transform_spectra(spectra, f; groupby=:fields)
    output = DataFrames.DataFrame()
    grouped = DataFrames.groupby(spectra, groupby)

    for spectrum in grouped
        transformed = f(spectrum)
        DataFrames.append!(output, transformed)
    end

    return output
end

include("analysis.jl")
include("plotting.jl")
include("molecule_specific.jl")

end # module
