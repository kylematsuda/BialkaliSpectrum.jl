"""
    calculate_induced_dipole_moments(
        spectra,
        hamiltonian_parts::HamiltonianParts;
        groupby=:fields
    )

Returns a new `DataFrame` containing `spectra` with an additional column `:d_ind`,
the induced dipole moment parallel to the applied electric field.
"""
function get_induced_dipole_moments(
    spectra,
    hamiltonian_parts::HamiltonianParts;
    groupby=:fields
)
    d_ind(state, fields) = tensor_dot(
        SphericalUnitVector(fields.E) |> T⁽¹⁾, # field orientation spherical tensor
        dressed_dipole_moment_vector(hamiltonian_parts, state, state)
    ) |> real

    add_d_ind(df) = DataFrames.transform(df,
        [:eigenstate, :fields] =>
            DataFrames.ByRow((s, f) -> d_ind(s, f)) => :d_ind,
    )

    return transform_spectra(spectra, add_d_ind; groupby=groupby)
end

"""
    with_ground_state(
        spectra,
        ground_basis_state::State;
        groupby=:fields
    )

Returns a new `DataFrame` with the contents of `spectra`, plus an additional
column `:ground_state` containing the `:index` of the state with the highest
overlap with `ground_basis_state`, in the same group determined by `groupby`.

Used internally but not necessary probably??
"""
function with_ground_state(
    spectra,
    ground_basis_state::State;
    groupby=:fields,
    tol=0.5
)
    output = DataFrames.DataFrame()

    function find_ground_state(df)
        closest = find_closest_eigenstate(df,
            ground_basis_state;
            tol=tol
        ) 
        DataFrames.transform(df,
            :index => DataFrames.ByRow(_ -> closest.index) => :ground_state
        )
    end

    return transform_spectra(spectra,
        find_ground_state;
        groupby=groupby
    )
end

"""
    get_transitions(
        spectra,
        hamiltonian_parts::HamiltonianParts,
        ground_basis_state::State,
        frequency_range::Union{Vector,Nothing} = nothing;
        groupby=:fields,
        tol=0.5,
        restrict_N=true,
        cutoff::Union{Float64,Nothing}=1e-3
    )

"""
function get_transitions(
    spectra,
    hamiltonian_parts::HamiltonianParts,
    ground_basis_state::State,
    frequency_range::Union{Vector,Nothing} = nothing;
    groupby=:fields,
    tol=0.5,
    restrict_N=true,
    cutoff::Union{Float64,Nothing}=1e-3
)

    function make_transitions(df)
        g_index = df.ground_state |> first
        g_row = filter(:index => i -> i == g_index, df) |> first

        g_energy = g_row.energy
        g_ev = g_row.eigenstate

        out = DataFrames.transform(df,
            :energy => 
                DataFrames.ByRow(en -> en - g_energy) => 
                    :transition_frequency
        )
        DataFrames.transform!(out,
            :eigenstate =>
                DataFrames.ByRow(ev ->
                    dressed_dipole_moment_vector(
                        hamiltonian_parts,
                        g_ev,
                        ev
                    )
                ) => [:d_minus, :d_0, :d_plus]
        )
        DataFrames.transform!(out,
            [:d_minus, :d_0, :d_plus] =>
                DataFrames.ByRow((dm, d0, dp) -> abs2.([dm, d0, dp]) |> sum)
                    => :transition_strength
        )
    end

    spectra_with_ground_state = with_ground_state(spectra,
        ground_basis_state;
        groupby=groupby,
        tol=tol
    )

    transformed = transform_spectra(
        spectra_with_ground_state,
        make_transitions;
        groupby
    )

    DataFrames.filter!(:transition_frequency => f -> f > 0, transformed)

    if restrict_N
        DataFrames.filter!(:N => n -> n == ground_basis_state.N + 1,
            transformed
        )
    end

    if cutoff !== nothing
        DataFrames.filter!(:transition_strength => s -> s > cutoff, transformed)
    end

    return transformed
end

# function _calculate_transition_strengths(
#     spectrum,
#     hamiltonian_parts::HamiltonianParts,
#     g::State,
#     frequency_range::Union{Vector,Nothing} = nothing;
#     tol = 0.5,
#     restrict_N = true,
#     use_cutoff = true,
#     cutoff = 1e-3,
# )
#     closest = find_closest_eigenstate(spectrum, g; tol = tol)
#     E_g, g_vec = closest.energy, closest.eigenstate

#     df = DataFrames.transform(
#         spectrum,
#         :energy => (es -> map(E -> E - E_g, es)) => :transition_frequency,
#     )
#     if frequency_range !== nothing
#         DataFrames.filter!(
#             :transition_frequency =>
#                 f -> f >= frequency_range[1] && f <= frequency_range[2],
#             df,
#         )
#     end

#     get_matrix_elements(eigenstate) = [
#         dressed_dipole_moment(hamiltonian_parts, g_vec, eigenstate, p) |> abs
#         for p = -1:1
#     ]
#     DataFrames.transform!(
#         df,
#         [:eigenstate] =>
#             (es -> map(ei -> get_matrix_elements(ei), es)) => [:d_minus, :d_0, :d_plus],
#     )

#     # The strength of |0> => |1> is D/sqrt(3).
#     get_strengths(dp, d0, dm) = sqrt(3) * sqrt.(abs2.(dp) + abs2.(d0) + abs2.(dm))
#     DataFrames.transform!(
#         df,
#         [:d_plus, :d_0, :d_minus] => get_strengths => [:transition_strength],
#     )

#     if use_cutoff
#         DataFrames.filter!(:transition_strength => ts -> ts > cutoff, df)
#     end

#     DataFrames.filter!(:transition_frequency => f -> f > 0, df)

#     if restrict_N
#         DataFrames.filter!(:N => n -> n == g.N + 1, df)
#     end

#     DataFrames.sort!(df, DataFrames.order(:transition_strength, rev = true))
#     return df
# end

# function calculate_transition_strengths(
#     spectra,
#     hamiltonian_parts::HamiltonianParts,
#     g::State,
#     frequency_range::Union{Vector,Nothing} = nothing;
#     tol = 0.5,
#     restrict_N = true,
#     use_cutoff = true,
#     cutoff = 1e-3,
# )
#     f(spectrum) = _calculate_transition_strengths(
#         spectrum,
#         hamiltonian_parts,
#         g,
#         frequency_range;
#         tol = tol,
#         restrict_N = restrict_N,
#         use_cutoff = use_cutoff,
#         cutoff = cutoff,
#     )
#     return transform_spectra(spectra, f)
# end

# function calculate_transitions_vs_E(
#     hamiltonian_parts::HamiltonianParts,
#     fields_scan::Vector{ExternalFields},
#     g::State,
#     frequency_range::Union{Vector,Nothing} = nothing;
#     restrict_N = true,
#     tol = 0.5,
#     use_cutoff = false,
#     cutoff = 1e-6,
# )
#     strengths(df) = calculate_transition_strengths(
#         df,
#         hamiltonian_parts,
#         g,
#         frequency_range;
#         tol = tol,
#         restrict_N = restrict_N,
#         use_cutoff = use_cutoff,
#         cutoff = cutoff,
#     )
#     add_E(df) =
#         DataFrames.transform(df, :fields => (fs -> map(f -> f.E.magnitude, fs)) => :E)
#     filter_N(df) = DataFrames.filter(:N => n -> n == g.N + 1, df)

#     spectra = calculate_spectra_vs_fields(hamiltonian_parts, fields_scan, add_E ∘ strengths)
#     DataFrames.sort!(spectra, [:index, :E])

#     return spectra
# end

# function calculate_induced_dipole_moments_vs_E(
#     hamiltonian_parts::HamiltonianParts,
#     fields_scan::Vector{ExternalFields},
#     hyperfine_manifold::Union{State,Nothing} = nothing,
# )
#     get_dipole_moments(df) = get_induced_dipole_moments(df, hamiltonian_parts)
#     add_E(df) =
#         DataFrames.transform(df, :fields => (fs -> map(f -> f.E.magnitude, fs)) => :E)

#     if hyperfine_manifold !== nothing
#         check_hyperfine(m_i1, m_i2) =
#             m_i1 == hyperfine_manifold.mᵢ[1] && m_i2 == hyperfine_manifold.mᵢ[2]

#         hyperfine_filter(df) =
#             DataFrames.filter([:m_i1, :m_i2] => (m1, m2) -> check_hyperfine(m1, m2), df)

#         transform = add_E ∘ get_dipole_moments ∘ hyperfine_filter
#     else
#         transform = add_E ∘ get_dipole_moments
#     end

#     return calculate_spectra_vs_fields(hamiltonian_parts, fields_scan, transform)
# end

# function calculate_chi_vs_E(
#     hamiltonian_parts::HamiltonianParts,
#     fields_scan::Vector{ExternalFields},
#     g::State;
#     p = 0,
#     tol = 0.5,
# )
#     @assert p == 0 || p == -1 || p == 1

#     spectra_with_induced =
#         calculate_induced_dipole_moments_vs_E(hamiltonian_parts, fields_scan, g)

#     grouped_by_E = DataFrames.groupby(spectra_with_induced, :E)

#     output = DataFrames.DataFrame()
#     for spectrum in grouped_by_E
#         gs = find_closest_eigenstate(spectrum, g; tol = tol)
#         d_g = gs.d_ind
#         spectrum = DataFrames.combine(spectrum, :, :d_ind => (x -> d_g) => :d_g)

#         strengths = calculate_transition_strengths(
#             spectrum,
#             hamiltonian_parts,
#             g;
#             restrict_N = false,
#             use_cutoff = false,
#             tol = tol,
#         )

#         # Keep only states that are connected by d^p
#         DataFrames.filter!(:m_n => m -> m == g.mₙ + p, strengths)
#         output = DataFrames.vcat(output, strengths)
#     end

#     dipole_factor = 0
#     if p == 0
#         col = :d_0
#         dipole_factor = 2
#     else
#         dipole_factor = -1
#         if p == -1
#             col = :d_minus
#         else
#             col = :d_plus
#         end
#     end
#     chi(μ_0, μ_1, μ_01) = -(μ_0 - μ_1)^2 + dipole_factor * μ_01^2

#     return DataFrames.combine(
#         output,
#         :,
#         [:d_g, :d_ind, col] => DataFrames.ByRow(chi) => :chi,
#     )
# end
