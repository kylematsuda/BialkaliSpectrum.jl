"""
    find_closest_basis_state(parts::HamiltonianParts, state)
"""
function find_closest_basis_state(parts::HamiltonianParts, state)
    weights = map(abs2, state)
    (weight, index) = findmax(weights)
    return (weight = weight, state = parts.basis[index])
end

"""
    dressed_dipole_moment(
        hamiltonian_parts::HamiltonianParts,
        g,
        e,
        polarization::Int,
    )::ComplexF64
"""
function dressed_dipole_moment(
    hamiltonian_parts::HamiltonianParts,
    g,
    e,
    polarization::Int,
)::ComplexF64
    @assert polarization <= 1 && polarization >= -1

    index = polarization + 2 # components are p = -1, 0, 1
    h_dipole = hamiltonian_parts.dipole_relative[index]
    return e' * h_dipole * g
end

"""
    dressed_dipole_moment(
        hamiltonian_parts::HamiltonianParts,
        g,
        e,
        polarization::Int,
    )::Vector{ComplexF64}
"""
function dressed_dipole_moment_vector(
    hamiltonian_parts::HamiltonianParts,
    g,
    e,
)::Vector{ComplexF64}
    return [dressed_dipole_moment(hamiltonian_parts, g, e, p) for p = -1:1]
end

"""
    get_induced_dipole_moments(
        spectra,
        hamiltonian_parts::HamiltonianParts;
        groupby=:fields
    )

Returns a new `DataFrame` containing `spectra` with an additional column `:d_ind`,
the induced dipole moment parallel to the applied electric field.
"""
function get_induced_dipole_moments(
    spectra,
    hamiltonian_parts::HamiltonianParts
)
    d_ind(state, fields) = tensor_dot(
        SphericalUnitVector(fields.E) |> T⁽¹⁾, # field orientation spherical tensor
        dressed_dipole_moment_vector(hamiltonian_parts, state, state)
    ) |> real

    return DataFrames.transform(spectra,
        [:eigenstate, :fields] =>
            DataFrames.ByRow((s, f) -> d_ind(s, f)) => :d_ind,
    )
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
    cutoff::Union{Float64,Nothing}=1e-3,
    normalization=1/3
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
                DataFrames.ByRow((dm, d0, dp) -> (1/normalization) * abs2.([dm, d0, dp]) |> sum)
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

    if frequency_range !== nothing
        DataFrames.filter!(
            :transition_frequency => f -> f >= frequency_range[1] && f <= frequency_range[2],
            transformed
        )
    end

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

"""
    connect_adiabatically(
        spectra;
        groupby=:fields,
        radius::Union{Int,Nothing}=nothing
    )

"""
function connect_adiabatically(
    spectra;
    groupby=:fields,
    radius::Union{Int,Nothing}=nothing
)
    function find_closest(prev_df, ev, e)
        if radius !== nothing
            prev_df = DataFrames.sort!(
                DataFrames.transform(
                    prev_df,
                    :energy => DataFrames.ByRow(en -> abs(en - e)) => :e_diff),
                :e_diff
            )
            prev_df = first(prev_df, radius)
        end
        
        overlap(e) = abs2(e' * ev)
        max_overlap_index = findmax(map(
            overlap,
            prev_df.eigenstate
        ))[2];

        return prev_df.adiabatic_index[max_overlap_index]
    end

    grouped = DataFrames.groupby(spectra, groupby)

    DataFrames.transform!(grouped[1], :index => (_ -> 1:DataFrames.nrow(grouped[1])) => :adiabatic_index)
    last = grouped[1]

    for i in 2:length(grouped)
        curr = grouped[i]
        DataFrames.transform!(
            curr, 
            [:eigenstate, :energy] => DataFrames.ByRow((es, en) -> find_closest(last, es, en)) => :adiabatic_index
        )
        last = curr
    end
    
    return DataFrames.combine(grouped, DataFrames.valuecols(grouped))
end