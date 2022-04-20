function find_closest_basis_state(parts::HamiltonianParts, state)
    weights = map(abs2, state)
    (weight, index) = findmax(weights)
    return (weight = weight, state = parts.basis[index])
end

function find_closest_eigenstate(spectrum, state::State; tol = 0.5)
    state_index = state_to_index(state)
    get_weight(e) = abs2(e[state_index])

    states =
        DataFrames.transform(spectrum, :eigenstate => (e -> map(get_weight, e)) => :weight)
    DataFrames.sort!(states, DataFrames.order(:weight, rev = true))
    out = DataFrames.first(states)

    if out.weight < tol
        @warn "The best overlap with your requested state is lower than $tol."
    end

    return out
end

function get_energy(spectrum, s::State; tol = 0.5)
    closest = find_closest_eigenstate(spectrum, s; tol = tol)
    return closest.energy
end

function get_energy_difference(spectrum, g::State, e::State; tol = 0.5)
    return find_closest_eigenstate(spectrum, e; tol = tol).energy -
           find_closest_eigenstate(spectrum, g; tol = tol).energy
end

function get_row_by_state(spectrum, s::State)
    index = state_to_index(s)
    return DataFrames.filter(:basis_index => bi -> bi == index, spectrum)
end


function transform_spectra(spectra, f; groupby = :fields)
    output = DataFrames.DataFrame()
    grouped = DataFrames.groupby(spectra, groupby)

    for spectrum in grouped
        transformed = f(spectrum)
        output = DataFrames.vcat(output, transformed)
    end

    return output
end

function _calculate_transition_strengths(
    spectrum,
    hamiltonian_parts::HamiltonianParts,
    g::State,
    frequency_range::Union{Vector,Nothing} = nothing;
    tol = 0.5,
    restrict_N = true,
    use_cutoff = true,
    cutoff = 1e-3,
)
    closest = find_closest_eigenstate(spectrum, g; tol = tol)
    E_g, g_vec = closest.energy, closest.eigenstate

    df = DataFrames.transform(
        spectrum,
        :energy => (es -> map(E -> E - E_g, es)) => :transition_frequency,
    )
    if frequency_range !== nothing
        DataFrames.filter!(
            :transition_frequency =>
                f -> f >= frequency_range[1] && f <= frequency_range[2],
            df,
        )
    end

    get_matrix_elements(eigenstate) = [
        calculate_dipole_matrix_element(hamiltonian_parts, g_vec, eigenstate, p) |> abs
        for p = -1:1
    ]
    DataFrames.transform!(
        df,
        [:eigenstate] =>
            (es -> map(ei -> get_matrix_elements(ei), es)) => [:d_minus, :d_0, :d_plus],
    )

    # The strength of |0> => |1> is D/sqrt(3).
    get_strengths(dp, d0, dm) = sqrt(3) * sqrt.(abs2.(dp) + abs2.(d0) + abs2.(dm))
    DataFrames.transform!(
        df,
        [:d_plus, :d_0, :d_minus] => get_strengths => [:transition_strength],
    )

    if use_cutoff
        DataFrames.filter!(:transition_strength => ts -> ts > cutoff, df)
    end

    DataFrames.filter!(:transition_frequency => f -> f > 0, df)

    if restrict_N
        DataFrames.filter!(:N => n -> n == g.N + 1, df)
    end

    DataFrames.sort!(df, DataFrames.order(:transition_strength, rev = true))
    return df
end

function calculate_transition_strengths(
    spectra,
    hamiltonian_parts::HamiltonianParts,
    g::State,
    frequency_range::Union{Vector,Nothing} = nothing;
    tol = 0.5,
    restrict_N = true,
    use_cutoff = true,
    cutoff = 1e-3,
)
    f(spectrum) = _calculate_transition_strengths(
        spectrum,
        hamiltonian_parts,
        g,
        frequency_range;
        tol = tol,
        restrict_N = restrict_N,
        use_cutoff = use_cutoff,
        cutoff = cutoff,
    )
    return transform_spectra(spectra, f)
end

function calculate_transitions_vs_E(
    hamiltonian_parts::HamiltonianParts,
    fields_scan::Vector{ExternalFields},
    g::State,
    frequency_range::Union{Vector,Nothing} = nothing;
    restrict_N = true,
    tol = 0.5,
    use_cutoff = false,
    cutoff = 1e-6,
)
    strengths(df) = calculate_transition_strengths(
        df,
        hamiltonian_parts,
        g,
        frequency_range;
        tol = tol,
        restrict_N = restrict_N,
        use_cutoff = use_cutoff,
        cutoff = cutoff,
    )
    add_E(df) =
        DataFrames.transform(df, :fields => (fs -> map(f -> f.E.magnitude, fs)) => :E)
    filter_N(df) = DataFrames.filter(:N => n -> n == g.N + 1, df)

    spectra = calculate_spectra_vs_fields(hamiltonian_parts, fields_scan, add_E ∘ strengths)
    DataFrames.sort!(spectra, [:index, :E])

    return spectra
end

function calculate_dipole_matrix_element(
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

function calculate_induced_dipole_moments(spectra, hamiltonian_parts::HamiltonianParts)
    output = DataFrames.DataFrame()
    grouped_by_fields = DataFrames.groupby(spectra, :fields)

    for spectrum in grouped_by_fields
        induced = _calculate_induced_dipole_moments(spectrum, hamiltonian_parts)
        output = DataFrames.vcat(output, induced)
    end

    return output
end

function _calculate_induced_dipole_moments(spectrum, hamiltonian_parts::HamiltonianParts)
    get_field_orientation(f::ExternalFields) = SphericalUnitVector(f.E) |> T⁽¹⁾
    get_matrix_elements(s) =
        [calculate_dipole_matrix_element(hamiltonian_parts, s, s, p) for p = -1:1]
    get_induced_dipole(s, f) = tensor_dot(get_field_orientation(f), get_matrix_elements(s))

    return DataFrames.combine(
        spectrum,
        :,
        [:eigenstate, :fields] =>
            DataFrames.ByRow((e, f) -> get_induced_dipole(e, f)) => :d_ind,
    )
end

function calculate_induced_dipole_moments_vs_E(
    hamiltonian_parts::HamiltonianParts,
    fields_scan::Vector{ExternalFields},
    hyperfine_manifold::Union{State,Nothing} = nothing,
)
    get_dipole_moments(df) = calculate_induced_dipole_moments(df, hamiltonian_parts)
    add_E(df) =
        DataFrames.transform(df, :fields => (fs -> map(f -> f.E.magnitude, fs)) => :E)

    if hyperfine_manifold !== nothing
        check_hyperfine(m_i1, m_i2) =
            m_i1 == hyperfine_manifold.mᵢ[1] && m_i2 == hyperfine_manifold.mᵢ[2]

        hyperfine_filter(df) =
            DataFrames.filter([:m_i1, :m_i2] => (m1, m2) -> check_hyperfine(m1, m2), df)

        transform = add_E ∘ get_dipole_moments ∘ hyperfine_filter
    else
        transform = add_E ∘ get_dipole_moments
    end

    return calculate_spectra_vs_fields(hamiltonian_parts, fields_scan, transform)
end


function calculate_chi_vs_E(
    hamiltonian_parts::HamiltonianParts,
    fields_scan::Vector{ExternalFields},
    g::State;
    p = 0,
    tol = 0.5,
)
    @assert p == 0 || p == -1 || p == 1

    spectra_with_induced =
        calculate_induced_dipole_moments_vs_E(hamiltonian_parts, fields_scan, g)

    grouped_by_E = DataFrames.groupby(spectra_with_induced, :E)

    output = DataFrames.DataFrame()
    for spectrum in grouped_by_E
        gs = find_closest_eigenstate(spectrum, g; tol = tol)
        d_g = gs.d_ind
        spectrum = DataFrames.combine(spectrum, :, :d_ind => (x -> d_g) => :d_g)

        strengths = calculate_transition_strengths(
            spectrum,
            hamiltonian_parts,
            g;
            restrict_N = false,
            use_cutoff = false,
            tol = tol,
        )

        # Keep only states that are connected by d^p
        DataFrames.filter!(:m_n => m -> m == g.mₙ + p, strengths)
        output = DataFrames.vcat(output, strengths)
    end

    dipole_factor = 0
    if p == 0
        col = :d_0
        dipole_factor = 2
    else
        dipole_factor = -1
        if p == -1
            col = :d_minus
        else
            col = :d_plus
        end
    end
    chi(μ_0, μ_1, μ_01) = -(μ_0 - μ_1)^2 + dipole_factor * μ_01^2

    return DataFrames.combine(
        output,
        :,
        [:d_g, :d_ind, col] => DataFrames.ByRow(chi) => :chi,
    )
end
