function find_closest_basis_state(parts::HamiltonianParts, state)
    weights = map(abs2, state)
    (weight, index) = findmax(weights)
    return (weight = weight, state = parts.basis[index])
end

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

function dressed_dipole_moment_vector(
    hamiltonian_parts::HamiltonianParts,
    g,
    e,
)::Vector{ComplexF64}
    return [dressed_dipole_moment(hamiltonian_parts, g, e, p) for p = -1:1]
end

function unpolarized_transition_strength(
    hamiltonian_parts::HamiltonianParts,
    g,
    e;
    normalization=1/3
)
    nonnormalized = abs2.(
        dressed_dipole_moment_vector(hamiltonian_parts, g, e)
    ) |> sum

    return nonnormalized / normalization
end
