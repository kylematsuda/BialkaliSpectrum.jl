n_hyperfine(I::HalfInt) = 2 * I + 1
n_hyperfine(s::State) = mapreduce(n_hyperfine, *, s.I)

function index_to_state(i::Int, I₁::HalfInt, I₂::HalfInt)::State
    N_Hyperfine = mapreduce(n_hyperfine, *, [I₁ I₂;])

    # Hyperfine part
    i_hyperfine = (i - 1) % N_Hyperfine
    m_1 = -I₁ + (i_hyperfine ÷ n_hyperfine(I₂))
    m_2 = -I₂ + (i_hyperfine % n_hyperfine(I₂))
    
    # Rotation part
    i_rotation = (i - 1) ÷ N_Hyperfine
    N::Int = floor(sqrt(i_rotation))
    mₙ = (i_rotation - N^2) - N
    return State(N, mₙ, I₁, m_1, I₂, m_2)
end

# TODO fix this!!!
# index_to_state(i::Int) = index_to_state(i, DEFAULT_MOLECULAR_PARAMETERS.I[1], DEFAULT_MOLECULAR_PARAMETERS.I[2])

# Todo: test for state_to_index(index_to_state(x)) == x
function state_to_index(s::State)::Int
    rotation = (s.N^2 + 1) + (s.N + s.mₙ)
    hyperfine = (s.I[1] + s.mᵢ[1]) * n_hyperfine(s.I[2]) + (s.I[2] + s.mᵢ[2])

    N_Hyperfine = n_hyperfine(s)
    return 1 + (rotation - 1) * N_Hyperfine + hyperfine
end

function order_by_overlap_with(spectrum::Spectrum, s::State)
    i = state_to_index(s)
    @assert i < size(spectrum.eigenstates, 1)
    return sortslices(spectrum.eigenstates, dims=2, lt=(x,y)->isless(abs2(x[i]), abs2(y[i])), rev=true)
end

# Returns tuple (overlap, index)
function max_overlap_with(spectrum::Spectrum, s::State)
    i = state_to_index(s)
    n_states = size(spectrum.eigenstates, 1)
    @assert i < n_states

    findmax(
        map(x -> abs2(x[i]), eachcol(spectrum.eigenstates))
    )
end

function get_energy(spectrum::Spectrum, s::State)
    return spectrum.energies[max_overlap_with(spectrum, s)[2]]
end

function get_energy_difference(spectrum::Spectrum, g::State, e::State)
    return mapreduce(x -> get_energy(spectrum, x), -, [e, g])
end

function find_closest_basis_state(spectrum::Spectrum, state::State)
    weights = map(abs2, state)
    index = findmax(weights)[2]
    return spectrum.hamiltonian_parts.basis[index]
end

function generate_basis(molecular_parameters::MolecularParameters, N_max::Int)
    n_elts::Int = (N_max + 1)^2 * mapreduce(n_hyperfine, *, molecular_parameters.I)
    return map(k -> index_to_state(k, molecular_parameters.I[1], molecular_parameters.I[2]), 1:n_elts)
end