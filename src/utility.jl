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

index_to_state(i::Int) = index_to_state(i, DEFAULT_MOLECULAR_PARAMETERS.I[1], DEFAULT_MOLECULAR_PARAMETERS.I[2])

# Todo: test for state_to_index(index_to_state(x)) == x
function state_to_index(s::State)::Int
    rotation = (s.N^2 + 1) + (s.N + s.mₙ)
    hyperfine = (s.I[1] + s.mᵢ[1]) * n_hyperfine(s.I[2]) + (s.I[2] + s.mᵢ[2])

    N_Hyperfine = n_hyperfine(s)
    return 1 + (rotation - 1) * N_Hyperfine + hyperfine
end

function order_by_overlap_with(s::State, eigenstates::Matrix)
    i = state_to_index(s)
    @assert i < size(eigenstates, 1)
    return sortslices(eigenstates, dims=2, lt=(x,y)->isless(abs2(x[i]), abs2(y[i])), rev=true)
end

# Returns tuple (overlap, index)
function max_overlap_with(s::State, eigenstates::Matrix)
    i = state_to_index(s)
    n_states = size(eigenstates, 1)
    @assert i < n_states

    findmax(
        map(x -> abs2(x[i]), eachcol(eigenstates))
    )
end

function get_energy(s::State, energies::Vector, eigenstates::Matrix)
    return energies[max_overlap_with(s, eigenstates)[2]]
end

function get_energy_difference(g::State, e::State, energies::Vector, eigenstates::Matrix)
    return mapreduce(x -> get_energy(x, energies, eigenstates), -, [e, g])
end

function generate_basis(molecular_parameters::MolecularParameters, N_max::Int)
    n_elts::Int = (N_max + 1)^2 * mapreduce(n_hyperfine, *, molecular_parameters.I)
    return map(index_to_state, 1:n_elts)
end