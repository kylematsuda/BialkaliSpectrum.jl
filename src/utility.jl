n_hyperfine(I) = 2 * I + 1
n_hyperfine(s::State) = mapreduce(n_hyperfine, *, s.I)

"""
    index_to_state(i, I₁, I₂)

Returns the `State` corresponding to the `i`th member of the basis.

The uncoupled basis ``|N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂⟩`` is ordered
with the quantum numbers on the left changing the slowest.

See also [`State`](@ref), [`state_to_index`](@ref).

# Examples
```jldoctest
julia> s = index_to_state(1, 4, 3/2)
State(0, 0, HalfIntegers.Half{Int64}[4, 3/2], HalfIntegers.Half{Int64}[-4, -3/2])
```

```jldoctest
julia> s = index_to_state(37, 4, 3/2)
State(1, -1, HalfIntegers.Half{Int64}[4, 3/2], HalfIntegers.Half{Int64}[-4, -3/2])
```
"""
function index_to_state(i::Int, I₁, I₂)::State
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

"""
    state_to_index(s::State)

Returns index of state `s` in the basis.

The uncoupled basis ``|N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂⟩`` is ordered
with the quantum numbers on the left changing the slowest.

See also [`State`](@ref), [`index_to_state`](@ref).

# Examples
```jldoctest
julia> state_to_index(KRbState(1, 1, -4, 1/2))
111
```

```jldoctest
julia> state_to_index(index_to_state(42, 4, 3/2))
42
``` 
"""
function state_to_index(s::State)::Int
    rotation = (s.N^2 + 1) + (s.N + s.mₙ)
    hyperfine = (s.I[1] + s.mᵢ[1]) * n_hyperfine(s.I[2]) + (s.I[2] + s.mᵢ[2])

    N_Hyperfine = n_hyperfine(s)
    return 1 + (rotation - 1) * N_Hyperfine + hyperfine
end

"""
    order_by_overlap_with(spectrum, target)

Orders the list of eigenstates from `spectrum` by their
wavefunction overlap with `target`.

See also [`calculate_spectrum`](@ref), [`max_overlap_with`](@ref).
"""
function order_by_overlap_with(spectrum::Spectrum, target::State)
    i = state_to_index(target)
    @assert i < size(spectrum.eigenstates, 1)
    return sortslices(
        spectrum.eigenstates,
        dims = 2,
        lt = (x, y) -> isless(abs2(x[i]), abs2(y[i])),
        rev = true,
    )
end

"""
    max_overlap_with(spectrum, target)

Find the eigenstate from `spectrum` with the most wavefunction overlap with `target`.

Returns a tuple `(overlap, index)`, where `overlap` is the wavefunction overlap with
`target` and `index` is the position of the eigenstates in `spectrum`.

See also [`calculate_spectrum`](@ref), [`order_by_overlap_with`](@ref).
"""
function max_overlap_with(spectrum::Spectrum, target::State)
    i = state_to_index(target)
    n_states = size(spectrum.eigenstates, 1)
    @assert i < n_states

    findmax(map(x -> abs2(x[i]), eachcol(spectrum.eigenstates)))
end

"""
    find_closest_basis_state(spectrum, index)

Find the basis state with the most wavefunction overlap with the eigenstate from `spectrum` at `index`.

This method can be thought of as an inverse to [`max_overlap_with`](@ref). 

See also [`calculate_spectrum`](@ref), [`max_overlap_with`](@ref).
"""
function find_closest_basis_state(spectrum::Spectrum, index)
    state = spectrum.eigenstates[:, index]
    weights = map(abs2, state)
    index = findmax(weights)[2]
    return spectrum.hamiltonian_parts.basis[index]
end

to_indices_and_weights(state) = sort!(
    collect(Tuple{Int, Float64}, enumerate(map(abs2, state)));
    by=k->k[2],
    rev = true
)

"""
    decompose_to_basis_states(spectrum::Spectrum, index::Int)
    decompose_to_basis_states(spectrum::Spectrum, state::Vector{ComplexF64})
    decompose_to_basis_states(state::Vector{ComplexF64}, I1, I2)

Decompose the given `state` into basis states in order of decreasing weight.

Outputs a `Vector{Tuple{State, Float64}}`, giving a list `[(basis_state, weight)]`. The first variant
provides the decomposition of the member of `spectrum.eigenstates` at position `index`. The second and third
variants provide the decomposition of the vector `state`. The third variant does not require `spectrum`
and recomputes the basis set according to `I1` and `I2`.

See also [`calculate_spectrum`](@ref), [`max_overlap_with`](@ref).
"""
function decompose_to_basis_states(spectrum::Spectrum, index::Int)
    return decompose_to_basis_states(spectrum, spectrum.eigenstates[:, index])
end

function decompose_to_basis_states(spectrum::Spectrum, state::Vector{ComplexF64})
    indices_and_weights = to_indices_and_weights(state)
    return map(
        ((index, weight),) -> (spectrum.hamiltonian_parts.basis[index], weight),
        indices_and_weights
    )
end

function decompose_to_basis_states(state::Vector{ComplexF64}, I1, I2)
    indices_and_weights = to_indices_and_weights(state)
    return map(
        ((index, weight),) -> (index_to_state(index, I1, I2), weight),
        indices_and_weights
    )
end


"""
    get_energy(spectrum, target)

Return the energy of the eigenstate with the most wavefunction overlap with `target`.

See also [`calculate_spectrum`](@ref), [`max_overlap_with`](@ref).
"""
function get_energy(spectrum::Spectrum, target::State)
    return spectrum.energies[max_overlap_with(spectrum, target)[2]]
end

"""
    get_energy_difference(spectrum, g, e)

Return the difference in energy between the eigenstates from `spectrum`
that have the most wavefunction overlap with `g` and `e`.

This method is used to calculate transition frequencies; `g` represents the ground
state, and `e` the excited state.

See also [`calculate_spectrum`](@ref), [`get_energy`](@ref).
"""
function get_energy_difference(spectrum::Spectrum, g::State, e::State)
    return mapreduce(x -> get_energy(spectrum, x), -, [e, g])
end

function generate_basis(molecular_parameters::MolecularParameters, N_max::Int)
    n_elts::Int = (N_max + 1)^2 * mapreduce(n_hyperfine, *, molecular_parameters.I)
    return map(
        k -> index_to_state(k, molecular_parameters.I[1], molecular_parameters.I[2]),
        1:n_elts,
    )
end
