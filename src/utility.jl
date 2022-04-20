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

function find_closest_eigenstate(
    spectrum,
    state::State;
    tol = 0.5
)
    state_index = state_to_index(state)
    get_weight(e) = abs2(e[state_index])

    states = DataFrames.transform(spectrum, :eigenstate => (e -> map(get_weight, e)) => :weight)
    DataFrames.sort!(states, DataFrames.order(:weight, rev=true))
    out = DataFrames.first(states)

    if out.weight < tol
        @warn "The best overlap with your requested state is lower than $tol."
    end

    return out
end

function get_energy(
    spectrum,
    s::State;
    tol = 0.5
)
    closest = find_closest_eigenstate(spectrum, s; tol=tol)
    return closest.energy
end

function get_energy_difference(
    spectrum,
    g::State,
    e::State;
    tol = 0.5
)
    return find_closest_eigenstate(spectrum, e; tol=tol).energy -
        find_closest_eigenstate(spectrum, g; tol=tol).energy
end

function get_row_by_state(
    spectrum,
    s::State
)
    index = state_to_index(s)
    return DataFrames.filter(:basis_index => bi -> bi == index, spectrum)
end
