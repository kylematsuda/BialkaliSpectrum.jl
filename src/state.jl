"""
    State

Represents a molecular state in the uncoupled basis.
"""
struct State
    N::Int
    mₙ::Int
    I::SVector{2,HalfInt} # [K, Rb]
    mᵢ::SVector{2,HalfInt}
end

"""
    State(N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂)

Creates a basis state ``|N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂⟩``.

"""
State(N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂) =
    State(N, mₙ, SVector(HalfInt(I₁), HalfInt(I₂)), SVector(HalfInt(mᵢ₁), HalfInt(mᵢ₂)))

"""
    KRbState(N, mₙ, mK, mRb)

Creates a basis state ``|N, m_n, m_{\\text{K}}, m_{\\text{Rb}}⟩`` for ``{}^{40}\\text{K}^{87}\\text{Rb}``.

This is a wrapper around [`State`](@ref) to avoid having to specify the nuclear spins ``I_k`` each time.

See also [`State`](@ref).
"""
KRbState(N, mₙ, mK, mRb) =
    State(N, mₙ, KRb_Parameters_Neyenhuis.I, [HalfInt(mK) HalfInt(mRb)])

"""
    state_to_named_tuple(s::State)

Creates a named tuple with fields `N`, `m_n`, `I_1`, `m_i1`, `I_2`, `m_i2` from `s`.

This is a utility function to simplify outputting to a `DataFrame`.
"""
state_to_named_tuple(s::State) =
    (N = s.N, m_n = s.mₙ, I_1 = s.I[1], m_i1 = s.mᵢ[1], I_2 = s.I[2], m_i2 = s.mᵢ[2])

"""
    state_to_string(s::State)

Pretty prints a [`State`](@ref).
"""
function state_to_string(s::State)
    N, m_n, m_i1, m_i2 = s.N, s.mₙ, s.mᵢ[1], s.mᵢ[2]
    return "|$N, $m_n, $m_i1, $m_i2⟩"
end

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
