"""
    δ(i, j)
    δ(i::State, j::State)

Returns the Kronecker δ between `i` and `j`.

Uses `==(i,j)` if possible. If `i` and `j` are both [`State`](@ref),
then this function computes the product of the Kronecker δs over all
quantum numbers.
"""
δ(i, j) = ==(i, j)
δ(i::State, j::State) = δ(i.N, j.N) * δ(i.mₙ, j.mₙ) * δ(i.I, j.I) * δ(i.mᵢ, j.mᵢ)

"""
    WignerJ1J(j, m)
    
Returns the same as `WignerSymbols.wigner3j(j, 1, j, -m, 0, m)`,
but computes it directly from the explicit formula.
"""
WignerJ1J(j, m) = (-1)^(j - m) * m / sqrt(j * (j + 1) * (2j + 1))

"""
    WignerJ2J(j, m)
    
Returns the same as `WignerSymbols.wigner3j(j, 2, j, -m, 0, m)`,
but computes it directly from the explicit formula.
"""
WignerJ2J(j, m) =
    (-1)^(j - m) * 2 * (3m^2 - j * (j + 1)) /
    sqrt((2j - 1) * (2j) * (2j + 1) * (2j + 2) * (2j + 3))

"""
    T⁽¹⁾N(p, bra, ket)

Compute the matrix elements of the `p`th component of the spherical tensor operator ``T⁽¹⁾(N)``.

Note: in this function, Kronecker deltas are only enforced between the rotational quantum numbers.
The nuclear quantum numbers are not required to match, because this tensor could be dotted with another operator,
e.g., N ⋅ Iₖ, which allows mᵢ to change.

Any Kronecker deltas on the nuclear quantum numbers must be enforced where the dot product is eventually taken.
"""
function T⁽¹⁾N(p::Int, bra::State, ket::State)::ComplexF64
    N, mₙ = bra.N, bra.mₙ
    N′, mₙ′ = ket.N, ket.mₙ

    if δ(N, N′) && (mₙ′ - mₙ + p == 0)
        return (-1)^(N - mₙ) * sqrt(N * (N + 1) * (2N + 1)) * wigner3j(N, 1, N, -mₙ, p, mₙ′)
    else
        return 0
    end
end

"""
    T⁽¹⁾Iₖ(p, k, bra, ket)

Compute the matrix elements of the `p`th component of the spherical tensor operator ``T⁽¹⁾(Iₖ)``,
which acts on the `k`th nucleus.

Note: in this function, Kronecker deltas are only enforced between the quantum numbers of the `k`th nucleus.
The other quantum numbers are not required to match, because this tensor could be dotted with another operator,
e.g., N ⋅ Iₖ, which allows mₙ to change.

Any Kronecker deltas on the remaining quantum numbers must be enforced where the dot product is eventually taken.
"""
function T⁽¹⁾Iₖ(p::Int, k::Int, bra::State, ket::State)::ComplexF64
    I, mᵢ = bra.I[k], bra.mᵢ[k]
    I′, mᵢ′ = ket.I[k], ket.mᵢ[k]

    if δ(I, I′) && (mᵢ′ - mᵢ + p == 0)
        return (-1)^(I - mᵢ) * sqrt(I * (I + 1) * (2I + 1)) * wigner3j(I, 1, I, -mᵢ, p, mᵢ′)
    else
        return 0
    end
end

"""
    rotation_matrix_element(bra::State, ket::State)::Float64
"""
rotation_matrix_element(bra::State, ket::State)::Float64 = ket.N * (ket.N + 1) * δ(bra, ket)

"""
    dipole_matrix_element(p::Int, bra::State, ket::State)::ComplexF64
"""
function dipole_matrix_element(p::Int, bra::State, ket::State)::ComplexF64
    N, mₙ, mᵢ = bra.N, bra.mₙ, bra.mᵢ
    N′, mₙ′, mᵢ′ = ket.N, ket.mₙ, ket.mᵢ

    if δ(mᵢ, mᵢ′) && (-mₙ + p + mₙ′ == 0)
        return (-1)^mₙ *
               sqrt((2N + 1) * (2N′ + 1)) *
               wigner3j(N, 1, N′, -mₙ, p, mₙ′) *
               wigner3j(N, 1, N′, 0, 0, 0)
    else
        return 0
    end
end

"""
    nuclear_quadrupole(k::Int, bra::State, ket::State)::ComplexF64
"""
function nuclear_quadrupole(k::Int, bra::State, ket::State)::ComplexF64
    N, mₙ, I, mᵢ = bra.N, bra.mₙ, bra.I[k], bra.mᵢ[k]
    N′, mₙ′, I′, mᵢ′ = ket.N, ket.mₙ, ket.I[k], ket.mᵢ[k]

    other = (k % 2) + 1 # States of other nucleus should Kronecker delta
    other_nucleus = δ(bra.I[other], ket.I[other]) * δ(bra.mᵢ[other], ket.mᵢ[other])

    if other_nucleus &&
       δ(I, I′) &&
       (δ(N, N′) || (abs(N - N′) == 2)) &&
       (mₙ′ - mₙ == mᵢ - mᵢ′) &&
       (abs(mₙ′ - mₙ) <= 2) &&
       (abs(mᵢ′ - mᵢ) <= 2)

        # Brown and Carrington pg. 477
        p_independent =
            (-1)^(I - mᵢ - mₙ) *
            sqrt(2 * N + 1) *
            sqrt(2 * N′ + 1) *
            wigner3j(N, 2, N′, 0, 0, 0) / WignerJ2J(I, I)
        p_dependent(p) =
            (-1)^(p) * wigner3j(N, 2, N′, -mₙ, p, mₙ′) * wigner3j(I, 2, I, -mᵢ, -p, mᵢ′)
        return p_independent * mapreduce(p_dependent, +, -2:2)
    else
        return 0.0
    end
end

"""
    nuclear_spin_spin(bra::State, ket::State)::ComplexF64
"""
function nuclear_spin_spin(bra::State, ket::State)::ComplexF64
    if δ(bra.N, ket.N) &&
       δ(bra.mₙ, ket.mₙ) &&
       (ket.mᵢ[1] - bra.mᵢ[1] == bra.mᵢ[2] - ket.mᵢ[2])
        return reduce(
            +,
            [(-1)^p * T⁽¹⁾Iₖ(p, 1, bra, ket) * T⁽¹⁾Iₖ(-p, 2, bra, ket) for p = -1:1],
        )
    else
        return 0.0
    end
end

"""
    nuclear_spin_rotation(k::Int, bra::State, ket::State)::ComplexF64
"""
function nuclear_spin_rotation(k::Int, bra::State, ket::State)::ComplexF64
    other = (k % 2) + 1 # States of other nucleus should Kronecker delta
    other_nucleus = δ(bra.I[other], ket.I[other]) * δ(bra.mᵢ[other], ket.mᵢ[other])

    if other_nucleus && (ket.mₙ - bra.mₙ == bra.mᵢ[k] - ket.mᵢ[k])
        return reduce(
            +,
            [(-1)^p * T⁽¹⁾N(p, bra, ket) * T⁽¹⁾Iₖ(-p, k, bra, ket) for p = -1:1],
        )
    else
        return 0
    end
end

"""
    zeeman_rotation(p::Int, bra::State, ket::State)::ComplexF64
"""
function zeeman_rotation(p::Int, bra::State, ket::State)::ComplexF64
    if δ(bra.I, ket.I) && δ(bra.mᵢ, ket.mᵢ)
        return T⁽¹⁾N(p, bra, ket)
    else
        return 0
    end
end

"""
    zeeman_nuclear(k::Int, p::Int, bra::State, ket::State)::ComplexF64
"""
function zeeman_nuclear(k::Int, p::Int, bra::State, ket::State)::ComplexF64
    other = 1 + (k % 2)
    rotation = δ(bra.N, ket.N) && δ(bra.mₙ, ket.mₙ)
    I_other = δ(bra.I[other], ket.I[other]) && δ(bra.mᵢ[other], ket.mᵢ[other])

    if rotation && I_other
        return T⁽¹⁾Iₖ(p, k, bra, ket)
    else
        return 0
    end
end

"""
    scalar_polarizability(bra::State, ket::State)
"""
scalar_polarizability(bra::State, ket::State) = δ(bra, ket)

"""
    tensor_polarizability(p::Int, bra::State, ket::State)
"""
function tensor_polarizability(p::Int, bra::State, ket::State)
    deltas = δ(bra.N, ket.N) * δ(bra.I, ket.I) * δ(bra.mᵢ, ket.mᵢ)
    N, mₙ = bra.N, bra.mₙ
    mₙ′ = ket.mₙ

    if deltas && (N > 0) && (mₙ′ - mₙ - p == 0)
        return sqrt(6) *
               (-1)^(mₙ) *
               (2 * N + 1) *
               WignerJ2J(N, 0) *
               wigner3j(N, 2, N, -mₙ, -p, mₙ′)
    else
        return 0
    end
end
