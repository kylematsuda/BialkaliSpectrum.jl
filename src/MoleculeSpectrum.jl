using Pkg
Pkg.add("WignerSymbols")
Pkg.add("HalfIntegers")
Pkg.add("SphericalHarmonics")

module MoleculeSpectrum

using WignerSymbols
using HalfIntegers
using SphericalHarmonics
using LinearAlgebra

greet() = print("Hello World!")

I_1 = HalfInt(4) # K
I_2 = HalfInt(3/2) # Rb

N_I_1 = 2 * I_1 + 1
N_I_2 = 2 * I_2 + 1
N_Hyperfine = N_I_1 * N_I_2

struct State
    N::Int
    m_N::Int
    m_1::HalfIntegers.HalfInt # K
    m_2::HalfIntegers.HalfInt # Rb
end

function index_to_state(i::Int)::State
    # Hyperfine part
    i_hyperfine = (i - 1) % N_Hyperfine
    m_1 = -I_1 + (i_hyperfine ÷ N_I_2)
    m_2 = -I_2 + (i_hyperfine % N_I_2)
    
    # Rotation part
    i_rotation = (i - 1) ÷ N_Hyperfine
    N::Int = floor(sqrt(i_rotation))
    m_N = (i_rotation - N^2) - N
    return State(N, m_N, m_1, m_2)
end

δ(i, j) = ==(i, j)

function dipole_matrix_element(p::Int, bra::State, ket::State)::Float64
    N, m_N, m_1, m_2 = bra.N, bra.m_N, bra.m_1, bra.m_2
    Np, m_Np, m_1p, m_2p = ket.N, ket.m_N, ket.m_1, ket.m_2

    hyperfine = δ(m_1, m_1p) * δ(m_2, m_2p)
    rotational = (-1)^m_N * sqrt((2*N + 1)*(2*Np + 1)) * WignerSymbols.wigner3j(N, 1, Np, -m_N, p, m_Np) * WignerSymbols.wigner3j(N, 1, Np, 0, 0, 0)
    return hyperfine * rotational
end

function rotation_matrix_element(bra::State, ket::State)::Float64
    N1, m_N1, m_11, m_21 = bra.N, bra.m_N, bra.m_1, bra.m_2
    N2, m_N2, m_12, m_22 = ket.N, ket.m_N, ket.m_1, ket.m_2

    return N1 * (N1 + 1) * δ(m_N1, m_N2) * δ(N1, N2) * δ(m_11, m_12) * δ(m_21, m_22) 
end

function h_rot(N_max::Int, ε::Float64)
    elts::Int = (N_max + 1)^2 * N_Hyperfine
    H = zeros(elts, elts)
    for i = 1:elts
        for j = 1:elts
            ket = index_to_state(i)
            bra = index_to_state(j)
            H[i, j] = rotation_matrix_element(bra, ket) + ε * dipole_matrix_element(0, bra, ket)
        end
    end
    return H
end

function I_plus(I::HalfIntegers.HalfInt, m::HalfIntegers.HalfInt)::Float64
    sqrt(I * (I + 1) - m * (m + 1))
end

function I_minus(I::HalfIntegers.HalfInt, m::HalfIntegers.HalfInt)::Float64
    sqrt(I * (I + 1) - m * (m - 1))
end

function I_p(p::Int, I::HalfIntegers.HalfInt, m::HalfIntegers.HalfInt)::Float64
    @assert p == 0 || p == 1 || p == -1

    if p == 0
        return m
    else
        if (p * m > I)
            return 0
        else 
            return -p * sqrt(I * (I + 1) - m * (m + p)) / sqrt(2)
        end
    end
end


# Spherical basis operator
function operator_I_1(p::Int, bra::State, ket::State)::Float64 
    N1, m_N1, m_11, m_21 = ket.N, ket.m_N, ket.m_1, ket.m_2
    N2, m_N2, m_12, m_22 = bra.N, bra.m_N, bra.m_1, bra.m_2

    deltas = δ(N1, N2) * δ(m_N1, m_N2) * δ(m_21, m_22) * δ(m_11 + p, m_12)
    return deltas * I_p(p, I_1, m_11)
end

# Spherical basis operator
function operator_I_2(p::Int, bra::State, ket::State)::Float64 
    N1, m_N1, m_11, m_21 = ket.N, ket.m_N, ket.m_1, ket.m_2
    N2, m_N2, m_12, m_22 = bra.N, bra.m_N, bra.m_1, bra.m_2

    deltas = δ(N1, N2) * δ(m_N1, m_N2) * δ(m_11, m_12)  * δ(m_21 + p, m_22)
    return deltas * I_p(p, I_2, m_21)
end


# Spherical basis operator I^(p_left) I^(p_right)
function operator_II_1(p_left::Int, p_right::Int, bra::State, ket::State)::Float64 
    N1, m_N1, m_11, m_21 = ket.N, ket.m_N, ket.m_1, ket.m_2
    N2, m_N2, m_12, m_22 = bra.N, bra.m_N, bra.m_1, bra.m_2

    deltas = δ(N1, N2) * δ(m_N1, m_N2) * δ(m_21, m_22) * δ(m_11 + p_left + p_right, m_12)
    return deltas * I_p(p_left, I_1, m_11 + p_right) * I_p(p_right, I_1, m_11)
end

# Spherical basis operator I^(p_left) I^(p_right)
function operator_II_2(p_left::Int, p_right::Int, bra::State, ket::State)::Float64 
    N1, m_N1, m_11, m_21 = ket.N, ket.m_N, ket.m_1, ket.m_2
    N2, m_N2, m_12, m_22 = bra.N, bra.m_N, bra.m_1, bra.m_2

    deltas = δ(N1, N2) * δ(m_N1, m_N2) * δ(m_11, m_12) * δ(m_21 + p_left + p_right, m_22)
    return deltas * I_p(p_left, I_2, m_21 + p_right) * I_p(p_right, I_2, m_21)
end

function T_1(α, β, bra::State, ket::State)::ComplexF64
    components = [
        operator_II_1(1, 1, bra, ket),
        (operator_II_1(1, 0, bra, ket) + operator_II_1(0, 1, bra, ket)) / sqrt(2),
        (operator_II_1(1, -1, bra, ket) + operator_II_1(-1, 1, bra, ket) + 2 * operator_II_1(0, 0, bra, ket)) / sqrt(6),
        (operator_II_1(-1, 0, bra, ket) + operator_II_1(0, -1, bra, ket)) / sqrt(2),
        operator_II_1(-1, -1, bra, ket)
    ]

    spherical_harmonics = [sqrt(4.0 * π / 5.0) * SphericalHarmonics.sphericalharmonic(α, β, l = 2, m = mm) for mm in -2:2]
    return LinearAlgebra.dot(components, spherical_harmonics)
end

function T_2(α, β, bra::State, ket::State)::ComplexF64
    components = [
        operator_II_2(1, 1, bra, ket),
        (operator_II_2(1, 0, bra, ket) + operator_II_2(0, 1, bra, ket)) / sqrt(2),
        (operator_II_2(1, -1, bra, ket) + operator_II_2(-1, 1, bra, ket) + 2 * operator_II_2(0, 0, bra, ket)) / sqrt(6),
        (operator_II_2(-1, 0, bra, ket) + operator_II_2(0, -1, bra, ket)) / sqrt(2),
        operator_II_2(-1, -1, bra, ket)
    ]

    spherical_harmonics = [sqrt(4.0 * π / 5.0) * SphericalHarmonics.sphericalharmonic(α, β, l = 2, m = mm) for mm in -2:2]
    return LinearAlgebra.dot(components, spherical_harmonics)
end

function h_quadrupole(N_max::Int, α, β)
    eqQ_1 = 0.452 # K, MHz
    eqQ_2 = -1.308 # Rb, MHz
    prefactors = [eqQ_1 / (I_1*(I_1-1)), eqQ_2 / (I_2*(I_2-1))]

    elts::Int = (N_max + 1)^2 * N_Hyperfine
    H = zeros(elts, elts)
    for i = 1:elts
        for j = i:elts
            ket = index_to_state(i)
            bra = index_to_state(j)
            H[i, j] = prefactors[1] * T_1(α, β, bra, ket) + prefactors[2] * T_2(α, β, bra, ket)
        end
    end
    return Hermitian(H)
end

# no angle dependence for now
function h_zeeman(N_max::Int, B_field::Float64)
    g_r = 0.014 # Aldegunde, PRA 78, 033434 (2008)
    μ_N = 7.622593285e-4 # MHz/G

    # Aldegunde, PRA 78, 033434 (2008)
    g_1 = -0.324 # g_K
    g_2 = 1.834 # g_Rb
    σ_1 = 1321e-6
    σ_2 = 3469e-6

    elts::Int = (N_max + 1)^2 * N_Hyperfine
    H = zeros(elts, elts)
    for i = 1:elts
        ket = index_to_state(i)
        H[i, i] = -g_r * μ_N * B_field * ket.m_N - g_1 * μ_N * B_field * ket.m_1 * (1 - σ_1) - g_2 * μ_N * B_field * ket.m_2 * (1 - σ_2)
    end
    return Hermitian(H)
end

function h(N_max::Int, α::Float64, β::Float64, ε::Float64, B_field::Float64)
    B_rot = 1113.9514
    return B_rot * h_rot(N_max, ε) + h_quadrupole(N_max, α, β) + h_zeeman(N_max, B_field)
end

end # module
