using Pkg
Pkg.add("WignerSymbols")
Pkg.add("HalfIntegers")
Pkg.add("SphericalHarmonics")

module MoleculeSpectrum

using WignerSymbols
using HalfIntegers
using SphericalHarmonics
using LinearAlgebra

I_K = HalfInt(4) # K
I_Rb = HalfInt(3/2) # Rb

struct State
    N::Int
    mₙ::Int
    I::Array{HalfIntegers.HalfInt, 2} # [K, Rb]
    mᵢ::Array{HalfIntegers.HalfInt, 2}
end

State(N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂) = State(N, mₙ, [I₁ I₂], [mᵢ₁ mᵢ₂])

n_hyperfine(I::HalfIntegers.HalfInt) = 2 * I + 1
n_hyperfine(s::State) = mapreduce(n_hyperfine, *, s.I)

N_Hyperfine = n_hyperfine(I_K) * n_hyperfine(I_Rb)

m_F(s::State) = reduce(+, s.mᵢ) + s.mₙ

function index_to_state(i::Int, I₁::HalfIntegers.HalfInt, I₂::HalfIntegers.HalfInt)::State
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

index_to_state(i::Int) = index_to_state(i, I_K, I_Rb)

δ(i, j) = ==(i, j)
δ(i::State, j::State) = δ(i.N, j.N) * δ(i.mₙ, j.mₙ) * δ(i.I, j.I) * δ(i.mᵢ, j.mᵢ)

function dipole_matrix_element(p::Int, bra::State, ket::State)::Float64
    N, mₙ, mᵢ = bra.N, bra.mₙ, bra.mᵢ
    N′, mₙ′, mᵢ′ = ket.N, ket.mₙ, ket.mᵢ

    hyperfine = δ(mᵢ, mᵢ′)
    rotational = (-1)^mₙ * sqrt((2*N + 1)*(2*N′ + 1)) * WignerSymbols.wigner3j(N, 1, N′, -mₙ, p, mₙ′) * WignerSymbols.wigner3j(N, 1, N′, 0, 0, 0)
    return hyperfine * rotational
end

rotation_matrix_element(bra::State, ket::State)::Float64 = ket.N * (ket.N + 1) * δ(bra, ket) 

function h_rot(N_max::Int, ε::Float64)
    elts::Int = (N_max + 1)^2 * N_Hyperfine
    H = zeros(elts, elts)
    for i = 1:elts
        for j = i:elts
            ket = index_to_state(i)
            bra = index_to_state(j)
            H[i, j] = rotation_matrix_element(bra, ket) + ε * dipole_matrix_element(0, bra, ket)
        end
    end
    return Hermitian(H)
end

function nuclear_quadrupole(k::Int, bra::State, ket::State)::Float64
    N, mₙ, I, mᵢ = bra.N, bra.mₙ, bra.I[k], bra.mᵢ[k]
    N′, mₙ′, I′, mᵢ′ = ket.N, ket.mₙ, ket.I[k], ket.mᵢ[k]

    # Brown and Carrington pg. 477
    p_independent = δ(I,I′) * (-1)^(I - mᵢ - mₙ) * WignerSymbols.wigner3j(N, 2, N′, 0, 0, 0) * sqrt(2*N + 1) * sqrt(2*N′ + 1) / WignerSymbols.wigner3j(I, 2, I, -I, 0, I)
    p_dependent(p) = (-1)^(p) * WignerSymbols.wigner3j(N, 2, N′, -mₙ, p, mₙ′) * WignerSymbols.wigner3j(I, 2, I, -mᵢ, -p, mᵢ′)

    return p_independent * mapreduce(p_dependent, +, -2:2)
end

function h_quadrupole(N_max::Int)
    eqQ_1 = 0.452 # K, MHz
    eqQ_2 = -1.308 # Rb, MHz
    prefactors = [eqQ_1 / 4, eqQ_2 / 4]

    elts::Int = (N_max + 1)^2 * N_Hyperfine
    H = zeros(elts, elts)
    for i = 1:elts
        for j = i:elts
            ket = index_to_state(i)
            bra = index_to_state(j)
            H[i, j] = dot(prefactors, [nuclear_quadrupole(k, bra, ket) for k in 1:2])
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
        H[i, i] = -g_r * μ_N * B_field * ket.mₙ - μ_N * B_field * dot([g_1*(1-σ_1), g_2*(1-σ_2)], ket.mᵢ)
    end
    return H
end

function h(N_max::Int, ε::Float64, B_field::Float64)
    B_rot = 1113.9514
    return B_rot * h_rot(N_max, ε) + h_quadrupole(N_max) + h_zeeman(N_max, B_field)
end

function fz(N_max::Int)
    elts::Int = (N_max + 1)^2 * N_Hyperfine
    H = zeros(elts, elts)
    for i = 1:elts
        ket = index_to_state(i)
        H[i, i] = m_F(ket)
    end
    return H
end

end # module
