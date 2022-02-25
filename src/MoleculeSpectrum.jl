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

# Todo: test for state_to_index(index_to_state(x)) == x
function state_to_index(s::State)::Int
    rotation = (s.N^2 + 1) + (s.N + s.mₙ)
    hyperfine = (s.I[1] + s.mᵢ[1]) * n_hyperfine(s.I[2]) + (s.I[2] + s.mᵢ[2])
    return 1 + (rotation - 1) * N_Hyperfine + hyperfine
end

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

    other = (k % 2) + 1 # States of other nucleus should Kronecker delta
    other_nucleus = δ(bra.I[other], ket.I[other]) * δ(bra.mᵢ[other], ket.mᵢ[other])

    # Brown and Carrington pg. 477
    p_independent = other_nucleus * δ(I,I′) * (-1)^(I - mᵢ - mₙ) * sqrt(2*N + 1) * sqrt(2*N′ + 1) * WignerSymbols.wigner3j(N, 2, N′, 0, 0, 0) / WignerSymbols.wigner3j(I, 2, I, -I, 0, I)
    p_dependent(p) = (-1)^(p) * WignerSymbols.wigner3j(N, 2, N′, -mₙ, p, mₙ′) * WignerSymbols.wigner3j(I, 2, I, -mᵢ, -p, mᵢ′)

    return p_independent * mapreduce(p_dependent, +, -2:2)
end

function nuclear_spin_spin(bra::State, ket::State)::Float64
    N, mₙ, (I_1, I_2), (mᵢ_1, mᵢ_2) = bra.N, bra.mₙ, bra.I, bra.mᵢ
    N′, mₙ′, (I_1′, I_2′), (mᵢ_1′, mᵢ_2′) = ket.N, ket.mₙ, ket.I, ket.mᵢ

    deltas = δ(N, N′) * δ(mₙ, mₙ′) * δ(I_1, I_1′) * δ(I_2, I_2′)
    p_independent = (-1)^(I_1 + I_2 - mᵢ_1 - mᵢ_2) * sqrt(I_1 * (I_1 + 1) * (2*I_1 + 1)) * sqrt(I_2 * (I_2 + 1) * (2*I_2 + 1))
    p_dependent(p) = (-1)^p * WignerSymbols.wigner3j(I_1, 1, I_1′, -mᵢ_1, p, mᵢ_1′) * WignerSymbols.wigner3j(I_2, 1, I_2′, -mᵢ_2, -p, mᵢ_2′)

    return deltas * p_independent * mapreduce(p_dependent, +, -1:1)
end

function nuclear_spin_rotation(k::Int, bra::State, ket::State)::Float64
    N, mₙ, I, mᵢ = bra.N, bra.mₙ, bra.I[k], bra.mᵢ[k]
    N′, mₙ′, I′, mᵢ′ = ket.N, ket.mₙ, ket.I[k], ket.mᵢ[k]

    other = (k % 2) + 1 # States of other nucleus should Kronecker delta
    other_nucleus = δ(bra.I[other], ket.I[other]) * δ(bra.mᵢ[other], ket.mᵢ[other])

    deltas = δ(N, N′) * δ(I, I′) * other_nucleus
    p_independent = (-1)^(N + I - mₙ - mᵢ) * sqrt(N*(N+1)*(2*N + 1)) * sqrt(I*(I+1)*(2*I + 1))
    p_dependent(p) = (-1)^p * WignerSymbols.wigner3j(N, 1, N′, -mₙ, p, mₙ′) * WignerSymbols.wigner3j(I, 1, I′, -mᵢ, -p, -mᵢ′)
    
    return deltas * p_independent * mapreduce(p_dependent, +, -1:1)
end

function h_quadrupole(N_max::Int)
    # Neyenhuis PRL
    eqQ_1 = 0.45 # K, MHz
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

function h_nuclear_spin_spin(N_max::Int)
    c4 = -2030.4e-6 # MHz, Aldegunde PRA 78, 033434 (2008)

    elts::Int = (N_max + 1)^2 * N_Hyperfine
    H = zeros(elts, elts)
    for i = 1:elts
        for j = i:elts
            ket = index_to_state(i)
            bra = index_to_state(j)
            H[i, j] = c4 * nuclear_spin_spin(bra, ket)
        end
    end
    return Hermitian(H)
end

function h_nuclear_spin_rotation(N_max::Int)
    # MHz, from Aldegunde PRA
    cK = -24.1e-6
    cRb = 420.1e-6
    prefactors = [cK, cRb]

    elts::Int = (N_max + 1)^2 * N_Hyperfine
    H = zeros(elts, elts)
    for i = 1:elts
        for j = i:elts
            ket = index_to_state(i)
            bra = index_to_state(j)
            H[i, j] = dot(prefactors,[nuclear_spin_rotation(k, bra, ket) for k in 1:2])
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
    B_rot = 1113.9514 # Neyenhuis PRL

    stark = B_rot * h_rot(N_max, ε)
    hf = h_quadrupole(N_max) + h_nuclear_spin_spin(N_max) + h_nuclear_spin_rotation(N_max)
    zeeman = h_zeeman(N_max, B_field)
    return stark + hf + zeeman
end

function order_by_overlap_with(s::State, eigenvectors::Matrix)
    i = state_to_index(s)
    @assert i < size(eigenvectors, 1)
    return sortslices(eigenvectors, dims=2, lt=(x,y)->isless(abs2(x[i]), abs2(y[i])), rev=true)
end

# Returns tuple (overlap, )
function max_overlap_with(s::State, eigenvectors::Matrix)
    i = state_to_index(s)
    n_states = size(eigenvectors, 1)
    @assert i < n_states

    findmax(
        map(x -> abs2(x[i]), eachcol(eigenvectors))
    )
end

end # module
