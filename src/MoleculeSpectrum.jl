module MoleculeSpectrum

using WignerSymbols
using HalfIntegers
using LinearAlgebra
using SparseArrays

# const I_K = HalfInt(4) # K
# const I_Rb = HalfInt(3/2) # Rb

struct ZeemanParameters
    "Rotational g factor"
    gᵣ::Float64
    "Nuclear g factor"
    gᵢ::Array{Float64, 2}
    "Nuclear shielding factor"
    σᵢ::Array{Float64, 2}
end

struct NuclearParameters
    "Nuclear electric quadrupole (MHz)"
    eqQᵢ::Array{Float64, 2}
    "Nuclear spin-rotation interaction (MHz)"
    cᵢ::Array{Float64, 2}
    "Nuclear spin-spin scalar interaction (MHz)"
    c₄::Float64
end

struct Polarizability
    "Parallel ac polarizability (MHz / (W / cm^2))"
    α_par::Float64
    "Perpendicular ac polarizability (MHz / (W / cm^2))"
    α_perp::Float64
end

struct MolecularParameters
    "Permanent dipole moment (Debye)"
    dₚ::Float64
    "Rotational constant (MHz)"
    Bᵣ::Float64
    "Nuclear angular momenta"
    I::Array{HalfInt, 2}
    "Zeeman parameters"
    zeeman::ZeemanParameters
    "Nuclear Parameters"
    nuclear::NuclearParameters
    "Molecular polarizability at the trapping wavelength"
    α::Polarizability
end

const KRb_Zeeman = ZeemanParameters(0.014, [-0.324 1.834], [1321e-6 3469e-6])
const KRb_Polarizability = Polarizability(10.0e-5, 3.3e-5)

const KRb_Nuclear_Neyenhuis = NuclearParameters([0.45 -1.308], [-24.1e-6 420.1e-6], -2030.4e-6)
const KRb_Nuclear_Ospelkaus = NuclearParameters([0.45 -1.41], [-24.1e-6 420.1e-6], -2030.4e-6)

const KRb_Parameters_Neyenhuis = MolecularParameters(0.574, 1113.9514, [HalfInt(4) HalfInt(3/2)], KRb_Zeeman, KRb_Nuclear_Neyenhuis, KRb_Polarizability)
const KRb_Parameters_Ospelkaus = MolecularParameters(0.574, 1113.9514, [HalfInt(4) HalfInt(3/2)], KRb_Zeeman, KRb_Nuclear_Ospelkaus, KRb_Polarizability)

const DEFAULT_MOLECULAR_PARAMETERS = KRb_Parameters_Neyenhuis

struct SphericalVector
    magnitude::Float64
    θ::Float64
    φ::Float64
end

struct ExternalFields
    B::SphericalVector
    E::SphericalVector
    Optical::Vector{SphericalVector}
end

struct State
    N::Int
    mₙ::Int
    I::Array{HalfIntegers.HalfInt, 2} # [K, Rb]
    mᵢ::Array{HalfIntegers.HalfInt, 2}
end

State(N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂) = State(N, mₙ, [I₁ I₂], [mᵢ₁ mᵢ₂])
State(N, mₙ, mᵢ₁::Number, mᵢ₂::Number) = State(N, mₙ, DEFAULT_MOLECULAR_PARAMETERS.I, [HalfIntegers.HalfInt(mᵢ₁) HalfIntegers.HalfInt(mᵢ₂)])

n_hyperfine(I::HalfIntegers.HalfInt) = 2 * I + 1
n_hyperfine(s::State) = mapreduce(n_hyperfine, *, s.I)

# const N_Hyperfine = n_hyperfine(I_K) * n_hyperfine(I_Rb)

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

index_to_state(i::Int) = index_to_state(i, DEFAULT_MOLECULAR_PARAMETERS.I[1], DEFAULT_MOLECULAR_PARAMETERS.I[2])

# Todo: test for state_to_index(index_to_state(x)) == x
function state_to_index(s::State)::Int
    rotation = (s.N^2 + 1) + (s.N + s.mₙ)
    hyperfine = (s.I[1] + s.mᵢ[1]) * n_hyperfine(s.I[2]) + (s.I[2] + s.mᵢ[2])

    N_Hyperfine = n_hyperfine(s)
    return 1 + (rotation - 1) * N_Hyperfine + hyperfine
end

function order_by_overlap_with(s::State, eigenvectors::Matrix)
    i = state_to_index(s)
    @assert i < size(eigenvectors, 1)
    return sortslices(eigenvectors, dims=2, lt=(x,y)->isless(abs2(x[i]), abs2(y[i])), rev=true)
end

# Returns tuple (overlap, index)
function max_overlap_with(s::State, eigenvectors::Matrix)
    i = state_to_index(s)
    n_states = size(eigenvectors, 1)
    @assert i < n_states

    findmax(
        map(x -> abs2(x[i]), eachcol(eigenvectors))
    )
end

function energy_difference(g::State, e::State, energies::Vector, eigenvectors::Matrix)
    indices = map(x -> max_overlap_with(x, eigenvectors)[2], [e, g])
    return mapreduce(x -> energies[x], -, indices)
end

function generate_basis(molecular_parameters::MolecularParameters, N_max::Int)
    n_elts::Int = (N_max + 1)^2 * mapreduce(n_hyperfine, *, molecular_parameters.I)
    return map(index_to_state, 1:n_elts)
end

δ(i, j) = ==(i, j)
δ(i::State, j::State) = δ(i.N, j.N) * δ(i.mₙ, j.mₙ) * δ(i.I, j.I) * δ(i.mᵢ, j.mᵢ)

# Same as WignerSymbols.wigner3j(j, 1, j, -m, 0, m)
WignerJ1J(j, m) = (-1)^(j-m) * m / sqrt(j*(j+1)*(2j+1))
# Same as WignerSymbols.wigner3j(j, 2, j, -m, 0, m)
WignerJ2J(j, m) = (-1)^(j-m) * 2 * (3m^2 - j*(j + 1)) / sqrt((2j - 1)*(2j)*(2j + 1)*(2j + 2)*(2j + 3))

function dipole_matrix_element(p::Int, bra::State, ket::State)::Float64
    N, mₙ, mᵢ = bra.N, bra.mₙ, bra.mᵢ
    N′, mₙ′, mᵢ′ = ket.N, ket.mₙ, ket.mᵢ

    if δ(mᵢ, mᵢ′) && (-mₙ + p + mₙ′ == 0)
        return (-1)^mₙ * sqrt((2N+1)*(2N′+1)) * WignerSymbols.wigner3j(N, 1, N′, -mₙ, p, mₙ′) * WignerSymbols.wigner3j(N, 1, N′, 0, 0, 0)
    else
        return 0
    end
end

rotation_matrix_element(bra::State, ket::State)::Float64 = ket.N * (ket.N + 1) * δ(bra, ket) 

function h_rot(basis::Vector{State}, ε::Float64)
    elts::Int = length(basis)
    H = spzeros(elts, elts)
    for i = 1:elts
        ket = basis[i]
        for j = i:elts
            bra = basis[j]
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

    if other_nucleus &&
        δ(I,I′) && (δ(N,N′) || (abs(N-N′) == 2)) &&
        (mₙ′- mₙ == mᵢ - mᵢ′) &&
        (abs(mₙ′-mₙ) <= 2) && (abs(mᵢ′-mᵢ) <= 2)

        # Brown and Carrington pg. 477
        p_independent =  (-1)^(I - mᵢ - mₙ) * sqrt(2*N + 1) * sqrt(2*N′ + 1) * WignerSymbols.wigner3j(N, 2, N′, 0, 0, 0) / WignerJ2J(I, I)
        p_dependent(p) = (-1)^(p) * WignerSymbols.wigner3j(N, 2, N′, -mₙ, p, mₙ′) * WignerSymbols.wigner3j(I, 2, I, -mᵢ, -p, mᵢ′)
        return p_independent * mapreduce(p_dependent, +, -2:2)
    else
        return 0.0
    end
end

function nuclear_spin_spin(bra::State, ket::State)::Float64
    N, mₙ, (I_1, I_2), (mᵢ_1, mᵢ_2) = bra.N, bra.mₙ, bra.I, bra.mᵢ
    N′, mₙ′, (I_1′, I_2′), (mᵢ_1′, mᵢ_2′) = ket.N, ket.mₙ, ket.I, ket.mᵢ

    deltas = δ(N, N′) * δ(mₙ, mₙ′) * δ(I_1, I_1′) * δ(I_2, I_2′)

    if deltas
        p_independent = (-1)^(I_1 + I_2 - mᵢ_1 - mᵢ_2) * sqrt(I_1 * (I_1 + 1) * (2*I_1 + 1)) * sqrt(I_2 * (I_2 + 1) * (2*I_2 + 1))
        p_dependent(p) = (-1)^p * WignerSymbols.wigner3j(I_1, 1, I_1, -mᵢ_1, p, mᵢ_1′) * WignerSymbols.wigner3j(I_2, 1, I_2, -mᵢ_2, -p, mᵢ_2′)
        return p_independent * mapreduce(p_dependent, +, -1:1)
    else
        return 0.0
    end
end

function nuclear_spin_rotation(k::Int, bra::State, ket::State)::Float64
    N, mₙ, I, mᵢ = bra.N, bra.mₙ, bra.I[k], bra.mᵢ[k]
    N′, mₙ′, I′, mᵢ′ = ket.N, ket.mₙ, ket.I[k], ket.mᵢ[k]

    other = (k % 2) + 1 # States of other nucleus should Kronecker delta
    other_nucleus = δ(bra.I[other], ket.I[other]) * δ(bra.mᵢ[other], ket.mᵢ[other])

    deltas = δ(N, N′) * δ(I, I′) * other_nucleus

    if deltas && (mₙ′ - mₙ == mᵢ - mᵢ′)
        p_independent = (-1)^(N + I - mₙ - mᵢ) * sqrt(N*(N+1)*(2*N + 1)) * sqrt(I*(I+1)*(2*I + 1))
        p_dependent(p) = (-1)^p * WignerSymbols.wigner3j(N, 1, N, -mₙ, p, mₙ′) * WignerSymbols.wigner3j(I, 1, I, -mᵢ, -p, mᵢ′)
        return p_independent * mapreduce(p_dependent, +, -1:1)
    else
        return 0
    end
end

function h_quadrupole(molecular_parameters::MolecularParameters, basis::Vector{State})
    (eqQ_1, eqQ_2) = molecular_parameters.nuclear.eqQᵢ

    # # Neyenhuis PRL
    # eqQ_1 = 0.45 # K, MHz
    # eqQ_2 = -1.308 # Rb, MHz

    # # Silke PRL
    # # eqQ_2 = -1.41 # Rb, MHz

    elts::Int = length(basis)
    H = spzeros(elts, elts)

    for i = 1:elts
        ket = basis[i]
        for j = i:elts
            bra = basis[j]
            H[i, j] = (1/4) * (eqQ_1 * nuclear_quadrupole(1, bra, ket) + eqQ_2 * nuclear_quadrupole(2, bra, ket))
        end
    end
    return Hermitian(H)
end

function h_nuclear_spin_spin(molecular_parameters::MolecularParameters, basis::Vector{State})
    c4 = molecular_parameters.nuclear.c₄

    # c4 = -2030.4e-6 # MHz, Aldegunde PRA 78, 033434 (2008)

    elts::Int = length(basis)
    H = spzeros(elts, elts)
    for i = 1:elts
        ket = basis[i]
        for j = i:elts
            bra = basis[j]
            H[i, j] = c4 * nuclear_spin_spin(bra, ket)
        end
    end
    return Hermitian(H)
end

function h_nuclear_spin_rotation(molecular_parameters::MolecularParameters, basis::Vector{State})
    # MHz, from Aldegunde PRA
    (c1, c2) = molecular_parameters.nuclear.cᵢ
    # cK = -24.1e-6
    # cRb = 420.1e-6

    elts::Int = length(basis)
    H = spzeros(elts, elts)
    for i = 1:elts
        ket = basis[i]
        for j = i:elts
            bra = basis[j]
            H[i, j] = c1 * nuclear_spin_rotation(1, bra, ket) + c2 * nuclear_spin_rotation(2, bra, ket)
        end
    end
    return Hermitian(H)
end

# no angle dependence for now
function h_zeeman(molecular_parameters::MolecularParameters, basis::Vector{State}, B_field::Float64)
    μ_N = 7.622593285e-4 # MHz/G

    g_r = molecular_parameters.zeeman.gᵣ
    (g_1, g_2) = molecular_parameters.zeeman.gᵢ
    (σ_1, σ_2) = molecular_parameters.zeeman.σᵢ

    # g_r = 0.014 # Aldegunde, PRA 78, 033434 (2008)
    # μ_N = 7.622593285e-4 # MHz/G

    # # Aldegunde, PRA 78, 033434 (2008)
    # g_1 = -0.324 # g_K
    # g_2 = 1.834 # g_Rb
    # σ_1 = 1321e-6
    # σ_2 = 3469e-6

    elts::Int = length(basis)
    H = spzeros(elts, elts)
    for i = 1:elts
        ket = basis[i]
        H[i, i] = -g_r * μ_N * B_field * ket.mₙ - μ_N * B_field * dot([g_1*(1-σ_1), g_2*(1-σ_2)], ket.mᵢ)
    end
    return H
end

function T2pol(θ, φ)
    x = sin(θ) * cos(φ)
    y = sin(θ) * sin(φ)
    z = cos(θ)

    T20 = (2*z^2 - x^2 - y^2) / sqrt(6)
    T21 = -(1/2)*(x*z + z*x + im * (y*z + z*y))
    T22 = (1/2)*(x*x - y*y + im*(x*y + y*x))

    return [conj(T22), -conj(T21), T20, T21, T22]
end

scalar_polarizability(bra::State, ket::State) = δ(bra, ket)
function tensor_polarizability(bra::State, ket::State, T2ϵϵ) 
    deltas = δ(bra.N, ket.N) * δ(bra.I, ket.I) * δ(bra.mᵢ, ket.mᵢ)
    N, mₙ = bra.N, bra.mₙ
    mₙ′ = ket.mₙ

    if deltas && (N > 0) && (abs(mₙ′-mₙ) <= 2)
        p_independent = sqrt(6) * (-1)^(mₙ) * (2*N + 1) * WignerJ2J(N, 0)
        p_dependent(p) = (-1)^p * T2ϵϵ[p+3] * WignerSymbols.wigner3j(N, 2, N, -mₙ, -p, mₙ′)
        return p_independent * mapreduce(p_dependent, +, -2:2)
    else
        return 0
    end
end

function h_ac(molecular_parameters::MolecularParameters, basis::Vector{State}, I_laser::Float64, θ_laser::Float64, φ_laser::Float64)
    α_par = molecular_parameters.α.α_par
    α_perp = molecular_parameters.α.α_perp

    # α_par = 10.0e-5 # MHz / (W / cm^2)
    # α_perp = 3.3e-5 # MHz / (W / cm^2)
    T2ϵϵ = T2pol(θ_laser, φ_laser)

    elts::Int = length(basis)
    H = spzeros(elts, elts)
    for i = 1:elts
        ket = basis[i]
        for j = i:elts
            bra = basis[j]
            scalar = ((α_par + 2 * α_perp) / 3) * scalar_polarizability(bra, ket)
            tensor = ((α_par - α_perp) / 3) * tensor_polarizability(bra, ket, T2ϵϵ)
            H[i, j] = -(scalar + tensor) * I_laser
        end
    end
    return Hermitian(H)
end


function h(molecular_parameters::MolecularParameters, N_max::Int, ε::Float64, B_field::Float64, I_laser::Float64, θ_laser::Float64, φ_laser::Float64)
    B_rot = molecular_parameters.Bᵣ

    # B_rot = 1113.9514 # Neyenhuis PRL
    # # B_rot = 1113.950 # Silke PRL

    basis = generate_basis(molecular_parameters, N_max)

    dc_stark = B_rot * h_rot(basis, ε)
    hf = h_quadrupole(molecular_parameters, basis) + h_nuclear_spin_spin(molecular_parameters, basis) + h_nuclear_spin_rotation(molecular_parameters, basis)
    zeeman = h_zeeman(molecular_parameters, basis, B_field)
    ac_stark = h_ac(molecular_parameters, basis, I_laser, θ_laser, φ_laser)
    return Array(Hermitian(dc_stark + hf + zeeman + ac_stark))
end



end # module
