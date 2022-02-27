module MoleculeSpectrum

using WignerSymbols
using HalfIntegers
using LinearAlgebra
using SparseArrays
using StaticArrays

using Test

module Constants
    const μN = 7.622593285e-4 # MHz/G
    const DToSI = 3.33564e-30
    const h = 6.62607004e-34
    const DVcm⁻¹ToMHz = (DToSI / h) * 1e-4
end # module

struct ZeemanParameters
    "Rotational g factor"
    gᵣ::Float64
    "Nuclear g factor"
    gᵢ::SVector{2, Float64}
    "Nuclear shielding factor"
    σᵢ::SVector{2, Float64}
end

struct NuclearParameters
    "Nuclear electric quadrupole (MHz)"
    eqQᵢ::SVector{2, Float64}
    "Nuclear spin-rotation interaction (MHz)"
    cᵢ::SVector{2, Float64}
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
    I::SVector{2, HalfInt}
    "Zeeman parameters"
    zeeman::ZeemanParameters
    "Nuclear Parameters"
    nuclear::NuclearParameters
    "Molecular polarizability at the trapping wavelength"
    α::Polarizability
end

const KRb_Zeeman = ZeemanParameters(0.014, [-0.324, 1.834], [1321e-6, 3469e-6])
const KRb_Polarizability = Polarizability(10.0e-5, 3.3e-5)

const KRb_Nuclear_Neyenhuis = NuclearParameters([0.45, -1.308], [-24.1e-6, 420.1e-6], -2030.4e-6)
const KRb_Nuclear_Ospelkaus = NuclearParameters([0.45, -1.41], [-24.1e-6, 420.1e-6], -2030.4e-6)

const KRb_Parameters_Neyenhuis = MolecularParameters(0.574, 1113.9514, [HalfInt(4), HalfInt(3/2)], KRb_Zeeman, KRb_Nuclear_Neyenhuis, KRb_Polarizability)
const KRb_Parameters_Ospelkaus = MolecularParameters(0.574, 1113.950, [HalfInt(4), HalfInt(3/2)], KRb_Zeeman, KRb_Nuclear_Ospelkaus, KRb_Polarizability)

const DEFAULT_MOLECULAR_PARAMETERS = KRb_Parameters_Neyenhuis

struct SphericalVector
    "Magnitude"
    magnitude::Float64
    "Polar angle (rad)"
    θ::Float64
    "Azimuthal angle (rad)"
    φ::Float64

    function SphericalVector(magnitude, θ, φ)
        if magnitude < 0
            error("Magnitude must be nonnegative")
        elseif θ < 0 || θ > π
            error("Polar angle must be in [0, π]")
        else
            φr = rem2pi(φ, RoundDown)
            if φr != φ
                @warn "φ was provided outside of [0, 2π)"
            end
            return new(magnitude, θ, φr)
        end
    end
end

function Base.:-(sv::SphericalVector)
    if sv.magnitude == 0
        return sv
    else
        return SphericalVector(sv.magnitude, π - sv.θ, sv.φ + π)
    end
end

VectorX(magnitude) = SphericalVector(magnitude, π/2, 0)
VectorY(magnitude) = SphericalVector(magnitude, π/2, π/2)
VectorZ(magnitude) = SphericalVector(magnitude, 0, 0)

struct SphericalUnitVector
    "Polar angle (rad)"
    θ::Float64
    "Azimuthal angle (rad)"
    φ::Float64

    function SphericalUnitVector(θ, φ)
        if θ < 0 || θ > π
            error("Polar angle must be in [0, π]")
        else
            φr = rem2pi(φ, RoundDown)
            if φr != φ
                @warn "φ was provided outside of [0, 2π)"
            end
            return new(θ, φr)
        end
    end

    SphericalUnitVector(v::SphericalVector) = new(v.θ, v.φ)
end

Base.:-(uv::SphericalUnitVector) = SphericalUnitVector(π - uv.θ, uv.φ + π)

UnitVectorX() = SphericalUnitVector(π/2, 0)
UnitVectorY() = SphericalUnitVector(π/2, π/2)
UnitVectorZ() = SphericalUnitVector(0, 0)

function T⁽¹⁾(v::SphericalUnitVector)::SVector{3, ComplexF64}
    θ = v.θ
    φ = v.φ
    
    x = sin(θ) * cos(φ)
    y = sin(θ) * sin(φ)
    z = cos(θ)

    T11 = -(1/sqrt(2)) * (x + im*y)
    T10 = z

    return SVector(-conj(T11), T10, T11)
end

function T⁽²⁾(v::SphericalUnitVector)::SVector{5, ComplexF64}
    θ = v.θ
    φ = v.φ

    x = sin(θ) * cos(φ)
    y = sin(θ) * sin(φ)
    z = cos(θ)

    T20 = (2*z^2 - x^2 - y^2) / sqrt(6)
    T21 = -(1/2)*(x*z + z*x + im * (y*z + z*y))
    T22 = (1/2)*(x*x - y*y + im*(x*y + y*x))

    return SVector(conj(T22), -conj(T21), T20, T21, T22)
end

function get_tensor_component(p::Int, tensor::Vector{ComplexF64})
    rank::Int = (length(tensor)-1) // 2 # Length should be 2*k + 1
    return tensor[1 + (p + rank)]
end

struct ExternalFields
    "Magnetic field (G)"
    B::SphericalVector
    "Electric field (V/cm)"
    E::SphericalVector
    "Laser fields (W/cm^2)"
    Optical::Vector{SphericalVector}

    ExternalFields(B::SphericalVector, E::SphericalVector, Optical) = new(B, E, Optical)
end

ExternalFields(B::Float64, E::Float64) = ExternalFields(VectorZ(B), VectorZ(E), [])
const DEFAULT_FIELDS = ExternalFields(545.9, 0.0)

struct State
    N::Int
    mₙ::Int
    I::SVector{2, HalfIntegers.HalfInt} # [K, Rb]
    mᵢ::SVector{2, HalfIntegers.HalfInt}
end

State(N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂) = State(N, mₙ, SVector(I₁, I₂), SVector(mᵢ₁, mᵢ₂))
State(N, mₙ, mᵢ₁::Number, mᵢ₂::Number) = State(N, mₙ, DEFAULT_MOLECULAR_PARAMETERS.I, [HalfIntegers.HalfInt(mᵢ₁) HalfIntegers.HalfInt(mᵢ₂)])

n_hyperfine(I::HalfIntegers.HalfInt) = 2 * I + 1
n_hyperfine(s::State) = mapreduce(n_hyperfine, *, s.I)

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

function get_energy(s::State, energies::Vector, eigenvectors::Matrix)
    return energies[max_overlap_with(s, eigenvectors)[2]]
end

function get_energy_difference(g::State, e::State, energies::Vector, eigenvectors::Matrix)
    return mapreduce(x -> get_energy(x, energies, eigenvectors), -, [e, g])
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

# Note: Need to enforce Kronecker deltas at the end!!!
# This is not checked inside of this function, because this tensor could be dotted with another operator,
# e.g., N ⋅ Iₖ, which allows mᵢ to change.
#
# Any Kronecker deltas need to be enforced at the scope where the dot product is taken.
function T⁽¹⁾N(p::Int, bra::State, ket::State)::ComplexF64
    N, mₙ = bra.N, bra.mₙ
    N′, mₙ′ = ket.N, ket.mₙ

    if δ(N, N′) && (mₙ′ - mₙ + p == 0)
        return (-1)^(N - mₙ) * sqrt(N*(N+1)*(2N+1)) * WignerSymbols.wigner3j(N, 1, N, -mₙ, p, mₙ′)
    else
        return 0
    end
end

# Note: Need to enforce Kronecker deltas at the end!!!
# This is not checked inside of this function, because this tensor could be dotted with another operator,
# e.g., N ⋅ Iₖ, which allows mₙ to change.
#
# Any Kronecker deltas need to be enforced at the scope where the dot product is taken.
function T⁽¹⁾Iₖ(p::Int, k::Int, bra::State, ket::State)::ComplexF64
    I, mᵢ = bra.I[k], bra.mᵢ[k]
    I′, mᵢ′ = ket.I[k], ket.mᵢ[k]

    if δ(I, I′) && (mᵢ′ - mᵢ + p == 0)
        return (-1)^(I - mᵢ) * sqrt(I*(I+1)*(2I+1)) * WignerSymbols.wigner3j(I, 1, I, -mᵢ, p, mᵢ′)
    else
        return 0
    end
end

rotation_matrix_element(bra::State, ket::State)::Float64 = ket.N * (ket.N + 1) * δ(bra, ket)

function h_rotation(basis::Vector{State})
    return Diagonal(map(x -> rotation_matrix_element(x, x), basis))
end

function dipole_matrix_element(p::Int, bra::State, ket::State)::ComplexF64
    N, mₙ, mᵢ = bra.N, bra.mₙ, bra.mᵢ
    N′, mₙ′, mᵢ′ = ket.N, ket.mₙ, ket.mᵢ

    if δ(mᵢ, mᵢ′) && (-mₙ + p + mₙ′ == 0)
        return (-1)^mₙ * sqrt((2N+1)*(2N′+1)) * WignerSymbols.wigner3j(N, 1, N′, -mₙ, p, mₙ′) * WignerSymbols.wigner3j(N, 1, N′, 0, 0, 0)
    else
        return 0
    end
end

function h_dipole_components(basis::Vector{State}, p::Int)
    elts::Int = length(basis)
    H = spzeros(ComplexF64, elts, elts)
    for i = 1:elts
        ket = basis[i]
        for j = i:elts
            bra = basis[j]
            H[i, j] = dipole_matrix_element(p, bra, ket)
        end
    end
    return Hermitian(H)
end

function h_dipole(basis::Vector{State})
    return [h_dipole_components(basis, p) for p in -1:1]
end

function nuclear_quadrupole(k::Int, bra::State, ket::State)::ComplexF64
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

function h_quadrupole(k::Int, basis::Vector{State})
    elts::Int = length(basis)
    H = spzeros(ComplexF64, elts, elts)

    for i = 1:elts
        ket = basis[i]
        for j = i:elts
            bra = basis[j]
            H[i, j] = nuclear_quadrupole(k, bra, ket)
        end
    end
    return Hermitian(H)
end

function nuclear_spin_spin(bra::State, ket::State)::ComplexF64
    if δ(bra.N, ket.N) && δ(bra.mₙ, ket.mₙ) && 
        (ket.mᵢ[1] - bra.mᵢ[1] == bra.mᵢ[2] - ket.mᵢ[2])

        return reduce(+, [(-1)^p * T⁽¹⁾Iₖ(p, 1, bra, ket) * T⁽¹⁾Iₖ(-p, 2, bra, ket) for p in -1:1])
    else
        return 0.0
    end
end

function h_nuclear_spin_spin(basis::Vector{State})
    elts::Int = length(basis)
    H = spzeros(ComplexF64, elts, elts)
    for i = 1:elts
        ket = basis[i]
        for j = i:elts
            bra = basis[j]
            H[i, j] = nuclear_spin_spin(bra, ket)
        end
    end
    return Hermitian(H)
end

function nuclear_spin_rotation(k::Int, bra::State, ket::State)::ComplexF64
    other = (k % 2) + 1 # States of other nucleus should Kronecker delta
    other_nucleus = δ(bra.I[other], ket.I[other]) * δ(bra.mᵢ[other], ket.mᵢ[other])

    if other_nucleus && (ket.mₙ - bra.mₙ == bra.mᵢ[k] - ket.mᵢ[k])
        return reduce(+, [(-1)^p * T⁽¹⁾N(p, bra, ket) * T⁽¹⁾Iₖ(-p, k, bra, ket) for p in -1:1])
    else
        return 0
    end
end

function h_nuclear_spin_rotation(k::Int, basis::Vector{State})
    elts::Int = length(basis)
    H = spzeros(ComplexF64, elts, elts)
    for i = 1:elts
        ket = basis[i]
        for j = i:elts
            bra = basis[j]
            H[i, j] = nuclear_spin_rotation(k, bra, ket)
        end
    end
    return Hermitian(H)
end

function h_zeeman_rotation_components(basis::Vector{State}, p::Int)
    elts::Int = length(basis)
    H = spzeros(ComplexF64, elts, elts)
    for i = 1:elts
        ket = basis[i]
        for j = i:elts
            bra = basis[j]

            if δ(bra.I, ket.I) && δ(bra.mᵢ, ket.mᵢ)
                H[i, j] = T⁽¹⁾N(p, bra, ket)
            end
        end
    end
    return Hermitian(H)
end

function h_zeeman_rotation(basis::Vector{State})
    return [h_zeeman_rotation_components(basis, p) for p in -1:1]
end

function h_zeeman_nuclear_components(basis::Vector{State}, k::Int, p::Int)
    elts::Int = length(basis)
    H = spzeros(ComplexF64, elts, elts)
    for i = 1:elts
        ket = basis[i]
        for j = i:elts
            bra = basis[j]

            other = 1 + (k % 2)
            rotation = δ(bra.N, ket.N) && δ(bra.mₙ, ket.mₙ)
            I_other = δ(bra.I[other], ket.I[other]) && δ(bra.mᵢ[other], ket.mᵢ[other])

            if rotation && I_other
                H[i, j] = T⁽¹⁾Iₖ(p, k, bra, ket)
            end
        end
    end
    return Hermitian(H)
end

function h_zeeman_nuclear(basis::Vector{State}, k::Int)
    return [h_zeeman_nuclear_components(basis, k, p) for p in -1:1]
end

scalar_polarizability(bra::State, ket::State) = δ(bra, ket)

function tensor_polarizability(bra::State, ket::State, p::Int) 
    deltas = δ(bra.N, ket.N) * δ(bra.I, ket.I) * δ(bra.mᵢ, ket.mᵢ)
    N, mₙ = bra.N, bra.mₙ
    mₙ′ = ket.mₙ

    if deltas && (N > 0) && (mₙ′- mₙ - p == 0)
        return sqrt(6) * (-1)^(mₙ) * (2*N + 1) * WignerJ2J(N, 0) * WignerSymbols.wigner3j(N, 2, N, -mₙ, -p, mₙ′)
    else
        return 0
    end
end

function h_ac_scalar(basis::Vector{State})
    elts::Int = length(basis)
    H = Diagonal(ones(elts))
    return Hermitian(H)
end

function h_ac_tensor_components(basis::Vector{State}, p::Int)
    elts::Int = length(basis)
    H = spzeros(ComplexF64, elts, elts)
    for i = 1:elts
        ket = basis[i]
        for j = i:elts
            bra = basis[j]
            H[i, j] = tensor_polarizability(bra, ket, p)
        end
    end
    return Hermitian(H)
end

function h_ac_tensor(basis::Vector{State})
    return [h_ac_tensor_components(basis, p) for p in -2:2]
end

function tensor_dot(a, b)
    @assert size(a, 1) == size(b, 1)
    @assert isodd(size(a, 1))

    mapreduce(p -> conj(a[p]) .* b[p], +, eachindex(a))
end

struct HamiltonianParts
    rotation
    dipole::SVector{3}
    hyperfine
    zeeman::SVector{3}
    ac_scalar
    ac_tensor::SVector{5}
end

function make_hamiltonian_parts(molecular_parameters::MolecularParameters, N_max::Int)::HamiltonianParts
    basis = generate_basis(molecular_parameters, N_max)

    rotation = molecular_parameters.Bᵣ * h_rotation(basis)
    dipole = (-1) * molecular_parameters.dₚ * Constants.DVcm⁻¹ToMHz * h_dipole(basis)
    hyperfine = make_hyperfine(molecular_parameters, basis)
    zeeman = make_zeeman(molecular_parameters, basis)
    (ac_scalar, ac_tensor) = make_ac(molecular_parameters, basis)

    return HamiltonianParts(rotation, dipole, hyperfine, zeeman, ac_scalar, ac_tensor)
end

function make_hyperfine(molecular_parameters::MolecularParameters, basis::Vector{State})
    eqQᵢ = molecular_parameters.nuclear.eqQᵢ
    quadrupole = mapreduce(k -> (1/4) * eqQᵢ[k] * h_quadrupole(k, basis), +, 1:2)

    cᵢ = molecular_parameters.nuclear.cᵢ
    spin_rotation = mapreduce(k -> cᵢ[k] * h_nuclear_spin_rotation(k, basis), +, 1:2)

    c₄ = molecular_parameters.nuclear.c₄
    spin_spin = c₄ * h_nuclear_spin_spin(basis)

    return quadrupole + spin_rotation + spin_spin
end

function make_zeeman(molecular_parameters::MolecularParameters, basis::Vector{State})
    μN = Constants.μN
    gᵣ = molecular_parameters.zeeman.gᵣ

    zeeman_n = gᵣ * μN * h_zeeman_rotation(basis)

    gᵢ = molecular_parameters.zeeman.gᵢ
    σᵢ = molecular_parameters.zeeman.σᵢ
    prefactors = [gᵢ[k] * μN * (1 - σᵢ[k]) for k in 1:2]
    zeeman_i = mapreduce(k -> prefactors[k] * h_zeeman_nuclear(basis, k), +, 1:2)

    return (-1) * (zeeman_n + zeeman_i)
end

function make_ac(molecular_parameters::MolecularParameters, basis::Vector{State})
    α_par = molecular_parameters.α.α_par
    α_perp = molecular_parameters.α.α_perp

    ac_scalar = (-1) * ((α_par + 2α_perp) / 3) * h_ac_scalar(basis)
    ac_tensor = (-1) * ((α_par - α_perp) / 3) * h_ac_tensor(basis)

    return (ac_scalar, ac_tensor)
end

function hamiltonian(parts::HamiltonianParts, external_fields::ExternalFields)
    h = parts.rotation + parts.hyperfine

    E = external_fields.E.magnitude
    E_n = SphericalUnitVector(external_fields.E)
    h += E * tensor_dot(T⁽¹⁾(E_n), parts.dipole)

    B = external_fields.B.magnitude
    B_n = SphericalUnitVector(external_fields.B)
    h += B * tensor_dot(T⁽¹⁾(B_n), parts.zeeman)

    for beam in external_fields.Optical
        I_laser = beam.magnitude
        ϵ = SphericalUnitVector(beam)
        h += parts.ac_scalar * I_laser
        h += tensor_dot(T⁽²⁾(ϵ), parts.ac_tensor) * I_laser
    end

    return Array(Hermitian(h))
end

@testset "Reproduces Ospelkaus et al., PRL 104, 030402 (2010)" begin
    N_max = 5
    tolerance = 0.0011 # MHz

    fields = ExternalFields(545.9, 0.0)
    parts = make_hamiltonian_parts(KRb_Parameters_Ospelkaus, N_max)
    h = hamiltonian(parts, fields)
    energies = eigvals(h)
    states = eigvecs(h)

    comparisons = [
        ((0, 0, -4, 1/2), (1, 1, -4, 1/2),  2227.835),
        ((0, 0, -4, 1/2), (1, 0, -4, 1/2),  2228.119),
        ((0, 0, -4, 1/2), (1, -1, -4, 1/2), 2227.776),
        ((0, 0, -4, 1/2), (1, 0, -4, 3/2),  2227.008),
        ((0, 0, -4, 1/2), (1, -1, -4, 3/2), 2227.128),
        ((0, 0, -4, 1/2), (1, 0, -3, 1/2),  2228.225),
        ((0, 0, -4, 1/2), (1, 1, -4, -1/2), 2228.593),
        ((0, 0, -4, 1/2), (1, 0, -4, -1/2), 2228.805),
        ((0, 0, -4, 3/2), (1, 0, -4, 3/2),  2227.761),
        ((0, 0, -3, 1/2), (1, 0, -3, 1/2),  2228.091),
    ]

    for c in comparisons
        (g, e) = map(i -> State(c[i]...), 1:2)
        transition = get_energy_difference(g, e, energies, states)
        expected = c[3]

        @test abs(transition - expected) < tolerance
    end
end

@testset verbose=true "Reproduces Neyenhuis et al., PRL 109, 230403 (2012)" begin
    N_max = 5
    parts = make_hamiltonian_parts(KRb_Parameters_Neyenhuis, N_max)
    B = 545.9

    fields = ExternalFields(B, 0.0)
    h = hamiltonian(parts, fields)
    energies = eigvals(h)
    states = eigvecs(h)

    @testset "No optical fields" begin
        tolerance = 0.005 # MHz

        comparisons = [
            ((0, 0, -4, 1/2), (1, 1, -4, 1/2),  2227.842),
            ((0, 0, -4, 1/2), (1, 0, -4, 1/2),  2228.110),
            ((0, 0, -4, 1/2), (1, -1, -4, 1/2), 2227.784),
        ]

        for c in comparisons
            (g, e) = map(i -> State(c[i]...), 1:2)
            transition = get_energy_difference(g, e, energies, states)
            expected = c[3]

            @test abs(transition - expected) < tolerance
        end
    end

    @testset "α(θ)" begin
        I_light = 2350.
        tolerance = 0.035 # Relative tolerance in α

        # These are generated by diagonalizing the simplified Hamiltonian H from the
        # main text and plugging in the experimental values for α_parallel and α_perpendicular.
        # A few % tolerance seems reasonable given that these values come from the simplified Hamiltonian.
        #
        # Note: The supplement shows the full calculation (which should match this, in principle).
        comparisons = [
            # θ, m =  -1,      +1,      0
            (0.0,   46.4e-6, 46.4e-6, 73.2e-6),
            (10.0,  46.9e-6, 46.9e-6, 72.2e-6),
            (20.0,  48.4e-6, 48.2e-6, 69.4e-6),
            (30.0,  50.9e-6, 49.8e-6, 65.3e-6),
            (40.0,  54.0e-6, 51.2e-6, 60.7e-6),
            (50.0,  57.5e-6, 52.4e-6, 56.2e-6),
            (60.0,  60.7e-6, 53.1e-6, 52.2e-6),
            (70.0,  63.3e-6, 53.6e-6, 49.0e-6),
            (80.0,  65.0e-6, 53.9e-6, 47.1e-6),
            (90.0,  65.6e-6, 54.0e-6, 46.4e-6),
        ]
        α00 = 55.3e-6

        states_to_check = [
            (0, 0, -4, 1/2),
            (1, -1, -4, 1/2),
            (1, 1, -4, 1/2),
            (1, 0, -4, 1/2),
        ]

        es = map(x -> get_energy(State(x...), energies, states), states_to_check)

        for c in comparisons
            θ = c[1] * π/180
            optical = SphericalVector(I_light, θ, 0.0)
            fields_with_light = ExternalFields(VectorZ(B), VectorZ(0.0), [optical])

            h_light = hamiltonian(parts, fields_with_light)
            energies_light = eigvals(h_light)
            states_light = eigvecs(h_light)

            es_light = map(
                x -> get_energy(State(x...), energies_light, states_light),
                states_to_check
            )

            αs = map(
                k -> -(es_light[k] - es[k]) / I_light,
                eachindex(es)
            )
            expected = [α00, c[2:end]...]
            
            errors = map(
                k -> abs(1 - (expected[k] / αs[k])),
                eachindex(expected)
            )

            # println(errors) # Uncomment to show the errors on every iteration

            @test all(errors .< tolerance)
        end
    end
end

@testset "Trivial angular dependence with one field" begin
    N_max = 5
    B = 545.9
    E = 1000.0
    parts = make_hamiltonian_parts(KRb_Parameters_Neyenhuis, N_max)

    b_z = ExternalFields(B, 0.0)
    b_x = ExternalFields(VectorX(B), VectorZ(0.0), [])
    b_y = ExternalFields(VectorY(B), VectorZ(0.0), [])
    b_xz = ExternalFields(SphericalVector(B, π/4, 0.0), VectorZ(0.0), [])
    b_xyz = ExternalFields(SphericalVector(B, π/3, π/3), VectorZ(0.0), [])


    e_z = ExternalFields(0.0, E)
    e_x = ExternalFields(VectorZ(0.0), VectorX(E), [])
    e_y = ExternalFields(VectorZ(0.0), VectorY(E), [])
    e_xz = ExternalFields(VectorZ(0.0), SphericalVector(E, π/4, 0.0), [])
    e_xyz = ExternalFields(VectorZ(0.0), SphericalVector(E, π/3, π/3), [])

    for fields in [(b_z, (b_x, b_y, b_xz, b_xyz)), (e_z, (e_x, e_y, e_xz, e_xyz))]
        h_z = hamiltonian(parts, fields[1])
        energies = eigvals(h_z)

        for f in fields[2]
            h = hamiltonian(parts, f)
            es = eigvals(h)
    
            @test es ≈ energies
        end
    end
end

end # module
