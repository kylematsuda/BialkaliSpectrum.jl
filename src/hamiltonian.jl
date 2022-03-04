"""
    h_diagonal(basis, matrix_element)

Create a diagonal matrix where the `(i, i)`th element is `matrix_element(basis[i], basis[i])`.
"""
function h_diagonal(basis::Vector{State}, matrix_element::Function)
    return Hermitian(Diagonal(map(x -> matrix_element(x, x), basis)))
end

"""
    h_rank_0(basis, matrix_element)

Create a matrix where the `(i, j)`th element is `matrix_element(basis[j], basis[i])`.

Note that `matrix_element` must have the signature (State, State) -> ComplexF64.
"""
function h_rank_0(basis::Vector{State}, matrix_element::Function)
    elts::Int = length(basis)
    H = spzeros(ComplexF64, elts, elts)
    for i = 1:elts
        ket = basis[i]
        for j = i:elts
            bra = basis[j]
            H[i, j] = matrix_element(bra, ket)
        end
    end
    return Hermitian(H)
end

"""
    h_tensor_component(basis, matrix_element, p)

Create a matrix where the `(i, j)`th element is `matrix_element(p, basis[j], basis[i])`.

Note that `matrix_element` must have the signature (Int, State, State) -> ComplexF64.
The first argument is the tensor component `p`.
"""
function h_tensor_component(basis::Vector{State}, matrix_element::Function, p::Int)
    elts::Int = length(basis)
    H = spzeros(ComplexF64, elts, elts)
    for i = 1:elts
        ket = basis[i]
        for j = 1:elts
            bra = basis[j]
            H[i, j] = matrix_element(p, bra, ket)
        end
    end
    return H
end

"""
    h_rank_1(basis, matrix_element)

Create a vector of 3 matrices to contract against a rank 1 spherical tensor.

Note that `matrix_element` must have the signature (Int, State, State) -> ComplexF64.
The first argument is the tensor component `p` in `-1:1`.
"""
function h_rank_1(basis::Vector{State}, matrix_element::Function)
    return [h_tensor_component(basis, matrix_element, p) for p = -1:1]
end

"""
    h_rank_2(basis, matrix_element)

Create a vector of 5 matrices to contract against a rank 2 spherical tensor.

Note that `matrix_element` must have the signature (Int, State, State) -> ComplexF64.
The first argument is the tensor component `p` in `-2:2`.
"""
function h_rank_2(basis::Vector{State}, matrix_element::Function)
    return [h_tensor_component(basis, matrix_element, p) for p = -2:2]
end

h_rotation(basis) = h_diagonal(basis, rotation_matrix_element)
h_dipole(basis) = h_rank_1(basis, dipole_matrix_element)

h_quadrupole(basis, k) = h_rank_0(basis, (bra, ket) -> nuclear_quadrupole(k, bra, ket))
h_nuclear_spin_spin(basis) = h_rank_0(basis, nuclear_spin_spin)
h_nuclear_spin_rotation(basis, k) =
    h_rank_0(basis, (bra, ket) -> nuclear_spin_rotation(k, bra, ket))

h_zeeman_rotation(basis) = h_rank_1(basis, zeeman_rotation)
h_zeeman_nuclear(basis, k) =
    h_rank_1(basis, (p, bra, ket) -> zeeman_nuclear(k, p, bra, ket))

h_ac_scalar(basis) = h_diagonal(basis, scalar_polarizability)
h_ac_tensor(basis) = h_rank_2(basis, tensor_polarizability)

function make_hyperfine(molecular_parameters::MolecularParameters, basis::Vector{State})
    eqQᵢ = molecular_parameters.nuclear.eqQᵢ
    quadrupole = mapreduce(k -> (1 / 4) * eqQᵢ[k] * h_quadrupole(basis, k), +, 1:2)

    cᵢ = molecular_parameters.nuclear.cᵢ
    spin_rotation = mapreduce(k -> cᵢ[k] * h_nuclear_spin_rotation(basis, k), +, 1:2)

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
    prefactors = [gᵢ[k] * μN * (1 - σᵢ[k]) for k = 1:2]
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

"""
    HamiltonianParts

Contains all parts of the Hamiltonian except external fields.

Zeeman, dc Stark, and ac Stark terms are stored as vectors of matrices
that can be contracted with the appropriate external field tensors.
This allows the full Hamiltonian to be re-constructed at various
field values and orientations without recalculating all of the matrix elements.

Should be created by [`make_hamiltonian_parts`](@ref) or [`make_krb_hamiltonian_parts`](@ref).
Used as an input to [`calculate_spectrum`](@ref) and [`hamiltonian`](@ref).

See also [`make_hamiltonian_parts`](@ref), [`make_krb_hamiltonian_parts`](@ref).
"""
struct HamiltonianParts
    basis::Any
    rotation::Any
    dipole::SVector{3}
    dipole_relative::SVector{3} # used for transition strengths
    hyperfine::Any
    zeeman::SVector{3}
    ac_scalar::Any
    ac_tensor::SVector{5}
end

"""
    make_hamiltonian_parts(molecular_parameters, N_max)

Construct all parts of the Hamiltonian that do not depend on external fields.

The size of the basis is determined by `molecular_parameters`, which contains the
nuclear spin quantum numbers `molecular_parameters.I`, and the rotational states
`0:N_max` to include.

See also [`make_krb_hamiltonian_parts`](@ref).
"""
function make_hamiltonian_parts(
    molecular_parameters::MolecularParameters,
    N_max::Int,
)::HamiltonianParts
    basis = generate_basis(molecular_parameters, N_max)

    rotation = molecular_parameters.Bᵣ * h_rotation(basis)
    dipole_relative = h_dipole(basis)
    dipole = (-1) * molecular_parameters.dₚ * Constants.DVcm⁻¹ToMHz * h_dipole(basis)
    hyperfine = make_hyperfine(molecular_parameters, basis)
    zeeman = make_zeeman(molecular_parameters, basis)
    (ac_scalar, ac_tensor) = make_ac(molecular_parameters, basis)

    return HamiltonianParts(
        basis,
        rotation,
        dipole,
        dipole_relative,
        hyperfine,
        zeeman,
        ac_scalar,
        ac_tensor,
    )
end

"""
    make_krb_hamiltonian_parts(N_max)

Construct all parts of the ``{}^{40}\\text{K}^{87}\\text{Rb}`` Hamiltonian
that do not depend on external fields.

The rotational states `0:N_max` are included. This is a shortcut method that
replaces `make_hamiltonian_parts` for KRb.

See also [`make_hamiltonian_parts`](@ref).
"""
function make_krb_hamiltonian_parts(N_max::Int)
    return make_hamiltonian_parts(KRb_Parameters_Neyenhuis, N_max)
end

"""
    hamiltonian(parts, external_fields)

Construct the full Hamiltonian including magnetic, electric, and optical fields.

The field-independent building blocks in `parts` can be reused over calls
to `hamiltonian` to avoid recalculating the matrix elements each time.

See also [`make_hamiltonian_parts`](@ref), [`make_krb_hamiltonian_parts`](@ref),
[`ExternalFields`](@ref).
"""
function hamiltonian(parts::HamiltonianParts, external_fields::ExternalFields)
    h = parts.rotation + parts.hyperfine

    E = external_fields.E.magnitude
    T1E = T⁽¹⁾(SphericalUnitVector(external_fields.E))
    # Manually writing these out seems to have slightly better performance
    # than using tensor_dot??? May not actually be faster
    h +=
        E * (
            (-1) * T1E[3] * parts.dipole[1] +
            T1E[2] * parts.dipole[2] +
            (-1) * T1E[1] * parts.dipole[3]
        )

    B = external_fields.B.magnitude
    T1B = T⁽¹⁾(SphericalUnitVector(external_fields.B))
    h +=
        B * (
            (-1) * T1B[3] * parts.zeeman[1] +
            T1B[2] * parts.zeeman[2] +
            (-1) * T1B[1] * parts.zeeman[3]
        )

    for beam in external_fields.Optical
        I_laser = beam.magnitude
        T2ϵ = T⁽²⁾(SphericalUnitVector(beam))
        h += parts.ac_scalar * I_laser
        h +=
            I_laser * (
                T2ϵ[5] * parts.ac_tensor[1] +
                (-1) * T2ϵ[4] * parts.ac_tensor[2] +
                T2ϵ[3] * parts.ac_tensor[3] +
                (-1) * T2ϵ[2] * parts.ac_tensor[4] +
                T2ϵ[1] * parts.ac_tensor[5]
            )
    end

    return Array(Hermitian(h))
end
