function h_rotation(basis::Vector{State})
    return Diagonal(map(x -> rotation_matrix_element(x, x), basis))
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

function h_ac_scalar(basis::Vector{State})
    elts::Int = length(basis)
    H = Diagonal(ones(elts))
    return Hermitian(H)
    # H = Diagonal(map(x -> scalar_polarizability(x, x), basis))
    # return Hermitian(H)
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
    T1E = T⁽¹⁾(SphericalUnitVector(external_fields.E))
    # Manually writing these out seems to have slightly better performance
    # than using tensor_dot??? May not actually be faster
    h += E * ((-1) * T1E[3] * parts.dipole[1] +
        T1E[2] * parts.dipole[2] +
        (-1) * T1E[1] * parts.dipole[3])

    B = external_fields.B.magnitude
    T1B = T⁽¹⁾(SphericalUnitVector(external_fields.B))
    h += B * ((-1) * T1B[3] * parts.zeeman[1] +
        T1B[2] * parts.zeeman[2] +
        (-1) * T1B[1] * parts.zeeman[3])

    for beam in external_fields.Optical
        I_laser = beam.magnitude
        T2ϵ = T⁽²⁾(SphericalUnitVector(beam))
        h += parts.ac_scalar * I_laser
        h += I_laser * (T2ϵ[5] * parts.ac_tensor[1] +
            (-1) * T2ϵ[4] * parts.ac_tensor[2] +
            T2ϵ[3] * parts.ac_tensor[3] +
            (-1) * T2ϵ[2] * parts.ac_tensor[4] +
            T2ϵ[1] * parts.ac_tensor[5])
    end

    return Array(Hermitian(h))
end