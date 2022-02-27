module MoleculeSpectrum

import WignerSymbols: wigner3j
import HalfIntegers: HalfInt
using LinearAlgebra, SparseArrays, StaticArrays, Test

export ZeemanParameters, NuclearParameters, Polarizability, MolecularParameters
export KRb_Parameters_Neyenhuis, KRb_Parameters_Ospelkaus, DEFAULT_MOLECULAR_PARAMETERS

export SphericalVector, VectorX, VectorY, VectorZ
export SphericalUnitVector, UnitVectorX, UnitVectorY, UnitVectorZ
export T⁽¹⁾, T⁽²⁾, get_tensor_component, tensor_dot
export ExternalFields, DEFAULT_FIELDS, TEST_FIELDS

export State, index_to_state, state_to_index
export order_by_overlap_with, max_overlap_with, get_energy, get_energy_difference

export HamiltonianParts, make_hamiltonian_parts, hamiltonian

module Constants
    const μN = 7.622593285e-4 # MHz/G
    const DToSI = 3.33564e-30
    const h = 6.62607004e-34
    const DVcm⁻¹ToMHz = (DToSI / h) * 1e-4
end # module

include("molecular_parameters.jl")
include("fields.jl")

struct State
    N::Int
    mₙ::Int
    I::SVector{2, HalfInt} # [K, Rb]
    mᵢ::SVector{2, HalfInt}
end

State(N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂) = State(N, mₙ, SVector(I₁, I₂), SVector(mᵢ₁, mᵢ₂))
State(N, mₙ, mᵢ₁::Number, mᵢ₂::Number) = State(N, mₙ, DEFAULT_MOLECULAR_PARAMETERS.I, [HalfInt(mᵢ₁) HalfInt(mᵢ₂)])

include("utility.jl")
include("matrix_elements.jl")
include("hamiltonian.jl")

end # module
