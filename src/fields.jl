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

function tensor_dot(a, b)
    @assert size(a, 1) == size(b, 1)
    @assert isodd(size(a, 1))

    mapreduce(p -> conj(a[p]) .* b[p], +, eachindex(a))
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
const TEST_FIELDS = ExternalFields(SphericalVector(545.9, π/4, π/4), SphericalVector(1000., 3π/4, 7π/4), [SphericalVector(2350., 0.0, 0.0)])