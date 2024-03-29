"""
    SphericalVector(magnitude, θ, φ)

Construct a vector with `magnitude`, polar angle `θ`, and azimuthal angle `φ`.

Represents an external field vector in spherical coordinates used
to construct [`ExternalFields`](@ref). Currently, `SphericalVector`s
may be negated but no other mathematical operations are implemented.

Vectors along ``x``, ``y``, or ``z`` can be quickly constructed
using [`VectorX`](@ref), [`VectorY`](@ref), and [`VectorZ`](@ref),
respectively.

See also [`SphericalUnitVector`](@ref), [`ExternalFields`](@ref).
"""
struct SphericalVector
    "Magnitude\n"
    magnitude::Float64
    "Polar angle (rad)\n"
    θ::Float64
    "Azimuthal angle (rad)\n"
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

"""
    VectorX(magnitude)

Construct a [`SphericalVector`](@ref) with `magnitude` along `x`.

# Examples
```jldoctest
julia> VectorX(5.0)
SphericalVector(5.0, 1.5707963267948966, 0.0)
```
"""
VectorX(magnitude) = SphericalVector(magnitude, π / 2, 0)

"""
    VectorY(magnitude)

Construct a [`SphericalVector`](@ref) with `magnitude` along `y`.

# Examples
```jldoctest
julia> VectorY(2.0)
SphericalVector(2.0, 1.5707963267948966, 1.5707963267948966)
```
"""
VectorY(magnitude) = SphericalVector(magnitude, π / 2, π / 2)

"""
    VectorZ(magnitude)

Construct a [`SphericalVector`](@ref) with `magnitude` along `z`.

# Examples
```jldoctest
julia> VectorZ(10.0)
SphericalVector(10.0, 0.0, 0.0)
```
"""
VectorZ(magnitude) = SphericalVector(magnitude, 0, 0)

"""
    SphericalUnitVector(magnitude, θ, φ)
    SphericalUnitVector(v::SphericalVector)

Construct a unit vector with polar angle `θ`, and azimuthal angle `φ`.

Represents the direction of an external field vector in spherical coordinates.
Currently, `SphericalUnitVector`s may be negated but no other
mathematical operations are implemented.

Vectors along ``x``, ``y``, or ``z`` can be quickly constructed
using [`UnitVectorX`](@ref), [`UnitVectorY`](@ref), and [`UnitVectorZ`](@ref),
respectively.

See also [`SphericalVector`](@ref), [`T⁽¹⁾`](@ref), [`T⁽²⁾`](@ref).
"""
struct SphericalUnitVector
    "Polar angle (rad)\n"
    θ::Float64
    "Azimuthal angle (rad)\n"
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

"""
    UnitVectorX()

Construct a [`SphericalUnitVector`](@ref) along `x`.
"""
UnitVectorX() = SphericalUnitVector(π / 2, 0)

"""
    UnitVectorY()

Construct a [`SphericalUnitVector`](@ref) along `y`.
"""
UnitVectorY() = SphericalUnitVector(π / 2, π / 2)

"""
    UnitVectorZ()

Construct a [`SphericalUnitVector`](@ref) along `z`.
"""
UnitVectorZ() = SphericalUnitVector(0, 0)

"""
    T⁽¹⁾(v)

Construct the components of the rank 1 spherical tensor ``T⁽¹⁾(v)``.
"""
function T⁽¹⁾(v::SphericalUnitVector)::SVector{3,ComplexF64}
    θ = v.θ
    φ = v.φ

    x = sin(θ) * cos(φ)
    y = sin(θ) * sin(φ)
    z = cos(θ)

    T11 = -(1 / sqrt(2)) * (x + im * y)
    T10 = z

    return SVector(-conj(T11), T10, T11)
end

"""
    T⁽²⁾(v)

Construct the components of the rank 2 spherical tensor ``T⁽²⁾(v, v)``.
"""
function T⁽²⁾(v::SphericalUnitVector)::SVector{5,ComplexF64}
    θ = v.θ
    φ = v.φ

    x = sin(θ) * cos(φ)
    y = sin(θ) * sin(φ)
    z = cos(θ)

    T20 = (2 * z^2 - x^2 - y^2) / sqrt(6)
    T21 = -(1 / 2) * (x * z + z * x + im * (y * z + z * y))
    T22 = (1 / 2) * (x * x - y * y + im * (x * y + y * x))

    return SVector(conj(T22), -conj(T21), T20, T21, T22)
end

function T⁽²⁾(T1a, T1b)::SVector{5,ComplexF64}
    (T1am1, T1a0, T1ap1) = T1a
    (T1bm1, T1b0, T1bp1) = T1b

    T22 = T1ap1 * T1bp1
    T21 = (1 / sqrt(2)) * (T1ap1 * T1b0 + T1a0 * T1bp1)
    T20 = (1 / sqrt(6)) * (2 * T1a0 * T1b0 + T1ap1 * T1bm1 + T1am1 * T1bp1)
    T2m1 = (1 / sqrt(2)) * (T1am1 * T1b0 + T1a0 * T1bm1)
    T2m2 = T1am1 * T1bm1
    return SVector(T2m2, T2m1, T20, T21, T22)
end

function get_tensor_component(p::Int, tensor)
    rank::Int = (length(tensor) - 1) // 2 # Length should be 2*k + 1
    return tensor[1+(p+rank)]
end

"""
    tensor_dot(a, b)

Contract two spherical tensors `a` and `b`.
"""
function tensor_dot(a, b)
    @assert size(a, 1) == size(b, 1)
    @assert isodd(size(a, 1))

    mapreduce(p -> conj(a[p]) .* b[p], +, eachindex(a))
end

"""
    ExternalFields(B::SphericalVector, E::SphericalVector, Optical::Vector{SphericalVector})
    ExternalFields(B::Float64, E::Float64)

External magnetic, electric, optical fields to use in constructing the Hamiltonian.

If `B` and `E` are provided as `Float64`s, then the fields are assumed to be along `z`.
The `Optical` argument can also be left as an empty vector `[]`.

See also [`get_spectrum`](@ref), [`hamiltonian`](@ref), [`SphericalVector`](@ref).

TODO: add the other signatures

# Examples
```jldoctest
julia> ExternalFields(VectorZ(545.9), VectorX(1020.0), [])
ExternalFields(SphericalVector(545.9, 0.0, 0.0), SphericalVector(1020.0, 1.5707963267948966, 0.0), SphericalVector[])
```

```jldoctest
julia> ExternalFields(545.9, 1020.0)
ExternalFields(SphericalVector(545.9, 0.0, 0.0), SphericalVector(1020.0, 0.0, 0.0), SphericalVector[])
```
"""
struct ExternalFields
    "Magnetic field (G)\n"
    B::SphericalVector
    "Electric field (V/cm)\n"
    E::SphericalVector
    "Laser fields (W/cm^2)\n"
    Optical::Vector{SphericalVector}

    ExternalFields(B::SphericalVector, E::SphericalVector, Optical) = new(B, E, Optical)
end

ExternalFields() = ExternalFields(0.0, 0.0, [])
ExternalFields(B, E) = ExternalFields(B, E, [])
ExternalFields(B::Float64, E, Optical) = ExternalFields(VectorZ(B), E, Optical)
ExternalFields(B, E::Float64, Optical) = ExternalFields(B, VectorZ(E), Optical)
ExternalFields(B::Float64, E::Float64, Optical) =
    ExternalFields(VectorZ(B), VectorZ(E), Optical)

const DEFAULT_FIELDS = ExternalFields(545.9, 0.0)
const TEST_FIELDS = ExternalFields(
    SphericalVector(545.9, π / 4, π / 4),
    SphericalVector(1000.0, 3π / 4, 7π / 4),
    [SphericalVector(2350.0, 0.0, 0.0)],
)

"""
    generate_fields_scan(Bs, Es, Opticals)

Produce a vector of [`ExternalFields`](@ref) for creating a scan of spectra as a function of fields.

TODO: explain the behavior as a product of the lists (not zipped)
"""
function generate_fields_scan(Bs, Es, Opticals)
    if length(Opticals) == 0
        Opticals = [[]]
    end

    return [ExternalFields(B, E, O) for B in Bs for E in Es for O in Opticals]
end
