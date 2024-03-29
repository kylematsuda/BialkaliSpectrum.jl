using LinearAlgebra, SparseArrays, StaticArrays, Test
using BialkaliSpectrum, BialkaliSpectrum.K40Rb87
using DataFrames

const N_max = 5
const parts_neyenhuis = make_krb_hamiltonian_parts(N_max)

@testset "Reproduces Ospelkaus et al., PRL 104, 030402 (2010)" begin
    tolerance = 0.0011 # MHz

    fields = ExternalFields(545.9, 0.0)
    parts = make_hamiltonian_parts(KRb_Parameters_Ospelkaus, N_max)
    spectrum = get_spectrum(parts, fields)

    # Table II in paper
    comparisons = [
        ((0, 0, -4, 1 / 2), (1, 1, -4, 1 / 2), 2227.835),
        ((0, 0, -4, 1 / 2), (1, 0, -4, 1 / 2), 2228.119),
        ((0, 0, -4, 1 / 2), (1, -1, -4, 1 / 2), 2227.776),
        ((0, 0, -4, 1 / 2), (1, 0, -4, 3 / 2), 2227.008),
        ((0, 0, -4, 1 / 2), (1, -1, -4, 3 / 2), 2227.128),
        ((0, 0, -4, 1 / 2), (1, 0, -3, 1 / 2), 2228.225),
        ((0, 0, -4, 1 / 2), (1, 1, -4, -1 / 2), 2228.593),
        ((0, 0, -4, 1 / 2), (1, 0, -4, -1 / 2), 2228.805),
        ((0, 0, -4, 3 / 2), (1, 0, -4, 3 / 2), 2227.761),
        # There is probably a typo in the paper in the last line of the table
        # Should be |1, 0, -3, 1/2> instead of |1, 0, -3, 3/2>?
        ((0, 0, -3, 1 / 2), (1, 0, -3, 1 / 2), 2228.091),
    ]

    for c in comparisons
        (g, e) = map(i -> KRbState(c[i]...), 1:2)
        transition = get_energy_difference(spectrum, g, e)
        expected = c[3]

        @test abs(transition - expected) < tolerance
    end
end

@testset verbose = true "Reproduces Neyenhuis et al., PRL 109, 230403 (2012)" begin
    B = 545.9

    fields = ExternalFields(B, 0.0)
    spectrum = get_spectrum(parts_neyenhuis, fields)

    @testset "No optical fields" begin
        tolerance = 0.005 # MHz

        # See supplement of the paper
        comparisons = [
            ((0, 0, -4, 1 / 2), (1, 1, -4, 1 / 2), 2227.842),
            ((0, 0, -4, 1 / 2), (1, 0, -4, 1 / 2), 2228.110),
            ((0, 0, -4, 1 / 2), (1, -1, -4, 1 / 2), 2227.784),
        ]

        for c in comparisons
            (g, e) = map(i -> KRbState(c[i]...), 1:2)
            transition = get_energy_difference(spectrum, g, e)
            expected = c[3]

            @test abs(transition - expected) < tolerance
        end
    end

    @testset "α(θ)" begin
        I_light = 2350.0
        tolerance = 0.035 # Relative tolerance in α

        # These are generated by diagonalizing the simplified Hamiltonian H from the
        # main text and plugging in the experimental values for ϵ_i, α_parallel and α_perpendicular.
        # A few % tolerance seems reasonable given that these values come from the simplified Hamiltonian.
        #
        # Note: The supplement shows the full calculation (which should match this, in principle).
        comparisons = [
            # θ, m =  -1,      +1,      0
            (0.0, 46.4e-6, 46.4e-6, 73.2e-6),
            (10.0, 46.9e-6, 46.9e-6, 72.2e-6),
            (20.0, 48.4e-6, 48.2e-6, 69.4e-6),
            (30.0, 50.9e-6, 49.8e-6, 65.3e-6),
            (40.0, 54.0e-6, 51.2e-6, 60.7e-6),
            (50.0, 57.5e-6, 52.4e-6, 56.2e-6),
            (60.0, 60.7e-6, 53.1e-6, 52.2e-6),
            (70.0, 63.3e-6, 53.6e-6, 49.0e-6),
            (80.0, 65.0e-6, 53.9e-6, 47.1e-6),
            (90.0, 65.6e-6, 54.0e-6, 46.4e-6),
        ]
        α00 = 55.3e-6

        states_to_check =
            [(0, 0, -4, 1 / 2), (1, -1, -4, 1 / 2), (1, 1, -4, 1 / 2), (1, 0, -4, 1 / 2)]

        es = map(x -> get_energy(spectrum, KRbState(x...)), states_to_check)

        for c in comparisons
            θ = c[1] * π / 180
            optical = SphericalVector(I_light, θ, 0.0)
            fields_with_light = ExternalFields(VectorZ(B), VectorZ(0.0), [optical])
            spectrum_light = get_spectrum(parts_neyenhuis, fields_with_light)

            es_light = map(x -> get_energy(spectrum_light, KRbState(x...)), states_to_check)

            αs = map(k -> -(es_light[k] - es[k]) / I_light, eachindex(es))
            expected = [α00, c[2:end]...]

            errors = map(k -> abs(1 - (expected[k] / αs[k])), eachindex(expected))

            # println(errors) # Uncomment to show the errors on every iteration
            @test all(errors .< tolerance)
        end
    end
end

@testset "No angular dependence of energies with one field" begin
    B = 545.9
    E = 1000.0

    b_z = ExternalFields(B, 0.0)
    b_x = ExternalFields(VectorX(B), VectorZ(0.0), [])
    b_y = ExternalFields(VectorY(B), VectorZ(0.0), [])
    b_xz = ExternalFields(SphericalVector(B, π / 4, 0.0), VectorZ(0.0), [])
    b_xyz = ExternalFields(SphericalVector(B, π / 3, π / 3), VectorZ(0.0), [])

    e_z = ExternalFields(0.0, E)
    e_x = ExternalFields(VectorZ(0.0), VectorX(E), [])
    e_y = ExternalFields(VectorZ(0.0), VectorY(E), [])
    e_xz = ExternalFields(VectorZ(0.0), SphericalVector(E, π / 4, 0.0), [])
    e_xyz = ExternalFields(VectorZ(0.0), SphericalVector(E, π / 3, π / 3), [])

    for fields in [(b_z, (b_x, b_y, b_xz, b_xyz)), (e_z, (e_x, e_y, e_xz, e_xyz))]
        spectrum_z = get_spectrum(parts_neyenhuis, fields[1])

        for f in fields[2]
            spectrum = get_spectrum(parts_neyenhuis, f)
            @test spectrum.energy ≈ spectrum_z.energy
        end
    end
end

@testset "No azimuthal dependence of energies with one optical field" begin
    B = 545.9
    E = 0.0
    I_light = 2350.0
    θ = π / 3 # Something random-ish but not on either pole

    fields_z = ExternalFields(VectorZ(B), VectorZ(E), [SphericalVector(I_light, θ, 0.0)])
    fields_test = [
        ExternalFields(VectorZ(B), VectorZ(E), [SphericalVector(I_light, θ, φ)]) for
        φ = 0:2π/7:2π
    ]

    sz = get_spectrum(parts_neyenhuis, fields_z)

    for f in fields_test
        s = get_spectrum(parts_neyenhuis, f)
        @test s.energy ≈ sz.energy
    end
end

@testset "No dependence of energies under rotation of all fields" begin
    B = 545.9
    E = 0. 
    Ilaser = 2000.
    
    # Here is the baseline: all fields oriented vertically
    fields_0 = ExternalFields(B, E, [VectorZ(Ilaser)])
    spectrum_0 = filter_rotational(
        get_spectrum(parts_neyenhuis, fields_0),
        [0, 1, 2, 3]
    );
    
    fields(θ, ϕ) = ExternalFields(
        SphericalVector(B, θ, ϕ),
        SphericalVector(E, θ, ϕ),
        [SphericalVector(Ilaser, θ, ϕ)]
    )
     
    spectra_theta = get_spectra(
        parts_neyenhuis,
        [fields(θ * π / 180, 0.0) for θ=0:5:180],
        df -> filter_rotational(df, [0, 1, 2, 3])
    );

    θ_0 = π / 4
    spectra_phi = get_spectra(
        parts_neyenhuis,
        [fields(θ_0, ϕ * π / 180) for ϕ=0:5:180],
        df -> filter_rotational(df, [0, 1, 2, 3])
    );
    
    tol = 1e-9 # in MHz
    compare((e1, e2)) = abs(e1 - e2) < tol
    
    test_func(grouped) = [
        all(map(compare,
            zip(spectrum_0.energy, grouped[i].energy))) 
            for i=1:length(grouped)
    ] |> all

    grouped_theta = groupby(spectra_theta, :fields)
    @test test_func(grouped_theta)

    grouped_phi = groupby(spectra_phi, :fields)
    @test test_func(grouped_phi)
end

# @testset "Transition strengths without hyperfine couplings" begin
#     parts = make_hamiltonian_parts(TOY_MOLECULE_PARAMETERS, N_max)
#     B = TOY_MOLECULE_PARAMETERS.Bᵣ

#     fields = ExternalFields(0.0, 0.0) # no fields
#     spectrum = spectrum(parts, fields)

#     frequency_range = [2B - 1, 2B + 1]

#     g = State(0, 0, 1, 0, 1, 0)
#     e10 = State(1, 0, 1, 0, 1, 0)
#     e1m1 = State(1, -1, 1, 0, 1, 0)
#     e1p1 = State(1, 1, 1, 0, 1, 0)

#     df = calculate_transition_strengths(spectrum, parts, g)

#     # π transitions
#     π_first = sort(df, order(:d_0, by=abs2, rev=true))
#     @test abs(first(π_first).d_0) .≈ 1/sqrt(3)
#     @test first(π_first).basis_index == state_to_index(e10)

#     # σ+ transitions
#     dp_first = sort(df, order(:d_plus, by=abs2, rev=true))
#     @test abs(first(dp_first).d_plus) .≈ 1/sqrt(3)
#     @test first(dp_first).basis_index == state_to_index(e1p1)

#     # σ- transitions
#     dm_first = sort(df, order(:d_minus, by=abs2, rev=true))
#     @test abs(first(dm_first).d_minus) .≈ 1/sqrt(3)
#     @test first(dm_first).basis_index == state_to_index(e1m1)
# end
