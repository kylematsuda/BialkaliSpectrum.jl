using LinearAlgebra, SparseArrays, StaticArrays, Test
using MoleculeSpectrum

@testset "Reproduces Ospelkaus et al., PRL 104, 030402 (2010)" begin
    N_max = 5
    tolerance = 0.0011 # MHz

    fields = ExternalFields(545.9, 0.0)
    parts = make_hamiltonian_parts(KRb_Parameters_Ospelkaus, N_max)
    spectrum = calculate_spectrum(parts, fields)

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
    N_max = 5
    parts = make_hamiltonian_parts(KRb_Parameters_Neyenhuis, N_max)
    B = 545.9

    fields = ExternalFields(B, 0.0)
    spectrum = calculate_spectrum(parts, fields)

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
            spectrum_light = calculate_spectrum(parts, fields_with_light)

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
    N_max = 5
    B = 545.9
    E = 1000.0
    parts = make_hamiltonian_parts(KRb_Parameters_Neyenhuis, N_max)

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
        spectrum_z = calculate_spectrum(parts, fields[1])
        # h_z = hamiltonian(parts, fields[1])
        # energies = eigvals(h_z)

        for f in fields[2]
            spectrum = calculate_spectrum(parts, f)
            # h = hamiltonian(parts, f)
            # es = eigvals(h)

            @test spectrum.energies ≈ spectrum_z.energies
        end
    end
end

@testset "No azimuthal dependence of energies with one optical field" begin
    N_max = 5
    parts = make_hamiltonian_parts(KRb_Parameters_Neyenhuis, N_max)

    B = 545.9
    E = 0.0
    I_light = 2350.0
    θ = π / 3 # Something random-ish but not on either pole

    fields_z = ExternalFields(VectorZ(B), VectorZ(E), [SphericalVector(I_light, θ, 0.0)])
    fields_test = [
        ExternalFields(VectorZ(B), VectorZ(E), [SphericalVector(I_light, θ, φ)]) for
        φ = 0:2π/7:2π
    ]

    sz = calculate_spectrum(parts, fields_z)

    for f in fields_test
        s = calculate_spectrum(parts, f)
        @test s.energies ≈ sz.energies
    end
end

@testset "Transition strengths without hyperfine couplings" begin
    N_max = 5
    parts = make_hamiltonian_parts(TOY_MOLECULE_PARAMETERS, N_max)
    B = TOY_MOLECULE_PARAMETERS.Bᵣ

    fields = ExternalFields(0.0, 0.0) # no fields
    spectrum = calculate_spectrum(parts, fields)

    g = State(0, 0, 1, 0, 1, 0)
    e10 = State(1, 0, 1, 0, 1, 0)
    e1m1 = State(1, -1, 1, 0, 1, 0)
    e1p1 = State(1, 1, 1, 0, 1, 0)

    π_transitions =
        transition_strengths(spectrum, g, 2B - 1, 2B + 1; polarization = UnitVectorZ())
    @test all(π_transitions[1][1:2] .≈ (2 * B, 1.0))
    @test π_transitions[1][3] == e10

    x_transitions =
        transition_strengths(spectrum, g, 2B - 1, 2B + 1; polarization = UnitVectorX())
    @test all(x_transitions[1][1:2] .≈ (2 * B, 1 / sqrt(2)))
    @test all(x_transitions[2][1:2] .≈ (2 * B, 1 / sqrt(2)))
    @test (x_transitions[1][3] == e1p1 && x_transitions[2][3] == e1m1) ||
          (x_transitions[1][3] == e1m1 && x_transitions[2][3] == e1p1)

    y_transitions =
        transition_strengths(spectrum, g, 2B - 1, 2B + 1; polarization = UnitVectorY())
    @test all(y_transitions[1][1:2] .≈ (2 * B, 1 / sqrt(2)))
    @test all(y_transitions[2][1:2] .≈ (2 * B, 1 / sqrt(2)))
    @test (y_transitions[1][3] == e1p1 && y_transitions[2][3] == e1m1) ||
          (y_transitions[1][3] == e1m1 && y_transitions[2][3] == e1p1)

    unpol = transition_strengths(spectrum, g, 2B - 1, 2B + 1)
    @test all(map(k -> ≈(k[2], 1 / sqrt(3)), unpol[1:3]))
end
