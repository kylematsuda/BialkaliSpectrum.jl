function plot_transition_strengths(
    spectrum,
    hamiltonian_parts::HamiltonianParts,
    g::State,
    frequency_range::Union{Vector, Nothing} = nothing;
    tol = 0.5,
    restrict_N = true,
    use_cutoff = true,
    cutoff = 1e-3,
)
    df = calculate_transition_strengths(spectrum, hamiltonian_parts, g, frequency_range; tol=tol, restrict_N=restrict_N, use_cutoff=use_cutoff, cutoff=cutoff)

    f = Figure(fontsize=18)
    ax = Axis(
        f[1,1],
        xlabel = "Frequency (MHz)",
        ylabel = "Relative transition dipole"
    )

    stem!(df.transition_frequency, df.transition_strength)
    ylims!(ax, 0, 1)
    return f
end

function plot_transitions_vs_E(
    hamiltonian_parts::HamiltonianParts,
    fields_scan::Vector{ExternalFields},
    g::State,
    frequency_range::Union{Vector, Nothing} = nothing;
    tol = 0.5,
    use_cutoff = false,
    cutoff = 1e-6,
)
    spectra = calculate_transitions_vs_E(hamiltonian_parts, fields_scan, g, frequency_range; tol=tol, use_cutoff=use_cutoff, cutoff=cutoff)
    return plot_transitions_vs_E(spectra)
end

function plot_transitions_vs_E(spectra)
    f = Figure(fontsize=18)
    ax = Axis(
        f[1,1],
        backgroundcolor=:gray75,
        xlabel = "E (V/cm)",
        ylabel = "Frequency (MHz)"
    )

    df = DataFrames.groupby(spectra, :index)

    cmap = :deep
    colorrange = (-4, 0)

    for group in df
        DataFrames.sort!(group, :E)
        lines!(
            group.E,
            group.transition_frequency,
            color=log10.(group.transition_strength),
            colormap=cmap,
            colorrange=colorrange
        )
    end

    Colorbar(
        f[1,2],
        limits=colorrange,
        colormap=cmap,
        label="log10(strength)"
    )

    return f
end

function plot_induced_dipole_vs_E(
    hamiltonian_parts::HamiltonianParts,
    fields_scan::Vector{ExternalFields},
    hyperfine_manifold::Union{State, Nothing} = nothing
)
    spectra = calculate_induced_dipole_moments_vs_E(hamiltonian_parts, fields_scan, hyperfine_manifold)
    return plot_induced_dipole_vs_E(spectra)
end

function plot_induced_dipole_vs_E(spectra)
    f = Figure(fontsize=18)
    ax = Axis(
        f[1,1],
        xlabel = "E (V/cm)",
        ylabel = "Induced dipole / permanent dipole"
    )

    df = DataFrames.groupby(spectra, :basis_index)

    for group in df
        DataFrames.sort!(group, :E)

        N = first(group.N)
        m_n = first(group.m_n)

        lines!(
            group.E,
            real.(group.d_ind),
            label = "$N, $m_n"
        )
    end

    f[1,2] = Legend(
        f,
        ax,
        "States",
        framevisible = true
    )

    return f
end

function plot_chi_vs_E(
    hamiltonian_parts::HamiltonianParts,
    fields_scan::Vector{ExternalFields},
    g::State;
    p = 0,
    tol = 0.5
)
    chis = calculate_chi_vs_E(hamiltonian_parts, fields_scan, g; p=p, tol=tol)
    plot_chi_vs_E(chis)
end

function plot_chi_vs_E(spectra)
    f = Figure(fontsize=18)
    ax = Axis(
        f[1,1],
        xlabel = "E (V/cm)",
        ylabel = L"$\chi / d_p^2$"
    )

    df = DataFrames.groupby(spectra, :basis_index)

    for group in df
        DataFrames.sort!(group, :E)
        N = first(group.N)
        m_n = first(group.m_n)

        lines!(
            group.E,
            real.(group.chi),
            label = "$N, $m_n",
        )
    end

    f[1,2] = Legend(
        f,
        ax,
        "States",
        framevisible = true
    )

    return f
end