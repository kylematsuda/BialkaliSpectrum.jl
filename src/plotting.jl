"""
    plot_transition_strengths(
        spectrum,
        hamiltonian_parts::HamiltonianParts,
        ground_basis_state::State,
        frequency_range::Union{Vector,Nothing} = nothing;
        tol=0.5,
        restrict_N=true,
        use_cutoff=true,
        cutoff::Union{Float64,Nothing}=1e-3,
    )

    plot_transition_strengths(
        spectra;
        groupby=:E,
    )

"""
function plot_transition_strengths(
    spectrum,
    hamiltonian_parts::HamiltonianParts,
    ground_basis_state::State,
    frequency_range::Union{Vector,Nothing} = nothing;
    groupby=:E,
    tol=0.5,
    restrict_N=true,
    use_cutoff=true,
    cutoff::Union{Float64,Nothing}=1e-3,
)
    df = get_transitions(
        spectrum,
        hamiltonian_parts,
        ground_basis_state,
        frequency_range;
        groupby=groupby,
        tol=tol,
        restrict_N=restrict_N,
        cutoff=cutoff,
    )
    return plot_transition_strengths(df; groupby=groupby)
end

function plot_transition_strengths(
    spectra;
    groupby=:E,
)
    f = Figure(fontsize = 18)
    ax = Axis(f[1,1], xlabel = "Frequency (MHz)", ylabel = "Relative transition dipole")

    grouped = DataFrames.groupby(spectra, groupby)
    for (g, k) in zip(grouped, keys(grouped))
        stem!(g.transition_frequency,
            g.transition_strength,
            label="$groupby = $(k[groupby])"
        )
    end

    ylims!(ax, 0, 1)
    Legend(f[1,2], ax)

    return f
end

"""
    plot_induced_dipole(
        spectra,
        hamiltonian_parts::HamiltonianParts;
        groupby=:fields
    )

    plot_induced_dipole(
        spectra,
        groupby=:fields
    )

"""
function plot_induced_dipole(
    spectra,
    hamiltonian_parts::HamiltonianParts;
    groupby=:E
)
    return plot_induced_dipole(
        get_induced_dipole_moments(spectra,
            hamiltonian_parts;
            groupby=groupby
        );
        groupby=groupby
    )
end

function plot_induced_dipole(spectra; groupby=:fields)
    f = Figure(fontsize = 18)

    if groupby == :E
        xlabel = "E (V/cm)"
    elseif groupby == :B
        xlabel = "B (G)"
    else
        xlabel = "$groupby"
    end

    ax = Axis(f[1, 1], xlabel=xlabel, ylabel="Induced dipole / permanent dipole")
    df = DataFrames.groupby(spectra, :basis_index)

    for group in df
        DataFrames.sort!(group, groupby)

        N = first(group.N)
        m_n = first(group.m_n)

        lines!(group.E, real.(group.d_ind), label = "$N, $m_n")
    end

    Legend(f[1, 2], ax, "States", framevisible = true)
    return f
end

function plot_states_adiabatic(
    spectra;
    groupby=:fields,
    radius::Union{Int,Nothing}=nothing 
)
    adiabatized = connect_adiabatic(spectra;
        groupby=groupby,
        radius=radius
    )

    if groupby == :E
        xlabel = "E (V/cm)"
    elseif groupby == :B
        xlabel = "B (G)"
    else
        xlabel = "$groupby"
    end

    f = Figure(
        fontsize = 14,
        font = "Helvetica",
        resolution = (300, 450),
    )

    ax = Axis(
        f[1,1],
        xlabel=xlabel,
        ylabel="Energy (MHz)",
    )

    groups = DataFrames.groupby(
        adiabatized,
        :adiabatic_index
    )

    for g in groups
        DataFrames.sort!(g, groupby)

        lines!(
            ax,
            g[!, groupby],
            g.energy,
            linewidth=2
        )
    end

    return f
end

# function plot_transitions_vs_E(
#     hamiltonian_parts::HamiltonianParts,
#     fields_scan::Vector{ExternalFields},
#     g::State,
#     frequency_range::Union{Vector,Nothing} = nothing;
#     tol = 0.5,
#     use_cutoff = false,
#     cutoff = 1e-6,
# )
#     spectra = calculate_transitions_vs_E(
#         hamiltonian_parts,
#         fields_scan,
#         g,
#         frequency_range;
#         tol = tol,
#         use_cutoff = use_cutoff,
#         cutoff = cutoff,
#     )
#     return plot_transitions_vs_E(spectra)
# end

# function plot_transitions_vs_E(spectra)
#     f = Figure(fontsize = 18)
#     ax = Axis(
#         f[1, 1],
#         backgroundcolor = :gray75,
#         xlabel = "E (V/cm)",
#         ylabel = "Frequency (MHz)",
#     )

#     df = DataFrames.groupby(spectra, :index)

#     cmap = :deep
#     colorrange = (-4, 0)

#     for group in df
#         DataFrames.sort!(group, :E)
#         lines!(
#             group.E,
#             group.transition_frequency,
#             color = log10.(group.transition_strength),
#             colormap = cmap,
#             colorrange = colorrange,
#         )
#     end

#     Colorbar(f[1, 2], limits = colorrange, colormap = cmap, label = "log10(strength)")

#     return f
# end



# function plot_chi_vs_E(
#     hamiltonian_parts::HamiltonianParts,
#     fields_scan::Vector{ExternalFields},
#     g::State;
#     p = 0,
#     tol = 0.5,
# )
#     chis = calculate_chi_vs_E(hamiltonian_parts, fields_scan, g; p = p, tol = tol)
#     plot_chi_vs_E(chis)
# end

# function plot_chi_vs_E(spectra)
#     f = Figure(fontsize = 18)
#     ax = Axis(f[1, 1], xlabel = "E (V/cm)", ylabel = L"$\chi / d_p^2$")

#     df = DataFrames.groupby(spectra, :basis_index)

#     for group in df
#         DataFrames.sort!(group, :E)
#         N = first(group.N)
#         m_n = first(group.m_n)

#         lines!(group.E, real.(group.chi), label = "$N, $m_n")
#     end

#     f[1, 2] = Legend(f, ax, "States", framevisible = true)

#     return f
# end

# function plot_states_vs_E(
#     hamiltonian_parts::HamiltonianParts,
#     fields_scan::Vector{ExternalFields},
#     states::Vector{State};
#     redraw_threshold = 0.2,
# )
#     Ns = [s.N for s in states]
#     filter_Ns(df) = DataFrames.filter(:N => n -> n in Ns, df)

#     addE(df) = DataFrames.transform(df, :fields => (f -> map(g -> g.E.magnitude, f)) => :E)

#     spectra = calculate_spectra_vs_fields(hamiltonian_parts, fields_scan, addE ∘ filter_Ns)
#     return plot_states_vs_E(spectra, states; redraw_threshold = redraw_threshold)
# end

# function plot_states_vs_E(spectra, states; redraw_threshold = 0.2)
#     f = Figure(fontsize = 18, resolution = (800, 600))
#     ax = Axis(f[1, 1], xlabel = "E (V/cm)", ylabel = "Frequency (MHz)")

#     cmap = :deep
#     colorrange = (-0.05, 1)

#     max_weight(ei) = maximum([abs2.(ei)[state_to_index(s)] for s in states])
#     df = DataFrames.transform(
#         spectra,
#         :eigenstate => (e -> map(max_weight, e)) => :max_weight,
#     )

#     grouped = DataFrames.groupby(df, :index)
#     redraw = []
#     for (i, group) in enumerate(grouped)
#         DataFrames.sort!(group, :E)

#         if maximum(group.max_weight) > redraw_threshold
#             append!(redraw, i)
#         end

#         lines!(
#             group.E,
#             group.energy,
#             color = group.max_weight,
#             colormap = cmap,
#             colorrange = colorrange,
#             transparency = false,
#             overdraw = false,
#         )
#     end

#     println("poop")

#     for i in redraw
#         group = grouped[i]
#         lines!(
#             group.E,
#             group.energy,
#             color = group.max_weight,
#             colormap = cmap,
#             colorrange = colorrange,
#             transparency = false,
#             overdraw = true,
#             fxaa = true,
#         )
#     end

#     return f
# end
