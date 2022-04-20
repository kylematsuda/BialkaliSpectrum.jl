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
    adiabatized = connect_adiabatically(spectra;
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

function plot_states_adiabatic_weighted(
    spectra,
    states::Vector{State};
    groupby=:E,
    radius::Union{Int,Nothing}=nothing 
)
    adiabatized = connect_adiabatically(spectra;
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
        backgroundcolor = :gray85,
    )

    cmap = :deep
    colorrange = (-4, 0)

    max_weight(ev) = maximum([
        ev[state_to_index(s)] |> abs2
        for s in states
    ])

    DataFrames.transform!(adiabatized,
        :eigenstate => DataFrames.ByRow(max_weight) => :max_weight
    )

    function add_sup_weight(df)
        DataFrames.sort!(df, :max_weight)
        sup_weight = last(df).max_weight
        DataFrames.transform!(df,
            :index => DataFrames.ByRow(x -> sup_weight) => :sup_weight
        )
        return df
    end

    by_sup_weight = transform_spectra(
        adiabatized,
        add_sup_weight;
        groupby=:adiabatic_index
    )


    groups = DataFrames.groupby(
        by_sup_weight,
        :sup_weight;
        sort=true
    )

    for g in groups
        DataFrames.sort!(g, groupby)

        lines!(
            ax,
            g[!, groupby],
            g.energy,
            linewidth=2,
            color=log10.(g.max_weight),
            colormap=cmap,
            colorrange=colorrange
        )
    end

    Colorbar(f[1, 2], limits=colorrange, colormap=cmap, label="log10(weight)")

    return f
end

function plot_transitions_adiabatic(
    spectra,
    hamiltonian_parts::HamiltonianParts,
    ground_basis_state::State,
    frequency_range::Union{Vector,Nothing} = nothing;
    groupby=:E,
    tol=0.5,
    restrict_N=true,
    radius::Union{Int,Nothing}=nothing 
)

    df = get_transitions(
        spectra,
        hamiltonian_parts,
        ground_basis_state,
        frequency_range;
        groupby=groupby,
        tol=tol,
        restrict_N=restrict_N,
        cutoff=nothing,
    )

    adiabatized = connect_adiabatically(df;
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
        backgroundcolor = :gray85,
    )

    function add_max_strength(df)
        DataFrames.sort!(df, :transition_strength)
        max_strength = last(df).transition_strength
        DataFrames.transform!(df,
            :index => DataFrames.ByRow(x -> max_strength) => :max_strength
        )
        return df
    end

    by_strength = transform_spectra(
        adiabatized,
        add_max_strength;
        groupby=:adiabatic_index
    )

    groups = DataFrames.groupby(
        adiabatized,
        :max_strength;
        sort=true
    )

    cmap = :deep
    colorrange = (-4, 0)

    for g in groups
        DataFrames.sort!(g, groupby)

        lines!(
            ax,
            g[!, groupby],
            g.transition_frequency,
            color = log10.(g.transition_strength),
            colormap=cmap,
            colorrange=colorrange,
            linewidth=2
        )
    end

    Colorbar(f[1, 2], limits=colorrange, colormap=cmap, label="log10(strength)")
    
    return f
end