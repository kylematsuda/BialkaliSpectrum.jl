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
    df = transitions(
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
    groupby=:E,
    adiabatic=true
)
    return plot_induced_dipole(
        induced_dipole_moments(spectra, hamiltonian_parts);
        groupby=groupby,
        adiabatic=adiabatic
    )
end

function plot_induced_dipole(spectra; groupby=:fields, adiabatic=adiabatic)
    f = Figure(fontsize = 18)

    if groupby == :E
        xlabel = "E (V/cm)"
    elseif groupby == :B
        xlabel = "B (G)"
    else
        xlabel = "$groupby"
    end

    ax = Axis(f[1, 1], xlabel=xlabel, ylabel="Induced dipole / permanent dipole")

    if adiabatic
        key = :adiabatic_index
    else
        key = :index
    end

    df = DataFrames.groupby(spectra, key)

    for group in df
        DataFrames.sort!(group, groupby)

        N = first(group.N)
        m_n = first(group.m_n)

        lines!(group.E, real.(group.d_ind), label = "$N, $m_n")
    end

    Legend(f[1, 2], ax, "States", framevisible = true)
    return f
end

"""
    plot_transition_dipole(
        spectra,
        hamiltonian_parts::HamiltonianParts,
        initial_state::State,
        p::Int;
        groupby=:fields
    )

    plot_transition_dipole(
        spectra,
        p::Int,
        groupby=:fields
    )

"""
function plot_transition_dipole(
    spectra,
    hamiltonian_parts::HamiltonianParts,
    initial_state::State,
    p::Int;
    groupby=:E,
    adiabatic=true,
)
    return plot_transition_dipole(
        transitions(
            spectra,
            hamiltonian_parts,
            initial_state;
            groupby=groupby,
            restrict_N=false,
            cutoff=nothing
        ),
        initial_state,
        p;
        groupby=groupby,
        adiabatic=adiabatic
    )
end

function plot_transition_dipole(
    spectra,
    initial_state::State,
    p::Int;
    groupby=:fields,
    adiabatic=adiabatic
)
    f = Figure(fontsize = 18)

    if groupby == :E
        xlabel = "E (V/cm)"
    elseif groupby == :B
        xlabel = "B (G)"
    else
        xlabel = "$groupby"
    end

    ax = Axis(f[1, 1], xlabel=xlabel, ylabel="Transition dipole / permanent dipole")

    if adiabatic
        key = :adiabatic_index
    else
        key = :index
    end

    if p == -1
        d = :d_minus
    elseif p == 0
        d = :d_0
    else
        d = :d_plus
    end

    df = DataFrames.groupby(spectra, key)

    for group in df
        DataFrames.sort!(group, groupby)

        firstrow = first(group)

        if !(firstrow.N == initial_state.N &&
            firstrow.m_n == initial_state.mₙ &&
            firstrow.m_i1 == initial_state.mᵢ[1] &&
            firstrow.m_i2 == initial_state.mᵢ[2])

            lines!(group.E, abs.(group[!, d]), label = "$(firstrow.N), $(firstrow.m_n)")
        end
    end

    Legend(f[1, 2], ax, "States", framevisible = true)
    return f
end

"""
    plot_states_adiabatic(
        spectra;
        groupby=:E,
        radius::Union{Int,Nothing}=nothing
    )
"""
function plot_states_adiabatic(
    spectra;
    groupby=:fields,
    radius::Union{Int,Nothing}=nothing 
)
    adiabatized = adiabatic(spectra;
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

"""
    plot_states_adiabatic_weighted(
        spectra,
        states::Vector{State};
        groupby=:E,
        radius::Union{Int,Nothing}=nothing
    )
"""
function plot_states_adiabatic_weighted(
    spectra,
    states::Vector{State};
    groupby=:E,
    radius::Union{Int,Nothing}=nothing 
)
    adiabatized = adiabatic(spectra;
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
        ev[basis_index(s)] |> abs2
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

"""
    plot_transitions_adiabatic(
        spectra,
        hamiltonian_parts::HamiltonianParts,
        ground_basis_state::State,
        frequency_range::Union{Vector,Nothing} = nothing;
        groupby=:E,
        tol=0.5,
        restrict_N=true,
        radius::Union{Int,Nothing}=nothing 
    )
"""
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

    df = transitions(
        spectra,
        hamiltonian_parts,
        ground_basis_state,
        frequency_range;
        groupby=groupby,
        tol=tol,
        restrict_N=restrict_N,
        cutoff=nothing,
    )

    adiabatized = adiabatic(df;
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