# Example: dipole moments vs E

```@contents
Pages = ["worked_example.md"]
```

This is a more complete worked example using `BialkaliSpectrum`.

Let's compute various dipole moments as a function of electric field!

For simplicity, we'll do the calculation without including hyperfine
structure. `BialkaliSpectrum` has a module `BialkaliSpectrum.Toy` that
defines molecular parameters for an idealized molecule without hyperfine
structure.

## Setup

```@setup dipoles
using BialkaliSpectrum, BialkaliSpectrum.Toy
import DataFrames, CairoMakie
```

We'll start by setting up a scan of the electric field:

```@repl dipoles
parts = make_toy_hamiltonian_parts(8);
fields = generate_fields_scan(0.0, 0.0:0.2:4.0, [[]]);
spectra = get_spectra(parts, fields, df -> filter_rotational(df, [0, 1, 2]))
```

This gives us the eigenstates and energies as a function of field.

For the toy molecule, the molecular parameters are such that the characteristic
field is ``E_c = B / d = 1 \text{ V/cm}``.
This is a bit of a weird situation, so let's try to make this somewhat 
more clear by adding a new column `:ϵ` to the
`DataFrame`, representing the unitless E-field parameter ``\epsilon = E/E_c``:

```@repl dipoles
DataFrames.transform!(spectra, :E => DataFrames.ByRow(identity) => :ϵ);
```

We'll also want to use [`adiabatic`](@ref) here.
This method adds a new column, `:adiabatic_index`. To obtain the values of
these columns, the function first groups the `DataFrame` by `:fields` and 
indexes the state in the first group in order of increasing energy; the value
of `:adiabatic_index` is assigned to these index.
The function next steps through the remaining groups, at each point giving each 
eigenstate's `:adiabatic_index` to match the value of the eigenstate in the previous
group with which it has the best overlap.

In other words, [`adiabatic`](@ref) connects states for us across different
fields based on the criterion of wavefunction adiabaticity as the field is changed.
This lets us track states through avoided crossings, and distinguishes between
real crossings and avoided ones.

```@repl dipoles
spectra = adiabatic(spectra);
```

## Plotting the induced dipoles

Having done this, let's now use [`induced_dipole_moments`](@ref)
to compute the dipole moments, which are added in a column called `:d_ind`:

```@repl dipoles
dips = induced_dipole_moments(spectra, parts);
dips.d_ind
```

We can now easily plot the dipole moments with [`plot_induced_dipole`](@ref):

```@repl dipoles
plot_induced_dipole(dips; groupby=:ϵ);
CairoMakie.save("induced_dipoles.svg", ans); nothing # hide
```

![](induced_dipoles.svg)

!!! note "Keyword arguments to the plotting function"

    A bit of explanation:
    - The `groupby` keyword argument sets the x-axis of the plot, and defaults to `:E`. In this case, we override the default since we want to plot against `:ϵ` instead.
    - The `use_adiabatic` kwarg sets how states are connected in the line plot across different fields. The default is `use_adiabatic=true`, which uses `:adiabatic_index` to connect the states. If `use_adiabatic=false`, then `:index` (counts energies at a given field from lowest to highest) is used instead.

## Plotting transition dipoles

Let's next plot the dipole matrix elements between different states.
To do this, we'll call [`transitions`](@ref):

```@repl dipoles
ts = transitions(spectra, parts, State(0, 0); restrict_N=false, cutoff=nothing);
```

The third argument, `State(0, 0)`, sets the "ground state" for this calculation.
The dipole moments between this state and each other state are computed and 
stored in the columns `:d_0`, `:d_plus`, and `:d_minus`.

Note that we also override some of the default keyword args.
The defaults have a nice behavior for some other calculations, but aren't
the best here. See the docs for [`transitions`](@ref) for more info.

We can plot these up. Here are the ``d^0`` matrix elements:

```@repl dipoles
plot_transition_dipole(ts, State(0, 0), 0; groupby=:ϵ);
CairoMakie.save("transition_dipoles_0.svg", ans); nothing # hide
```

![](transition_dipoles_0.svg)

The third argument sets which spherical tensor component of the dipole operator is plotted. We can plot the ``d^{+1}`` matrix elements instead:

```@repl dipoles
plot_transition_dipole(ts, State(0, 0), 1; groupby=:ϵ);
CairoMakie.save("transition_dipoles_1.svg", ans); nothing # hide
```

![](transition_dipoles_1.svg)

## Dumping the outputs to a file

For later use, we often want to output the results to a file.
In this example, we'll put induced dipole moments of ``|0,0\rangle`` and
``|1,0\rangle`` and the transition dipole moment between them into a CSV file.

We've already calculated all the necessary dipole moments, but the `DataFrame` 
also has a bunch of extra information. We can remind ourselves by running the following command:

```@repl dipoles
names(ts)
```

Let's `select` what we want into a new `DataFrame`:

```@repl dipoles
n_0_and_1 = filter_rotational(dips, [0, 1], 0);
d_inds = DataFrames.select(n_0_and_1, [:ϵ, :N, :m_n, :d_ind])
```

The data is currently in "long format", with separate rows for ``|0,0\rangle`` and
``|1,0\rangle`` at each field. For this example, let's put it in "wide format"
instead, so we have a single row at each field with the induced dipoles (`:d00` and `:d11`) as the columns. For this, we can use the [`wide_format`](@ref) function:

```@repl dipoles
formatted = wide_format(d_inds, :d_ind, row -> "d$(row.N)$(row.N)"; groupby=:ϵ)
```

To get the transition dipole, we can filter our original `DataFrame` and join it
with `formatted`:

```@repl dipoles
d_transition = DataFrames.select(filter_rotational(ts, 1, 0), [:ϵ, :d_0]);
DataFrames.transform!(d_transition, :d_0 => DataFrames.ByRow(abs) => :d01);
DataFrames.select!(d_transition, [:ϵ, :d01]);

formatted = DataFrames.DataAPI.rightjoin(formatted, d_transition; on=:ϵ)
```

Now that we have it in the right format, we can easily [save it to a CSV](https://dataframes.juliadata.org/stable/man/importing_and_exporting/).

## Custom plotting with CairoMakie.jl

Plan: Put a custom plot here using the contents of `formatted`.

## Plotting shortcuts

Plan: Show examples of the other signatures of the above plotting methods (the ones that take `hamiltonian_parts`).