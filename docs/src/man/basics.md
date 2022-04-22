# Getting started

```@contents
Pages = ["basics.md"]
```

## Installation

This package is not yet registered with the Julia package manager.

From the REPL, do
```julia
julia> import Pkg
julia> Pkg.add(url="https://github.com/kylematsuda/BialkaliSpectrum.jl")
```

Now bring the package contents into scope:

```@repl bialkali
using BialkaliSpectrum
```

You can optionally also pull in molecule-specific definitions, which are kept
in different submodules (`BialkaliSpectrum.K40Rb87`, `BialkaliSpectrum.Toy`, etc.):

```@repl bialkali
using BialkaliSpectrum.K40Rb87 # optional, brings all KRb-specific stuff into scope
```

## First calculation

### Setting up

Before doing the calculation, we need to do a little bit of setup: 
construct the Hamiltonian and choose the external fields to use in the calculation:

```@repl bialkali
parts = make_krb_hamiltonian_parts(5);
fields = ExternalFields(545.9, 0.0, []) # B = 545.9 G, E = 0, no optical fields
```

### Finding the energy levels

Next, we call [`get_spectrum`](@ref), which diagonalizes the Hamiltonian
(constructed from `parts`) at `fields`.

```@repl bialkali
spectrum = get_spectrum(parts, fields)
```

The output is a [`DataFrame`](https://dataframes.juliadata.org/stable/).
We can inspect its columns:

```@repl bialkali
names(spectrum)
```

Descriptions of these fields can be found in the documentation for [`get_spectrum`](@ref).

### Working with the output DataFrame

[`DataFrames.jl`](https://dataframes.juliadata.org/stable/) provides a nice API for manipulating the data that comes out
of [`get_spectrum`](@ref).

`BialkaliSpectrum.jl` provides a few convenience methods for working with `DataFrame`s.
Documentation of these methods can be found here: [DataFrame helpers](@ref).
Let's try a few of them:

```@repl bialkali
filter_rotational(spectrum, 0, 0) # get only the N = 0, m_N = 0 states
filter_rotational(spectrum, [1, 2, 3]) # get everything in N = 1, 2, 3
filter_hyperfine(spectrum, -4, [1/2, 3/2]) # get m_K = -4, m_Rb = 1/2 or 3/2 states
filter_basis_state(spectrum, KRbState(2, 2, -3, -1/2)) # get the states whose nearest basis state is |2, 2, -3, 1/2>
```

All of these methods copy the results into a new `DataFrame` to avoid mutating the
original. Typically, you'll want to calculate the spectrum first with 
all of the levels you need (including for the dipole moments, etc. to converge),
then use these methods to filter the states you actually care about into a new `DataFrame`.

!!! tip "Use DataFrames.jl"

    For more general transformations on `DataFrame`s, see the docs at
    [`DataFrames.jl`](https://dataframes.juliadata.org/stable/).
    Here are a few methods from `DataFrames.jl` that I've found especially useful:

    - [`filter`](https://dataframes.juliadata.org/stable/lib/functions/#Base.filter): 
        keep only the rows that match some predicate.
    - [`select`](https://dataframes.juliadata.org/stable/lib/functions/#DataFrames.select): 
        keep only some of the columns of the `DataFrame`.
    - [`transform`](https://dataframes.juliadata.org/stable/lib/functions/#DataFrames.transform): 
        apply a function to the values in a column
    - [`sort`](https://dataframes.juliadata.org/stable/lib/functions/#Base.sort): 
        sort the rows according to the value of some column
    - [`groupby`](https://dataframes.juliadata.org/stable/lib/functions/#DataFrames.groupby): 
        group the rows according to the value of some column.
        This is especially useful when working with a big `DataFrame` that's the
        aggregate of several spectra at different fields, and you want to split it
        back up into the individual spectrum at each field.

    I've only listed the non-mutating versions, but all of these (except `groupby`) have 
    a mutating version (e.g., `filter!`, `select!`, etc.) that may be more efficient when
    working with large `DataFrame`s, since you can avoid copying the entire `DataFrame`.

### Plotting the results

Let's start off by using our calculated `spectrum` to plot the rotational transitions
out of some state.
We'll use the state we normally populate in KRb, ``|0,0,-4,1/2\rangle``.

If you're doing this in the REPL, as opposed to an Jupyter (IJulia) notebook,
you'll need to add `ElectronDisplay.jl`,

```@repl
using ElectronDisplay
```

Then, run the following line,

```@setup plot
using BialkaliSpectrum, BialkaliSpectrum.K40Rb87
import CairoMakie

parts = make_krb_hamiltonian_parts(5);
fields = ExternalFields(545.9, 0.0, []);
spectrum = get_spectrum(parts, fields);
```

```@repl plot
plot_transition_strengths(spectrum, parts, KRbState(0, 0, -4, 1/2));
CairoMakie.save("transition_strengths.svg", ans); nothing # hide
```

This may be very slow on the first iteration, since Julia will need to compile
`CairoMakie` if this is your first plot in the session.
After a while, a new window should appear with an image that looks like this:

![](transition_strengths.svg)

These are the transitions coming out of `KRbState(0, 0, -4, 1/2)`
as a function of the transition frequency.
The height of the lines is the strength of their transition dipole, 
normalized to ``d_\text{perm} / \sqrt{3}``.
Transitions with all three polarizations, ``\pi``, ``\sigma^+`` and ``\sigma^-``, are shown.

By default, [`plot_transition_strengths`](@ref) only plots transitions with a 
relative strength (keyword arg `cutoff`) greater than `1e-3`.
To plot all of the transitions with ``N' = N \pm 1``, pass `cutoff=nothing`,

```@repl plot
plot_transition_strengths(spectrum, parts, KRbState(0, 0, -4, 1/2); cutoff=nothing);
CairoMakie.save("transition_strengths_nocutoff.svg", ans); nothing # hide
```
which gives a plot like this:

![](transition_strengths_nocutoff.svg)

Check out the docs for [`plot_transition_strengths`](@ref) for a full list of the
allowed keyword arguments.

What if we want to plot the transitions out of ``N = 1`` instead?

```@repl plot
plot_transition_strengths(spectrum, parts, KRbState(1, 0, -4, 1/2));
CairoMakie.save("transitions_N=1.svg", ans); nothing # hide
```

This produces a weird-looking plot:

![](transitions_N=1.svg)

[`plot_transition_strengths`](@ref) plots transitions to both higher and lower states
by default, so we are seeing transition to both the ``N = 0`` and ``N = 2`` manifolds.
To make things clearer, we can supply the `frequency_range` parameter,
which we have omitted so far. For example, if we want to plot the ``N = 2`` transitions,

```@repl plot
plot_transition_strengths(spectrum, parts, KRbState(1, 0, -4, 1/2), [4000, 5000]);
CairoMakie.save("transitions_N=1_to0.svg", ans); nothing # hide
```

This produces the following image:

![](transitions_N=1_to0.svg)


!!! tip "hamiltonian_parts"

    Note that [`plot_transition_strengths`](@ref) takes `parts` as a parameter 
    (written as `hamiltonian_parts` in the function signature).
    This is a common requirement for many of the analysis functions -- in this case,
    the method needs the dipole matrix elements contained within `parts` 
    to calculate the transition strengths.

    There is another signature of [`plot_transition_strengths`](@ref) that doesn't take
    `hamiltonian_parts` as an argument. This variant assumes that you are passing in a
    `DataFrame` with the requisite fields already defined (coming from e.g., the output
    of [`transitions`](@ref)).

## Scanning the external fields

So far, we have shown a few examples of what we can do with the spectrum calculated
for some particular value of the external fields.
But we often want to calculate some property of the molecules as the external fields
are changed.

To do this, we will use a vector of `ExternalField`s instead of a single one,

```@repl plot
fields = generate_fields_scan(545.9, 0.0:1000.0:10000.0, [[]]);
```

This produces a `Vector{ExternalFields}` with a constant magnetic field of
545.9 G, a constant optical intensity of zero, and electric field strength
increasing from 0 V/cm to 10 kV/cm in steps of 1 kV/cm.
(This is a very coarse scan for the sake of this tutorial.)

Next, we call [`get_spectra`](@ref) to calculate the spectrum at each field configuration,

```@repl plot
spectra = get_spectra(parts, fields, df -> filter_rotational(df, [0, 1]))
```

Notice that we passed a third parameter to [`get_spectra`](@ref), an anonymous function
`df -> filter_rotational(df, [0, 1])`.
The third argument of [`get_spectra`](@ref) is a closure that's evaluated on
the spectrum obtained at each field configuration before it's appended to the output
`DataFrame`.
In this case, we just applied a filter that removes all of the rows from the output
except those corresponding to states with ``N \leq 1``.

!!! tip "Removing unneeded states"

    For field scans containing hundreds of points and many rotational levels,
    the amount of rows stored in the `DataFrame` can quickly become very large.
    In these cases, it can sometimes speed up the calculation a bit to
    throw away the unneeded states as soon as possible, as is done in this
    example.

As an example, let's plot the energies as a function of the electric field:

```@repl plot
plot_states_adiabatic(spectra; groupby=:E);
CairoMakie.save("energies.svg", ans); nothing # hide
```

The following plot should pop up on your screen:

![](energies.svg)

The "fuzziness" of the lines is because there are 36 hyperfine states per rotational state,
which all have the same Stark shift.

Hopefully this is enough to get started!
See the [Public API](@ref) docs for more details.
We'll add a more complete example calculation at some point here: [Example: dipole moments vs E](@ref).