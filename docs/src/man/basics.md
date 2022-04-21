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

```jldoctest bialkali
julia> using BialkaliSpectrum
```

You can optionally also pull in molecule-specific definitions, which are kept
in different submodules (`BialkaliSpectrum.K40Rb87`, `BialkaliSpectrum.Toy`, etc.):

```jldoctest bialkali
julia> using BialkaliSpectrum.K40Rb87 # optional, brings all KRb-specific stuff into scope
```

## First calculation

### Setting up

At the start of any calculation, we need to do a little bit of setup, namely,
construct the Hamiltonian and choose the external fields to use in the calculation:

```jldoctest bialkali
julia> parts = make_krb_hamiltonian_parts(5);

julia> fields = ExternalFields(545.9, 0.0, []) # B = 545.9 G, E = 0, no optical fields
ExternalFields(SphericalVector(545.9, 0.0, 0.0), SphericalVector(0.0, 0.0, 0.0), SphericalVector[])
```

### Finding the energy levels

Next, we call [`get_spectrum`](@ref), which diagonalizes the Hamiltonian
(constructed from `parts`) at `fields`.

```jldoctest bialkali
julia> spectrum = get_spectrum(parts, fields)
1296×13 DataFrame
  Row │ fields                             index  energy        eigenstate     ⋯
      │ External…                          Int64  Float64       SubArray…      ⋯
──────┼─────────────────────────────────────────────────────────────────────────
    1 │ ExternalFields(SphericalVector(5…      1     -1.6672    ComplexF64[3.9 ⋯
    2 │ ExternalFields(SphericalVector(5…      2     -1.53563   ComplexF64[-8.
    3 │ ExternalFields(SphericalVector(5…      3     -1.40404   ComplexF64[4.0
    4 │ ExternalFields(SphericalVector(5…      4     -1.27245   ComplexF64[-6.
    5 │ ExternalFields(SphericalVector(5…      5     -1.14086   ComplexF64[0.0 ⋯
    6 │ ExternalFields(SphericalVector(5…      6     -1.00925   ComplexF64[0.0
    7 │ ExternalFields(SphericalVector(5…      7     -0.914827  ComplexF64[4.8
    8 │ ExternalFields(SphericalVector(5…      8     -0.877637  ComplexF64[0.0
  ⋮   │                 ⋮                    ⋮         ⋮                       ⋱
 1290 │ ExternalFields(SphericalVector(5…   1290  33420.2       ComplexF64[0.0 ⋯
 1291 │ ExternalFields(SphericalVector(5…   1291  33420.2       ComplexF64[2.5
 1292 │ ExternalFields(SphericalVector(5…   1292  33420.3       ComplexF64[0.0
 1293 │ ExternalFields(SphericalVector(5…   1293  33420.3       ComplexF64[0.0
 1294 │ ExternalFields(SphericalVector(5…   1294  33420.3       ComplexF64[0.0 ⋯
 1295 │ ExternalFields(SphericalVector(5…   1295  33420.3       ComplexF64[0.0
 1296 │ ExternalFields(SphericalVector(5…   1296  33420.4       ComplexF64[6.4
                                                10 columns and 1281 rows omitted

```

The output is a [`DataFrame`](https://dataframes.juliadata.org/stable/).
We can inspect its columns,

```jldoctest bialkali
julia> names(spectrum)
13-element Vector{String}:
 "fields"
 "index"
 "energy"
 "eigenstate"
 "basis_index"
 "N"
 "m_n"
 "I_1"
 "m_i1"
 "I_2"
 "m_i2"
 "B"
 "E"
```

Descriptions of these fields can be found in the documentation for [`get_spectrum`](@ref).

### Working with the output DataFrame

We put the output of `get_spectrum` into a `DataFrame` because
[`DataFrames.jl`](https://dataframes.juliadata.org/stable/) provides a nice API for manipulating all of the data that comes out.
Even the simple calculation we just did, which only involves a single field
and ``N \leq 5`` rotational states, gives a decent amount of data.

`BialkaliSpectrum.jl` provides a few convenience methods for working with `DataFrame`s in the API docs at [ DataFrame helpers](@ref).
Here, we'll demonstrate a few of them:

```jldoctest bialkali
julia> filter_rotational(spectrum, 0, 0) # get only the N = 0, m_N = 0 states
36×13 DataFrame
 Row │ fields                             index  energy     eigenstate         ⋯
     │ External…                          Int64  Float64    SubArray…          ⋯
─────┼──────────────────────────────────────────────────────────────────────────
   1 │ ExternalFields(SphericalVector(5…      1  -1.6672    ComplexF64[3.97483 ⋯
   2 │ ExternalFields(SphericalVector(5…      2  -1.53563   ComplexF64[-8.1453
   3 │ ExternalFields(SphericalVector(5…      3  -1.40404   ComplexF64[4.08372
   4 │ ExternalFields(SphericalVector(5…      4  -1.27245   ComplexF64[-6.9833
   5 │ ExternalFields(SphericalVector(5…      5  -1.14086   ComplexF64[0.0-0.0 ⋯
   6 │ ExternalFields(SphericalVector(5…      6  -1.00925   ComplexF64[0.0-0.0
   7 │ ExternalFields(SphericalVector(5…      7  -0.914827  ComplexF64[4.87699
   8 │ ExternalFields(SphericalVector(5…      8  -0.877637  ComplexF64[0.0-0.0
  ⋮  │                 ⋮                    ⋮        ⋮                      ⋮  ⋱
  30 │ ExternalFields(SphericalVector(5…     30   0.92291   ComplexF64[0.0-0.0 ⋯
  31 │ ExternalFields(SphericalVector(5…     31   1.00312   ComplexF64[1.6655e
  32 │ ExternalFields(SphericalVector(5…     32   1.14082   ComplexF64[-4.7147
  33 │ ExternalFields(SphericalVector(5…     33   1.27851   ComplexF64[-2.2568
  34 │ ExternalFields(SphericalVector(5…     34   1.41619   ComplexF64[0.0-0.0 ⋯
  35 │ ExternalFields(SphericalVector(5…     35   1.55387   ComplexF64[0.0-0.0
  36 │ ExternalFields(SphericalVector(5…     36   1.69154   ComplexF64[0.0-0.0
                                                  10 columns and 21 rows omitted

julia> filter_rotational(spectrum, [1, 2, 3]) # get everything in N = 1, 2, 3
540×13 DataFrame
 Row │ fields                             index  energy    eigenstate          ⋯
     │ External…                          Int64  Float64   SubArray…           ⋯
─────┼──────────────────────────────────────────────────────────────────────────
   1 │ ExternalFields(SphericalVector(5…     37   2226.11  ComplexF64[0.0-0.0i ⋯
   2 │ ExternalFields(SphericalVector(5…     38   2226.2   ComplexF64[0.0-0.0i
   3 │ ExternalFields(SphericalVector(5…     39   2226.22  ComplexF64[0.0-0.0i
   4 │ ExternalFields(SphericalVector(5…     40   2226.28  ComplexF64[0.0-0.0i
   5 │ ExternalFields(SphericalVector(5…     41   2226.31  ComplexF64[0.0-0.0i ⋯
   6 │ ExternalFields(SphericalVector(5…     42   2226.37  ComplexF64[0.0-0.0i
   7 │ ExternalFields(SphericalVector(5…     43   2226.42  ComplexF64[0.0-0.0i
   8 │ ExternalFields(SphericalVector(5…     44   2226.44  ComplexF64[0.0-0.0i
  ⋮  │                 ⋮                    ⋮       ⋮                      ⋮   ⋱
 534 │ ExternalFields(SphericalVector(5…    570  13369.1   ComplexF64[0.0-0.0i ⋯
 535 │ ExternalFields(SphericalVector(5…    571  13369.1   ComplexF64[0.0-0.0i
 536 │ ExternalFields(SphericalVector(5…    572  13369.1   ComplexF64[0.0-0.0i
 537 │ ExternalFields(SphericalVector(5…    573  13369.1   ComplexF64[0.0-0.0i
 538 │ ExternalFields(SphericalVector(5…    574  13369.2   ComplexF64[0.0-0.0i ⋯
 539 │ ExternalFields(SphericalVector(5…    575  13369.2   ComplexF64[0.0-0.0i
 540 │ ExternalFields(SphericalVector(5…    576  13369.2   ComplexF64[0.0-0.0i
                                                 10 columns and 525 rows omitted

julia> filter_hyperfine(spectrum, -4, [1/2, 3/2]) # get m_K = -4, m_Rb = 1/2 or 3/2 states
72×13 DataFrame
 Row │ fields                             index  energy        eigenstate      ⋯
     │ External…                          Int64  Float64       SubArray…       ⋯
─────┼──────────────────────────────────────────────────────────────────────────
   1 │ ExternalFields(SphericalVector(5…      1     -1.6672    ComplexF64[3.97 ⋯
   2 │ ExternalFields(SphericalVector(5…      7     -0.914827  ComplexF64[4.87
   3 │ ExternalFields(SphericalVector(5…     37   2226.11      ComplexF64[0.0-
   4 │ ExternalFields(SphericalVector(5…     39   2226.22      ComplexF64[0.0-
   5 │ ExternalFields(SphericalVector(5…     40   2226.28      ComplexF64[0.0- ⋯
   6 │ ExternalFields(SphericalVector(5…     54   2226.87      ComplexF64[0.0-
   7 │ ExternalFields(SphericalVector(5…     56   2226.92      ComplexF64[0.0-
   8 │ ExternalFields(SphericalVector(5…     66   2227.19      ComplexF64[0.0-
  ⋮  │                 ⋮                    ⋮         ⋮                        ⋱
  66 │ ExternalFields(SphericalVector(5…    976  33417.6       ComplexF64[1.27 ⋯
  67 │ ExternalFields(SphericalVector(5…    979  33417.7       ComplexF64[-1.5
  68 │ ExternalFields(SphericalVector(5…    983  33417.7       ComplexF64[-1.1
  69 │ ExternalFields(SphericalVector(5…    984  33417.7       ComplexF64[1.47
  70 │ ExternalFields(SphericalVector(5…    988  33417.7       ComplexF64[3.47 ⋯
  71 │ ExternalFields(SphericalVector(5…    992  33417.7       ComplexF64[9.27
  72 │ ExternalFields(SphericalVector(5…    993  33417.7       ComplexF64[1.16
                                                  10 columns and 57 rows omitted

julia> filter_basis_state(spectrum, KRbState(2, 2, -3, -1/2)) # get the states whose nearest basis state is |2, 2, -3, 1/2>
1×13 DataFrame
 Row │ fields                             index  energy   eigenstate           ⋯
     │ External…                          Int64  Float64  SubArray…            ⋯
─────┼──────────────────────────────────────────────────────────────────────────
   1 │ ExternalFields(SphericalVector(5…    225  6683.58  ComplexF64[-5.56021e ⋯
                                                              10 columns omitted

```

These methods all copy the results into a new `DataFrame` to avoid mutating the
original. Typically, you'll want to calculate the spectrum once initially with 
all of the levels you need (including for the dipole moments, etc. to converge).
Then, before further using the data, you can use these methods to filter the
states you want into a new `DataFrame`.

!!! tip "Use DataFrames.jl"

    For more general transformations on `DataFrame`s, see the docs at
    [`DataFrames.jl`](https://dataframes.juliadata.org/stable/).
    Here are a few methods from `DataFrames.jl` that I found especially useful
    when writing `BialkaliSpectrum`:

    - [`filter`](https://dataframes.juliadata.org/stable/lib/functions/#Base.filter): 
        retain only the rows that match some predicate.
    - [`select`](https://dataframes.juliadata.org/stable/lib/functions/#DataFrames.select): 
        retain only some of the columns of the `DataFrame`.
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
out of some state, say ``|0,0,-4,1/2\rangle``, which is the state we normally
populate in KRb.

If you're doing this in the REPL, as opposed to an Jupyter (IJulia) notebook,
you'll need to add `ElectronDisplay.jl`,

```julia
julia> using ElectronDisplay
```

Then, run the following method,

```julia
julia> plot_transition_strengths(spectrum, parts, KRbState(0, 0, -4, 1/2))
```

This may be very slow on the first iteration, since Julia will need to compile
`CairoMakie` if this is your first plot in the session.
After a while, a new window should appear with an image that looks like this:

![](basics_firstplot.png)

These are the transitions coming out of our chosen state `KRbState(0, 0, -4, 1/2)`
as a function of the transition frequency.
The height of the lines is the strength of their transition dipole, 
normalized to ``d_\text{perm} / \sqrt{3}``.

By default, [`plot_transition_strengths`](@ref) only plots transitions with a 
relative strength (`cutoff`) greater than `1e-3`.
To plot all of the transitions with ``N' = N \pm 1``, pass `cutoff=nothing`,

```julia
julia> plot_transition_strengths(spectrum, parts, KRbState(0, 0, -4, 1/2); cutoff=nothing)
```
which gives a plot like this:

![](basics_secondplot.png)

Check out the docs for [`plot_transition_strengths`](@ref) for a full list of the
allowed keyword arguments.

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