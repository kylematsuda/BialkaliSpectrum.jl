"""
    filter_rotational(df, N, m_n)

Returns a new `DataFrame` containing the rows of `df` where `(row.N in N) && (row.m_n in m_n)`.
"""
filter_rotational(df, N, m_n) = DataFrames.filter(
    [:N, :m_n] =>
        (n, m) -> (n in N) && (m in m_n),
    df
)

"""
    filter_rotational!(df, N, m_n)

Retains the rows of `df` where `(row.N in N) && (row.m_n in m_n)`, removing the rest.
This mutates `df`.
"""
filter_rotational!(df, N, m_n) = DataFrames.filter!(
    [:N, :m_n] =>
        (n, m) -> (n in N) && (m in m_n),
    df
)

"""
    filter_hyperfine(df, m_i1, m_i2)

Returns a new `DataFrame` containing the rows of `df` where `(row.m_i1 in m_i1) && (row.m_i2 in m_i2)`.
"""
filter_hyperfine(df, m_i1, m_i2) = DataFrames.filter(
    [:m_i1, :m_i2] =>
        (m1, m2) -> (m1 in m_i1) && (m2 in m_i2),
    df
)

"""
    filter_hyperfine!(df, m_i1, m_i2)

Retains the rows of `df` where `(row.m_i1 in m_i1) && (row.m_i2 in m_i2)`, removing the rest.
This mutates `df`.
"""
filter_hyperfine!(df, m_i1, m_i2) = DataFrames.filter(
    [:m_i1, :m_i2] =>
        (m1, m2) -> (m1 in m_i1) && (m2 in m_i2),
    df
)

"""
    filter_basis_state(df, basis_state::State)

Returns a new `DataFrame` containing the rows of `df` whose `:eigenstate` is
most overlapping with `basis_state`.
"""
filter_basis_state(df, basis_state::State) = 
    DataFrames.filter(:basis_index => bi -> bi == state_to_index(basis_state), df)

"""
    filter_basis_state!(df, basis_state::State)

Retains the rows of `df` whose `:eigenstate` is most overlapping with `basis_state`.
This mutates `df`.
"""
filter_basis_state(df, basis_state::State) = 
    DataFrames.filter!(:basis_index => bi -> bi == state_to_index(basis_state), df)

"""
    expand_fields!(df; magnitude_only=true)
    expand_fields!(df, col_names::Vector{Symbol}, expander::Function)

Expands the `:fields` column of `df` to new columns with the individual field
parameters.

In the first variant, the optical fields are ignored. 
If `magnitude_only=true` (default), then the `:fields`
column is expanded to new columns `:B` and `:E` contaning `fields.B.magnitude`
and `fields.E.magnitude`, respectively. This is convenient for plotting, so the
values of E can be accessed as `df.E` (and similarly for `df.B`).
If `magnitude_only=false`, then new columns for the field angles `:θ_B`, `:ϕ_B`,
`:θ_E`, and `:ϕ_E` are added, in addition to the field magnitudes.

The second variant allows full control of the behavior. `col_names` is the
list of the columns to add to `df`. `expander` must be a function
`fields -> Vector`, which defines how a value in the `:fields` column is 
transformed to the quantities named in `col_names`.
This variant should be used for extracting the optical field parameters.

Add examples!
"""
expand_fields!(df, field_names::Vector{Symbol}, expander::Function) =
    transform!(df, :fields => ByRow(expander) => fields_names)

function expand_fields!(df; magnitude_only=true)
    if magnitude_only
        field_names = [:B, :E]
        expander = fields -> [fields.B.magnitude, fields.E.magnitude]
    else
        field_names = [:B, :θ_B, :ϕ_B, :E, :θ_E, :ϕ_E]
        expander = fields -> [
            fields.B.magnitude,
            fields.B.θ,
            fields.B.ϕ,
            fields.E.magnitude,
            fields.E.θ,
            fields.E.ϕ
        ]
    end

    return expand_fields!(df, fields_names, expander)
end

