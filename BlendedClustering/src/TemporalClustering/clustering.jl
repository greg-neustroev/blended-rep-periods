
"""
    combine_periods!(df)

Modifies a dataframe `df` by combining the columns `timestep` and `period`
into a single column `timestep` of global time steps. The period duration is
inferred automatically from the maximum time step value, assuming that
periods start with time step 1.

# Examples

```
julia> df = DataFrame([:period => [1, 1, 2], :timestep => [1, 2, 1], :value => 1:3])
3×3 DataFrame
 Row │ period  timestep  value
     │ Int64   Int64      Int64
─────┼──────────────────────────
   1 │      1          1      1
   2 │      1          2      2
   3 │      2          1      3

julia> combine_periods!(df)
3×2 DataFrame
 Row │ timestep  value
     │ Int64      Int64
─────┼──────────────────
   1 │         1      1
   2 │         2      2
   3 │         3      3
```
"""
function combine_periods!(df::AbstractDataFrame)
  # First check that df contains a timestep column
  if columnindex(df, :timestep) == 0
    throw(DomainError(df, "DataFrame does not contain a column `timestep`"))
  end
  if columnindex(df, :period) == 0
    return  # if there is no column df.period, leave df as is
  end
  max_t = maximum(df.timestep)
  df.timestep .= (df.period .- 1) .* max_t .+ df.timestep
  select!(df, Not(:period))
end

"""
    split_into_periods!(df; period_duration=nothing)

Modifies a dataframe `df` by separating the column `timestep` into periods of
length `period_duration`. The new data is written into two columns:

  - `period`: the period ID;
  - `timestep`: the time step within the current period.

If `period_duration` is `nothing`, then all of the time steps are within the
same period with index 1.

# Examples

```
julia> df = DataFrame([:timestep => 1:4, :value => 5:8])
4×2 DataFrame
 Row │ timestep  value
     │ Int64      Int64
─────┼──────────────────
   1 │         1      5
   2 │         2      6
   3 │         3      7
   4 │         4      8

julia> split_into_periods!(df; period_duration=2)
4×3 DataFrame
 Row │ period  timestep  value
     │ Int64   Int64      Int64
─────┼──────────────────────────
   1 │      1          1      5
   2 │      1          2      6
   3 │      2          1      7
   4 │      2          2      8

julia> df = DataFrame([:period => [1, 1, 2], :timestep => [1, 2, 1], :value => 1:3])
3×3 DataFrame
 Row │ period  timestep  value
     │ Int64   Int64      Int64
─────┼──────────────────────────
   1 │      1          1      1
   2 │      1          2      2
   3 │      2          1      3

julia> split_into_periods!(df; period_duration=1)
3×3 DataFrame
 Row │ period  timestep  value
     │ Int64   Int64      Int64
─────┼──────────────────────────
   1 │      1          1      1
   2 │      2          1      2
   3 │      3          1      3

julia> split_into_periods!(df)
3×3 DataFrame
 Row │ period  timestep  value
     │ Int64   Int64      Int64
─────┼──────────────────────────
   1 │      1          1      1
   2 │      1          2      2
   3 │      1          3      3
```
"""
function split_into_periods!(df::AbstractDataFrame; period_duration::Union{Int,Nothing}=nothing)
  # If the periods already exist, combine them into the time steps if necessary
  combine_periods!(df)

  if isnothing(period_duration)
    # If period_duration is nothing, then leave the time steps as is and
    # everything is just the same period with index 1.
    insertcols!(df, :period => 1)
  else
    # Otherwise, split the time step index using 1-based modular arithmetic
    indices = fldmod1.(df.timestep, period_duration)  # find the new indices
    indices = reinterpret(reshape, Int, indices)       # change to an array for slicing

    df.period = indices[1, :]     # first row is the floor quotients, i.e., the period indices
    df.timestep = indices[2, :]  # second row is the remainders, i.e., the new time steps
  end
  select!(df, :period, :timestep, :)  # move the time-related columns to the front
end

"""
    validate_df_and_find_key_columns(df)

Checks that dataframe `df` contains the necessary columns and returns a list of
columns that act as keys (i.e., unique data identifiers within different periods).

# Examples

```
julia> df = DataFrame([:period => [1, 1, 2], :timestep => [1, 2, 1], :a .=> "a", :value => 1:3])
3×4 DataFrame
 Row │ period  timestep  a       value
     │ Int64   Int64      String  Int64
─────┼──────────────────────────────────
   1 │      1          1  a           1
   2 │      1          2  a           2
   3 │      2          1  a           3

julia> validate_df_and_find_key_columns(df)
2-element Vector{Symbol}:
 :timestep
 :a

julia> df = DataFrame([:value => 1])
1×1 DataFrame
 Row │ value
     │ Int64
─────┼───────
   1 │     1

julia> validate_df_and_find_key_columns(df)
ERROR: DomainError with 1×1 DataFrame
 Row │ value
     │ Int64
─────┼───────
   1 │     1:
DataFrame must contain columns `timestep` and `value`
```
"""
function validate_df_and_find_key_columns(df::AbstractDataFrame)::Vector{Symbol}
  columns = propertynames(df)
  if :timestep ∉ columns || :value ∉ columns
    throw(DomainError(df, "DataFrame must contain columns `timestep` and `value`"))
  end
  if :period ∉ columns
    throw(
      DomainError(
        df,
        "DataFrame must contain column `period`; call split_into_periods! to split it into periods.",
      ),
    )
  end
  non_key_columns = [:period, :value]
  key_columns = filter!(col -> col ∉ non_key_columns, columns)
  return key_columns
end

"""
    find_auxiliary_data(clustering_data)

Calculates auxiliary data associated with the `clustering_data`. These include:

  - `key_columns_demand`: key columns in the demand dataframe
  - `key_columns_generation_availability`: key columns in the generation availability dataframe
  - `period_duration`: duration of time periods (in time steps)
  - `last_period_duration`: duration of the last period
  - `n_periods`: total number of periods
"""
function find_auxiliary_data(clustering_data::AbstractDataFrame)
  key_columns = validate_df_and_find_key_columns(clustering_data)
  n_periods = maximum(clustering_data.period)
  period_duration = maximum(clustering_data.timestep)
  last_period_duration = maximum(clustering_data[clustering_data.period.==n_periods, :timestep])

  return AuxiliaryClusteringData(key_columns, period_duration, last_period_duration, n_periods)
end

"""
    find_period_weights(period_duration, last_period_duration, n_periods, drop_incomplete_periods)

Finds weights of two different types of periods in the clustering data:

  - complete periods: these are all of the periods with length equal to `period_duration`.
  - incomplete last period: if last period duration is less than `period_duration`, it is incomplete.
"""
function find_period_weights(
  period_duration::Int,
  last_period_duration::Int,
  n_periods::Int,
  drop_incomplete_periods::Bool,
)::Tuple{Float64,Union{Float64,Nothing}}
  if last_period_duration == period_duration
    complete_period_weight = 1.0
    incomplete_period_weight = nothing
  elseif drop_incomplete_periods
    full_period_timesteps = period_duration * (n_periods - 1)
    total_timesteps = full_period_timesteps + last_period_duration
    complete_period_weight = total_timesteps / full_period_timesteps
    incomplete_period_weight = nothing
  else
    complete_period_weight = 1.0
    incomplete_period_weight = 1.0
  end
  return complete_period_weight, incomplete_period_weight
end

"""
    df_to_matrix_and_keys(df, key_columns)

Converts a dataframe `df` (in a long format) to a matrix, ignoring the columns
specified as `key_columns`. The key columns are converted from long to wide
format and returned alongside the matrix.

# Examples

```
julia> df = DataFrame([:period => [1, 1, 2, 2], :timestep => [1, 2, 1, 2], :a .=> "a", :value => 1:4])
4×4 DataFrame
 Row │ period  timestep  a       value
     │ Int64   Int64      String  Int64
─────┼──────────────────────────────────
   1 │      1          1  a           1
   2 │      1          2  a           2
   3 │      2          1  a           3
   4 │      2          2  a           4

julia> m, k = df_to_matrix_and_keys(df, [:timestep, :a]); m
2×2 Matrix{Float64}:
 1.0  3.0
 2.0  4.0

julia> k
2×2 DataFrame
 Row │ timestep  a
     │ Int64      String
─────┼───────────────────
   1 │         1  a
   2 │         2  a
```
"""
function df_to_matrix_and_keys(df::AbstractDataFrame, key_columns::Vector{Symbol})
  wide_df = unstack(df, key_columns, :period, :value)
  matrix = select(wide_df, Not(key_columns)) |> dropmissing |> Matrix{Float64}
  keys = select(wide_df, key_columns)
  return matrix, keys
end

"""
    matrix_and_keys_to_df(matrix, keys)

Converts a a matrix `matrix` to a dataframe, appending the key columns given by
`keys`.

# Examples

```
julia> m = [1.0 3.0; 2.0 4.0]
2×2 Matrix{Float64}:
 1.0  3.0
 2.0  4.0

julia> k = DataFrame([:timestep => 1:2, :a .=> "a"])
2×2 DataFrame
 Row │ timestep  a
     │ Int64      String
─────┼───────────────────
   1 │         1  a
   2 │         2  a

julia> matrix_and_keys_to_df(m, k)
4×4 DataFrame
 Row │ rep_period  timestep  a       value
     │ Int64       Int64      String  Float64
─────┼────────────────────────────────────────
   1 │          1          1  a           1.0
   2 │          1          2  a           2.0
   3 │          2          1  a           3.0
   4 │          2          2  a           4.0
```
"""
function matrix_and_keys_to_df(matrix::Matrix{Float64}, keys::AbstractDataFrame)
  n_columns = size(matrix, 2)
  result = DataFrame(matrix, string.(1:n_columns))
  result = hcat(keys, result)            # prepend the previously deleted columns
  result = stack(result; variable_name=:rep_period) |> dropmissing # convert from wide to long format
  result.rep_period = parse.(Int, result.rep_period)  # change the type of rep_period column to Int
  select!(result, :rep_period, :timestep, :)         # move the rep_period column to the front

  return result
end

"""
    append_period_from_source_df_as_rp!(df; source_df, period, rp, key_columns)

Extracts a period with index `period` from `source_df` and appends it as a
representative period with index `rp` to `df`, using `key_columns` as keys.

# Examples

```
julia> source_df = DataFrame([:period => [1, 1, 2, 2], :timestep => [1, 2, 1, 2], :a .=> "b", :value => 5:8])
4×4 DataFrame
 Row │ period  timestep  a       value
     │ Int64   Int64      String  Int64
─────┼──────────────────────────────────
   1 │      1          1  b           5
   2 │      1          2  b           6
   3 │      2          1  b           7
   4 │      2          2  b           8

julia> df = DataFrame([:rep_period => [1, 1, 2, 2], :timestep => [1, 2, 1, 2], :a .=> "a", :value => 1:4])
4×4 DataFrame
 Row │ rep_period  timestep  a       value
     │ Int64       Int64      String  Int64
─────┼──────────────────────────────────────
   1 │          1          1  a           1
   2 │          1          2  a           2
   3 │          2          1  a           3
   4 │          2          2  a           4

julia> append_period_from_source_df_as_rp!(df; source_df, period = 2, rp = 3, key_columns = [:timestep, :a])
6×4 DataFrame
 Row │ rep_period  timestep  a       value
     │ Int64       Int64      String  Int64
─────┼──────────────────────────────────────
   1 │          1          1  a           1
   2 │          1          2  a           2
   3 │          2          1  a           3
   4 │          2          2  a           4
   5 │          3          1  b           7
   6 │          3          2  b           8
```
"""
function append_period_from_source_df_as_rp!(
  df::AbstractDataFrame;
  source_df::AbstractDataFrame,
  period::Int,
  rp::Int,
  key_columns::Vector{Symbol},
)
  period_df = source_df[source_df.period.==period, :]
  period_df.period .= rp
  select!(period_df, :period => :rep_period, key_columns..., :value)
  append!(df, period_df)
end

"""
    greedy_convex_hull(matrix; n_points, initial_indices=nothing, mean_vector=nothing, cache=true)

Greedily selects `n_points` columns of `matrix` whose convex hull approximates the
data. Starting from the column furthest from the column mean (or from
`initial_indices`, if given), each step adds the column whose Euclidean distance
to the current hull is largest, where that distance is the residual of the
projection onto the hull. Returns the selected column indices.

All distances are Euclidean. When `cache` is `true`, a column's projection
distance is reused across iterations whenever the obtuse-angle certificate
guarantees the cached projection is still exact; `cache=false`
recomputes every projection and yields identical results.

# Examples

```
julia> M = [0.0 1.0 0.0 0.3; 0.0 0.0 1.0 0.3]  # three hull vertices and one interior point
2×4 Matrix{Float64}:
 0.0  1.0  0.0  0.3
 0.0  0.0  1.0  0.3

julia> greedy_convex_hull(M; n_points=3)  # interior point (column 4) is never selected
3-element Vector{Int64}:
 2
 3
 1
```
"""
function greedy_convex_hull(
  matrix::AbstractMatrix{Float64};
  n_points::Int,
  initial_indices::Union{Vector{Int},Nothing}=nothing,
  mean_vector::Union{Vector{Float64},Nothing}=nothing,
  cache::Bool=true,
  tol::Float64=1e-2,
  stats::Union{Nothing,Dict{Symbol,Any}}=nothing,
)
  if initial_indices ≡ nothing
    if mean_vector ≡ nothing
      mean_vector = vec(mean(matrix, dims=2))
    end
    distances_from_mean = [norm(mean_vector - matrix[:, j]) for j ∈ axes(matrix, 2)]
    initial_indices = [argmax(distances_from_mean)]
  end
  if length(initial_indices) ≥ n_points
    return initial_indices[1:n_points]
  end
  hull_indices = initial_indices
  # Cache maps a candidate index d to (q_d, δ_d): its Euclidean projection onto
  # the current hull and the associated distance. After a new representative
  # c_new is added, the cached projection q_d stays exact for the enlarged hull
  # iff the obtuse-angle certificate (c_d − q_d)ᵀ(c_new − q_d) ≤ 0 holds;
  # otherwise it is recomputed.
  projection_cache = Dict{Int,Tuple{Vector{Float64},Float64}}()
  starting_index = length(initial_indices) + 1
  # Cache accounting: per-outer-iteration hit/miss counts (so the hit-rate-vs-
  # iteration curve and the total can be reported) when `stats` is provided.
  hits_per_iter = Int[]
  misses_per_iter = Int[]
  for _ ∈ starting_index:n_points
    max_distance = -Inf
    furthest_vector_index = nothing
    hull_matrix = matrix[:, hull_indices]
    projection_matrix = pinv(hull_matrix)
    # Principled PGD step size: 1 / L with L = σ_max(hull_matrix)^2, the
    # Lipschitz constant of the projection objective's gradient.
    step_size = 1 / opnorm(hull_matrix, 2)^2
    last_added_vector = matrix[:, last(hull_indices)]
    iter_hits = 0
    iter_misses = 0
    # Soundness requires re-testing EVERY candidate against the newest
    # representative on EVERY outer iteration; do not skip candidates here.
    for column_index ∈ axes(matrix, 2)
      if column_index ∈ hull_indices
        continue
      end
      target_vector = matrix[:, column_index]
      cached = get(projection_cache, column_index, nothing)
      if cache && cached !== nothing &&
         dot(target_vector - cached[1], last_added_vector - cached[1]) ≤ 0
        # Cached projection is still exact for the enlarged hull.
        d = cached[2]
        iter_hits += 1
      else
        gradient = x -> hull_matrix' * (hull_matrix * x - target_vector)
        x = projection_matrix * target_vector
        x, _ = projected_gradient_descent!(x; gradient, projection=project_onto_simplex, learning_rate=step_size, tol)
        projected_target = hull_matrix * x
        d = norm(projected_target - target_vector)
        projection_cache[column_index] = (projected_target, d)
        iter_misses += 1
      end
      if d > max_distance
        max_distance = d
        furthest_vector_index = column_index
      end
    end
    push!(hits_per_iter, iter_hits)
    push!(misses_per_iter, iter_misses)
    if furthest_vector_index ≡ nothing
      throw(ArgumentError("Point not found"))
    end
    push!(hull_indices, furthest_vector_index)
  end
  if stats ≢ nothing
    stats[:cache_hits] = get(stats, :cache_hits, 0) + sum(hits_per_iter)
    stats[:cache_misses] = get(stats, :cache_misses, 0) + sum(misses_per_iter)
    stats[:cache_hits_per_iter] = hits_per_iter
    stats[:cache_misses_per_iter] = misses_per_iter
  end
  return hull_indices
end

"""
    select_representatives(matrix, n_rp, n_periods, method; tol, args...)

Select `n_rp` representative columns of `matrix` using `method` and assign every
one of the `n_periods` base periods to a representative. Returns the tuple
`(rp_matrix, assignments, selected_indices)`:

  - `rp_matrix`: the representative columns (or k-means centroids) in the space
    of `matrix`;
  - `assignments`: for each base period, the index of its nearest representative;
  - `selected_indices`: the column indices chosen as representatives, or `nothing`
    for `:k_means` (whose representatives are synthetic centroids, not columns).

This is the method dispatch shared by every clustering path. `matrix` is the
space in which selection and nearest-representative assignment happen; callers
that cluster in a transformed feature space (e.g. economic normalization) pass
that space here and reconstruct the output profiles separately.
"""
function select_representatives(
  matrix::AbstractMatrix{Float64},
  n_rp::Int,
  n_periods::Int,
  method::Symbol;
  tol::Float64=1e-2,
  cache::Bool=true,
  stats::Union{Nothing,Dict{Symbol,Any}}=nothing,
  args...,
)
  if method ≡ :k_means
    kmeans_result = kmeans(matrix, n_rp; distance=Euclidean(), args...)
    rp_matrix = kmeans_result.centers
    assignments = kmeans_result.assignments
    selected_indices = nothing
  elseif method ≡ :k_medoids
    # k-medoids uses distance matrix instead of clustering matrix
    distance_matrix = pairwise(Euclidean(), matrix; dims=2)
    kmedoids_result = kmedoids(distance_matrix, n_rp; args...)
    rp_matrix = matrix[:, kmedoids_result.medoids]
    assignments = kmedoids_result.assignments
    selected_indices = kmedoids_result.medoids
  elseif method ≡ :hierarchical
    # Agglomerative (Ward) clustering on the same Euclidean distances as
    # k-medoids; each cluster's representative is its medoid (the member nearest
    # the rest), so the representatives are real periods like the k-medoids/hull ones.
    distance_matrix = pairwise(Euclidean(), matrix; dims=2)
    assignments = cutree(hclust(distance_matrix; linkage=:ward); k=n_rp)
    medoids = Vector{Int}(undef, n_rp)
    for c ∈ 1:n_rp
      members = findall(==(c), assignments)
      medoids[c] = members[argmin([sum(@view distance_matrix[m, members]) for m ∈ members])]
    end
    rp_matrix = matrix[:, medoids]
    selected_indices = medoids
  elseif method ≡ :convex_hull
    hull_indices = greedy_convex_hull(matrix; n_points=n_rp, tol, cache, stats)
    rp_matrix = matrix[:, hull_indices]
    assignments = [argmin([norm(matrix[:, h] - matrix[:, p]) for h ∈ hull_indices]) for p ∈ 1:n_periods]
    selected_indices = hull_indices
  elseif method ≡ :convex_hull_with_null
    # Add a null vector as a column, run the convex-hull method initialized at
    # it, then drop it; the remaining positive weights sum to at most one.
    augmented = [zeros(size(matrix, 1), 1) matrix]
    hull_indices = greedy_convex_hull(augmented; n_points=n_rp + 1, initial_indices=[1], tol, cache, stats)
    popfirst!(hull_indices)
    hull_indices .-= 1
    rp_matrix = matrix[:, hull_indices]
    assignments = [argmin([norm(matrix[:, h] - matrix[:, p]) for h ∈ hull_indices]) for p ∈ 1:n_periods]
    selected_indices = hull_indices
  elseif method ≡ :conical_hull
    normal_vector = vec(mean(matrix, dims=2))
    normalize!(normal_vector)
    projection_coefficients = [1.0 / dot(normal_vector, matrix[:, j]) for j ∈ axes(matrix, 2)]
    projected_matrix = [matrix[i, j] * projection_coefficients[j] for i ∈ axes(matrix, 1), j ∈ axes(matrix, 2)]
    hull_indices = greedy_convex_hull(projected_matrix; n_points=n_rp, mean_vector=normal_vector, tol, cache, stats)
    rp_matrix = matrix[:, hull_indices]
    assignments = [argmin([norm(matrix[:, h] - matrix[:, p]) for h ∈ hull_indices]) for p ∈ 1:n_periods]
    selected_indices = hull_indices
  else
    throw(ArgumentError("Clustering method is not supported"))
  end
  return rp_matrix, assignments, selected_indices
end

"""
    original_space_centroids(matrix, assignments, n_rp)

Recompute the `n_rp` cluster centroids of `matrix` from `assignments` (the column →
cluster map). Used by the normalized paths to express k-means representatives —
which are centroids of the *transformed* feature space — back in the original units
the model needs: the centroid of the original columns of a cluster equals the
original-space image of that cluster's transformed centroid. Throws on an empty
cluster rather than emitting an all-zero representative (a fabricated period of zero
everything), which the unscaled path's `kmeans.centers` would never produce.
"""
function original_space_centroids(matrix::AbstractMatrix{Float64}, assignments::AbstractVector{<:Integer}, n_rp::Int)
  centroids = zeros(size(matrix, 1), n_rp)
  for k ∈ 1:n_rp
    members = findall(==(k), assignments)
    isempty(members) && throw(ArgumentError(
      "k-means produced an empty cluster ($k of $n_rp) in the normalized selection " *
      "space, so its representative profile cannot be reconstructed in original units. " *
      "Use fewer representative periods or a hull/k-medoids method."))
    centroids[:, k] = vec(mean(view(matrix, :, members), dims=2))
  end
  return centroids
end

"""
    minmax_rescale(matrix) -> (scaled, keep)

Min-max normalize each row of `matrix` to `[0,1]` (row minimum → 0, row maximum →
1), dropping rows that are constant across periods (zero range — they carry no
inter-period information and would give a 0/0). Returns the `scaled` matrix over the
kept rows and the `keep` mask. Unlike the economic scaling, this is an affine per-row
transform: it *centers* (subtracts the row minimum), so it is used only for
selection/weight fitting, never for the profiles returned to the model.
"""
function minmax_rescale(matrix::AbstractMatrix{Float64})
  row_min = vec(minimum(matrix, dims=2))
  row_range = vec(maximum(matrix, dims=2)) .- row_min
  keep = row_range .> 0.0
  scaled = (matrix[keep, :] .- row_min[keep]) ./ row_range[keep]
  return scaled, keep
end

"""
  find_representative_periods(
    clustering_data;
    n_rp = 10,
    drop_incomplete_last_period = false,
    method = :k_means,
    args...,
  )

Finds representative periods via data clustering. All distances are Euclidean.

  - `clustering_data`: the data to perform clustering on.
  - `n_rp`: number of representative periods to find.
  - `drop_incomplete_last_period`: controls how the last period is treated if it
    is not complete: if this parameter is set to `true`, the incomplete period
    is dropped and the weights are rescaled accordingly; otherwise, clustering
    is done for `n_rp - 1` periods, and the last period is added as a special
    shorter representative period
  - `method`: clustering method to use, either `:k_means` and `:k_medoids`
  - `tol`: the projected gradient descent tolerance `ε` used by the hull methods
    to rank candidate periods by their distance to the current hull (ignored by
    `:k_means`/`:k_medoids`).
  - `feature_scale`: optional economic normalization. When `nothing` (the default),
    clustering runs on the dimensionless profiles as-is. When a `Dict` mapping
    `(asset id, profile_type) => scale`, every feature row is multiplied by its
    scale and rows that are constant across periods (zero range) are dropped *for
    selection and weight fitting only* — the representative-period profiles returned
    for the model are always reconstructed in the original, unscaled units (see
    [`build_economic_feature_scale`](@ref)).
  - `minmax`: optional per-feature-row min-max normalization. When `true`, each
    feature row is rescaled to `[0,1]` (its minimum across periods → 0, its maximum
    → 1) for selection and weight fitting only; the returned profiles stay in the
    original units. Mutually exclusive with `feature_scale`.
  - other named arguments can be provided; they are passed to the clustering method.
"""
function find_representative_periods(
  clustering_data::AbstractDataFrame,
  n_rp::Int;
  drop_incomplete_last_period::Bool=false,
  method::Symbol=:k_means,
  tol::Float64=1e-2,
  feature_scale::Union{Nothing,AbstractDict}=nothing,
  minmax::Bool=false,
  cache::Bool=true,
  args...,
)
  feature_scale ≡ nothing || !minmax ||
    throw(ArgumentError("feature_scale and minmax normalization are mutually exclusive"))

  # Find auxiliary data and pre-compute additional constants that are used multiple times alter
  aux = find_auxiliary_data(clustering_data)
  has_incomplete_last_period = aux.last_period_duration ≠ aux.period_duration
  is_last_period_excluded = has_incomplete_last_period && !drop_incomplete_last_period
  n_periods = aux.n_periods
  n_complete_periods = has_incomplete_last_period ? n_periods - 1 : n_periods

  # 2. Find the weights of the two types of periods and pre-build the weight matrix.
  # We assume that the only period that can be incomplete (i.e., has a duration
  # that is less than aux.period_duration) is the very last one. All other periods
  # are complete periods.
  complete_period_weight, incomplete_period_weight = find_period_weights(
    aux.period_duration,
    aux.last_period_duration,
    n_periods,
    drop_incomplete_last_period,
  )

  # In both cases, the weights of the complete periods will be found after clustering.
  if is_last_period_excluded
    weight_matrix = sparse([n_periods], [n_rp], [incomplete_period_weight])
    n_rp -= 1  # incomplete last period becomes its own representative, exclude it from clustering
  else
    weight_matrix = spzeros(n_complete_periods, n_rp)
  end

  # 3. Build the clustering matrix

  # First, find the demand matrix and rescale it if needed
  clustering_matrix, keys = df_to_matrix_and_keys(
    clustering_data[clustering_data.period.≤n_complete_periods, :],
    aux.key_columns,
  )

  # 4. Build the selection space — the features clustering and weight-fitting run
  # on. `:unscaled` (the default) clusters on the profiles as stored. `:economic`
  # (a `feature_scale` dict) and `:minmax` instead cluster on a transformed space,
  # but always reconstruct the representative-period profiles in the original units
  # the model expects — only the *selection* and the fitted *weights* live in the
  # transformed space, never the profiles handed to the model.
  # Selection diagnostics (greedy-hull cache hit/miss counts) accumulate here and
  # are attached to the returned ClusteringResult.
  selection_stats = Dict{Symbol,Any}()
  if feature_scale ≡ nothing && !minmax
    rp_matrix, assignments, selected_indices = select_representatives(clustering_matrix, n_rp, n_periods, method; tol, cache, stats=selection_stats, args...)
    representative_profiles = rp_matrix
    selection_matrix = clustering_matrix
  else
    if feature_scale ≢ nothing
      # Economic: multiply each row by its (id, profile_type) scale; absolute, never
      # centered (zero stays physically zero). Drop rows constant across periods —
      # they carry no inter-period information and bias the origin-anchored conic
      # hulls; the test uses the *original* variation, so it is scale-independent.
      scale_vector = economic_scale_vector(feature_scale, keys)
      row_range = vec(maximum(clustering_matrix, dims=2) .- minimum(clustering_matrix, dims=2))
      keep = row_range .> 0.0
      selection_matrix = (scale_vector .* clustering_matrix)[keep, :]
    else
      # Min-max: rescale each row to [0,1] (its minimum across periods → 0, its
      # maximum → 1). Rows constant across periods have no range to rescale and are
      # dropped (which also avoids a 0/0).
      selection_matrix, keep = minmax_rescale(clustering_matrix)
    end
    rp_matrix, assignments, selected_indices =
      select_representatives(selection_matrix, n_rp, n_periods, method; tol, cache, stats=selection_stats, args...)
    # Original-unit representative profiles: the selected period columns for the
    # hull / k-medoids methods, or the original-space centroids of the transformed
    # clusters for k-means (whose representatives are synthetic, not columns).
    representative_profiles = selected_indices ≡ nothing ?
                              original_space_centroids(clustering_matrix, assignments, n_rp) :
                              clustering_matrix[:, selected_indices]
    log_scaled_selection_diagnostics(keys.profile_type, keep, selection_matrix)
  end

  # Record the nearest-representative assignment and the selected base-period
  # indices (hull / k-medoids; `nothing` for synthetic k-means centroids) as
  # clustering artifacts so they can be dumped without re-running the selection.
  selection_stats[:assignments] = assignments
  selection_stats[:selected_indices] = selected_indices

  # Fill in the weight matrix using the assignments

  for (p, rp) ∈ enumerate(assignments)
    weight_matrix[p, rp] = complete_period_weight
  end

  # 5. Reinterpret the clustering results into a format we need

  # First, convert the matrix data back to dataframes using the previously saved key columns
  rp_df = matrix_and_keys_to_df(representative_profiles, keys)

  # Next, re-append the last period if it was excluded from clustering
  if is_last_period_excluded
    n_rp += 1
    append_period_from_source_df_as_rp!(
      rp_df;
      source_df=clustering_data,
      period=n_periods,
      rp=n_rp,
      key_columns=aux.key_columns,
    )
  end

  # `selection_matrix`/`rp_matrix` are the space weights are fitted in: the original
  # profiles for `:unscaled` (where `selection_matrix == clustering_matrix`), the
  # scaled, pruned features for `:economic`.
  result = ClusteringResult(rp_df, weight_matrix, selection_matrix, rp_matrix)
  merge!(result.diagnostics, selection_stats)
  return result
end

function cluster_using_experiment_data(experiment_data, connection)
  # Extract parameters from the experiment data into local variables
  n_rep_periods = experiment_data.n_rep_periods
  clustering_type = experiment_data.clustering_type
  weight_type = experiment_data.weight_type
  normalization = experiment_data.normalization

  # Fast path: when a single representative period is requested and the data
  # already consists of a single period (i.e. the period length covers the
  # whole horizon), clustering is a no-op — the representative period *is* the
  # full profile. Skip the profiles → matrix → k-means → DataFrame round-trip.
  if n_rep_periods == 1 && count_profile_periods(connection) == 1
    return single_period_clustering_result(connection)
  end

  # Collect the clustering data into a data frame
  clustering_df = get_clustering_profiles(connection)
  # Determine the clustering method
  clustering_method = clustering_type_to_method(
    clustering_type,
    weight_type
  )
  # Select the normalization. `:unscaled` (default) clusters on the profiles as
  # stored; `:minmax` rescales each feature row to [0,1]; `:economic` builds the
  # per-feature physical/economic scale and logs the dataset-level diagnostics.
  feature_scale = nothing
  minmax = false
  if normalization ≡ :economic
    feature_scale, info = build_economic_feature_scale(connection)
    log_economic_scale_info(info)
  elseif normalization ≡ :minmax
    minmax = true
  elseif normalization ≢ :unscaled
    throw(ArgumentError("Unsupported normalization $(normalization); expected :unscaled, :minmax, or :economic"))
  end
  # Perform the clustering
  find_representative_periods(
    clustering_df,
    n_rep_periods;
    method=clustering_method,
    tol=experiment_data.tol,
    feature_scale,
    minmax,
    cache=experiment_data.cache,
    init=:kmcen
  )
end

"""
    single_period_clustering_result(connection)

Build the trivial `ClusteringResult` for the case of a single representative
period that spans the entire horizon. The representative profile is the full
profile data relabelled with `rep_period = 1`, and the single period represents
itself with weight `1.0` (matching `find_period_weights` for one complete
period). The `1×1` `rp_matrix`/`clustering_matrix` make the projection error
zero and keep `fit_rep_period_weights!` well-defined for every weight type.
"""
function single_period_clustering_result(connection)
  rp_df = get_single_period_profiles(connection)
  rp_df.rep_period = Int.(rp_df.rep_period)
  weight_matrix = ones(1, 1)
  clustering_matrix = ones(1, 1)
  rp_matrix = ones(1, 1)
  return ClusteringResult(rp_df, weight_matrix, clustering_matrix, rp_matrix)
end

function clustering_type_to_method(clustering_type, weight_type)
    if clustering_type ≡ :hull
        if weight_type ≡ :conical
            :conical_hull
        elseif weight_type ≡ :conical_bounded
            :convex_hull_with_null
        else
            :convex_hull
        end
    else
        clustering_type
    end
end
