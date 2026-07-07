
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
    # Precompute the Gram matrix once per outer iteration and reuse it for every
    # candidate's projection gradient: G w − Rᵀc_d costs O(n_rp²) per PGD step
    # instead of two matmuls against the tall hull matrix (O(|features|·n_rp)).
    # Principled PGD step size: eigmax(G) = σ_max²(hull_matrix) = L, so α = 1/L
    # exactly as before.
    gram = Symmetric(hull_matrix' * hull_matrix)
    step_size = 1 / eigmax(gram)
    grad_buf = Vector{Float64}(undef, length(hull_indices))
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
        b = hull_matrix' * target_vector
        gradient = let gram = gram, b = b, grad_buf = grad_buf
          w -> (mul!(grad_buf, gram, w); grad_buf .-= b; grad_buf)
        end
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
  elseif method ≡ :chronological
    # Sequential linked time-period blocks: partition the base periods into `n_rp`
    # contiguous, (near-)equal-length segments in chronological order; each segment's
    # representative is its medoid (the member nearest the rest, like k-medoids), and
    # every period is assigned to its own segment. An order-preserving, non-clustering
    # baseline (no feature-space grouping — grouping is by time alone).
    n_cols = size(matrix, 2)
    distance_matrix = pairwise(Euclidean(), matrix; dims=2)
    edges = round.(Int, range(0, n_cols; length=n_rp + 1))
    assignments = Vector{Int}(undef, n_cols)
    medoids = Vector{Int}(undef, n_rp)
    for k ∈ 1:n_rp
      members = (edges[k]+1):edges[k+1]
      assignments[members] .= k
      medoids[k] = members[argmin([sum(@view distance_matrix[m, members]) for m ∈ members])]
    end
    rp_matrix = matrix[:, medoids]
    selected_indices = medoids
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
    build_inflow_integral_rows(clustering_matrix, keys; energy_scale=nothing, timestep_duration=1.0) -> Matrix{Float64}

Build the integrated-inflow feature rows `E_int` appended to the clustering matrix when
an inflow-integral weight λ is set. For every inflow asset `s` (a `keys` row with
`profile_type == "inflows"`) one row is produced from the per-period sum, over the
period's timesteps `h`, of its dimensionless inflow profile `Σ_h E_{s,d,h}` (the
`inflows` rows of `clustering_matrix`). Summing over timesteps collapses the time
dimension, so each inflow asset contributes a single seasonal-energy row. How that sum
is scaled depends on the selection space:

  - **Unscaled / min-max** (`energy_scale === nothing`): normalize by the per-period
    peak integral `H` (the asset's number of timesteps), giving the dimensionless
    period-mean profile `E_int[s,d] = (Σ_h E_{s,d,h}) / H ∈ [0,1]`. The physical `τ`
    and `E^max_s` cancel against this per-period-peak normalizer, so the row sits on
    the same `1 = peak` scale as the hourly profile block.
  - **Economic** (`energy_scale` given): keep the physical, cost-weighted energy
    `E_int[s,d] = τ · (V̄·E^max_s) · Σ_h E_{s,d,h}`, where `energy_scale[s] = V̄·E^max_s`
    is the asset's economic inflow scale and `τ` is `timestep_duration`. Here `E^max_s`
    does its inter-asset weighting natively, alongside the economically-scaled blocks.

Rows are ordered by first appearance of the asset id in `keys`. Returns a `0 × n_periods`
matrix when the data carries no inflow profile, and (economic only) throws if an inflow
asset has no `energy_scale` entry.

# Examples

```
julia> C = [1.0 2.0; 0.5 0.5; 0.2 0.4]  # timestep rows; last two are asset "h" inflows
3×2 Matrix{Float64}:
 1.0  2.0
 0.5  0.5
 0.2  0.4

julia> keys = DataFrame(timestep=[1, 1, 2], id=["d", "h", "h"],
                        profile_type=["demand", "inflows", "inflows"]);

julia> build_inflow_integral_rows(C, keys)  # unscaled period mean: (0.5+0.2)/2, (0.5+0.4)/2
1×2 Matrix{Float64}:
 0.35  0.45

julia> build_inflow_integral_rows(C, keys; energy_scale=Dict("h" => 10.0), timestep_duration=2.0)
1×2 Matrix{Float64}:
 14.0  18.0
```
"""
function build_inflow_integral_rows(
  clustering_matrix::AbstractMatrix{Float64},
  keys::AbstractDataFrame;
  energy_scale::Union{Nothing,AbstractDict}=nothing,
  timestep_duration::Float64=1.0,
)
  n_periods = size(clustering_matrix, 2)
  (hasproperty(keys, :profile_type) && hasproperty(keys, :id)) ||
    return Matrix{Float64}(undef, 0, n_periods)
  inflow_rows = findall(==("inflows"), string.(keys.profile_type))
  isempty(inflow_rows) && return Matrix{Float64}(undef, 0, n_periods)
  inflow_ids = string.(keys.id[inflow_rows])
  assets = unique(inflow_ids)  # one integrated row per inflow asset, in first-seen order
  integral = Matrix{Float64}(undef, length(assets), n_periods)
  for (i, asset) ∈ enumerate(assets)
    asset_rows = inflow_rows[inflow_ids.==asset]
    period_sum = vec(sum(view(clustering_matrix, asset_rows, :); dims=1))
    if energy_scale ≡ nothing
      # Unscaled / min-max: divide by H (the asset's timestep count) → period mean in
      # [0,1]; τ and E^max cancel against this per-period-peak normalizer.
      integral[i, :] = period_sum ./ length(asset_rows)
    else
      # Economic: physical, cost-weighted energy τ·(V̄·E^max_s)·Σ_h.
      haskey(energy_scale, asset) || error(
        "No economic inflow scale (V̄·E^max) for inflow asset $asset; cannot build " *
        "the integrated-inflow clustering rows.")
      integral[i, :] = (timestep_duration * energy_scale[asset]) .* period_sum
    end
  end
  return integral
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
  - `method`: clustering method to use; one of `:k_means`, `:k_medoids`,
    `:hierarchical`, `:chronological`, `:convex_hull`, `:convex_hull_with_null`,
    `:conical_hull` (selection is chosen independently of the weight class)
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
  - `inflow_integral_weight`: the weight λ on the integrated-inflow rows. When it is
    non-zero (and the data has inflow profiles), the clustering matrix `C = [D; A; E]`
    is augmented to `C = [D; A; E; λ·E_int]`, where each `E_int` row is an inflow
    asset's per-period inflow integral (see [`build_inflow_integral_rows`](@ref)). The
    row is the dimensionless period-mean profile under `:unscaled`/`:minmax`, and the
    physical cost-weighted energy `τ·V̄·E^max_s·Σ_h` under `:economic`. The integrated
    rows live only in the selection/weight-fitting space — appended *after* any
    `feature_scale`/`minmax` transform — so the representative-period profiles handed to
    the model stay in the original units. Defaults to `0.0` (no augmentation), which
    leaves every existing result unchanged and makes the row a clean ablation toggle.
  - `timestep_duration`: the timestep duration `τ`. It enters only the `:economic`
    integrated-inflow rows (under `:unscaled`/`:minmax` it cancels against the
    per-period-peak normalizer); defaults to `1.0` and is otherwise unused.
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
  inflow_integral_weight::Float64=0.0,
  timestep_duration::Float64=1.0,
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

  # 4a. Base selection space (profiles only), per normalization mode.
  scaled = feature_scale ≢ nothing || minmax
  if !scaled
    base_selection = clustering_matrix
  elseif feature_scale ≢ nothing
    # Economic: multiply each row by its (id, profile_type) scale; absolute, never
    # centered (zero stays physically zero). Drop rows constant across periods —
    # they carry no inter-period information and bias the origin-anchored conic
    # hulls; the test uses the *original* variation, so it is scale-independent.
    scale_vector = economic_scale_vector(feature_scale, keys)
    row_range = vec(maximum(clustering_matrix, dims=2) .- minimum(clustering_matrix, dims=2))
    keep = row_range .> 0.0
    base_selection = (scale_vector .* clustering_matrix)[keep, :]
  else
    # Min-max: rescale each row to [0,1] (its minimum across periods → 0, its
    # maximum → 1). Rows constant across periods have no range to rescale and are
    # dropped (which also avoids a 0/0).
    base_selection, keep = minmax_rescale(clustering_matrix)
  end
  scaled && log_scaled_selection_diagnostics(keys.profile_type, keep, base_selection)

  # 4b. Optional integrated-inflow augmentation: append the λ·E_int rows so each inflow
  # asset's per-period inflow integral (summed over timesteps) drives selection and
  # weight fitting alongside the time-resolved profiles. Under `:unscaled`/`:minmax` the
  # row is the dimensionless period-mean profile (τ, E^max cancel); under `:economic` it
  # is the physical cost-weighted energy τ·V̄·E^max·Σ_h, where `V̄·E^max` is read back
  # from `feature_scale`. The rows are appended *after* any economic/minmax transform and
  # live only in the selection/weight space — the profiles are reconstructed in original
  # units below. `λ = 0` (or no inflow data) leaves the base selection untouched, so the
  # weight is a clean on/off ablation toggle.
  if inflow_integral_weight ≠ 0.0
    energy_scale = feature_scale ≡ nothing ? nothing :
                   Dict(id => sc for ((id, pt), sc) ∈ feature_scale if pt == "inflows")
    inflow_integral = build_inflow_integral_rows(clustering_matrix, keys; energy_scale, timestep_duration)
  else
    inflow_integral = Matrix{Float64}(undef, 0, size(clustering_matrix, 2))
  end
  augmented = size(inflow_integral, 1) > 0
  selection_matrix = augmented ?
                     vcat(base_selection, inflow_integral_weight .* inflow_integral) :
                     base_selection
  if augmented
    selection_stats[:inflow_integral_weight] = inflow_integral_weight
    selection_stats[:n_inflow_integral_rows] = size(inflow_integral, 1)
  end

  # 4c. Select representatives in the (possibly scaled, possibly augmented) space.
  rp_matrix, assignments, selected_indices =
    select_representatives(selection_matrix, n_rp, n_periods, method; tol, cache, stats=selection_stats, args...)

  # 4d. Reconstruct the representative-period profiles in the original units the model
  # needs. The unscaled, un-augmented path returns the selection-space representatives
  # directly (they already are the original profiles); every transformed selection
  # space (economic, minmax, or the inflow-integral augmentation) instead maps back to
  # the original columns — selected base periods for the hull / k-medoids methods, or
  # the original-space centroids of the clusters for the synthetic k-means centers.
  if !scaled && !augmented
    representative_profiles = rp_matrix
  else
    representative_profiles = selected_indices ≡ nothing ?
                              original_space_centroids(clustering_matrix, assignments, n_rp) :
                              clustering_matrix[:, selected_indices]
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
  clustering_method = clustering_type_to_method(clustering_type)
  # Select the normalization. `:unscaled` (default) clusters on the profiles as
  # stored; `:minmax` rescales each feature row to [0,1]; `:economic` builds the
  # per-feature physical/economic scale and logs the dataset-level diagnostics.
  feature_scale = nothing
  minmax = false
  # τ enters only the `:economic` integrated-inflow rows; the economic scale builder
  # already reads and validates it, so reuse it from `info` rather than re-querying.
  timestep_duration = 1.0
  if normalization ≡ :economic
    feature_scale, info = build_economic_feature_scale(connection)
    log_economic_scale_info(info)
    timestep_duration = info.timestep_duration
  elseif normalization ≡ :minmax
    minmax = true
  elseif normalization ≢ :unscaled
    throw(ArgumentError("Unsupported normalization $(normalization); expected :unscaled, :minmax, or :economic"))
  end
  # Perform the clustering. A non-zero `inflow_integral_weight` λ augments the
  # clustering matrix with the integrated-inflow rows; `0.0` (the default) is the
  # un-augmented baseline arm of the ablation.
  find_representative_periods(
    clustering_df,
    n_rep_periods;
    method=clustering_method,
    tol=experiment_data.tol,
    feature_scale,
    minmax,
    cache=experiment_data.cache,
    inflow_integral_weight=experiment_data.inflow_integral_weight,
    timestep_duration,
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

function clustering_type_to_method(clustering_type)
    supported = (:k_means, :k_medoids, :hierarchical, :chronological,
        :convex_hull, :convex_hull_with_null, :conical_hull)
    clustering_type ∈ supported || throw(ArgumentError(
        "Unsupported clustering_type $(clustering_type); expected one of $(supported). " *
        "Selection is now specified independently of the weight class: name :convex_hull " *
        "or :conical_hull directly (the old :hull alias, whose hull geometry was inferred " *
        "from the weight type, has been removed)."))
    return clustering_type
end
