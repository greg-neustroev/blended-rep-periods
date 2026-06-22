export fit_rep_period_weights!, projected_gradient_descent!, project_onto_simplex

"""
  project_onto_simplex(vector)

Projects `vector` onto a unit simplex using Michelot's algorithm in
Condat's accelerated implementation (2017). See Figure 2 of
[Condat, L. _Fast projection onto the simplex and the  ball._ Math. Program. 158,
575–585 (2016).](https://doi.org/10.1007/s10107-015-0946-6). For the details on
the meanings of v, ṽ, ρ and other variables, see the original paper.

# Examples

```
julia> project_onto_simplex([0.5, 0.5, 0.5])
3-element Vector{Float64}:
 0.33333333333333337
 0.33333333333333337
 0.33333333333333337

julia> project_onto_simplex([2.0, 0.0, 0.0])
3-element Vector{Float64}:
 1.0
 0.0
 0.0
```
"""
function project_onto_simplex(vector::AbstractVector{Float64})
  # There is a trivial solution when it's a one-element vector
  if length(vector) == 1
    return [1.0]
  end
  # step 1
  v = [vector[1]]
  ṽ = Vector{Float64}()
  ρ = vector[1] - 1.0
  # step 2
  for y ∈ vector[2:end]
    if y > ρ
      ρ += (y - ρ) / (length(v) + 1)
      if ρ > y - 1.0
        push!(v, y)
      else
        append!(ṽ, v)
        v = [y]
        ρ = y - 1.0
      end
    end
  end
  # step 3
  for y ∈ ṽ
    if y > ρ
      push!(v, y)
      ρ += (y - ρ) / length(v)
    end
  end
  # step 4
  while true
    length_v = length(v)
    to_be_removed = Vector{Int}()
    for (i, y) ∈ enumerate(v)
      if y ≤ ρ
        push!(to_be_removed, i)
        length_v -= 1
        ρ += (ρ - y) / length_v
      end
    end
    if length(to_be_removed) == 0
      break
    end
    deleteat!(v, to_be_removed)
  end
  # step 5 is skipped because it computes the data that we do not need
  # step 6
  return max.(vector .- ρ, 0.0)
end

"""
  project_onto_nonnegative_orthant(vector)

Projects `vector` onto the nonnegative_orthant. This projection is trivial:
replace negative components of the vector with zeros.

# Examples

```
julia> project_onto_nonnegative_orthant([-1.0, 2.0, -3.0])
3-element Vector{Float64}:
 0.0
 2.0
 0.0
```
"""
function project_onto_nonnegative_orthant(vector::AbstractVector{Float64})
  return max.(vector, 0.0)
end

"""
  projected_gradient_descent!(x; gradient, projection, niters, rtol, learning_rate)

Fits `x` using the projected gradient descent scheme.

The arguments:

  - `x`: the value to fit
  - `gradient`: the gradient operator, that is, a function that takes
    vectors of the same shape as `x` as inputs and returns a gradient of the
    loss at that point; the fitting is done to minimize the corresponding
    implicit loss
  - `projection`: the projection operator, that is, a function that, given a
    vector `x`, finds a point within some subspace that is closest to `x`
  - `niters`: maximum number of projected gradient descent iterations
  - `tol`: tolerance; when no components of `x` improve by more than `tol`, the
    algorithm stops
  - `learning_rate`: learning rate of the algorithm
"""
function projected_gradient_descent!(
  x::AbstractVector{Float64};
  gradient::Function,
  projection::Function,
  niters::Int=1000,
  tol::Float64=1e-5,
  learning_rate::Float64=1e-3,
)
  # It is possible that the initial guess is not in the required subspace;
  # project it first.
  x = projection(x)

  for _ ∈ 1:niters
    g = gradient(x)  # find the gradient
    α = learning_rate
    y = x .- α .* g            # gradient step, may leave the domain
    x_new = projection(y)      # projection step, return to the domain

    diff = maximum(abs.(x_new .- x))  # ‖x_prev − x‖_∞: how much the vector moved
    if diff ≤ tol / niters
      break
    end

    x = x_new
  end
  return x
end

"""
  fit_rep_period_weights!(weight_matrix, clustering_matrix, rp_matrix; weight_type, tol, args...)

Given the initial weight guesses, finds better weights for convex or conical
combinations of representative periods. For conical weights, it is possible to
bound the total weight by one.

The arguments:

  - `weight_matrix`: the initial guess for weights; the weights are adjusted
    using a projected gradient descent method
  - `clustering_matrix`: the matrix of raw clustering data
  - `rp_matrix`: the matrix of raw representative period data
  - `weight_type`: the type of weights to find; possible values are:
      - `:convex`: each period is represented as a convex sum of the
        representative periods (a sum with nonnegative weights adding into one)
      - `:conical`: each period is represented as a conical sum of the
        representative periods (a sum with nonnegative weights)
      - `:conical_bounded`: each period is represented as a conical sum of the
        representative periods (a sum with nonnegative weights) with the total
        weight bounded from above by one.
  - `tol`: algorithm's tolerance; when the weights are adjusted by a value less
    then or equal to `tol`, they stop being fitted further.
  - other arguments control the projected gradient method; they are passed
    through to `projected_gradient_descent!`.
"""
function fit_rep_period_weights!(
  weight_matrix::Union{SparseMatrixCSC{Float64,Int64},Matrix{Float64}},
  clustering_matrix::Matrix{Float64},
  rp_matrix::Matrix{Float64};
  weight_type::Symbol=:dirac,
  tol::Float64=1e-2,
  args...,
)
  # Determine the appropriate projection method
  if weight_type ≡ :dirac
    return weight_matrix
  elseif weight_type ≡ :convex
    projection = project_onto_simplex
  elseif weight_type ≡ :conical
    projection = project_onto_nonnegative_orthant
  elseif weight_type ≡ :conical_bounded
    # Conic bounded method does convex fitting, but adds a zero component.
    # The weight of a zero vector is then discarded without affecting the
    # total, and the resulting weights will always have sums between zero and
    # one.
    projection = project_onto_simplex
    n_data_points = size(rp_matrix, 1)
    rp_matrix = hcat(repeat([0.0], n_data_points), rp_matrix)
    #weight_matrix = hcat(repeat([0.0], size(weight_matrix, 1)), weight_matrix)
  else
    throw(ArgumentError("Unsupported weight type."))
  end

  n_periods = size(clustering_matrix, 2)
  n_rp = size(rp_matrix, 2)

  is_sparse = issparse(weight_matrix)

  default_weight_matrix = weight_matrix' |> Matrix{Float64}
  if weight_type ≡ :conical_bounded
    # A zero column needs to be added to the initial weight matrix
    default_weight_matrix = vcat(zeros(1, size(default_weight_matrix, 2)), default_weight_matrix)
  end
  for col in eachcol(default_weight_matrix)
    col .= projection(col)
  end
  moore_penrose_weight_matrix = rp_matrix \ clustering_matrix
  for col in eachcol(moore_penrose_weight_matrix)
    col .= projection(col)
  end
  # check which is the better initial guess and use it
  default_weight_projection_error = sum((rp_matrix * default_weight_matrix - clustering_matrix) .^ 2)
  moore_penrose_projection_error = sum((rp_matrix * moore_penrose_weight_matrix - clustering_matrix) .^ 2)
  initial_weight_matrix = default_weight_projection_error <= moore_penrose_projection_error ?
                          default_weight_matrix :
                          moore_penrose_weight_matrix

  # Principled PGD step size: 1 / L with L = σ_max(rp_matrix)^2, the Lipschitz
  # constant of the projection objective's gradient (Proposition 2).
  step_size = 1 / opnorm(rp_matrix, 2)^2

  for period ∈ 1:n_periods
    target_vector = clustering_matrix[:, period]
    x = projection(initial_weight_matrix[:, period])
    gradient = x -> rp_matrix' * (rp_matrix * x - target_vector)
    initial_projection_eror = norm(rp_matrix * x - target_vector)
    if initial_projection_eror ≤ tol
      continue
    end
    x = projected_gradient_descent!(x; gradient, projection, tol=tol * 0.01, learning_rate=step_size, args...)
    fitted_projection_error = norm(rp_matrix * x - target_vector)
    if fitted_projection_error > initial_projection_eror
      @warn "Projection error after fitting is larger than before fitting. Using the initial guess instead."
      @info "Initial projection error: $initial_projection_eror"
      @info "Fitted projection error: $fitted_projection_error"
      x = initial_weight_matrix[:, period]
    end
    x[x.<tol] .= 0.0  # replace insignificant small values with zeros
    if weight_type ≡ :convex || weight_type ≡ :conical_bounded
      # Because some values might have been removed, convexity can be lost.
      # In the upper-bounded case, sometimes the sum can be slightly more than one
      # due to floating-point arithmetic and rounding.
      # To account for these cases, the weights are re-normalized.
      sum_x = sum(x)
      if weight_type ≡ :convex || sum_x > 1.0
        x = x ./ sum_x
      end
    end
    if weight_type ≡ :conical_bounded
      popfirst!(x)
    end
    if is_sparse
      x = sparse(x)
    end
    weight_matrix[period, 1:length(x)] = x
  end
  return weight_matrix
end

"""
  fit_rep_period_weights!(weight_matrix, clustering_matrix, rp_matrix; weight_type, tol, args...)

  Given the initial weight guesses, finds better weights for convex or conical
combinations of representative periods. For conical weights, it is possible to
bound the total weight by one.

The arguments:

  - `clustering_result`: the result of running
    `find_representative_periods`
  - `weight_type`: the type of weights to find; possible values are:
      - `:convex`: each period is represented as a convex sum of the
        representative periods (a sum with nonnegative weights adding into one)
      - `:conical`: each period is represented as a conical sum of the
        representative periods (a sum with nonnegative weights)
      - `:conical_bounded`: each period is represented as a conical sum of the
        representative periods (a sum with nonnegative weights) with the total
        weight bounded from above by one.
  - `tol`: algorithm's tolerance; when the weights are adjusted by a value less
    then or equal to `tol`, they stop being fitted further.
  - other arguments control the projected gradient method; they are passed
    through to `projected_gradient_descent!`.
"""
function fit_rep_period_weights!(
  clustering_result::ClusteringResult;
  weight_type::Symbol=:dirac,
  tol::Float64=1e-2,
  args...,
)
  fit_rep_period_weights!(
    clustering_result.weight_matrix,
    clustering_result.clustering_matrix,
    clustering_result.rp_matrix;
    weight_type,
    tol,
    args...,
  )
end
