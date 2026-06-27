
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

# Numerical backstop on the PGD iteration count. In practice the resolution-based
# stopping test in `projected_gradient_descent!` fires first; this constant only
# rules out pathological non-termination from floating-point plateaus.
const PGD_MAX_ITERS = 100_000

"""
  projected_gradient_descent!(x; gradient, projection, tol, learning_rate)

Fits `x` using the projected gradient descent scheme.

`tol` is the resolution to which the components of `x` are determined: the
algorithm iterates until no component moves by more than `tol` between steps, i.e.
the solution has settled to within the precision that is kept downstream (anything
finer is discarded). No iteration count is imposed; the loop is bounded only by the
numerical backstop `PGD_MAX_ITERS`.

The arguments:

  - `x`: the value to fit
  - `gradient`: the gradient operator, that is, a function that takes
    vectors of the same shape as `x` as inputs and returns a gradient of the
    loss at that point; the fitting is done to minimize the corresponding
    implicit loss
  - `projection`: the projection operator, that is, a function that, given a
    vector `x`, finds a point within some subspace that is closest to `x`
  - `tol`: the resolution; the algorithm stops once no component of `x` changes
    by more than `tol` in an iteration
  - `learning_rate`: the step size `α` (the callers set it to `1/L` per instance)

Returns the tuple `(x, iterations)`: the fitted point and the number of iterations
the loop actually ran (so callers can record the observed `N(ε)` without imposing it).
"""
function projected_gradient_descent!(
  x::AbstractVector{Float64};
  gradient::Function,
  projection::Function,
  tol::Float64=1e-2,
  learning_rate::Float64=1e-3,
)
  # It is possible that the initial guess is not in the required subspace;
  # project it first.
  x = projection(x)

  # Certified stop (paper's Algorithm 2): the projected-gradient mapping residual
  # ‖x − x_new‖₂ equals α·‖G(x)‖₂, so stopping at ‖x − x_new‖₂ ≤ tol·α certifies
  # ‖G(x)‖₂ ≤ tol. Scaling by the step size α keeps the projections accurate
  # enough for the greedy-hull cache certificate (Lemma 1) to stay sound — a flat
  # ∞-norm move threshold is too coarse at the resolutions we use.
  threshold = tol * learning_rate
  iterations = 0
  for _ ∈ 1:PGD_MAX_ITERS
    iterations += 1
    g = gradient(x)              # find the gradient
    y = x .- learning_rate .* g  # gradient step, may leave the domain
    x_new = projection(y)        # projection step, return to the domain

    converged = norm(x_new .- x) ≤ threshold  # ‖projected-gradient mapping‖₂ ≤ tol
    x = x_new
    if converged
      break
    end
  end
  # Return the realized iteration count alongside the fitted point so callers can
  # record the observed N(ε) without imposing it.
  return x, iterations
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
  - `tol`: the resolution to which weights are determined; the single tolerance
    used throughout fitting. The projected gradient descent stops once no weight
    moves by more than `tol`, base periods already reconstructed within `tol` skip
    fitting, and fitted weights below `tol` are dropped (e.g. with `tol = 1e-2`
    weights are blending percentages and contributions under 1% are noise).
"""
function fit_rep_period_weights!(
  weight_matrix::Union{SparseMatrixCSC{Float64,Int64},Matrix{Float64}},
  clustering_matrix::Matrix{Float64},
  rp_matrix::Matrix{Float64};
  weight_type::Symbol=:dirac,
  tol::Float64=1e-2,
  diagnostics::Union{Nothing,Dict{Symbol,Any}}=nothing,
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

  is_sparse = issparse(weight_matrix)

  # Initial guess: the projected default weights (the seed passed in). The
  # Moore-Penrose least-squares alternative (and the auto-selector that chose
  # between them) was removed — it left the clustering-space projection error
  # essentially unchanged but produced markedly worse storage-regret trajectories,
  # and the selector always kept this default anyway.
  initial_weight_matrix = weight_matrix' |> Matrix{Float64}
  if weight_type ≡ :conical_bounded
    # A zero column needs to be added to the initial weight matrix.
    initial_weight_matrix = vcat(zeros(1, size(initial_weight_matrix, 2)), initial_weight_matrix)
  end
  for col in eachcol(initial_weight_matrix)
    col .= projection(col)
  end

  # Principled PGD step size: 1 / L with L = σ_max(rp_matrix)^2, the Lipschitz
  # constant of the projection objective's gradient.
  step_size = 1 / opnorm(rp_matrix, 2)^2

  # Record per-fit PGD iteration counts (one entry per base period that is
  # actually fitted, i.e. not already within `tol`), so the observed N(ε) range
  # can be reported without re-running.
  pgd_iters = Int[]

  for period ∈ 1:n_periods
    target_vector = clustering_matrix[:, period]
    x = projection(initial_weight_matrix[:, period])
    gradient = x -> rp_matrix' * (rp_matrix * x - target_vector)
    initial_projection_eror = norm(rp_matrix * x - target_vector)
    if initial_projection_eror ≤ tol
      continue  # already reconstructed within the resolution; keep the initial weights
    end
    x, n_iter = projected_gradient_descent!(x; gradient, projection, tol, learning_rate=step_size)
    push!(pgd_iters, n_iter)
    fitted_projection_error = norm(rp_matrix * x - target_vector)
    if fitted_projection_error > initial_projection_eror
      @warn "Projection error after fitting is larger than before fitting. Using the initial guess instead."
      @info "Initial projection error: $initial_projection_eror"
      @info "Fitted projection error: $fitted_projection_error"
      x = initial_weight_matrix[:, period]
    end
    # Drop weights below the resolution. The projections keep `x` non-negative, so
    # `x .< tol` is equivalent to `abs.(x) .< tol` (no `abs` needed) and captures
    # exactly the contributions too small to matter (e.g. < 1% for tol = 1e-2).
    x[x.<tol] .= 0.0
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
  if diagnostics ≢ nothing
    diagnostics[:pgd_iters] = pgd_iters
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
  - `tol`: the resolution to which weights are determined; the single tolerance
    used throughout fitting (see the primary method).
"""
function fit_rep_period_weights!(
  clustering_result::ClusteringResult;
  weight_type::Symbol=:dirac,
  tol::Float64=1e-2,
)
  fit_rep_period_weights!(
    clustering_result.weight_matrix,
    clustering_result.clustering_matrix,
    clustering_result.rp_matrix;
    weight_type,
    tol,
    diagnostics=clustering_result.diagnostics,
  )
end

"""
  fit_chain_weights(increment_matrix, selected_indices) -> Matrix{Float64}

Fit the signed "chain" weights `W^ch` used by the inter-period storage chain
(the prolongation role of the weights). Unlike the operational weights, which
solve the constrained projection (4) per period with projected gradient descent,
`W^ch` solves the *same* least-squares projection on the storage-increment data
`increment_matrix` over the **unconstrained (signed)** domain, for which the
minimiser is the closed-form pseudoinverse — no PGD.

`increment_matrix` `G` has one row per seasonal-storage asset and one column per
base period; column `g_d` is period `d`'s net storage-increment proxy (per-period
net inflow energy). `selected_indices` are the base-period indices chosen as
representatives (`G_R = G[:, selected_indices]`), in representative-period order.

Each asset row of `G` is first centred to net zero over the year. That places `G`
on the cyclically closed manifold `1ᵀG = 0`, on which the unconstrained optimum
automatically satisfies operator-level closure — the column sums of `W^ch` vanish,
`1ᵀW^ch = 0` — so the reconstructed inter-period trajectory returns to its start
(constraint (2k)) for any intra-period solution, with no constrained solve needed.

Returns `W^ch` as a dense `n_base_periods × n_rp` matrix (the same shape and RP
ordering as the operational `weight_matrix`), so it drops straight into the chain
constraints. The fit is exact for period `d` when `g_d ∈ col(G_R)`; otherwise the
pseudoinverse returns the minimum-norm least-squares blend.
"""
function fit_chain_weights(
  increment_matrix::AbstractMatrix{Float64},
  selected_indices::AbstractVector{<:Integer},
)
  # Centre each asset row to net-zero over the year (cyclic closure): 1ᵀG = 0.
  n_periods = size(increment_matrix, 2)
  G = increment_matrix .- sum(increment_matrix; dims=2) ./ n_periods
  G_R = G[:, selected_indices]
  # Signed, unconstrained least squares per base period: w^ch_d = G_R⁺ g_d, i.e.
  # W^ch = (G_R⁺ G)ᵀ. `pinv` gives the minimum-norm solution when G_R is rank
  # deficient, so the fit degrades gracefully rather than failing.
  Wch = Matrix{Float64}((pinv(G_R) * G)')
  # Enforce *exact* operator closure 1ᵀW^ch = 0 (zero column sums over base periods).
  # De-meaning G only drives the column sums to O(ε); the storage chain pins both
  # σ^inter_0 = S^0 and the cyclic σ^inter_D = σ^inter_0, and the increment dynamics
  # give σ^inter_D = σ^inter_0 + (column sums)·y, so a residual O(ε) column sum would
  # ε-conflict the two endpoints into spurious infeasibility. Projecting each column to
  # mean-zero over d removes it exactly, an O(ε) change to the fit itself.
  Wch .-= sum(Wch; dims=1) ./ size(Wch, 1)
  return Wch
end
