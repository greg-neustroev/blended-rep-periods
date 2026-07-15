
"""
  project_onto_simplex(vector)

Projects `vector` onto a unit simplex using Michelot's algorithm in
Condat's accelerated implementation (2017). See Figure 2 of
[Condat, L. _Fast projection onto the simplex and the  ball._ Math. Program. 158,
575–585 (2016).](https://doi.org/10.1007/s10107-015-0946-6). For the details on
the meanings of v, ṽ, ρ and other variables, see the original paper.

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
  ṽ = Vector{Float64}()
  ρ = vector[1] - 1.0
  # step 2
  for y ∈ vector[2:end]
    if y > ρ
      ρ += (y - ρ) / (length(v) + 1)
      if ρ > y - 1.0
        push!(v, y)
      else
        append!(ṽ, v)
        v = [y]
        ρ = y - 1.0
      end
    end
  end
  # step 3
  for y ∈ ṽ
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

Fits `x` by accelerated projected gradient descent (FISTA with adaptive restart).

Minimizes a smooth convex objective over the feasible set defined by `projection`,
starting from `x`. The iteration uses Nesterov extrapolation, which reaches a target
accuracy in `O(√κ·log(1/ε))` iterations — with `κ = σ_max²/σ_min²` the condition
number of the objective's Hessian — versus `O(κ·log(1/ε))` for un-accelerated PGD.
That is a quadratic-to-linear improvement in the conditioning, which is exactly the
factor that blows up as the number of representative periods (and hence the near
collinearity of the RP columns) grows. O'Donoghue–Candès gradient restart recovers
the linear rate on the locally active face without needing the strong-convexity
constant `μ`, so the method stays parameter-free: the callers still pass only `α=1/L`.

`tol` is the resolution to which `x` is determined: the loop stops once the
projected-gradient mapping residual satisfies `‖G(x)‖₂ ≤ tol` at a *feasible* iterate
(a genuine first-order optimality certificate, checked as
`‖x − proj(x − α·g(x))‖₂ ≤ tol·α`). A cheap consecutive-move test gates that check so
the one extra gradient it costs is only paid near convergence. Certifying at a
feasible point — rather than at the extrapolated point — keeps the returned
projection accurate enough for the greedy-hull cache certificate (Lemma 1) to stay
sound. No iteration count is imposed; the loop is bounded only by the numerical
backstop `PGD_MAX_ITERS`.

The arguments:

  - `x`: the value to fit
  - `gradient`: the gradient operator, a function taking vectors of the same shape as
    `x` and returning the gradient of the (implicit) loss at that point
  - `projection`: the Euclidean projection onto the feasible set
  - `tol`: the resolution / first-order optimality tolerance (see above)
  - `learning_rate`: the step size `α` (the callers set it to `1/L` per instance)
  - `backend`: `:fista` (default, the accelerated method described above) or `:plain`
    — un-accelerated projected gradient descent (no momentum/extrapolation, no
    restart), kept ONLY to let the acceleration ablation study isolate FISTA's
    contribution from the Gram-matrix reformulation's (see `fit_rep_period_weights!`'s
    `pgd_variant` argument). `:plain` uses the identical certified stopping test, so
    both backends converge to the same fixed point — only the iteration count and
    wall-clock differ.

Returns the tuple `(x, iterations)`: the fitted point and the number of iterations
the loop actually ran (so callers can record the observed `N(ε)` without imposing it).
"""
function projected_gradient_descent!(
  x::AbstractVector{Float64};
  gradient::Function,
  projection::Function,
  tol::Float64=1e-2,
  learning_rate::Float64=1e-3,
  backend::Symbol=:fista,
)
  backend ∈ (:fista, :plain) || throw(ArgumentError("Unsupported PGD backend $(backend); expected :fista or :plain."))
  α = learning_rate
  # The initial guess may lie outside the feasible set; project it first. `x` always
  # holds the latest feasible iterate; `y` is the (possibly infeasible) extrapolation
  # point at which the gradient is evaluated (for backend=:plain, y is always x, i.e.
  # plain un-accelerated PGD with no momentum).
  x = projection(x)
  y = copy(x)
  t = 1.0
  # Certified stop: at a feasible iterate the projected-gradient mapping residual is
  # ‖x − proj(x − α·g(x))‖₂ = α·‖G(x)‖₂, so ‖…‖₂ ≤ tol·α certifies ‖G(x)‖₂ ≤ tol.
  threshold = tol * α
  iterations = 0
  for _ ∈ 1:PGD_MAX_ITERS
    iterations += 1
    g = gradient(y)
    x_new = projection(y .- α .* g)          # (accelerated, if backend=:fista) gradient + projection step
    # Adaptive restart (O'Donoghue–Candès): if the momentum (x_new − x) points against
    # the projected descent step (y − x_new), the extrapolation is overshooting. Drop
    # the momentum and retake the step from the feasible iterate x. Not applicable to
    # the :plain backend, which never builds momentum in the first place.
    if backend ≡ :fista && dot(y .- x_new, x_new .- x) > 0
      y = x
      g = gradient(y)
      x_new = projection(y .- α .* g)
      t = 1.0
    end
    # First-order optimality certificate at the feasible iterate x_new, gated by the
    # cheap consecutive-move test so the extra gradient is only paid near convergence.
    converged = false
    if norm(x_new .- x) ≤ threshold
      g_cert = gradient(x_new)
      converged = norm(projection(x_new .- α .* g_cert) .- x_new) ≤ threshold
    end
    if backend ≡ :fista
      # Nesterov extrapolation for the next iterate.
      t_next = (1.0 + sqrt(1.0 + 4.0 * t^2)) / 2.0
      y = x_new .+ ((t - 1.0) / t_next) .* (x_new .- x)
      t = t_next
    else
      y = x_new   # :plain — next gradient is evaluated at the feasible iterate itself
    end
    x = x_new
    converged && break
  end
  # Return the realized iteration count alongside the fitted point so callers can
  # record the observed N(ε) without imposing it.
  return x, iterations
end

"""
  fit_rep_period_weights!(weight_matrix, clustering_matrix, rp_matrix; weight_type, tol, args...)

Given the initial weight guesses, finds better weights for convex or conical
combinations of representative periods.

The arguments:

  - `weight_matrix`: the initial guess for weights; the weights are adjusted
    using a projected gradient descent method
  - `clustering_matrix`: the matrix of raw clustering data
  - `rp_matrix`: the matrix of raw representative period data
  - `weight_type`: the type of weights to find; possible values are:
      - `:dirac`: each period is represented by its single nearest representative
        (no blending; the seeded nearest-representative weights are kept as-is)
      - `:convex`: each period is represented as a convex sum of the
        representative periods (a sum with nonnegative weights adding into one)
      - `:conical`: each period is represented as a conical sum of the
        representative periods (a sum with nonnegative weights)
      - `:conical_bounded`: a conical sum (nonnegative weights) whose total weight is
        bounded above by one (the sub-unit conic class)
  - `tol`: the resolution to which weights are determined; the single tolerance
    used throughout fitting. The projected gradient descent stops once no weight
    moves by more than `tol`, base periods already reconstructed within `tol` skip
    fitting, and fitted weights below `tol` are dropped (e.g. with `tol = 1e-2`
    weights are blending percentages and contributions under 1% are noise).
  - `pgd_variant`: which of the three acceleration-ablation arms to run (exists so the
    contribution of the Gram-matrix reformulation and of FISTA-with-restart can be
    isolated experimentally; NOT a knob for production use — always leave this at its
    default outside the ablation study):
      - `:proposed` (default): Gram-matrix gradient (precomputed `G=RᵀR`, `B=RᵀC`) +
        FISTA with adaptive restart — the method described above and used everywhere
        else in the package.
      - `:gram_only`: the same Gram-matrix gradient, but plain (un-accelerated) PGD —
        isolates FISTA's contribution alone (same gradient cost per iteration as
        `:proposed`, more iterations to converge).
      - `:plain`: the "tall-R" gradient `Rᵀ(Rw−c_d)`, recomputed every iteration
        without the Gram-matrix precomputation, and plain PGD — the pre-acceleration
        baseline both changes were compared against.
    All three converge to the identical fixed point (same objective, same certified
    stopping test) — only iteration count and wall-clock time differ.
"""
function fit_rep_period_weights!(
  weight_matrix::Union{SparseMatrixCSC{Float64,Int64},Matrix{Float64}},
  clustering_matrix::Matrix{Float64},
  rp_matrix::Matrix{Float64};
  weight_type::Symbol=:dirac,
  tol::Float64=1e-2,
  diagnostics::Union{Nothing,Dict{Symbol,Any}}=nothing,
  pgd_variant::Symbol=:proposed,
)
  pgd_variant ∈ (:proposed, :gram_only, :plain) ||
    throw(ArgumentError("Unsupported pgd_variant $(pgd_variant); expected :proposed, :gram_only, or :plain."))
  use_gram = pgd_variant ≡ :proposed || pgd_variant ≡ :gram_only
  pgd_backend = pgd_variant ≡ :proposed ? :fista : :plain
  # Determine the appropriate projection method
  if weight_type ≡ :dirac
    return weight_matrix
  elseif weight_type ≡ :convex
    projection = project_onto_simplex
  elseif weight_type ≡ :conical
    projection = project_onto_nonnegative_orthant
  elseif weight_type ≡ :conical_bounded
    # Sub-unit conic: convex fitting with an extra zero component appended. The
    # weight of that zero vector is then discarded without affecting the total, so
    # the remaining weights are nonnegative with a sum between zero and one.
    projection = project_onto_simplex
    n_data_points = size(rp_matrix, 1)
    rp_matrix = hcat(repeat([0.0], n_data_points), rp_matrix)
  else
    throw(ArgumentError("Unsupported weight type $(weight_type); expected :dirac, :convex, :conical, or :conical_bounded."))
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
    # A zero row matches the zero column appended to rp_matrix above.
    initial_weight_matrix = vcat(zeros(1, size(initial_weight_matrix, 2)), initial_weight_matrix)
  end
  for col in eachcol(initial_weight_matrix)
    col .= projection(col)
  end

  # Precompute the Gram matrix G = RᵀR and the projected targets B = RᵀC once and
  # reuse them for every base period. The gradient R ᵀ(R w − c_d) = G w − B[:,d] is
  # then a single n_rp×n_rp matvec (O(n_rp²)) instead of two matmuls against the tall
  # R (O(|features|·n_rp) each); since |features| ≫ n_rp this is the dominant
  # per-iteration saving, and the same G/B are shared across all base periods.
  # Principled PGD step size: eigmax(G) = σ_max²(R) = L is the Lipschitz constant of
  # the objective's gradient, so α = 1/L exactly as before (now read off the small
  # n_rp×n_rp Gram matrix rather than an SVD of the tall R). `L` itself is needed for
  # the step size regardless of `pgd_variant` (even the :plain ablation arm uses the
  # same principled α = 1/L — only the per-iteration gradient cost differs).
  gram = Symmetric(rp_matrix' * rp_matrix)
  step_size = 1 / eigmax(gram)
  # `projected_targets`/`grad_buf` are only needed by the Gram-based gradient path.
  projected_targets = use_gram ? rp_matrix' * clustering_matrix : nothing
  grad_buf = use_gram ? Vector{Float64}(undef, size(gram, 1)) : nothing

  # Record per-fit PGD iteration counts (one entry per base period that is
  # actually fitted, i.e. not already within `tol`), so the observed N(ε) range
  # can be reported without re-running.
  pgd_iters = Int[]

  for period ∈ 1:n_periods
    target_vector = clustering_matrix[:, period]
    x = projection(initial_weight_matrix[:, period])
    gradient = if use_gram
      b = view(projected_targets, :, period)
      let gram = gram, b = b, grad_buf = grad_buf
        w -> (mul!(grad_buf, gram, w); grad_buf .-= b; grad_buf)
      end
    else
      # Ablation-only baseline: recompute the "tall-R" gradient Rᵀ(Rw−c_d) every call,
      # without the Gram-matrix precomputation — the pre-7c80509 cost profile.
      let rp_matrix = rp_matrix, target_vector = target_vector
        w -> rp_matrix' * (rp_matrix * w - target_vector)
      end
    end
    initial_projection_eror = norm(rp_matrix * x - target_vector)
    if initial_projection_eror ≤ tol
      continue  # already reconstructed within the resolution; keep the initial weights
    end
    x, n_iter = projected_gradient_descent!(x; gradient, projection, tol, learning_rate=step_size, backend=pgd_backend)
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
      # Because some values might have been removed, convexity can be lost. In the
      # upper-bounded case the sum can also drift slightly above one from rounding.
      # To account for these cases, the weights are re-normalized.
      sum_x = sum(x)
      if weight_type ≡ :convex || sum_x > 1.0
        x = x ./ sum_x
      end
    end
    if weight_type ≡ :conical_bounded
      # Drop the null component; the remaining weights sum to at most one.
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
  fit_rep_period_weights!(clustering_result; weight_type, tol)

Given the initial weight guesses, finds better weights for convex or conical
combinations of representative periods.

The arguments:

  - `clustering_result`: the result of running
    `find_representative_periods`
  - `weight_type`: the type of weights to find; possible values are `:dirac`,
    `:convex`, `:conical`, and `:conical_bounded` (see the primary method).
  - `tol`: the resolution to which weights are determined; the single tolerance
    used throughout fitting (see the primary method).
  - `pgd_variant`: acceleration-ablation arm; see the primary method (leave at its
    `:proposed` default outside the ablation study).
"""
function fit_rep_period_weights!(
  clustering_result::ClusteringResult;
  weight_type::Symbol=:dirac,
  tol::Float64=1e-2,
  pgd_variant::Symbol=:proposed,
)
  fit_rep_period_weights!(
    clustering_result.weight_matrix,
    clustering_result.clustering_matrix,
    clustering_result.rp_matrix;
    weight_type,
    tol,
    diagnostics=clustering_result.diagnostics,
    pgd_variant,
  )
end

