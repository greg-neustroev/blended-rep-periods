
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
  project_box_sum(v; lo=0.0, hi=1.0, s=1.0)

Euclidean projection of `v` onto `{ w : sum(w)=s, lo ≤ w_i ≤ hi }` — the capped-simplex
projection generalized to an arbitrary box. One projector serves all three sum-constrained
chain-weight classes, differing only in `lo`/`hi`:

  - `lo=0,  hi=1`  → the probability simplex = **convex** weights (the `hi=1` clip is
    vacuous since `sum=1, w≥0 ⇒ w≤1`);
  - `lo=-1, hi=1`  → sum-1 with `|w|≤1` = **clipped affine** weights (convex with bounded
    negative weights allowed);
  - `lo=-Inf` or `hi=Inf` → the closed-form hyperplane projection = **affine** weights.

The KKT optimum has the form `w_i = clip(v_i + θ, lo, hi)` for a single scalar `θ` fixed by
`sum(w)=s`. This is NOT project-then-clip (projecting onto the hyperplane then clamping
breaks the sum); the correct object finds the `θ` such that clipping *already* sums to `s`.
`θ` is found exactly by walking the `2n` breakpoints of the piecewise-linear, nondecreasing
`g(θ)=Σ clip(v_i+θ,lo,hi)`.
"""
function project_box_sum(v::AbstractVector{T}; lo::Real=0.0, hi::Real=1.0, s::Real=1.0) where {T<:Real}
  n = length(v)
  lo = T(lo); hi = T(hi); s = T(s)
  lo < hi || throw(ArgumentError("need lo < hi"))
  # Unbounded (affine) case: closed-form hyperplane projection.
  if isinf(lo) || isinf(hi)
    return v .+ (s - sum(v)) / n
  end
  g(θ) = (acc = zero(T); @inbounds for x in v; acc += clamp(x + θ, lo, hi); end; acc)
  bps = Vector{T}(undef, 2n)
  @inbounds for i in 1:n
    bps[2i-1] = lo - v[i]
    bps[2i] = hi - v[i]
  end
  sort!(bps)
  a = bps[1] - one(T); ga = n * lo  # left of all breakpoints everything clamps to lo
  @inbounds for k in 1:2n
    b = bps[k]; gb = g(b)
    if gb >= s
      # root in [a,b]; g linear here with slope = #interior coords
      mid = (a + b) / 2; slope = zero(T)
      @inbounds for x in v
        y = x + mid
        slope += (lo < y < hi) ? one(T) : zero(T)
      end
      θ = slope > 0 ? a + (s - ga) / slope : b
      return clamp.(v .+ θ, lo, hi)
    end
    a, ga = b, gb
  end
  return clamp.(v .+ (bps[end] + one(T)), lo, hi)  # right of all breakpoints: clamp to hi
end

"""
  project_subunit_conic(v)

Euclidean projection of `v` onto `{ w : w ≥ 0, sum(w) ≤ 1 }` — the sub-unit conic class.
Like convex it is a contraction in the ℓ∞-induced norm (`‖w‖₁ = sum(w) ≤ 1`, so no expansion
blowup), but unlike convex it does NOT force the sum to 1, so it is *balance-relaxed*. If the
nonnegative part already sums to ≤ 1 the sum constraint is inactive (`w = max(v,0)`); otherwise
it binds and the projection lands on the simplex.
"""
function project_subunit_conic(v::AbstractVector{Float64})
  w = max.(v, 0.0)
  return sum(w) <= 1.0 ? w : project_box_sum(v; lo=0.0, hi=1.0, s=1.0)
end

"""
  project_l1_ball(v)

Euclidean projection of `v` onto the ℓ1 ball `{ w : Σ_r|w_r| ≤ 1 }` — the *signed*
contraction class. Like sub-unit conic it is a contraction in the ℓ∞-induced norm
(`‖W‖_{∞→∞} = max row-ℓ1 ≤ 1`, so no expansion blow-up), but signs are free, so it has the
reconstruction reach convex/sub-unit-conic lack while staying contractive — the
norm-correct version of clipped affine (which bounded ℓ∞ and let ℓ1 run to R). Standard
soft-threshold projection: if `‖v‖₁ ≤ 1` keep `v`; else project `|v|` onto the simplex and
restore signs (`w = sign(v) ⊙ max(|v|-θ, 0)`).
"""
function project_l1_ball(v::AbstractVector{Float64})
  sum(abs, v) <= 1.0 && return copy(v)
  return sign.(v) .* project_box_sum(abs.(v); lo=0.0, hi=1.0, s=1.0)
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
  fit_chain_weights(increment_matrix, selected_indices; weight_type, tol) -> (W^ch, residual)

Fit the "chain" weights `W^ch` used by the inter-period storage chain (the prolongation
role of the weights), as a *separate* fit from the operational weights, on the
storage-increment data `increment_matrix`. `increment_matrix` `G` has one row per
seasonal-storage asset and one column per base period; column `g_d` is period `d`'s net
storage-increment proxy (per-period net inflow energy). `selected_indices` are the
base-period indices chosen as representatives (`G_R = G[:, selected_indices]`), in RP order.

Each asset row of `G` is first centred to net zero over the year (`1ᵀG = 0`). `weight_type`
selects the class:

  - `:signed` — unconstrained least squares, closed-form pseudoinverse `W^ch=(G_R⁺G)ᵀ`,
    then projected to *exact* column-sum-zero so the closure `1ᵀW^ch=0` holds to machine
    precision (the reconstructed trajectory returns to its start for any dispatch). Exact
    when `g_d ∈ col(G_R)`; min-norm least-squares otherwise. Its weights are unbounded, so
    it blows up when `G_R` is rank-deficient.
  - `:convex` / `:clipped_affine` / `:affine` / `:conical` — per-period projected gradient
    descent on the *same* target, differing ONLY in the projection (`project_box_sum` with
    `lo=0` / `lo=-1` / `lo=-Inf`, or the nonneg orthant). The sum-1 classes have positive
    column sums (no closure gauge); annual balance is earned by dispatch via the cyclic
    constraint, as the operational weights do. The `lo`-only difference between convex and
    clipped-affine makes any downstream regret difference attributable to the lower bound.

Returns `(W^ch, residual)`: `W^ch` is a dense `n_base_periods × n_rp` matrix (same shape and
RP ordering as the operational `weight_matrix`); `residual = ‖G_R W^chᵀ − G‖_F/‖G‖_F` is the
increment-reconstruction error (~0 for full-rank `:signed`; the hull-containment error for
the bounded classes — small in-hull, large in the sign-changing seasonal regime).
"""
function fit_chain_weights(
  increment_matrix::AbstractMatrix{Float64},
  selected_indices::AbstractVector{<:Integer};
  weight_type::Symbol=:signed,
  tol::Float64=1e-2,
)
  # Centre each asset row to net-zero over the year (cyclic closure): 1ᵀG = 0. All
  # classes fit the *same* centred-increment target, so their reconstruction residuals
  # (returned alongside) are directly comparable as a hull-containment diagnostic.
  n_periods = size(increment_matrix, 2)
  G = increment_matrix .- sum(increment_matrix; dims=2) ./ n_periods
  G_R = G[:, selected_indices]
  R = length(selected_indices)
  if weight_type ≡ :signed
    # Signed, unconstrained least squares per base period: w^ch_d = G_R⁺ g_d, i.e.
    # W^ch = (G_R⁺ G)ᵀ. `pinv` gives the minimum-norm solution when G_R is rank
    # deficient, so the fit degrades gracefully rather than failing. Then enforce
    # *exact* operator closure 1ᵀW^ch = 0: de-meaning G only drives the column sums to
    # O(ε), but the chain pins both σ^inter_0 = S^0 and the cyclic σ^inter_D = σ^inter_0
    # while the dynamics give σ^inter_D = σ^inter_0 + (column sums)·y, so a residual
    # O(ε) column sum would ε-conflict the endpoints into spurious infeasibility.
    # Projecting each column to mean-zero over d removes it exactly.
    Wch = Matrix{Float64}((pinv(G_R) * G)')
    Wch .-= sum(Wch; dims=1) ./ size(Wch, 1)
  else
    # Sum-1 / nonneg classes via per-base-period projected gradient descent, all sharing
    # the SAME PGD and target and differing ONLY in the projection — so any difference
    # between, e.g., convex and clipped_affine is attributable to the projection's lower
    # bound (`lo`: 0 vs -1) alone. None of these carry the column-sum-zero closure gauge
    # (the sum-1 classes have positive column sums; conical nonneg) — annual balance is
    # earned by dispatch through the cyclic constraint, exactly as the operational weights.
    projection =
      weight_type ≡ :convex         ? (w -> project_box_sum(w; lo=0.0, hi=1.0, s=1.0)) :
      weight_type ≡ :clipped_affine ? (w -> project_box_sum(w; lo=-1.0, hi=1.0, s=1.0)) :
      weight_type ≡ :affine         ? (w -> w .+ (1.0 - sum(w)) / length(w)) :
      weight_type ≡ :conical        ? project_onto_nonnegative_orthant :
      weight_type ≡ :conical_bounded ? project_subunit_conic :
      weight_type ≡ :l1_ball        ? project_l1_ball :
      throw(ArgumentError("Unsupported chain weight_type $(weight_type); expected :signed, " *
                          ":convex, :clipped_affine, :affine, :conical, :conical_bounded, or :l1_ball"))
    Wch = Matrix{Float64}(undef, n_periods, R)
    step_size = 1.0 / opnorm(G_R, 2)^2
    for d in 1:n_periods
      target = G[:, d]
      gradient = w -> G_R' * (G_R * w - target)
      x = projection(zeros(R))
      x, _ = projected_gradient_descent!(x; gradient, projection, tol, learning_rate=step_size)
      Wch[d, :] = x
    end
  end
  # Relative reconstruction residual ‖G_R W^chᵀ − G‖_F / ‖G‖_F: ~0 for :signed (exact
  # at full row rank), and the hull-containment error for the bounded classes (small
  # when the base increments are in-hull, large in the sign-changing seasonal regime).
  residual = norm(G_R * Wch' - G) / max(norm(G), eps())
  return Wch, residual
end
