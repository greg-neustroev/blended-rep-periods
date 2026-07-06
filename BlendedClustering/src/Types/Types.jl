"""
Core data types shared across the package: experiment configuration, clustering
results, and the per-experiment result record.
"""
module Types

using CSV
using DataFrames
using JuMP
using LinearAlgebra
using SparseArrays
import Tables

export ExperimentData, ExperimentResult, read_run_data
export ClusteringResult, AuxiliaryClusteringData

# Default weight resolution tolerance used when a configuration CSV omits the
# optional `tol` column: weights are fit and reported only down to this precision.
const DEFAULT_PGD_TOL = 1e-2

# Default clustering-feature normalization used when a configuration CSV omits the
# optional `normalization` column. Supported values:
#   - `:unscaled` (default): the historical behavior — cluster directly on the
#     dimensionless profiles as stored, where each value is a fraction of the asset's
#     peak (0 means physically zero, 1 means the peak), with no scaling of its own.
#   - `:minmax`: rescale each feature row to [0,1] (its minimum across periods → 0,
#     its maximum → 1). The classic per-feature min-max baseline; it centers (removes
#     the absolute level), so it is applied only to the clustering features.
#   - `:economic`: rescale the features into common physical/economic units (peak ×
#     profile, per-block economic weights), never centered, so zero keeps meaning zero.
# `:minmax` and `:economic` transform only the selection/weight-fitting space; the
# representative-period profiles handed to the model stay in the original units.
const DEFAULT_NORMALIZATION = :unscaled

# Default inflow-integral weight λ used when a configuration CSV omits the optional
# `inflow_integral_weight` column. `0.0` means the integrated-inflow rows are not
# added, so the clustering matrix stays `[D; A; E]` and every existing result is
# unchanged; a non-zero λ augments it to `[D; A; E; λ·E_int]`, stacking each inflow
# asset's per-period total inflow energy (see `find_representative_periods`).
const DEFAULT_INFLOW_INTEGRAL_WEIGHT = 0.0

# Default storage-regret fixing cadence used when a configuration CSV omits the
# optional `fix_every` column. `1` pins every base-period boundary state of charge
# in the full-horizon evaluation model (the original, strictest grading); a coarser
# `k > 1` pins only every k-th boundary so the full model re-optimizes the storage
# trajectory between checkpoints (see `evaluate_solution!`). Ignored for any
# `evaluation_type` other than `:storage_regret`.
const DEFAULT_FIX_EVERY = 1

# Default for the optional `inject_inflow` column. `true` (the standard model) switches on
# exact-inflow injection: the inter-period storage chain blends only the dispatch response
# (Δσ_r − Ē_r), adds each base period's *real* inflow E_d exactly, closes the per-base-day
# ledger with conservation slacks (spill/borrow at the model's own prices), and enforces the
# seasonal capacity band as a hard [min,max]·cap constraint (defaulting to [0,1] where a
# dataset lacks reservoir-level profiles — a data gap the model assumes filled). This keeps
# every day-resolution term of the storage conservation law at day resolution; only the
# compressed dispatch response stays blended. `false` recovers the historical chain (blend
# the whole increment Δσ_r, inflow included, with the soft VOLL band) — retained as the
# `-injection` ablation arm. See the `storage_inter` / `soc_cap_inter` blocks in
# `create_optimization_model!`.
const DEFAULT_INJECT_INFLOW = true


"""
    read_run_data(path) -> DataFrame

Read the experiment configuration at `path` (a CSV with one experiment per row)
into a `DataFrame`, checking that the required columns are present and converting
the categorical columns (`clustering_type`, `weight_type`, `evaluation_type`, and
the optional `normalization`) to `Symbol`s. The projected gradient descent
tolerance column `tol` is optional and defaults to `DEFAULT_PGD_TOL` when absent;
the `normalization` column is optional and defaults to `DEFAULT_NORMALIZATION`; the
`inflow_integral_weight` column is optional and defaults to
`DEFAULT_INFLOW_INTEGRAL_WEIGHT`; the `fix_every` column is optional and defaults to
`DEFAULT_FIX_EVERY`.
"""
function read_run_data(path)
    df = CSV.read(path, DataFrame)
    required_columns = Set([
        "n_rep_periods",
        "period_length",
        "clustering_type",
        "weight_type",
        "evaluation_type",
    ])
    if !issubset(required_columns, Set(names(df)))
        error("Input CSV file is missing required columns: $required_columns")
    end
    symbol_columns = [:clustering_type, :weight_type, :evaluation_type]
    # `normalization` is optional; convert it to a Symbol only when present.
    if "normalization" in names(df)
        push!(symbol_columns, :normalization)
    end
    transform!(df, symbol_columns .=> ByRow(Symbol) .=> symbol_columns)
    return df
end


"""
Data needed to run a single experiment (i.e., a single optimization model)
"""
struct ExperimentData
    name::String

    # Clustering
    n_rep_periods::Int
    period_length::Int
    clustering_type::Symbol
    weight_type::Symbol
    tol::Float64
    normalization::Symbol
    cache::Bool
    inflow_integral_weight::Float64
    fix_every::Int
    inject_inflow::Bool
    evaluation_type::Symbol

    function ExperimentData(run_data_row::DataFrameRow{DataFrame,DataFrames.Index}, base_name::String)
        # `tol` (the PGD tolerance ε) is optional; fall back to the default when
        # the configuration CSV does not provide the column.
        tol = hasproperty(run_data_row, :tol) ? Float64(run_data_row.tol) : DEFAULT_PGD_TOL
        # `normalization` is optional and defaults to the historical `:minmax`.
        normalization = hasproperty(run_data_row, :normalization) ?
                        Symbol(run_data_row.normalization) : DEFAULT_NORMALIZATION
        # `cache` toggles the greedy-hull projection cache; optional and on by
        # default. Setting it off is only for the cached-vs-uncached benchmark.
        cache = hasproperty(run_data_row, :cache) ? Bool(run_data_row.cache) : true
        # `inflow_integral_weight` (the weight λ on the integrated-inflow rows) is
        # optional and defaults to off; only a non-zero value augments the clustering
        # matrix with the per-period total-inflow rows.
        inflow_integral_weight = hasproperty(run_data_row, :inflow_integral_weight) ?
                                 Float64(run_data_row.inflow_integral_weight) :
                                 DEFAULT_INFLOW_INTEGRAL_WEIGHT
        # `fix_every` (the storage-regret fixing cadence k) is optional and defaults
        # to pinning every boundary; only a coarser cadence relaxes the grading.
        fix_every = hasproperty(run_data_row, :fix_every) ?
                    Int(run_data_row.fix_every) : DEFAULT_FIX_EVERY
        # `inject_inflow` toggles exact-inflow injection in the inter-period storage chain;
        # optional and on by default (the standard model). Set it false for the ablation arm.
        inject_inflow = hasproperty(run_data_row, :inject_inflow) ?
                        Bool(run_data_row.inject_inflow) : DEFAULT_INJECT_INFLOW
        # Keep the experiment name (and thus output paths) unchanged for the
        # default normalization/cache; only the non-default arms get a suffix, so the
        # same RP grid can be run several ways without colliding.
        name_parts = [
            base_name,
            run_data_row.n_rep_periods,
            run_data_row.period_length,
            string(run_data_row.clustering_type),
            string(run_data_row.weight_type),
            tol,
        ]
        if normalization ≠ DEFAULT_NORMALIZATION
            push!(name_parts, string(normalization))
        end
        if inflow_integral_weight ≠ DEFAULT_INFLOW_INTEGRAL_WEIGHT
            push!(name_parts, "inflowint$(inflow_integral_weight)")
        end
        if fix_every ≠ DEFAULT_FIX_EVERY
            push!(name_parts, "fixevery$(fix_every)")
        end
        if inject_inflow ≠ DEFAULT_INJECT_INFLOW
            push!(name_parts, inject_inflow ? "inject" : "noinject")
        end
        if !cache
            push!(name_parts, "nocache")
        end
        name = join(name_parts, "_")
        return new(
            name,
            run_data_row.n_rep_periods,
            run_data_row.period_length,
            run_data_row.clustering_type,
            run_data_row.weight_type,
            tol,
            normalization,
            cache,
            inflow_integral_weight,
            fix_every,
            inject_inflow,
            run_data_row.evaluation_type
        )
    end
end

"""
Structure to hold the time series used in clustering together with some
summary statistics on the data.
"""
mutable struct AuxiliaryClusteringData
    key_columns::AbstractVector{Symbol}
    period_duration::Int
    last_period_duration::Int
    n_periods::Int
end

"""
Structure to hold the clustering result.
"""
mutable struct ClusteringResult
    profiles::AbstractDataFrame
    weight_matrix::Union{SparseMatrixCSC{Float64,Int64},Matrix{Float64}}
    clustering_matrix::Union{Matrix{Float64},Nothing}
    rp_matrix::Union{Matrix{Float64},Nothing}
    # Free-form diagnostics captured while selecting representatives and fitting
    # weights (greedy-hull cache hit/miss counts, per-fit PGD iteration counts).
    # Populated in place so the experiment layer can record them without re-running
    # the clustering.
    diagnostics::Dict{Symbol,Any}
end

# Convenience constructors: the clustering/representative-period matrices are not
# always available (e.g. the single-period fast path), so default them to nothing;
# diagnostics default to an empty dict.
ClusteringResult(profiles, weight_matrix, clustering_matrix, rp_matrix) =
    ClusteringResult(profiles, weight_matrix, clustering_matrix, rp_matrix, Dict{Symbol,Any}())
ClusteringResult(profiles, weight_matrix) =
    ClusteringResult(profiles, weight_matrix, nothing, nothing, Dict{Symbol,Any}())

# --- Helpers for populating the per-run summary record (capture-once stats) ---

# Value of a named model expression/variable after solving, or `missing` when the
# model is absent, unsolved/infeasible, or does not register that name.
_solved_value(m, sym::Symbol) =
    (m !== nothing && JuMP.is_solved_and_feasible(m) && haskey(m, sym)) ? JuMP.value(m[sym]) : missing

# Singular-value spectrum summary of the representative-period matrix R: the
# extreme singular values, the 2-norm condition number, the stable rank
# (‖R‖_F²/σ_max²) and the numerical rank at a relative tolerance. All `missing`
# when R is unavailable. The full spectrum itself is dumped as an artifact, so
# any other spectral quantity stays derivable offline.
function _spectrum_stats(R)
    (R === nothing || isempty(R)) && return (missing, missing, missing, missing, missing)
    s = svdvals(Matrix{Float64}(R))
    isempty(s) && return (missing, missing, missing, missing, missing)
    smax = maximum(s)
    smin = minimum(s)
    kappa = smin > 0 ? smax / smin : Inf
    stable_rank = smax > 0 ? sum(abs2, s) / smax^2 : missing
    rtol = smax * eps(Float64) * maximum(size(R))
    numerical_rank = count(>(rtol), s)
    return (smax, smin, kappa, stable_rank, numerical_rank)
end

# Per-base-period weight-sum (ΣW) summary: mean and max blending mass, the
# fractions of base periods whose weights sum above/below one (conic over-shoot
# vs sub-unit/convex shortfall), and the weight-matrix density.
function _weight_stats(W)
    (W === nothing || isempty(W)) && return (missing, missing, missing, missing, missing)
    row_sums = vec(sum(W, dims=2))
    n = length(row_sums)
    mean_sum = sum(row_sums) / n
    max_sum = maximum(row_sums)
    frac_gt1 = count(>(1 + 1e-9), row_sums) / n
    frac_lt1 = count(<(1 - 1e-9), row_sums) / n
    density = count(!iszero, W) / length(W)
    return (mean_sum, max_sum, frac_gt1, frac_lt1, density)
end

# Total variables / constraints of a JuMP model (0 when the model is absent).
_n_variables(m) = m === nothing ? 0 : JuMP.num_variables(m)
_n_constraints(m) = m === nothing ? 0 : JuMP.num_constraints(m; count_variable_in_set_constraints=true)

struct ExperimentResult
    name::String
    seed::Int
    n_rep_periods::Int
    period_length::Int
    clustering_type::Symbol
    weight_type::Symbol
    normalization::Symbol
    tol::Float64
    fix_every::Int
    # Clustering geometry / quality.
    projection_error::Float64
    sigma_max::Union{Float64,Missing}
    sigma_min::Union{Float64,Missing}
    kappa::Union{Float64,Missing}
    stable_rank::Union{Float64,Missing}
    numerical_rank::Union{Int,Missing}
    # Blended-weight structure.
    mean_weight_sum::Union{Float64,Missing}
    max_weight_sum::Union{Float64,Missing}
    frac_weight_sum_gt1::Union{Float64,Missing}
    frac_weight_sum_lt1::Union{Float64,Missing}
    weight_density::Union{Float64,Missing}
    # PGD / greedy-hull cache diagnostics.
    pgd_total_iters::Int
    pgd_max_iters_per_fit::Int
    cache_hits::Int
    cache_misses::Int
    # Reduced-model solve, objective and its decomposition, problem size.
    termination_status::String
    objective_value::Union{Float64,Missing}
    cost_of_investment::Union{Float64,Missing}
    cost_of_operations::Union{Float64,Missing}
    cost_of_spillage::Union{Float64,Missing}
    cost_of_borrow::Union{Float64,Missing}
    cost_of_soc_band::Union{Float64,Missing}
    n_variables::Int
    n_constraints::Int
    # Evaluation (full-horizon) model: regret numerator, cost decomposition, size.
    evaluation_termination_status::String
    evaluated_objective_value::Union{Float64,Missing}
    eval_cost_of_investment::Union{Float64,Missing}
    eval_cost_of_operations::Union{Float64,Missing}
    eval_cost_of_spillage::Union{Float64,Missing}
    eval_cost_of_borrow::Union{Float64,Missing}
    eval_cost_of_soc_band::Union{Float64,Missing}
    eval_n_variables::Int
    eval_n_constraints::Int
    total_spillage::Float64
    total_borrow::Float64
    total_soc_band::Float64
    # Runtime breakdown.
    time_to_preprocess::Float64
    time_to_cluster::Float64
    time_to_fit_weights::Float64
    time_to_formulate_model::Float64
    time_to_solve::Float64

    function ExperimentResult(
        data::ExperimentData,
        seed::Int,
        solved_model::JuMP.AbstractModel,
        eval_model::Union{JuMP.AbstractModel,Nothing},
        projection_error::Float64,
        clustering_result::ClusteringResult,
        time_to_preprocess::Float64,
        time_to_cluster::Float64,
        time_to_fit_weights::Float64,
        time_to_formulate_model::Float64,
        time_to_solve::Float64,
    )
        name = data.name
        n_rep_periods = data.n_rep_periods
        period_length = data.period_length
        clustering_type = data.clustering_type
        weight_type = data.weight_type

        # Clustering geometry and blended-weight structure.
        sigma_max, sigma_min, kappa, stable_rank, numerical_rank =
            _spectrum_stats(clustering_result.rp_matrix)
        mean_weight_sum, max_weight_sum, frac_weight_sum_gt1, frac_weight_sum_lt1, weight_density =
            _weight_stats(clustering_result.weight_matrix)

        # PGD / greedy-hull cache diagnostics (populated in place during clustering
        # and weight fitting; absent for the single-period fast path or Dirac weights).
        diag = clustering_result.diagnostics
        pgd_iters = get(diag, :pgd_iters, Int[])
        pgd_total_iters = isempty(pgd_iters) ? 0 : sum(pgd_iters)
        pgd_max_iters_per_fit = isempty(pgd_iters) ? 0 : maximum(pgd_iters)
        cache_hits = get(diag, :cache_hits, 0)
        cache_misses = get(diag, :cache_misses, 0)

        termination_status = solved_model |> JuMP.termination_status |> string
        objective_value = if JuMP.is_solved_and_feasible(solved_model)
            solved_model |> JuMP.objective_value
        else
            missing
        end
        cost_of_investment = _solved_value(solved_model, :cost_of_investment)
        cost_of_operations = _solved_value(solved_model, :cost_of_operations)
        cost_of_spillage = _solved_value(solved_model, :cost_of_spillage)
        cost_of_borrow = _solved_value(solved_model, :cost_of_borrow)
        cost_of_soc_band = _solved_value(solved_model, :cost_of_soc_band)
        n_variables = _n_variables(solved_model)
        n_constraints = _n_constraints(solved_model)

        evaluation_termination_status = if !isnothing(eval_model)
            eval_model |> JuMP.termination_status |> string
        else
            "N/A"
        end
        evaluated_objective_value = if !isnothing(eval_model) && JuMP.is_solved_and_feasible(eval_model)
            eval_model |> JuMP.objective_value
        else
            missing
        end
        eval_cost_of_investment = _solved_value(eval_model, :cost_of_investment)
        eval_cost_of_operations = _solved_value(eval_model, :cost_of_operations)
        eval_cost_of_spillage = _solved_value(eval_model, :cost_of_spillage)
        eval_cost_of_borrow = _solved_value(eval_model, :cost_of_borrow)
        eval_cost_of_soc_band = _solved_value(eval_model, :cost_of_soc_band)
        eval_n_variables = _n_variables(eval_model)
        eval_n_constraints = _n_constraints(eval_model)
        total_spillage = if !isnothing(eval_model) && haskey(eval_model, :spillage) && !isempty(eval_model[:spillage])
            value.(eval_model[:spillage]) |> sum
        else
            0.0
        end
        total_borrow = if !isnothing(eval_model) && haskey(eval_model, :borrow) && !isempty(eval_model[:borrow])
            value.(eval_model[:borrow]) |> sum
        else
            0.0
        end
        # Total inter-period band violation (over- plus under-shoot) in the eval
        # model: the physical seasonal-band breach the soft relaxation absorbed.
        total_soc_band = if !isnothing(eval_model) &&
                            haskey(eval_model, :soc_band_over) && !isempty(eval_model[:soc_band_over])
            (value.(eval_model[:soc_band_over]) |> sum) + (value.(eval_model[:soc_band_under]) |> sum)
        else
            0.0
        end
        return new(
            name,
            seed,
            n_rep_periods,
            period_length,
            clustering_type,
            weight_type,
            data.normalization,
            data.tol,
            data.fix_every,
            projection_error,
            sigma_max,
            sigma_min,
            kappa,
            stable_rank,
            numerical_rank,
            mean_weight_sum,
            max_weight_sum,
            frac_weight_sum_gt1,
            frac_weight_sum_lt1,
            weight_density,
            pgd_total_iters,
            pgd_max_iters_per_fit,
            cache_hits,
            cache_misses,
            termination_status,
            objective_value,
            cost_of_investment,
            cost_of_operations,
            cost_of_spillage,
            cost_of_borrow,
            cost_of_soc_band,
            n_variables,
            n_constraints,
            evaluation_termination_status,
            evaluated_objective_value,
            eval_cost_of_investment,
            eval_cost_of_operations,
            eval_cost_of_spillage,
            eval_cost_of_borrow,
            eval_cost_of_soc_band,
            eval_n_variables,
            eval_n_constraints,
            total_spillage,
            total_borrow,
            total_soc_band,
            time_to_preprocess,
            time_to_cluster,
            time_to_fit_weights,
            time_to_formulate_model,
            time_to_solve,
        )
    end
end

Tables.columns(res::ExperimentResult) = (;
    name=[res.name],
    seed=[res.seed],
    n_rep_periods=[res.n_rep_periods],
    period_length=[res.period_length],
    clustering_type=[string(res.clustering_type)],
    weight_type=[string(res.weight_type)],
    normalization=[string(res.normalization)],
    tol=[res.tol],
    fix_every=[res.fix_every],
    projection_error=[res.projection_error],
    sigma_max=[res.sigma_max],
    sigma_min=[res.sigma_min],
    kappa=[res.kappa],
    stable_rank=[res.stable_rank],
    numerical_rank=[res.numerical_rank],
    mean_weight_sum=[res.mean_weight_sum],
    max_weight_sum=[res.max_weight_sum],
    frac_weight_sum_gt1=[res.frac_weight_sum_gt1],
    frac_weight_sum_lt1=[res.frac_weight_sum_lt1],
    weight_density=[res.weight_density],
    pgd_total_iters=[res.pgd_total_iters],
    pgd_max_iters_per_fit=[res.pgd_max_iters_per_fit],
    cache_hits=[res.cache_hits],
    cache_misses=[res.cache_misses],
    termination_status=[res.termination_status],
    objective_value=[res.objective_value],
    cost_of_investment=[res.cost_of_investment],
    cost_of_operations=[res.cost_of_operations],
    cost_of_spillage=[res.cost_of_spillage],
    cost_of_borrow=[res.cost_of_borrow],
    cost_of_soc_band=[res.cost_of_soc_band],
    n_variables=[res.n_variables],
    n_constraints=[res.n_constraints],
    evaluation_termination_status=[res.evaluation_termination_status],
    evaluated_objective_value=[res.evaluated_objective_value],
    eval_cost_of_investment=[res.eval_cost_of_investment],
    eval_cost_of_operations=[res.eval_cost_of_operations],
    eval_cost_of_spillage=[res.eval_cost_of_spillage],
    eval_cost_of_borrow=[res.eval_cost_of_borrow],
    eval_cost_of_soc_band=[res.eval_cost_of_soc_band],
    eval_n_variables=[res.eval_n_variables],
    eval_n_constraints=[res.eval_n_constraints],
    total_spillage=[res.total_spillage],
    total_borrow=[res.total_borrow],
    total_soc_band=[res.total_soc_band],
    time_to_preprocess=[res.time_to_preprocess],
    time_to_cluster=[res.time_to_cluster],
    time_to_fit_weights=[res.time_to_fit_weights],
    time_to_formulate_model=[res.time_to_formulate_model],
    time_to_solve=[res.time_to_solve],
)

end # module Types
