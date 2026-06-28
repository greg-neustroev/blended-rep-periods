"""
    save_result_to_csv(path, result, time_to_read)

Append `result` (with `time_to_read` filled in) as a row to the CSV at `path`,
writing the header only when the file does not yet exist.
"""
function save_result_to_csv(path::String, result::ExperimentResult, time_to_read::Float64)
    row = result |> DataFrame
    row.time_to_read .= time_to_read
    write_header = !isfile(path)
    CSV.write(path, row; append=true, writeheader=write_header)
end

# Every primal variable container worth dumping in full, with the names of its
# index dimensions (the first axis is always the asset/line id). Saving these
# (plus the duals below) means downstream metrics — curtailment, generation mix,
# storage dispatch, ENS, congestion — are derivable offline without re-solving.
const SOLUTION_VARIABLES = [
    (:invested_units, Symbol[]),
    (:power_out, [:rep_period, :timestep]),
    (:power_in, [:rep_period, :timestep]),
    (:state_of_charge_intra, [:rep_period, :timestep]),
    (:state_of_charge_intra_0, [:rep_period]),
    (:state_of_charge_inter, [:period]),
    (:state_of_charge_inter_0, Symbol[]),
    (:spillage, [:rep_period, :timestep]),
    (:borrow, [:rep_period, :timestep]),
    (:soc_band_over, [:period]),
    (:soc_band_under, [:period]),
    (:flow, [:rep_period, :timestep]),
]

# Tidy DataFrame of a solved variable container: one row per index tuple, columns
# `[:seed, :id, index_names..., :value]`. Returns `nothing` when the model lacks
# the variable or the container is empty (e.g. `invested_units` on a dispatch-only
# dataset).
function _variable_dataframe(model, varname::Symbol, index_names::Vector{Symbol}, seed::Int)
    haskey(model, varname) || return nothing
    header = [:id, index_names..., :variable]
    df = Containers.rowtable(model[varname]; header=header) |> DataFrame
    isempty(df) && return nothing
    df.value = value.(df.variable)
    select!(df, Not(:variable))
    df.seed .= seed
    select!(df, Cols(:seed, Not(:seed)))
    return df
end

"""
    save_solution_to_arrow(model, outputs_dir, result_name, seed; kind)

Dump the full primal solution of `model` — every variable in
[`SOLUTION_VARIABLES`](@ref) — and, when available, the nodal balance duals
(marginal/clearing prices) to Arrow files under `<outputs_dir>/<result_name>/`,
each prefixed with `kind` (`"reduced"` or `"eval"`). Arrow keeps the full-horizon
arrays compact. Does nothing if `model` is not solved to a feasible point.
"""
function save_solution_to_arrow(
    model,
    outputs_dir::AbstractString,
    result_name::AbstractString,
    seed::Int;
    kind::AbstractString,
)
    JuMP.is_solved_and_feasible(model) || return
    subdir = joinpath(outputs_dir, result_name)
    mkpath(subdir)
    for (varname, index_names) in SOLUTION_VARIABLES
        df = _variable_dataframe(model, varname, index_names, seed)
        df === nothing && continue
        Arrow.write(joinpath(subdir, "$(kind)_$(varname).arrow"), df)
    end
    # Nodal/marginal prices: duals of the retained balance constraints (LMPs).
    # Collect keys and constraints into aligned vectors before reading duals.
    if JuMP.has_duals(model) && haskey(model.ext, :balance_constraints)
        bc = model.ext[:balance_constraints]
        if !isempty(bc)
            ks = collect(keys(bc))
            price_df = DataFrame(
                location=[string(k[1]) for k in ks],
                carrier=[string(k[2]) for k in ks],
                rep_period=[k[3] for k in ks],
                timestep=[k[4] for k in ks],
                price=[JuMP.dual(bc[k]) for k in ks],
            )
            price_df.seed .= seed
            Arrow.write(joinpath(subdir, "$(kind)_nodal_prices.arrow"), price_df)
        end
    end
    return
end

# Tidy (period, rep_period, weight) table of the non-zero blending weights.
function _weight_table(W)
    periods = Int[]
    reps = Int[]
    weights = Float64[]
    nrows, ncols = size(W)
    for p in 1:nrows, r in 1:ncols
        w = W[p, r]
        iszero(w) && continue
        push!(periods, p)
        push!(reps, r)
        push!(weights, w)
    end
    return DataFrame(period=periods, rep_period=reps, weight=weights)
end

"""
    save_clustering_artifacts(clustering_result, outputs_dir, result_name, seed)

Dump the clustering-side artifacts to Arrow under `<outputs_dir>/<result_name>/`:
the blended-weight matrix, the nearest-representative assignment, the selected
representative base-period indices, the singular-value spectrum of the RP matrix,
the per-base-period projection-error residuals, the per-fit PGD iteration counts,
and the per-greedy-iteration cache hit/miss counts. These make every weight,
geometry, and convergence analysis reproducible offline. Missing pieces (e.g. the
single-period fast path, or `selected_indices` for k-means centroids) are skipped.
"""
function save_clustering_artifacts(
    clustering_result::ClusteringResult,
    outputs_dir::AbstractString,
    result_name::AbstractString,
    seed::Int,
)
    subdir = joinpath(outputs_dir, result_name)
    mkpath(subdir)
    diag = clustering_result.diagnostics
    function write_df(name, df)
        df[!, :seed] = fill(seed, nrow(df))
        Arrow.write(joinpath(subdir, name), df)
    end

    W = clustering_result.weight_matrix
    W === nothing || write_df("weights.arrow", _weight_table(W))

    # Storage-chain weights W^ch (when the chain split was fit), dumped in the same tidy
    # (period, rep_period, weight) form as the operational weights, plus a one-row
    # `chain_fit.arrow` recording the chain class and the increment-reconstruction
    # residual (the in-hull / out-of-hull diagnostic).
    Wch = clustering_result.chain_weight_matrix
    if Wch !== nothing
        write_df("chain_weights.arrow", _weight_table(Wch))
        write_df("chain_fit.arrow", DataFrame(
            chain_weight_type=[string(get(diag, :chain_weight_type, missing))],
            chain_fit_residual=[get(diag, :chain_fit_residual, missing)],
            chain_max_abs_weight=[get(diag, :chain_max_abs_weight, missing)],
            chain_max_row_l1=[get(diag, :chain_max_row_l1, missing)],
        ))
    end

    R = clustering_result.rp_matrix
    if R !== nothing && !isempty(R)
        svals = svdvals(Matrix{Float64}(R))
        write_df("rp_spectrum.arrow", DataFrame(index=1:length(svals), singular_value=svals))
        C = clustering_result.clustering_matrix
        if W !== nothing && C !== nothing && size(C, 2) == size(W, 1)
            residual = R * Matrix(W)' - C                       # features × base periods
            per_period = vec(sqrt.(sum(abs2, residual; dims=1)))
            write_df("projection_error_per_period.arrow",
                DataFrame(period=1:length(per_period), residual=per_period))
        end
    end

    if haskey(diag, :assignments) && diag[:assignments] !== nothing
        a = diag[:assignments]
        write_df("rp_assignment.arrow", DataFrame(period=1:length(a), rep_period=a))
    end
    if haskey(diag, :selected_indices) && diag[:selected_indices] !== nothing
        idx = diag[:selected_indices]
        write_df("selected_rp_indices.arrow", DataFrame(rep_period=1:length(idx), base_period=idx))
    end
    if haskey(diag, :pgd_iters) && !isempty(diag[:pgd_iters])
        it = diag[:pgd_iters]
        write_df("pgd_iterations.arrow", DataFrame(fit_index=1:length(it), iterations=it))
    end
    if haskey(diag, :cache_hits_per_iter)
        hits = diag[:cache_hits_per_iter]
        misses = get(diag, :cache_misses_per_iter, zeros(Int, length(hits)))
        write_df("cache_per_iteration.arrow",
            DataFrame(greedy_iteration=1:length(hits), hits=hits, misses=misses))
    end
    return
end
