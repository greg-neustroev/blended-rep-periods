#!/usr/bin/env julia
#
# Fig. central_regret_overview.pdf (analysis/R/regret.R) + case_studies.tex §Results prose
# (investment case studies §sec:results and storage-dispatch §sec:results subsections).
#
# regret_by_seed.csv: long, per-seed rows for every (case_study, clustering, weight, normalization,
# n_rep_periods, seed) — feeds regret.R's error bars and the Pareto-front time-vs-regret plot.
# regret_summary.csv: one row per (case_study, clustering, weight, normalization, n_rep_periods),
# mean/SD over the 5 seeds plus the Pareto-frontier and "reasonable in group" flags computed with
# analysis/common.jl's own functions, so the figure and the console diagnostic (summarize.jl) can
# never silently disagree.
#
# Usage: julia --project=. analysis/paper/regret.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

function export_regret()
    by_seed = DataFrame()
    summary = DataFrame()
    for (path, meta) in CASE_META
        loaded = load_case(path)
        loaded === nothing && (@warn "no reference optimum for $path — skipping"; continue)
        arms = canonical(loaded.arms)   # canonical tol=0.01, exclude nocache twins

        sub = select(arms, :n_rep_periods, :clustering_type, :weight_type, :normalization, :seed,
                     :regret_pct, :mismatch_pct, :total_time => :total_time_s,
                     :time_to_preprocess, :time_to_cluster, :time_to_fit_weights,
                     :time_to_formulate_model, :time_to_solve, :time_to_read,
                     :projection_error, :max_weight_sum, :cache_hit_rate, :pgd_total_iters,
                     :termination_status, :evaluation_termination_status,
                     :n_variables, :n_constraints, :eval_n_variables, :eval_n_constraints)
        sub.case_study .= meta.case_study; sub.problem_class .= meta.problem_class; sub.source .= meta.source
        by_seed = vcat(by_seed, sub; cols=:union)

        # Pareto-frontier and per-group "reasonable" flags, computed with common.jl's own
        # functions so the figures can never disagree with the console diagnostic. Computed
        # PER NORMALIZATION (not once over the pooled economic+unscaled+minmax rows): otherwise a
        # cell could be excluded from the frontier because some OTHER normalization's cell
        # dominates it, even though it is itself non-dominated within its own normalization -- the
        # figure only ever plots one normalization at a time, so cross-normalization domination
        # would mark points as "not on the frontier" for a comparison the plot never actually makes.
        nondom = union((nondominated_keys(arms[arms.normalization .== nrm, :]) for nrm in ordered_norms(arms))...)
        for nrm in ordered_norms(arms), n in sort(unique(arms.n_rep_periods))
            g = arms[(arms.normalization .== nrm) .& (arms.n_rep_periods .== n), :]
            isempty(g) && continue
            clusts = unique(g.clustering_type); weights = unique(g.weight_type)
            keys_ = [(cl, w) for cl in clusts for w in weights]
            samps = [samples(g[(g.clustering_type .== cl) .& (g.weight_type .== w), :]) for (cl, w) in keys_]
            tied = Set(keys_[i] for i in tied_with_best(samps))
            for (cl, w) in keys_
                rows = g[(g.clustering_type .== cl) .& (g.weight_type .== w), :]
                isempty(rows) && continue
                (m_r, s_r) = meanstd(rows, :regret_pct)
                (m_t, s_t) = meanstd(rows, :total_time)
                (m_m, s_m) = meanstd(rows, :mismatch_pct)
                sdr = isempty(s_r) ? 0.0 : parse(Float64, replace(s_r, "±" => ""))
                push!(summary, (
                    case_study = meta.case_study, problem_class = meta.problem_class, source = meta.source,
                    n_rep_periods = n, clustering_type = cl, weight_type = w, normalization = nrm,
                    n_seeds = nrow(rows), regret_mean_pct = something(m_r, NaN), regret_sd_pct = sdr,
                    mismatch_mean_pct = something(m_m, NaN),
                    time_mean_s = something(m_t, NaN),
                    is_deterministic = sdr < 1e-9,
                    pareto_frontier = (nrm, n, cl, w) in nondom,
                    reasonable_in_group = (cl, w) in tied,
                ); cols=:union)
            end
        end
    end
    write_csv("regret_by_seed.csv", by_seed)
    write_csv("regret_summary.csv", summary)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_regret()
