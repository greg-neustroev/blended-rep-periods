#!/usr/bin/env julia
#
# PROPOSED vs. single-component knockouts (economic normalization, conic selection, convex
# weights) on the synth dev system — the same arms `ablation_label`/`ablation_table` in
# summarize.jl print to the console. NOTE: unlike the other scripts in this directory, no table or
# prose paragraph in the currently-compiled paper cites ablation_synth.csv directly (tab:knockout
# covers the cache/FISTA/Gram knockouts instead, via knockout_ablation.jl) — kept for the console
# diagnostic and in case a future draft adds a component-ablation table; verify against the .tex
# before relying on it for a specific claim.
#
# Usage: julia --project=. analysis/paper/ablation_synth.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

function export_ablation()
    isfile(joinpath(REPO_ROOT, SYNTH_FILE)) || return
    df = CSV.read(joinpath(REPO_ROOT, SYNTH_FILE), DataFrame)
    for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
    arms = df[df.n_rep_periods .!= 1, :]
    refrows = df[df.n_rep_periods .== 1, :]
    isempty(refrows) && return
    ref_opt = mean(skipmissing(refrows.objective_value))
    arms.regret_pct = map(eachrow(arms)) do r
        (ismissing(r.evaluated_objective_value) || r.evaluated_objective_value <= 0) && return missing
        100 * (r.evaluated_objective_value / ref_opt - 1)
    end
    arms.total_time = reduce(.+, (coalesce.(getcol(arms, s), 0.0) for s in TIME_STAGES))
    arms.ablabel = [ablation_label(r.clustering_type, r.weight_type, r.normalization) for r in eachrow(arms)]
    sub = arms[.!ismissing.(arms.ablabel), :]
    isempty(sub) && return
    out = DataFrame()
    order = ["PROPOSED", "ablate: -economic", "ablate: -conic-selection", "ablate: -convex-weights"]
    for label in order, n in sort(unique(sub.n_rep_periods))
        rows = sub[(sub.ablabel .=== label) .& (sub.n_rep_periods .== n), :]
        isempty(rows) && continue
        (m_r, s_r) = meanstd(rows, :regret_pct); (m_t, s_t) = meanstd(rows, :total_time)
        sdr = isempty(s_r) ? 0.0 : parse(Float64, replace(s_r, "±" => ""))
        push!(out, (arm_label = label, n_rep_periods = n, n_seeds = nrow(rows),
                    regret_mean_pct = something(m_r, NaN), regret_sd_pct = sdr,
                    time_mean_s = something(m_t, NaN)); cols=:union)
    end
    write_csv("ablation_synth.csv", out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_ablation()
