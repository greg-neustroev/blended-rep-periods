#!/usr/bin/env julia
#
# case_studies.tex §Results prose inline regret callouts (both subsections: "Investment Case
# Studies: Clustering Dominates" and "Economic Dispatch with Seasonal Storage"), e.g. "At n_rep=20
# with Dirac weights, 5-bus's regret ranges from 1.2% (conical hull) to 662.4% (k-means)". One row
# per cited range, labeled to match the sentence, so each can be checked against the .tex at a
# glance. Reads regret_summary.csv (already written by regret.jl) rather than recomputing regret
# from raw outputs, so this can never disagree with the figure it accompanies.
#
# Must run AFTER regret.jl (see run_all.jl's ordering).
#
# Usage: julia --project=. analysis/paper/prose_callouts_regret.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

const FOCUS_NRP = 20

# min/max regret_mean_pct over `varying` (:clustering_type or :weight_type), holding the other
# fixed to `fixed_value`, at n_rep=FOCUS_NRP, economic normalization — plus which method achieves
# each extreme, matching how the prose phrases each range ("from X% (method) to Y% (method)").
function regret_range(df, case_study, varying, fixed_col, fixed_value)
    g = df[(df.case_study .== case_study) .& (df.normalization .== "economic") .&
           (df.n_rep_periods .== FOCUS_NRP) .& (getproperty(df, fixed_col) .== fixed_value), :]
    isempty(g) && return nothing
    g = g[.!isnan.(g.regret_mean_pct), :]
    isempty(g) && return nothing
    lo = argmin(g.regret_mean_pct); hi = argmax(g.regret_mean_pct)
    (min_pct = g.regret_mean_pct[lo], min_method = string(getproperty(g, varying)[lo]),
     max_pct = g.regret_mean_pct[hi], max_method = string(getproperty(g, varying)[hi]))
end

function export_prose_callouts_regret()
    path = joinpath(OUTDIR, "regret_summary.csv")
    isfile(path) || (@warn "no regret_summary.csv in $OUTDIR — run regret.jl first"; return)
    df = CSV.read(path, DataFrame)
    for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end

    out = DataFrame()
    add(case_study, callout, r) = r === nothing ? nothing : push!(out, (
        case_study = case_study, callout = callout, min_pct = r.min_pct, min_method = r.min_method,
        max_pct = r.max_pct, max_method = r.max_method); cols=:union)

    for cs in ("5bus", "gep")
        add(cs, "clustering_varies_dirac_weight_n_rp20",
            regret_range(df, cs, :clustering_type, :weight_type, "dirac"))
        add(cs, "weight_varies_conical_hull_clustering_n_rp20",
            regret_range(df, cs, :weight_type, :clustering_type, "conical_hull"))
    end
    for cs in ("p2x", "118bus")
        add(cs, "weight_varies_conical_hull_clustering_n_rp20",
            regret_range(df, cs, :weight_type, :clustering_type, "conical_hull"))
        add(cs, "clustering_varies_dirac_weight_n_rp20",
            regret_range(df, cs, :clustering_type, :weight_type, "dirac"))
    end
    write_csv("prose_regret_callouts.csv", out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_prose_callouts_regret()
