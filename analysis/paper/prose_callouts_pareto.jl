#!/usr/bin/env julia
#
# case_studies.tex §Joint Time--Regret Trade-off prose (the paragraphs following Figure~fig:pareto):
# "on GEP, conical hull with Dirac weights reaches 14.2% regret at n_rep=10 in 1.9s, while the best
# conventional method... is still worse at four times as many RPs (16.0% at n_rep=40, 13.0s)..."
# etc. For every (case_study, n_rep_periods), the best (lowest mean regret) hull-family
# configuration and the best conventional configuration, with regret/time mean+-SD, so each cited
# pair can be checked against the .tex.
#
# Reads regret_summary.csv (already written by regret.jl) rather than recomputing from raw
# outputs, so this can never disagree with Figure~fig:pareto, which is built from the same file.
#
# Must run AFTER regret.jl (see run_all.jl's ordering).
#
# Usage: julia --project=. analysis/paper/prose_callouts_pareto.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

# best (lowest mean regret) row in `g` restricted to `group`'s clustering types, or nothing if none
# of `group`'s methods have a result at this cell.
function best_in_group(g, group)
    sub = g[in.(g.clustering_type, Ref(group)), :]
    sub = sub[.!isnan.(sub.regret_mean_pct), :]
    isempty(sub) && return nothing
    i = argmin(sub.regret_mean_pct)
    (regret_pct = sub.regret_mean_pct[i], regret_sd_pct = sub.regret_sd_pct[i],
     time_s = sub.time_mean_s[i], time_sd_s = sub.time_sd_s[i],
     clustering_type = sub.clustering_type[i], weight_type = sub.weight_type[i])
end

function export_prose_callouts_pareto()
    path = joinpath(OUTDIR, "regret_summary.csv")
    isfile(path) || (@warn "no regret_summary.csv in $OUTDIR — run regret.jl first"; return)
    df = CSV.read(path, DataFrame)
    for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
    df = df[df.normalization .== "economic", :]

    out = DataFrame()
    add(case_study, n, group_label, r) = r === nothing ? nothing : push!(out, (
        case_study = case_study, n_rep_periods = n, group = group_label,
        clustering_type = r.clustering_type, weight_type = r.weight_type,
        regret_pct = r.regret_pct, regret_sd_pct = r.regret_sd_pct,
        time_s = r.time_s, time_sd_s = r.time_sd_s); cols=:union)

    for cs in ("5bus", "gep", "p2x", "118bus"), n in sort(unique(df.n_rep_periods))
        g = df[(df.case_study .== cs) .& (df.n_rep_periods .== n), :]
        isempty(g) && continue
        add(cs, n, "hull", best_in_group(g, HULL))
        add(cs, n, "conventional", best_in_group(g, CONVENTIONAL))
    end
    write_csv("prose_pareto_callouts.csv", out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_prose_callouts_pareto()
