#!/usr/bin/env julia
#
# case_studies.tex §Runtime Breakdown and Scaling prose: "5-bus is a close call, with its n_rep=80
# time (2.34s) staying just under its own reference (2.42s); for GEP... the reduced model's total
# time overtakes the full-horizon reference between n_rep=40 (27.0s) and n_rep=80 (116.2s, vs. a
# 72.8s reference)... P2X and 118-bus... (45.2s vs 65.4s, and 43.2s vs 73.7s, at n_rep=80)". One
# row per case study x n_rep, PROPOSED's (conical hull + convex weights) total time (summed over
# every pipeline stage) next to the n_rep=1 full-horizon reference solve time, so each cited pair
# can be checked against the .tex.
#
# Computed directly from the raw per-case-study CSVs (NOT from runtime_breakdown.csv): that CSV's
# clustering grid is chosen for the figure's baseline spread and does not include conical_hull, so
# it can no longer answer a PROPOSED-specific prose callout.
#
# Usage: julia --project=. analysis/paper/prose_callouts_runtime.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

function export_prose_callouts_runtime()
    out = DataFrame()
    for (path, meta) in CASE_META
        loaded = load_case(path)
        loaded === nothing && (@warn "no reference optimum for $path — skipping"; continue)
        arms = canonical(loaded.arms)
        df = CSV.read(joinpath(REPO_ROOT, path), DataFrame)
        for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
        refrows = df[df.n_rep_periods .== 1, :]
        rt = collect(skipmissing(getcol(refrows, :time_to_solve)))
        ref_s = isempty(rt) ? missing : mean(rt)
        ref_sd_s = length(rt) > 1 ? std(rt) : 0.0

        prop = arms[(arms.normalization .== "economic") .& (arms.clustering_type .== PROPOSED[1]) .&
                    (arms.weight_type .== PROPOSED[2]), :]
        for n in sort(unique(prop.n_rep_periods))
            sub = prop[prop.n_rep_periods .== n, :]
            isempty(sub) && continue
            total_s = mean(sub.total_time)
            total_sd_s = nrow(sub) > 1 ? std(sub.total_time) : 0.0
            push!(out, (case_study = meta.case_study, n_rep_periods = n, proposed_total_time_s = total_s,
                        proposed_total_time_sd_s = total_sd_s,
                        full_reference_time_s = ref_s, full_reference_time_sd_s = ref_sd_s,
                        proposed_faster_than_reference = ismissing(ref_s) ? missing : total_s < ref_s
                       ); cols=:union)
        end
    end
    write_csv("prose_runtime_callouts.csv", out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_prose_callouts_runtime()
