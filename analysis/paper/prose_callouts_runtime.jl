#!/usr/bin/env julia
#
# case_studies.tex §Runtime Breakdown and Scaling prose: "5-bus is a close call, with its n_rep=80
# time (2.34s) staying just under its own reference (2.42s); for GEP... the reduced model's total
# time overtakes the full-horizon reference between n_rep=40 (27.0s) and n_rep=80 (116.2s, vs. a
# 72.8s reference)... P2X and 118-bus... (45.2s vs 65.4s, and 43.2s vs 73.7s, at n_rep=80)". One
# row per case study x n_rep, PROPOSED's total time (summed over every pipeline stage) next to the
# n_rep=1 full-horizon reference solve time, so each cited pair can be checked against the .tex.
# Reads runtime_breakdown.csv (already written by runtime_breakdown.jl) rather than recomputing,
# so this can never disagree with the figure it accompanies.
#
# Must run AFTER runtime_breakdown.jl (see run_all.jl's ordering).
#
# Usage: julia --project=. analysis/paper/prose_callouts_runtime.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

function export_prose_callouts_runtime()
    path = joinpath(OUTDIR, "runtime_breakdown.csv")
    isfile(path) || (@warn "no runtime_breakdown.csv in $OUTDIR — run runtime_breakdown.jl first"; return)
    df = CSV.read(path, DataFrame)

    out = DataFrame()
    for cs in unique(df.case_study)
        d = df[df.case_study .== cs, :]
        ref_rows = d[d.method_label .== "full_reference", :]
        ref_s = isempty(ref_rows) ? missing : sum(ref_rows.time_mean_s)
        for n in sort(unique(d.n_rep_periods[d.n_rep_periods .!= 1]))
            prop = d[(d.method_label .== "Conic, Convex (proposed)") .& (d.n_rep_periods .== n), :]
            isempty(prop) && continue
            total_s = sum(prop.time_mean_s)
            push!(out, (case_study = cs, n_rep_periods = n, proposed_total_time_s = total_s,
                        full_reference_time_s = ref_s,
                        proposed_faster_than_reference = ismissing(ref_s) ? missing : total_s < ref_s
                       ); cols=:union)
        end
    end
    write_csv("prose_runtime_callouts.csv", out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_prose_callouts_runtime()
