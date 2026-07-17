#!/usr/bin/env julia
#
# Fig. runtime_breakdown.pdf (analysis/R/runtime_breakdown.R) + case_studies.tex §Runtime
# Breakdown and Scaling prose. All 5 timing stages (keeps model formulation as its own segment),
# PROPOSED vs a k_means/k_medoids/hierarchical-dirac baseline, per case study, vs n_rp; plus the
# n_rep=1 full-horizon reference solve time (the prose's "2.42s"/"72.8s"/etc. reference numbers).
#
# Usage: julia --project=. analysis/paper/runtime_breakdown.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

const STAGE_COLS = [:time_to_preprocess, :time_to_cluster, :time_to_fit_weights,
                    :time_to_formulate_model, :time_to_solve]
const STAGE_NAMES = ["read_preprocess", "cluster", "fit_weights", "formulate_model", "solve"]

function export_runtime_breakdown()
    out = DataFrame()
    for (path, meta) in CASE_META
        df = CSV.read(joinpath(REPO_ROOT, path), DataFrame)
        for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
        refrows = df[df.n_rep_periods .== 1, :]
        methods = [("PROPOSED", PROPOSED[1], PROPOSED[2]), ("k_means", "k_means", "dirac"),
                   ("k_medoids", "k_medoids", "dirac"), ("hierarchical", "hierarchical", "dirac")]
        for (label, cl, w) in methods
            for n in sort(unique(df.n_rep_periods))
                n == 1 && continue
                sub = df[(df.normalization .== "economic") .& isapprox.(df.tol, 0.01; atol=1e-12) .&
                         .!occursin.("nocache", df.name) .&
                         (df.clustering_type .== cl) .& (df.weight_type .== w) .& (df.n_rep_periods .== n), :]
                isempty(sub) && continue
                for (col, stage) in zip(STAGE_COLS, STAGE_NAMES)
                    v = collect(skipmissing(getcol(sub, col))); isempty(v) && continue
                    push!(out, (case_study = meta.case_study, method_label = label, n_rep_periods = n,
                                stage = stage, time_mean_s = mean(v),
                                time_sd_s = length(v) > 1 ? std(v) : 0.0); cols=:union)
                end
            end
        end
        # n_rep=1 full-model reference solve time (read_preprocess + solve; no clustering/fitting)
        if !isempty(refrows)
            rt = collect(skipmissing(getcol(refrows, :time_to_solve)))
            isempty(rt) || push!(out, (case_study = meta.case_study, method_label = "full_reference",
                                       n_rep_periods = 1, stage = "solve", time_mean_s = mean(rt),
                                       time_sd_s = length(rt) > 1 ? std(rt) : 0.0); cols=:union)
        end
    end
    write_csv("runtime_breakdown.csv", out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_runtime_breakdown()
