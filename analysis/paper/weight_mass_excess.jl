#!/usr/bin/env julia
#
# methodology.tex §Sub-unit conic weights / Conic weights prose: "every conic-weight configuration
# we ran has at least one base period with sum_r W_{d,r} > 1, up to 3.85x on 5-bus and 1.27-1.66x
# on the other three" — the weight-mass excess of the unbounded conic weight class (Sigma w <= 1
# is what the sub-unit class enforces; plain conic drops that constraint entirely, Section
# sec:conical). `max_weight_sum` and `frac_weight_sum_gt1` are already tracked per run (raw output
# CSV columns); this just pulls the extremes across every clustering x n_rep x seed configuration
# run at weight_type=="conical", per case study.
#
# Usage: julia --project=. analysis/paper/weight_mass_excess.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

function export_weight_mass_excess()
    out = DataFrame()
    for (path, meta) in CASE_META
        loaded = load_case(path)
        loaded === nothing && continue
        arms = canonical(loaded.arms)
        conic = arms[(arms.weight_type .== "conical") .& (arms.normalization .== "economic"), :]
        isempty(conic) && continue
        maxsum = collect(skipmissing(getcol(conic, :max_weight_sum)))
        fracgt1 = collect(skipmissing(getcol(conic, :frac_weight_sum_gt1)))
        isempty(maxsum) && continue
        push!(out, (
            case_study = meta.case_study, n_configs = nrow(conic),
            max_weight_sum_excess = maximum(maxsum),
            n_configs_with_any_period_over_1 = count(>(0), fracgt1),
            n_configs_total_with_frac_data = length(fracgt1),
        ); cols=:union)
    end
    write_csv("weight_mass_excess.csv", out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_weight_mass_excess()
