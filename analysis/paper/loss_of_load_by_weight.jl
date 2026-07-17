#!/usr/bin/env julia
#
# methodology.tex prose: "maximum conic loss-of-load is 0.0014-3.09% of annual demand, against
# 0.0009-3.27% for the other three weight types, with Dirac the worst offender by up to 24x."
#
# loss_of_load_by_weight.csv: for every case study x weight type, the WORST (maximum)
# loss-of-load-as-%-of-demand across every clustering method and n_rep in the economic-
# normalization canonical sweep — backs the "0.0014-3.09%" / "0.0009-3.27%" ranges (the min/max
# across case studies of this per-weight-type max).
# loss_of_load_dirac_ratio.csv: the "Dirac the worst offender by up to 24x" claim specifically —
# read literally as a PAIRED comparison (same case study, clustering method, and n_rep), not
# independent grid-wide maxima: for every (case_study, clustering, n_rep) cell, the ratio of
# Dirac's loss-of-load to the smallest of the other three weight types' loss-of-load AT THAT SAME
# cell, then the worst (largest) such ratio per case study.
#
# Reads eval_power_out.arrow for every (case_study, weight_type, clustering_type, n_rep) arm
# directory — the full grid, not just PROPOSED/baselines like feasibility_slack.jl and
# curtailment_renewable.jl — so this is the most expensive script in analysis/paper/ to run
# (roughly 4 case studies x 4 weight types x 7 clustering methods x 4 n_rep Arrow reads).
#
# Usage: julia --project=. analysis/paper/loss_of_load_by_weight.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

const FEASIBILITY_DATASETS = [
    ("sienna/5bus", "5bus"), ("tyndp/gep", "gep"), ("tyndp/p2x", "p2x"), ("nrel/118bus", "118bus"),
]
const NRPS = (10, 20, 40, 80)

function loss_of_load_pct(ds, dsname, cl, w, n, ens_ids, demand)
    (isempty(ens_ids) || demand <= 0) && return missing
    dir = joinpath(REPO_ROOT, "outputs", split(ds, "/")[1])
    d = arm_dir(dir, dsname, cl, w, n)
    d === nothing && return missing
    f = joinpath(d, "eval_power_out.arrow")
    isfile(f) || return missing
    df = DataFrame(Arrow.Table(f))
    df = df[in.(df.id, Ref(ens_ids)), :]
    isempty(df) && return 0.0
    by_seed = combine(groupby(df, :seed), :value => (v -> sum(skipmissing(v))) => :ens)
    100 * mean(by_seed.ens) / demand
end

function export_loss_of_load_by_weight()
    out = DataFrame()
    ratios = DataFrame()
    for (ds, case_study) in FEASIBILITY_DATASETS
        dsname = split(ds, "/")[2]
        ens_ids = asset_ids_by_technology(ds, "ENS")
        demand = total_demand_mwh(ds)

        # cache every (clustering, n_rep, weight) cell once, reused by both outputs below.
        cell = Dict{Tuple{String,Int,String},Float64}()
        for cl in METHOD_ORDER, n in NRPS, w in WEIGHT_ORDER
            pct = loss_of_load_pct(ds, dsname, cl, w, n, ens_ids, demand)
            ismissing(pct) || (cell[(cl, n, w)] = pct)
        end

        for w in WEIGHT_ORDER
            vals = [v for ((cl, n, ww), v) in cell if ww == w]
            push!(out, (case_study = case_study, weight_type = w,
                        max_loss_of_load_pct_of_demand = isempty(vals) ? missing : maximum(vals)); cols=:union)
        end

        worst_ratio = missing; worst_combo = ""
        for cl in METHOD_ORDER, n in NRPS
            haskey(cell, (cl, n, "dirac")) || continue
            others = [cell[(cl, n, w)] for w in ("convex", "conical_bounded", "conical") if haskey(cell, (cl, n, w))]
            isempty(others) && continue
            other_min = minimum(others)
            other_min <= 0 && continue
            ratio = cell[(cl, n, "dirac")] / other_min
            (ismissing(worst_ratio) || ratio > worst_ratio) && (worst_ratio = ratio; worst_combo = "$cl/n_rp=$n")
        end
        push!(ratios, (case_study = case_study, dirac_worst_ratio = worst_ratio,
                       worst_combo = worst_combo); cols=:union)
    end
    write_csv("loss_of_load_by_weight.csv", out)
    write_csv("loss_of_load_dirac_ratio.csv", ratios)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_loss_of_load_by_weight()
