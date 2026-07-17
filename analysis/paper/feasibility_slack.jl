#!/usr/bin/env julia
#
# case_studies.tex's closing paragraph (§Weight-mass excess / feasibility check, end of the
# section): feasibility-slack usage (loss-of-load, spillage, borrow) at the paper's recommended
# configuration (conical hull, convex weights, n_rep=20), for all four real case studies — loss-of-
# load in MWh and as a % of annual demand, spillage in MWh, and borrow in MWh + its cost.
#
# Loss-of-load is dispatch from the virtual "ENS" generator (assets.csv technology=="ENS", priced
# at VOLL) — summed from eval_power_out.arrow, the same Arrow-dump pattern
# curtailment_renewable.jl uses for true renewable curtailment. Spillage/borrow come straight off
# the raw CSV's total_spillage/total_borrow/eval_cost_of_borrow columns (already averaged over
# seeds), since those are already tracked at the model level (no Arrow read needed).
#
# Usage: julia --project=. analysis/paper/feasibility_slack.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

const FEASIBILITY_DATASETS = [
    ("sienna/5bus", "5bus"), ("tyndp/gep", "gep"), ("tyndp/p2x", "p2x"), ("nrel/118bus", "118bus"),
]
const FOCUS_NRP = 20

function loss_of_load_mwh(ds, dsname, n)
    dir = joinpath(REPO_ROOT, "outputs", split(ds, "/")[1])
    d = arm_dir(dir, dsname, PROPOSED[1], PROPOSED[2], n)
    d === nothing && return missing
    f = joinpath(d, "eval_power_out.arrow")
    isfile(f) || return missing
    ens_ids = asset_ids_by_technology(ds, "ENS")
    isempty(ens_ids) && return 0.0
    df = DataFrame(Arrow.Table(f))
    df = df[in.(df.id, Ref(ens_ids)), :]
    isempty(df) && return 0.0
    by_seed = combine(groupby(df, :seed), :value => (v -> sum(skipmissing(v))) => :ens)
    mean(by_seed.ens)
end

function export_feasibility_slack()
    out = DataFrame()
    for (ds, case_study) in FEASIBILITY_DATASETS
        path = "outputs/$(ds).csv"
        # dataset id used in output directory/arm names is the LAST path component (gep, 5bus,
        # p2x, 118bus), matching CASE_META's file-vs-name convention.
        dsname = split(ds, "/")[2]
        loaded = load_case(path)
        loaded === nothing && (@warn "no reference optimum for $path — skipping feasibility slack"; continue)
        arms = canonical(loaded.arms)
        sub = arms[(arms.normalization .== "economic") .& (arms.clustering_type .== PROPOSED[1]) .&
                   (arms.weight_type .== PROPOSED[2]) .& (arms.n_rep_periods .== FOCUS_NRP), :]
        isempty(sub) && (@warn "no PROPOSED n_rp=$FOCUS_NRP rows for $ds"; continue)

        spillage = collect(skipmissing(getcol(sub, :total_spillage)))
        borrow = collect(skipmissing(getcol(sub, :total_borrow)))
        borrow_cost = collect(skipmissing(getcol(sub, :eval_cost_of_borrow)))
        ens = loss_of_load_mwh(ds, dsname, FOCUS_NRP)
        demand = total_demand_mwh(ds)

        push!(out, (
            case_study = case_study, n_rep_periods = FOCUS_NRP,
            loss_of_load_mwh = ens, loss_of_load_pct_of_demand = (ismissing(ens) || demand <= 0) ? missing : 100 * ens / demand,
            spillage_mwh = isempty(spillage) ? missing : mean(spillage),
            borrow_mwh = isempty(borrow) ? missing : mean(borrow),
            borrow_cost = isempty(borrow_cost) ? missing : mean(borrow_cost),
            annual_demand_mwh = demand,
        ); cols=:union)
    end
    write_csv("feasibility_slack.csv", out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_feasibility_slack()
