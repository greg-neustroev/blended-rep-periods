#!/usr/bin/env julia
#
# Table tab:knockout (case_studies.tex): cache / FISTA / Gram-matrix knockouts, crossed over four
# hull x weight combinations (configs/synth_knockout.toml). Four arms per combination: proposed
# (cache=true, pgd=proposed), no_cache (cache=false, pgd=proposed), gram_only (cache=true,
# pgd=gram_only), plain (cache=true, pgd=plain). Regret + mismatch + total runtime, mean (± sd
# over 5 seeds; every arm here is deterministic so sd is ~0, per the table's own caption).
#
# knockout_table.csv: the printed table itself — one row per (hull/weight, arm), columns
# r_10rp/m_10rp/t_10rp/r_80rp/m_80rp/t_80rp, mirroring tab:knockout's "10 RPs"/"80 RPs" column
# blocks directly (the table only shows these two n_rp values, though the sweep also ran
# n_rp=20,40 — see knockout_ablation.csv for the full grid).
# knockout_ablation.csv: the general (every n_rep_periods in the sweep) long-format view.
#
# Usage: julia --project=. analysis/paper/knockout_ablation.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

const KNOCKOUT_COMBOS = [("convex_hull", "convex"), ("convex_hull", "conical"),
                          ("conical_hull", "convex"), ("conical_hull", "conical")]
const KNOCKOUT_ARMS = ["proposed", "no_cache", "gram_only", "plain"]
const TABLE_NRPS = (10, 80)   # the two n_rp values tab:knockout actually prints

# Name-suffix order is base_n_periodlen_clustering_weight_tol[_norm][_pgd_variant][_nocache]
# (see ExperimentData's constructor in Types.jl) — "nocache" is always the last suffix when
# present, so it cannot collide with a pgd_variant suffix (our sweep never combines cache=false
# with a non-proposed pgd_variant, but the check order below is robust even if it did).
function knockout_arm_of(name)
    occursin("_nocache", name) && return "no_cache"
    occursin("_gram_only", name) && return "gram_only"
    occursin("_plain", name) && return "plain"
    return "proposed"
end

function knockout_arms_df()
    isfile(joinpath(REPO_ROOT, SYNTH_FILE)) || return nothing
    df = CSV.read(joinpath(REPO_ROOT, SYNTH_FILE), DataFrame)
    for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
    refrows = df[df.n_rep_periods .== 1, :]
    isempty(refrows) && return nothing
    ref_opt = mean(skipmissing(refrows.objective_value))

    is_knockout_combo = [(cl, w) in KNOCKOUT_COMBOS for (cl, w) in zip(df.clustering_type, df.weight_type)]
    arms = df[(df.n_rep_periods .!= 1) .& (df.normalization .== "economic") .&
              isapprox.(df.tol, 0.01; atol=1e-12) .& is_knockout_combo, :]
    arms = copy(arms)
    arms.arm = knockout_arm_of.(arms.name)
    isempty(arms) && (@warn "no knockout-ablation rows found in $SYNTH_FILE — run configs/synth_knockout.toml first"; return nothing)
    arms.regret_pct = map(eachrow(arms)) do r
        (ismissing(r.evaluated_objective_value) || r.evaluated_objective_value <= 0) && return missing
        100 * (r.evaluated_objective_value / ref_opt - 1)
    end
    arms.mismatch_pct = map(eachrow(arms)) do r
        (ismissing(r.objective_value) || ismissing(r.evaluated_objective_value) ||
         r.evaluated_objective_value == 0) && return missing
        100 * abs(r.objective_value - r.evaluated_objective_value) / abs(r.evaluated_objective_value)
    end
    arms.total_time = reduce(.+, (coalesce.(getcol(arms, s), 0.0) for s in TIME_STAGES))
    return arms
end

function export_knockout_ablation()
    arms = knockout_arms_df()
    arms === nothing && return

    out = DataFrame()
    for (cl, w) in KNOCKOUT_COMBOS, arm in KNOCKOUT_ARMS, n in sort(unique(arms.n_rep_periods))
        sub = arms[(arms.clustering_type .== cl) .& (arms.weight_type .== w) .&
                   (arms.arm .== arm) .& (arms.n_rep_periods .== n), :]
        isempty(sub) && continue
        regret = collect(skipmissing(sub.regret_pct))
        mismatch = collect(skipmissing(sub.mismatch_pct))
        times = collect(skipmissing(sub.total_time))
        iters = collect(skipmissing(getcol(sub, :pgd_total_iters)))
        push!(out, (clustering_type = cl, weight_type = w, arm = arm, n_rep_periods = n,
                    n_seeds = nrow(sub),
                    regret_mean_pct = isempty(regret) ? missing : mean(regret),
                    regret_sd_pct = length(regret) > 1 ? std(regret) : 0.0,
                    mismatch_mean_pct = isempty(mismatch) ? missing : mean(mismatch),
                    mismatch_sd_pct = length(mismatch) > 1 ? std(mismatch) : 0.0,
                    time_mean_s = isempty(times) ? missing : mean(times),
                    time_sd_s = length(times) > 1 ? std(times) : 0.0,
                    iters_mean = isempty(iters) ? missing : mean(iters),
                    iters_sd = length(iters) > 1 ? std(iters) : 0.0); cols=:union)
    end
    write_csv("knockout_ablation.csv", out)

    # --- knockout_table.csv: the printed table, wide by n_rp in {10, 80} ---
    table_out = DataFrame()
    for (cl, w) in KNOCKOUT_COMBOS, arm in KNOCKOUT_ARMS
        row = Dict{Symbol,Any}(:clustering_type => cl, :weight_type => w, :arm => arm)
        for n in TABLE_NRPS
            cell = out[(out.clustering_type .== cl) .& (out.weight_type .== w) .&
                       (out.arm .== arm) .& (out.n_rep_periods .== n), :]
            row[Symbol("r_$(n)rp")] = isempty(cell) ? missing : first(cell.regret_mean_pct)
            row[Symbol("m_$(n)rp")] = isempty(cell) ? missing : first(cell.mismatch_mean_pct)
            row[Symbol("t_$(n)rp")] = isempty(cell) ? missing : first(cell.time_mean_s)
        end
        push!(table_out, row; cols=:union)
    end
    write_csv("knockout_table.csv", table_out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_knockout_ablation()
