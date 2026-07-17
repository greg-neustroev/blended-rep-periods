#!/usr/bin/env julia
#
# Table tab:sensitivity (case_studies.tex) + its prose, and the "backstop reached in exactly 20"
# claim right below it. Synth dev system (outputs/synth/hydro.csv), n_rep=10.
#
# sensitivity_table.csv: the table itself — one row per tolerance (the PROPOSED combo, conical
# hull + convex weights, at n_rep=10), columns iters, then {economic,peak,minmax} x {regret_pct,
# mismatch_pct, time_s} — mirrors tab:sensitivity's column blocks directly ("peak" = normalization
# "unscaled" in the raw data; the paper calls it "peak-scaled").
# sensitivity_tol.csv: the general (every clustering x weight x normalization) tol-sensitivity
# view, min/median/max/range — supporting data, not itself a printed table.
# pgd_iters_by_tol.csv: total PGD iterations vs tol, pooled over methods/normalizations (the
# "1.7x, 42,163 to 71,366" prose claim).
# sensitivity_by_seed.csv: per-seed regret/mismatch for PROPOSED, tol x normalization — lets a
# paired t-test be run downstream (kept long/tidy; not a printed table).
# sensitivity_stats.csv: the paired t-tests themselves (tol: loosest vs tightest; normalization:
# every present pair) so the "moves regret by at most 0.29pp" / "order of magnitude smaller than
# the gap between normalization schemes" prose numbers are never hand-copied from a console.
# sensitivity_backstop.csv: the "1,505 configurations... backstop reached in exactly 20" claim.
#
# Usage: julia --project=. analysis/paper/sensitivity.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

const NORM_COLNAME = Dict("economic" => "economic", "unscaled" => "peak", "minmax" => "minmax")

function export_sensitivity()
    isfile(joinpath(REPO_ROOT, SYNTH_FILE)) || return
    df = CSV.read(joinpath(REPO_ROOT, SYNTH_FILE), DataFrame)
    for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
    refrows = df[df.n_rep_periods .== 1, :]
    isempty(refrows) && return
    ref_opt = mean(skipmissing(refrows.objective_value))
    # Exclude both cache and pgd_variant knockout-ablation arms (configs/synth_knockout.toml shares
    # this same output file): this sweep is about tolerance/normalization sensitivity for the
    # PROPOSED method only, not acceleration variants, which need far more iterations and would
    # otherwise skew the pooled iteration-count statistics below.
    sens = df[(df.n_rep_periods .== 10) .& .!occursin.("nocache", df.name) .&
              .!occursin.("_gram_only", df.name) .& .!occursin.("_plain", df.name), :]
    (isempty(sens) || length(unique(sens.tol)) < 2) && return
    sens.regret_pct = map(eachrow(sens)) do r
        (ismissing(r.evaluated_objective_value) || r.evaluated_objective_value <= 0) && return missing
        100 * (r.evaluated_objective_value / ref_opt - 1)
    end
    sens.mismatch_pct = map(eachrow(sens)) do r
        (ismissing(r.objective_value) || ismissing(r.evaluated_objective_value) ||
         r.evaluated_objective_value == 0) && return missing
        100 * abs(r.objective_value - r.evaluated_objective_value) / abs(r.evaluated_objective_value)
    end
    sens.total_time = reduce(.+, (coalesce.(getcol(sens, s), 0.0) for s in TIME_STAGES))

    # --- sensitivity_table.csv: the printed table, PROPOSED only, wide by normalization ---
    prop_sens = sens[(sens.clustering_type .== PROPOSED[1]) .& (sens.weight_type .== PROPOSED[2]), :]
    table_out = DataFrame()
    for t in sort(unique(prop_sens.tol); rev=true)
        gt = prop_sens[isapprox.(prop_sens.tol, t; atol=1e-12), :]
        isempty(gt) && continue
        row = Dict{Symbol,Any}(:tol => t)
        it = collect(skipmissing(gt[gt.normalization .== "economic", :pgd_total_iters]))
        row[:iters_economic] = isempty(it) ? missing : round(Int, mean(it))
        for nrm in ordered_norms(gt)
            colname = NORM_COLNAME[nrm]
            g = gt[gt.normalization .== nrm, :]
            r = collect(skipmissing(g.regret_pct)); m = collect(skipmissing(g.mismatch_pct))
            tt = collect(skipmissing(g.total_time))
            row[Symbol("$(colname)_regret_pct")] = isempty(r) ? missing : mean(r)
            row[Symbol("$(colname)_mismatch_pct")] = isempty(m) ? missing : mean(m)
            row[Symbol("$(colname)_time_s")] = isempty(tt) ? missing : mean(tt)
        end
        push!(table_out, row; cols=:union)
    end
    write_csv("sensitivity_table.csv", table_out)

    # --- sensitivity_tol.csv: general (every clustering x weight) tol-sensitivity view ---
    tol_out = DataFrame()
    for nrm in ordered_norms(sens), cl in METHOD_ORDER, w in WEIGHT_ORDER
        mask = (sens.normalization .== nrm) .& (sens.clustering_type .== cl) .& (sens.weight_type .== w)
        g = sens[mask, :]; length(unique(g.tol)) < 2 && continue
        for t in sort(unique(g.tol))
            gt = g[isapprox.(g.tol, t; atol=1e-12), :]
            v = collect(skipmissing(gt.regret_pct)); isempty(v) && continue
            tt = collect(skipmissing(gt.total_time))
            push!(tol_out, (normalization = nrm, clustering_type = cl, weight_type = w, tol = t,
                            regret_min_pct = minimum(v), regret_median_pct = median(v),
                            regret_max_pct = maximum(v), regret_range_pct = maximum(v) - minimum(v),
                            median_time_s = isempty(tt) ? missing : median(tt)); cols=:union)
        end
    end
    write_csv("sensitivity_tol.csv", tol_out)

    iters_out = DataFrame()
    for t in sort(unique(sens.tol); rev=true)
        it = collect(skipmissing(sens[isapprox.(sens.tol, t; atol=1e-12), :pgd_total_iters]))
        isempty(it) && continue
        push!(iters_out, (tol = t, median_pgd_iters = median(it)); cols=:union)
    end
    write_csv("pgd_iters_by_tol.csv", iters_out)

    # Per-seed tol × normalization grid for PROPOSED (conical_hull/convex), with mismatch --
    # keeps the seed column (unlike sensitivity_tol.csv's min/median/max) so a paired t-test across
    # tol or across normalization is possible downstream.
    prop = prop_sens
    by_seed = select(prop, :normalization, :tol, :seed, :regret_pct, :mismatch_pct)
    write_csv("sensitivity_by_seed.csv", by_seed)

    stats = DataFrame()
    tols = sort(unique(prop.tol))
    if length(tols) >= 2
        tlo, thi = minimum(tols), maximum(tols)
        for nrm in ordered_norms(prop)
            g = prop[prop.normalization .== nrm, :]
            glo = Dict(r.seed => r.regret_pct for r in eachrow(g[isapprox.(g.tol, tlo; atol=1e-15), :]) if !ismissing(r.regret_pct))
            ghi = Dict(r.seed => r.regret_pct for r in eachrow(g[isapprox.(g.tol, thi; atol=1e-15), :]) if !ismissing(r.regret_pct))
            ks = sort(collect(intersect(keys(glo), keys(ghi))))
            length(ks) < 2 && continue
            d = [ghi[k] - glo[k] for k in ks]   # tightest − loosest (negative ⇒ tighter tol helps)
            p = all(iszero, d) ? 1.0 : pvalue(OneSampleTTest(d))
            push!(stats, (comparison = "tol", normalization = nrm, a = thi, b = tlo,
                          mean_diff = mean(d), p_value = p, n_paired = length(ks)); cols=:union)
        end
    end
    norms_present = ordered_norms(prop)
    for i in eachindex(norms_present), j in (i+1):length(norms_present)
        r = normalization_compare_pair(prop, norms_present[i], norms_present[j])
        (r === nothing || r.n_total == 0) && continue
        push!(stats, (comparison = "normalization", normalization = "$(r.normA) vs $(r.normB)",
                      a = missing, b = missing, mean_diff = r.mean_diff, p_value = r.p_value,
                      n_paired = r.n_total, n_a_better = r.n_a_better, n_b_better = r.n_b_better,
                      n_comparable = r.n_comparable); cols=:union)
    end
    write_csv("sensitivity_stats.csv", stats)

    # --- sensitivity_backstop.csv: "1,505 configurations... backstop reached in exactly 20" ---
    # Over the FULL synth/hydro.csv (every tested tolerance, normalization, clustering method,
    # weight type, and n_rep — the entire proposed-method sensitivity+ablation sweep), not just the
    # n_rep=10 PROPOSED slice above. The backstop (Algorithm PGD's N_max=1e5) caps a SINGLE PGD
    # call; `pgd_total_iters` sums iterations over every weight-fitting/hull-distance call in the
    # whole run (routinely hundreds of thousands even when no individual call is capped), so the
    # right column to threshold is `pgd_max_iters_per_fit` — the worst single call in that run.
    all_arms = df[df.n_rep_periods .!= 1, :]
    proposed_arms = all_arms[.!occursin.("nocache", all_arms.name) .& .!occursin.("_gram_only", all_arms.name) .&
                              .!occursin.("_plain", all_arms.name), :]
    backstop = proposed_arms[coalesce.(proposed_arms.pgd_max_iters_per_fit, 0) .>= 100_000, :]
    combo_key(r) = (r.clustering_type, r.weight_type, r.normalization, r.n_rep_periods, r.tol)
    distinct_combos = unique(combo_key(r) for r in eachrow(backstop))
    backstop_out = DataFrame(
        total_configurations = nrow(proposed_arms),
        n_seeds_at_backstop = nrow(backstop),
        n_distinct_configurations_at_backstop = length(distinct_combos),
        distinct_configurations = join([join(string.(k), "/") for k in distinct_combos], "; "),
    )
    write_csv("sensitivity_backstop.csv", backstop_out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_sensitivity()
