#!/usr/bin/env julia
#
# Summarize the revision-campaign results. For each result CSV it prints:
#   1. Comparison  — regret (%) vs n_rp, traditional (k-means/k-medoids/hierarchical/
#                    chronological) vs the PROPOSED method, mean ± std over seeds.
#   2. Ablation    — PROPOSED vs each single-component knockout, regret (%) vs n_rp.
#   3. Secondary   — at a representative n_rp: curtailment (spillage), infeasibility (borrow),
#                    upper-bound feasibility (max weight-sum), and the true cost breakdown.
#   4. Capacity    — (investment datasets) ‖Δ invested_units‖₁ vs the full-horizon optimum.
#   5. Runtime     — per-stage timing (cluster / weight-fit / solve) vs n_rp, + realized PGD
#                    iterations N(ε) and greedy-hull cache hit-rate.
#   6. Sensitivity — PROPOSED at n_rp=10 vs tol × normalization, with N(ε); + a cache-on/off
#                    validation (identical regret, cache hit-rate, speed-up).
#
# Regret = evaluated_objective_value / reference_optimum − 1 (reference = the n_rep=1 row).
# Only k-means is seed-sensitive; the hull/k-medoids/hierarchical/chronological arms are
# deterministic (std ≈ 0). Everything is derived from the recorded result columns, except
# the capacity diff, which reads the saved reduced_invested_units.arrow dumps.
#
# Usage: julia --project=BlendedClustering analysis/summarize.jl [result.csv ...]
# Default: outputs/tyndp/gep.csv, outputs/tyndp/p2x.csv, outputs/sienna/5bus.csv,
#          outputs/sienna/118bus.csv, outputs/gridmod/rts.csv (whichever exist).

using CSV, DataFrames, Printf, Statistics
using Arrow

# ---- arm classification ---------------------------------------------------------------
function arm_label(ct, wt, nrm, has_chain, is_storage)
    trad = Dict("k_means" => "k-means", "k_medoids" => "k-medoids",
        "hierarchical" => "hierarchical", "chronological" => "chronological")
    if wt == "dirac" && haskey(trad, ct) && nrm == "unscaled" && !has_chain
        return "$(trad[ct]) (traditional)"
    elseif ct == "conical_hull" && wt == "convex" && nrm == "economic" && has_chain == is_storage
        return "PROPOSED"
    elseif ct == "conical_hull" && wt == "convex" && nrm == "unscaled" && has_chain == is_storage
        return "ablate: -economic"
    elseif ct == "convex_hull" && wt == "convex" && nrm == "economic" && has_chain == is_storage
        return "ablate: -conic-selection"
    elseif ct == "conical_hull" && wt == "dirac" && nrm == "economic" && has_chain == is_storage
        return "ablate: -convex-weights"
    elseif is_storage && ct == "conical_hull" && wt == "convex" && nrm == "economic" && !has_chain
        return "ablate: -chain-split"
    end
    return "other: $ct/$wt/$nrm" * (has_chain ? "/chain" : "")
end

COMPARISON = ["k-means (traditional)", "k-medoids (traditional)", "hierarchical (traditional)",
    "chronological (traditional)", "PROPOSED"]
ABLATION = ["PROPOSED", "ablate: -economic", "ablate: -conic-selection",
    "ablate: -convex-weights", "ablate: -chain-split"]
METHOD_ORDER = ["k_means", "k_medoids", "hierarchical", "chronological", "convex_hull", "conical_hull"]

# mean ± std of `col` over the rows in `sub`, formatted as a percentage cell.
function cell(sub, col)
    vals = collect(skipmissing(sub[!, col]))
    isempty(vals) && return "-"
    m = mean(vals)
    length(vals) > 1 && std(vals) > 1e-9 ? @sprintf("%.1f±%.1f", m, std(vals)) : @sprintf("%.1f", m)
end

# label × n_rp table for `col`, rows in `order`. Skipped unless ≥2 of the family's arms
# are present (so a comparison-only file shows no lonely ablation table, and vice versa).
function grid_table(title, df, order, col; pct=true)
    count(l -> any(df.label .== l), order) < 2 && return
    println("\n", title)
    grid = sort(unique(df.n_rep_periods))
    @printf("  %-30s", "arm \\ n_rp")
    for n in grid; @printf("%12d", n); end
    println()
    for label in order
        rows = df[df.label .== label, :]
        isempty(rows) && continue
        @printf("  %-30s", label)
        for n in grid
            r = rows[rows.n_rep_periods .== n, :]
            s = isempty(r) ? "-" : cell(r, col)
            @printf("%11s%s", s, pct && s != "-" ? "%" : (pct ? " " : " "))
        end
        println()
    end
end

# ‖arm − reference‖₁ / ‖reference‖₁ of invested_units, read from the arrow dumps.
function capacity_l1_diff(csv_path, arm_name, ref_name)
    dir = dirname(csv_path)
    af = joinpath(dir, basename(arm_name), "reduced_invested_units.arrow")
    rf = joinpath(dir, basename(ref_name), "reduced_invested_units.arrow")
    (isfile(af) && isfile(rf)) || return missing
    a = combine(groupby(DataFrame(Arrow.Table(af)), :id), :value => mean => :a)
    r = combine(groupby(DataFrame(Arrow.Table(rf)), :id), :value => mean => :r)
    j = outerjoin(a, r; on=:id)
    j.a = coalesce.(j.a, 0.0); j.r = coalesce.(j.r, 0.0)
    tot = sum(j.r)
    tot <= 0 && return missing
    return 100 * sum(abs, j.a .- j.r) / tot
end

function summarize_file(path)
    println("\n", "="^80, "\nRESULTS: ", path, "\n", "="^80)
    df = CSV.read(path, DataFrame)
    for c in ("clustering_type", "weight_type", "normalization")
        df[!, c] = string.(df[!, c])
    end
    refrows = df[df.n_rep_periods .== 1, :]
    (isempty(refrows) || all(ismissing, refrows.objective_value)) &&
        (println("  !! no n_rep=1 reference optimum — skipping"); return)
    ref_opt = mean(skipmissing(refrows.objective_value))
    ref_name = first(refrows.name)
    @printf("  reference optimum (n_rep=1): %.6g\n", ref_opt)

    arms = df[df.n_rep_periods .!= 1, :]
    is_storage = any(skipmissing(arms.fix_every) .> 1)
    arms.has_chain = occursin.("chain", arms.name)
    arms.label = [arm_label(r.clustering_type, r.weight_type, r.normalization, r.has_chain, is_storage)
                  for r in eachrow(arms)]
    arms.regret_pct = map(eachrow(arms)) do r
        (ismissing(r.evaluated_objective_value) || r.evaluated_objective_value <= 0) && return missing
        100 * (r.evaluated_objective_value / ref_opt - 1)
    end
    getcol(r, c) = hasproperty(r, c) ? r[!, c] : fill(missing, nrow(r))

    # Main sweep = the default tolerance (0.01) and cached; the sensitivity rows carry other tols.
    main = arms[isapprox.(arms.tol, 0.01; atol=1e-12) .& .!occursin.("nocache", arms.name), :]

    # 1–2: comparison & ablation (regret, mean ± std over seeds).
    grid_table("1. COMPARISON — regret (%) vs n_rp (mean ± std over seeds)", main, COMPARISON, :regret_pct)
    grid_table("2. ABLATION — regret (%) vs n_rp (each row knocks out one component)", main, ABLATION, :regret_pct)

    # 3: secondary metrics at a representative n_rp (the largest in the main grid).
    if !isempty(main)
        nfocus = maximum(main.n_rep_periods)
        sub = main[main.n_rep_periods .== nfocus, :]
        println("\n3. SECONDARY METRICS at n_rp=$nfocus  (curtailment, infeasibility, feasibility, true cost split)")
        @printf("  %-30s %12s %12s %10s %12s %12s\n", "arm", "curtail", "borrow", "maxΣW", "ops-cost", "ens/borrow-cost")
        for label in COMPARISON
            rows = sub[sub.label .== label, :]
            isempty(rows) && continue
            spill = mean(skipmissing(getcol(rows, :total_spillage)))
            borrow = mean(skipmissing(getcol(rows, :total_borrow)))
            maxw = mean(skipmissing(getcol(rows, :max_weight_sum)))
            opsc = mean(skipmissing(getcol(rows, :eval_cost_of_operations)))
            ensc = mean(skipmissing(getcol(rows, :eval_cost_of_borrow)))
            f(x) = ismissing(x) || isnan(x) ? "-" : @sprintf("%.4g", x)
            @printf("  %-30s %12s %12s %10s %12s %12s\n", label, f(spill), f(borrow), f(maxw), f(opsc), f(ensc))
        end
    end

    # 4: capacity-decision differences (investment datasets only).
    if !is_storage && !isempty(main)
        println("\n4. CAPACITY DECISIONS — ‖Δ invested_units‖₁ vs the full-horizon optimum (%)")
        grid = sort(unique(main.n_rep_periods))
        for label in COMPARISON
            names = unique(main[main.label .== label, :name])
            isempty(names) && continue
            @printf("  %-30s", label)
            for n in grid
                rr = main[(main.label .== label) .& (main.n_rep_periods .== n), :]
                d = isempty(rr) ? missing : capacity_l1_diff(path, first(rr.name), ref_name)
                @printf("%11s%s", ismissing(d) ? "-" : @sprintf("%.1f", d), ismissing(d) ? " " : "%")
            end
            println()
        end
    end

    # 5: runtime breakdown + N(ε) + cache hit-rate, for the PROPOSED arm across n_rp.
    prop = main[main.label .== "PROPOSED", :]
    if !isempty(prop)
        println("\n5. RUNTIME (PROPOSED) vs n_rp — seconds per stage, PGD iters, cache hit-rate")
        @printf("  %-8s %10s %12s %10s %10s %10s\n", "n_rp", "cluster", "weight-fit", "solve", "PGDiters", "cacheHit%")
        for n in sort(unique(prop.n_rep_periods))
            r = prop[prop.n_rep_periods .== n, :]
            tc = mean(skipmissing(getcol(r, :time_to_cluster)))
            tf = mean(skipmissing(getcol(r, :time_to_fit_weights)))
            ts = mean(skipmissing(getcol(r, :time_to_solve)))
            pit = mean(skipmissing(getcol(r, :pgd_total_iters)))
            ch = sum(skipmissing(getcol(r, :cache_hits))); cm = sum(skipmissing(getcol(r, :cache_misses)))
            hit = (ch + cm) > 0 ? 100 * ch / (ch + cm) : NaN
            @printf("  %-8d %10.2f %12.2f %10.2f %10.0f %9.1f%%\n", n, tc, tf, ts, pit, hit)
        end
    end

    # 6: sensitivity (development system) — EVERY clustering type × normalization × tol at
    # n_rp=10 (convex weights, chain off). Shows which selection is most hyperparameter-robust,
    # so the method is not presupposed. Chain is off here so k-means (no base-period anchors)
    # is comparable to the real-period methods.
    sens = arms[(arms.n_rep_periods .== 10) .& (arms.weight_type .== "convex") .&
                (arms.has_chain .== false) .& .!occursin.("nocache", arms.name), :]
    if nrow(sens) > 1 && length(unique(sens.clustering_type)) > 1 && length(unique(sens.tol)) > 1
        tols = sort(unique(sens.tol); rev=true)
        methods = [m for m in METHOD_ORDER if m in sens.clustering_type]
        println("\n6. SENSITIVITY (development system) — regret (%) vs tol, per clustering type, n_rp=10")
        for nrm in sort(unique(sens.normalization))
            println("  --- normalization: $nrm ---")
            @printf("    %-16s", "method \\ tol")
            for t in tols; @printf("%12s", t); end
            println()
            for m in methods
                @printf("    %-16s", m)
                for t in tols
                    r = sens[(sens.clustering_type .== m) .& (sens.normalization .== nrm) .&
                             isapprox.(sens.tol, t; atol=1e-12), :]
                    s = (isempty(r) || all(ismissing, r.regret_pct)) ? "-" :
                        @sprintf("%.1f", mean(skipmissing(r.regret_pct)))
                    @printf("%11s%s", s, s == "-" ? " " : "%")
                end
                println()
            end
        end
        # robustness: regret spread across every (tol × normalization) cell, per method.
        println("  robustness — regret spread over tol × normalization (smaller range ⇒ more robust)")
        @printf("    %-16s %10s %10s %10s\n", "method", "min", "max", "range")
        for m in methods
            v = collect(skipmissing(sens[sens.clustering_type .== m, :regret_pct]))
            isempty(v) && continue
            @printf("    %-16s %9.1f%% %9.1f%% %9.1f%%\n", m, minimum(v), maximum(v), maximum(v) - minimum(v))
        end
        # realized N(ε): median PGD iterations per tol (grows as ε shrinks; α, N are fixed).
        println("  realized PGD iterations N(ε) — median over methods × normalizations")
        @printf("    %-16s", "tol")
        for t in tols; @printf("%12g", t); end
        println()
        @printf("    %-16s", "median N(ε)")
        for t in tols
            it = collect(skipmissing(sens[isapprox.(sens.tol, t; atol=1e-12), :pgd_total_iters]))
            @printf("%12s", isempty(it) ? "-" : string(round(Int, median(it))))
        end
        println()
        # cache validation: uncached hull twins — identical regret, cluster-time comparison.
        nocache = arms[occursin.("nocache", arms.name), :]
        if !isempty(nocache)
            println("  cache validation (hull methods): cluster time cached vs uncached (regret identical)")
            for r in eachrow(nocache)
                cached = main[main.name .== replace(r.name, "_nocache" => ""), :]
                isempty(cached) && continue
                cr = first(eachrow(cached))
                @printf("    %-14s tol=%g : %.3fs (cache) vs %.3fs (no cache); regret %.2f%% vs %.2f%%\n",
                    r.clustering_type, r.tol, coalesce(cr.time_to_cluster, NaN),
                    coalesce(r.time_to_cluster, NaN), coalesce(cr.regret_pct, NaN), coalesce(r.regret_pct, NaN))
            end
        end
    end
    return
end

function main()
    files = isempty(ARGS) ? ["outputs/tyndp/gep.csv", "outputs/tyndp/p2x.csv",
        "outputs/sienna/5bus.csv", "outputs/sienna/118bus.csv", "outputs/gridmod/rts.csv"] : ARGS
    found = false
    for f in files
        if isfile(f); found = true; summarize_file(f)
        else; println("(skipping $f — not found; run its case study first)"); end
    end
    found || println("\nNo result CSVs. Run e.g. run_case_studies(\"configs/5bus.toml\").")
end

main()
