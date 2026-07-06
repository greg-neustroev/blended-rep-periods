#!/usr/bin/env julia
#
# Summarize the revision-campaign results. Two kinds of result file:
#   * development (5bus): sensitivity (clustering × weight × normalization × tol) + ablation.
#   * held-out (GEP/P2X/118bus/RTS): comparison over the clustering × weight matrix.
# The file is classified by whether it carries multiple tolerances at n_rp=10 (sensitivity).
#
# Sections printed as applicable:
#   COMPARISON  — regret (%) vs n_rp, one row per clustering × weight combo (mean ± std/seeds).
#   ABLATION    — PROPOSED vs each single-component knockout, regret (%) vs n_rp.
#   SECONDARY   — curtailment, infeasibility (borrow), feasibility (max ΣW), cost split.
#   CAPACITY    — (investment) ‖Δ invested_units‖₁ vs the full-horizon optimum.
#   RUNTIME     — per-stage timing + realized PGD iterations N(ε) + cache hit-rate (PROPOSED).
#   SENSITIVITY — per clustering × weight, regret spread over tol × normalization (robustness).
#
# Regret = evaluated_objective_value / reference_optimum − 1 (reference = the n_rep=1 row).
# Storage is graded day-exact (fix_every=1); the convex chain reconstructs every boundary, so
# every real-period method carries it (k-means, no anchors, is the no-chain baseline).
#
# Usage: julia --project=BlendedClustering analysis/summarize.jl [result.csv ...]

using CSV, DataFrames, Printf, Statistics
using Arrow

const METHOD_ORDER = ["k_means", "k_medoids", "hierarchical", "chronological", "convex_hull", "conical_hull"]
const WEIGHT_ORDER = ["dirac", "convex", "conical", "conical_bounded"]
const PROPOSED = ("conical_hull", "convex")   # marked with * in the tables

# mean ± std of `col` over rows, as a fixed-width percentage cell.
function cell(sub, col)
    vals = collect(skipmissing(sub[!, col]))
    isempty(vals) && return "-"
    m = mean(vals)
    length(vals) > 1 && std(vals) > 1e-9 ? @sprintf("%.1f±%.1f", m, std(vals)) : @sprintf("%.1f", m)
end
getcol(r, c) = hasproperty(r, c) ? r[!, c] : fill(missing, nrow(r))

# regret (%) matrix: rows = clustering × weight combos present, cols = n_rp.
function comparison_matrix(df)
    grid = sort(unique(df.n_rep_periods))
    combos = [(cl, w) for cl in METHOD_ORDER for w in WEIGHT_ORDER
              if any((df.clustering_type .== cl) .& (df.weight_type .== w))]
    isempty(combos) && return
    println("\nCOMPARISON — regret (%) vs n_rp, per clustering × weight (mean ± std over seeds; * = PROPOSED)")
    @printf("  %-26s", "clustering / weight")
    for n in grid; @printf("%12d", n); end
    println()
    for (cl, w) in combos
        star = (cl, w) == PROPOSED ? " *" : ""
        @printf("  %-26s", "$cl / $w$star")
        for n in grid
            r = df[(df.clustering_type .== cl) .& (df.weight_type .== w) .& (df.n_rep_periods .== n), :]
            s = isempty(r) ? "-" : cell(r, :regret_pct)
            @printf("%11s%s", s, s == "-" ? " " : "%")
        end
        println()
    end
end

# Ablation arm labels (5bus): PROPOSED with one component knocked out.
function ablation_label(ct, wt, nrm, has_chain)
    if ct == "conical_hull" && wt == "convex" && nrm == "economic" && has_chain; return "PROPOSED"
    elseif ct == "conical_hull" && wt == "convex" && nrm == "unscaled" && has_chain; return "ablate: -economic"
    elseif ct == "convex_hull" && wt == "convex" && nrm == "economic" && has_chain; return "ablate: -conic-selection"
    elseif ct == "conical_hull" && wt == "dirac" && nrm == "economic" && has_chain; return "ablate: -convex-weights"
    elseif ct == "conical_hull" && wt == "convex" && nrm == "economic" && !has_chain; return "ablate: -chain-split"
    end
    return missing
end

function ablation_table(df)
    order = ["PROPOSED", "ablate: -economic", "ablate: -conic-selection", "ablate: -convex-weights", "ablate: -chain-split"]
    sub = df[.!ismissing.(df.ablabel), :]
    count(l -> any(sub.ablabel .=== l), order) < 2 && return
    grid = sort(unique(sub.n_rep_periods))
    println("\nABLATION (5bus) — regret (%) vs n_rp (each row knocks out one component of PROPOSED)")
    @printf("  %-28s", "arm \\ n_rp")
    for n in grid; @printf("%12d", n); end
    println()
    for label in order
        rows = sub[sub.ablabel .=== label, :]
        isempty(rows) && continue
        @printf("  %-28s", label)
        for n in grid
            r = rows[rows.n_rep_periods .== n, :]
            s = isempty(r) ? "-" : cell(r, :regret_pct)
            @printf("%11s%s", s, s == "-" ? " " : "%")
        end
        println()
    end
end

# ‖arm − reference‖₁ / ‖reference‖₁ of invested_units, from the arrow dumps.
function capacity_l1_diff(csv_path, arm_name, ref_name)
    dir = dirname(csv_path)
    af = joinpath(dir, basename(arm_name), "reduced_invested_units.arrow")
    rf = joinpath(dir, basename(ref_name), "reduced_invested_units.arrow")
    (isfile(af) && isfile(rf)) || return missing
    a = combine(groupby(DataFrame(Arrow.Table(af)), :id), :value => mean => :a)
    r = combine(groupby(DataFrame(Arrow.Table(rf)), :id), :value => mean => :r)
    j = outerjoin(a, r; on=:id); j.a = coalesce.(j.a, 0.0); j.r = coalesce.(j.r, 0.0)
    tot = sum(j.r); tot <= 0 && return missing
    return 100 * sum(abs, j.a .- j.r) / tot
end

# per clustering × weight, regret spread over the tol × normalization grid (robustness).
function sensitivity_tables(df)
    # keep the per-method-default chain (convex for real-period, none for k-means); this drops
    # the -chain ablation (conical_hull with chain off) so it does not pollute the grouping.
    sens = df[(df.n_rep_periods .== 10) .& .!occursin.("nocache", df.name) .&
              (df.has_chain .== (df.clustering_type .!= "k_means")), :]
    (nrow(sens) < 2 || length(unique(sens.tol)) < 2) && return
    combos = [(cl, w) for cl in METHOD_ORDER for w in WEIGHT_ORDER
              if any((sens.clustering_type .== cl) .& (sens.weight_type .== w))]
    # one table per normalization so a single fragile normalization (e.g. minmax) does not
    # inflate the range of an otherwise-robust method; here range = spread over tol alone.
    norm_order = ["economic", "unscaled", "minmax"]
    norms = [nrm for nrm in norm_order if any(sens.normalization .== nrm)]
    append!(norms, sort(unique(sens.normalization[.!in.(sens.normalization, Ref(norm_order))])))
    println("\nSENSITIVITY (5bus dev) — regret (%) spread over tol, n_rp=10, per normalization (smaller range ⇒ more robust; * = PROPOSED)")
    for nrm in norms
        ns = sens[sens.normalization .== nrm, :]
        length(unique(ns.tol)) < 2 && continue
        println("  normalization = $nrm")
        @printf("  %-26s %10s %10s %10s %10s\n", "clustering / weight", "min", "median", "max", "range")
        for (cl, w) in combos
            v = collect(skipmissing(ns[(ns.clustering_type .== cl) .& (ns.weight_type .== w), :regret_pct]))
            isempty(v) && continue
            star = (cl, w) == PROPOSED ? " *" : ""
            @printf("  %-26s %9.1f%% %9.1f%% %9.1f%% %9.1f%%\n", "$cl / $w$star",
                minimum(v), median(v), maximum(v), maximum(v) - minimum(v))
        end
    end
    # realized N(ε): median PGD iterations per tol (grows as ε shrinks; α, N fixed).
    tols = sort(unique(sens.tol); rev=true)
    println("  realized PGD iterations N(ε) — median over methods × normalizations")
    @printf("    %-14s", "tol"); for t in tols; @printf("%12g", t); end; println()
    @printf("    %-14s", "median N(ε)")
    for t in tols
        it = collect(skipmissing(sens[isapprox.(sens.tol, t; atol=1e-12), :pgd_total_iters]))
        @printf("%12s", isempty(it) ? "-" : string(round(Int, median(it))))
    end
    println()
    # cache validation: uncached hull twins — identical regret, cluster-time comparison.
    nocache = df[occursin.("nocache", df.name), :]
    if !isempty(nocache)
        println("  cache validation (hull): cluster time cached vs uncached (regret identical)")
        for r in eachrow(nocache)
            cached = df[df.name .== replace(r.name, "_nocache" => ""), :]
            isempty(cached) && continue
            cr = first(eachrow(cached))
            @printf("    %-14s tol=%g : %.3fs (cache) vs %.3fs (no cache); regret %.2f%% vs %.2f%%\n",
                r.clustering_type, r.tol, coalesce(cr.time_to_cluster, NaN),
                coalesce(r.time_to_cluster, NaN), coalesce(cr.regret_pct, NaN), coalesce(r.regret_pct, NaN))
        end
    end
end

function summarize_file(path)
    println("\n", "="^84, "\nRESULTS: ", path, "\n", "="^84)
    df = CSV.read(path, DataFrame)
    for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
    refrows = df[df.n_rep_periods .== 1, :]
    (isempty(refrows) || all(ismissing, refrows.objective_value)) &&
        (println("  !! no n_rep=1 reference optimum — skipping"); return)
    ref_opt = mean(skipmissing(refrows.objective_value)); ref_name = first(refrows.name)
    @printf("  reference optimum (n_rep=1): %.6g\n", ref_opt)

    arms = df[df.n_rep_periods .!= 1, :]
    arms.has_chain = occursin.("chain", arms.name)
    arms.regret_pct = map(eachrow(arms)) do r
        (ismissing(r.evaluated_objective_value) || r.evaluated_objective_value <= 0) && return missing
        100 * (r.evaluated_objective_value / ref_opt - 1)
    end
    arms.ablabel = [ablation_label(r.clustering_type, r.weight_type, r.normalization, r.has_chain) for r in eachrow(arms)]
    is_investment = !any(skipmissing(getcol(arms, :fix_every)) .>= 1) || all(ismissing, getcol(arms, :total_borrow))
    # development file iff it carries several tolerances at n_rp=10 (the sensitivity sweep).
    at10 = arms[arms.n_rep_periods .== 10, :]
    is_dev = nrow(at10) > 0 && length(unique(at10.tol)) > 1

    if is_dev
        ablation_table(arms)
        sensitivity_tables(arms)
    else
        comparison_matrix(arms)
    end

    # SECONDARY — at the largest n_rp, PROPOSED vs the k-means / k-medoids dirac baselines.
    grid = sort(unique(arms.n_rep_periods)); nfocus = maximum(grid)
    show_combos = [PROPOSED, ("k_means", "dirac"), ("k_medoids", "dirac")]
    # pin to the canonical config so ablation twins (unscaled, chain-off, other tol) that share a
    # clustering×weight do not pool into these metrics; a no-op on comparison files (economic/0.01 only).
    storage_graded = any(skipmissing(getcol(arms, :fix_every)) .>= 1)
    chain_ok = storage_graded ? (arms.has_chain .== (arms.clustering_type .!= "k_means")) : trues(nrow(arms))
    canon = arms[chain_ok .& .!occursin.("nocache", arms.name) .&
                 (arms.normalization .== "economic") .& isapprox.(arms.tol, 0.01; atol=1e-12), :]
    sub = canon[canon.n_rep_periods .== nfocus, :]
    haveany = any(any((sub.clustering_type .== cl) .& (sub.weight_type .== w)) for (cl, w) in show_combos)
    if haveany
        println("\nSECONDARY METRICS at n_rp=$nfocus  (curtailment, borrow-infeasibility, max ΣW, true cost split)")
        @printf("  %-26s %11s %11s %8s %11s %11s\n", "clustering / weight", "curtail", "borrow", "maxΣW", "ops-cost", "ens/borrow")
        for (cl, w) in show_combos
            r = sub[(sub.clustering_type .== cl) .& (sub.weight_type .== w), :]
            isempty(r) && continue
            f(c) = (v = collect(skipmissing(getcol(r, c))); isempty(v) ? "-" : @sprintf("%.4g", mean(v)))
            @printf("  %-26s %11s %11s %8s %11s %11s\n", "$cl / $w",
                f(:total_spillage), f(:total_borrow), f(:max_weight_sum), f(:eval_cost_of_operations), f(:eval_cost_of_borrow))
        end
    end

    # CAPACITY (investment only) — ‖Δ invested_units‖₁ vs the optimum, PROPOSED + baselines.
    if is_investment
        println("\nCAPACITY DECISIONS — ‖Δ invested_units‖₁ vs the full-horizon optimum (%), at n_rp=$nfocus")
        for (cl, w) in show_combos
            r = sub[(sub.clustering_type .== cl) .& (sub.weight_type .== w), :]
            isempty(r) && continue
            d = capacity_l1_diff(path, first(r.name), ref_name)
            @printf("  %-26s %s\n", "$cl / $w", ismissing(d) ? "-" : @sprintf("%.1f%%", d))
        end
    end

    # RUNTIME (PROPOSED) vs n_rp — per stage + PGD iterations + cache hit-rate.
    prop = arms[(arms.clustering_type .== PROPOSED[1]) .& (arms.weight_type .== PROPOSED[2]) .&
                (arms.normalization .== "economic") .& (arms.has_chain .| .!any(skipmissing(getcol(arms,:fix_every)).>=1)) .&
                isapprox.(arms.tol, 0.01; atol=1e-12), :]
    if !isempty(prop)
        println("\nRUNTIME (PROPOSED) vs n_rp — seconds per stage, PGD iters, cache hit-rate")
        @printf("  %-8s %10s %12s %10s %10s %10s\n", "n_rp", "cluster", "weight-fit", "solve", "PGDiters", "cacheHit%")
        for n in sort(unique(prop.n_rep_periods))
            r = prop[prop.n_rep_periods .== n, :]
            tc = mean(skipmissing(getcol(r, :time_to_cluster))); tf = mean(skipmissing(getcol(r, :time_to_fit_weights)))
            ts = mean(skipmissing(getcol(r, :time_to_solve))); pit = mean(skipmissing(getcol(r, :pgd_total_iters)))
            ch = sum(skipmissing(getcol(r, :cache_hits))); cm = sum(skipmissing(getcol(r, :cache_misses)))
            hit = (ch + cm) > 0 ? 100 * ch / (ch + cm) : NaN
            @printf("  %-8d %10.2f %12.2f %10.2f %10.0f %9.1f%%\n", n, tc, tf, ts, pit, hit)
        end
    end
    return
end

function main()
    files = isempty(ARGS) ? ["outputs/tyndp/gep.csv", "outputs/tyndp/p2x.csv",
        "outputs/sienna/118bus.csv", "outputs/sienna/5bus.csv", "outputs/gridmod/rts.csv"] : ARGS
    found = false
    for f in files
        if isfile(f); found = true; summarize_file(f)
        else; println("(skipping $f — not found; run its case study first)"); end
    end
    found || println("\nNo result CSVs. Run e.g. run_case_studies(\"configs/5bus.toml\").")
end

main()
