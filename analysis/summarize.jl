#!/usr/bin/env julia
#
# Summarize the revision-campaign results. Two kinds of result file:
#   * development (5bus): sensitivity (clustering × weight × normalization × tol) + ablation.
#   * held-out (GEP/P2X/118bus/RTS): comparison over the clustering × weight matrix.
# The file is classified by whether it carries multiple tolerances at n_rp=10 (sensitivity).
#
# Sections printed as applicable:
#   COMPARISON  — regret (%) + total time (s), clustering (rows) × weight (cols), one table per
#                 (normalization, n_rp) group (mean ± std over seeds). Within each group the
#                 lowest-mean-regret cell is bold, along with every cell that is reasonable to use —
#                 i.e. NOT both meaningfully worse (mean regret beyond the equivalence margin) AND
#                 statistically worse (paired t-test p<ALPHA). The margin keeps reliable near-ties;
#                 the t-test drops noisy methods whose big mean gap survives their spread.
#                 CYAN marks cells on the global Pareto frontier (non-dominated across clustering ×
#                 weight × n_rp × normalization): no other combo is no-slower-on-average (within a
#                 small time margin) with either lower mean regret, or comparable mean regret and
#                 lower regret variance.
#   NORMALIZATION — economic vs unscaled, paired cell-by-cell: how many combos differ meaningfully,
#                 plus a pooled two-sided paired t-test on Δregret. Prints only when both are present.
#   ABLATION    — PROPOSED vs each single-component knockout, regret (%) + total time (s) vs n_rp (bold = reasonable per n_rp).
#   SECONDARY   — curtailment, infeasibility (borrow), feasibility (max ΣW), cost split.
#   CAPACITY    — (investment) ‖Δ invested_units‖₁ vs the full-horizon optimum.
#   RUNTIME     — per-stage timing + realized PGD iterations N(ε) + cache hit-rate (PROPOSED).
#   SENSITIVITY — per clustering × weight, regret spread over tol × normalization (robustness).
#
# Regret = evaluated_objective_value / reference_optimum − 1 (reference = the n_rep=1 row).
# Storage is graded day-exact; the inter-period storage chain (using the operational weights,
# with exact-inflow injection) reconstructs every boundary for every method.
#
# Usage: julia --project=BlendedClustering analysis/summarize.jl [result.csv ...]

isdefined(Main, :METHOD_ORDER) || include(joinpath(@__DIR__, "common.jl"))

# Fixed field widths so wide cells never break alignment. Each n_rp splits into a regret and a
# time sub-column; the label holds the longest combo "convex_hull_0 / subunit" (23 chars).
const LABELW = 28
const REGW = 14    # regret sub-column; widest cell "441.1±138.5%"
const TIMEW = 13   # time sub-column;   widest cell "146.7±20.0s"

# regret / total-time cells, each "mean ± std over seeds"; "-" when there is no matching row.
regret_cell(sub) = ((m, s) = meanstd(sub, :regret_pct); m === nothing ? "-" : @sprintf("%.1f%s%%", m, s))
time_cell(sub)   = ((m, s) = meanstd(sub, :total_time); m === nothing ? "-" : @sprintf("%.1f%ss", m, s))

center(s, w) = (p = max(0, w - length(s)); l = p ÷ 2; " "^l * s * " "^(p - l))

# header for the regret+time matrices: each column label spans a (regret | time) pair.
function matrix_header(cols, corner)
    @printf("  %-*s", LABELW, corner)
    for c in cols; print(center(c, REGW + TIMEW)); end
    @printf("\n  %-*s", LABELW, "")
    for _ in cols; @printf("%*s%*s", REGW, "regret", TIMEW, "time"); end
    println()
end

# one (regret | time) pair for the rows matching `r`, or a "-" when empty. ANSI styling is applied
# AFTER width-padding so the escape codes never disturb alignment: `bold` = reasonable in group,
# `cyan` = non-dominated on the global Pareto frontier (the two compose).
function matrix_pair(r; bold=false, cyan=false)
    s = @sprintf("%*s%*s", REGW, isempty(r) ? "-" : regret_cell(r), TIMEW, isempty(r) ? "" : time_cell(r))
    cyan && (s = "\e[36m$s\e[39m")
    bold && (s = "\e[1m$s\e[22m")
    print(s)
end

# regret (%) matrix: rows = clustering, cols = weight, one table per (normalization, n_rp) group.
# Grouping by normalization keeps economic / unscaled / minmax out of the same mean ± std cell;
# grouping by n_rp lets each table's "reasonable" set be judged among directly comparable runs.
function comparison_matrix(df)
    grid = sort(unique(df.n_rep_periods))
    isempty(grid) && return
    weights = [w for w in WEIGHT_ORDER if any(df.weight_type .== w)]
    nondom = nondominated_keys(df)   # Pareto frontier across ALL (normalization, n_rp, clustering, weight)
    println("\nCOMPARISON — regret (%) and total time (s): clustering (rows) × weight (cols), grouped by normalization × n_rp (mean ± std over seeds; BOLD = reasonable in that group: not both >margin worse in regret AND paired-t-test worse, p<$ALPHA; CYAN = non-dominated across clustering×weight×n_rp×normalization: no other combo is no-slower-on-average with lower mean regret, or with comparable mean regret and lower regret variance)")
    for nrm in ordered_norms(df)
        nd = df[df.normalization .== nrm, :]
        clusts = [cl for cl in METHOD_ORDER if any(nd.clustering_type .== cl)]
        isempty(clusts) && continue
        for n in grid
            g = nd[nd.n_rep_periods .== n, :]
            isempty(g) && continue
            cellrows(cl, w) = g[(g.clustering_type .== cl) .& (g.weight_type .== w), :]
            # "reasonable" set judged over every clustering × weight cell in this group.
            keys_ = [(cl, w) for cl in clusts for w in weights]
            tied = Set(keys_[i] for i in tied_with_best([samples(cellrows(k...)) for k in keys_]))
            println("\n  normalization = $nrm,  n_rp = $n")
            matrix_header([disp(w) for w in weights], "clustering \\ weight")
            for cl in clusts
                @printf("  %-*s", LABELW, disp(cl))
                for w in weights
                    matrix_pair(cellrows(cl, w); bold = (cl, w) in tied,
                                cyan = (nrm, n, cl, w) in nondom)
                end
                println()
            end
        end
    end
end

# NORMALIZATION — is economic meaningfully different from unscaled? If no cell differs
# meaningfully, the choice is immaterial (drop one). Thin printing wrapper around
# `normalization_compare_pair`; output text unchanged from before the refactor.
function normalization_compare(arms)
    r = normalization_compare_pair(arms, "economic", "unscaled")
    (r === nothing || r.n_total == 0) && return
    println("\nNORMALIZATION — economic vs unscaled (paired over seeds, canonical tol; same margin+t-test rule as the tables): does the choice matter?")
    @printf("  %d combos (clustering × weight × n_rp): comparable %d, economic better %d, unscaled better %d\n",
            r.n_total, r.n_comparable, r.n_a_better, r.n_b_better)
    @printf("  pooled Δregret (economic − unscaled): mean %+.3f pp, largest |mean| %.3f pp (%s); two-sided paired t-test p=%.2g\n",
            r.mean_diff, abs(r.worst_diff), isempty(r.worst_combo) ? "-" : r.worst_combo, r.p_value)
    if r.n_a_better + r.n_b_better == 0
        println("  ⇒ economic and unscaled are statistically indistinguishable here — the normalization choice does not matter.")
    else
        println("  ⇒ they differ meaningfully in $(r.n_a_better + r.n_b_better)/$(r.n_total) combos — the normalization choice matters.")
    end
end

# The "few RPs beat many" story the grouped tables hide: the fewest representative periods at
# which the best blended-hull method already matches or beats the best CONVENTIONAL method run at
# the largest n_rp. Uses the canonical economic runs; compares best-in-family so the claim is
# method-agnostic on both sides (strongest baseline vs strongest hull).
function few_rp_headline(arms)
    df = arms[.!occursin.("nocache", arms.name) .& (arms.normalization .== "economic") .&
              isapprox.(arms.tol, 0.01; atol=1e-12), :]
    grid = sort(unique(df.n_rep_periods))
    length(grid) < 2 && return
    nmax = maximum(grid)
    meanof(sub, c) = (v = collect(skipmissing(getcol(sub, c))); isempty(v) ? nothing : mean(v))
    # lowest-mean-regret (clustering, weight) within `family` at a given n_rp.
    function best(family, n)
        b = nothing
        for cl in family, w in WEIGHT_ORDER
            sub = df[(df.clustering_type .== cl) .& (df.weight_type .== w) .& (df.n_rep_periods .== n), :]
            r = meanof(sub, :regret_pct); r === nothing && continue
            (b === nothing || r < b.r) && (b = (r = r, t = meanof(sub, :total_time), name = combo_label(cl, w)))
        end
        b
    end
    bar = best(CONVENTIONAL, nmax); bar === nothing && return
    for n in grid
        n >= nmax && break
        h = best(HULL, n); h === nothing && continue
        if h.r <= bar.r
            println("\nHEADLINE — few RPs beat many (economic): a blended-hull method needs far fewer periods than the best conventional baseline for the same regret")
            @printf("  best blended-hull @ n_rp=%-4d : %-24s regret %5.1f%%, %6.1fs\n", n, h.name, h.r, something(h.t, NaN))
            @printf("  best conventional @ n_rp=%-4d : %-24s regret %5.1f%%, %6.1fs\n", nmax, bar.name, bar.r, something(bar.t, NaN))
            if h.t !== nothing && bar.t !== nothing && h.t > 0
                @printf("  ⇒ %.0f× fewer representative periods and %.0f× faster, at equal-or-better regret.\n", nmax / n, bar.t / h.t)
            else
                @printf("  ⇒ %.0f× fewer representative periods, at equal-or-better regret.\n", nmax / n)
            end
            return
        end
    end
end

function ablation_table(df)
    order = ["PROPOSED", "ablate: -economic", "ablate: -conic-selection", "ablate: -convex-weights"]
    sub = df[.!ismissing.(df.ablabel), :]
    count(l -> any(sub.ablabel .=== l), order) < 2 && return
    grid = sort(unique(sub.n_rep_periods))
    println("\nABLATION (5bus) — regret (%) and total time (s) vs n_rp (each row knocks out one component of PROPOSED; BOLD = reasonable at that n_rp: not both >margin worse AND paired-t-test worse, p<$ALPHA)")
    matrix_header(["n_rp=$n" for n in grid], "arm \\ n_rp")
    labels = [l for l in order if any(sub.ablabel .=== l)]
    cellrows(l, n) = (rows = sub[sub.ablabel .=== l, :]; rows[rows.n_rep_periods .== n, :])
    # per n_rp column: the arms whose regret is tied with the column's best.
    tied = Dict(n => Set(labels[i] for i in tied_with_best([samples(cellrows(l, n)) for l in labels]))
                for n in grid)
    for label in labels
        @printf("  %-*s", LABELW, label)
        for n in grid
            matrix_pair(cellrows(label, n); bold = label in tied[n])
        end
        println()
    end
end

# per clustering × weight, regret spread over the tol × normalization grid (robustness).
function sensitivity_tables(df)
    sens = df[(df.n_rep_periods .== 10) .& .!occursin.("nocache", df.name), :]
    (nrow(sens) < 2 || length(unique(sens.tol)) < 2) && return
    combos = [(cl, w) for cl in METHOD_ORDER for w in WEIGHT_ORDER
              if any((sens.clustering_type .== cl) .& (sens.weight_type .== w))]
    # one table per normalization so a single fragile normalization (e.g. minmax) does not
    # inflate the range of an otherwise-robust method; here range = spread over tol alone.
    println("\nSENSITIVITY (5bus dev) — regret (%) spread over tol + median total time (s), n_rp=10, per normalization (smaller range ⇒ more robust)")
    for nrm in ordered_norms(sens)
        ns = sens[sens.normalization .== nrm, :]
        length(unique(ns.tol)) < 2 && continue
        println("  normalization = $nrm")
        @printf("  %-*s %10s %10s %10s %10s %10s\n", LABELW, "clustering / weight", "min", "median", "max", "range", "med time")
        for (cl, w) in combos
            mask = (ns.clustering_type .== cl) .& (ns.weight_type .== w)
            v = collect(skipmissing(ns[mask, :regret_pct]))
            isempty(v) && continue
            t = collect(skipmissing(ns[mask, :total_time]))
            tstr = isempty(t) ? "-" : @sprintf("%.2fs", median(t))
            @printf("  %-*s %9.1f%% %9.1f%% %9.1f%% %9.1f%% %10s\n", LABELW, combo_label(cl, w),
                minimum(v), median(v), maximum(v), maximum(v) - minimum(v), tstr)
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
    arms.regret_pct = map(eachrow(arms)) do r
        (ismissing(r.evaluated_objective_value) || r.evaluated_objective_value <= 0) && return missing
        100 * (r.evaluated_objective_value / ref_opt - 1)
    end
    arms.ablabel = [ablation_label(r.clustering_type, r.weight_type, r.normalization) for r in eachrow(arms)]
    # total time = cluster + weight-fit + solve (the same stages broken out under RUNTIME).
    arms.total_time = reduce(.+, (coalesce.(getcol(arms, s), 0.0) for s in TIME_STAGES))
    # investment (GEP) has no seasonal storage, so total_borrow is entirely missing/NA there.
    is_investment = all(ismissing, getcol(arms, :total_borrow))
    # development file iff it carries several tolerances at n_rp=10 (the sensitivity sweep).
    at10 = arms[arms.n_rep_periods .== 10, :]
    is_dev = nrow(at10) > 0 && length(unique(at10.tol)) > 1

    if is_dev
        ablation_table(arms)
        sensitivity_tables(arms)
    else
        few_rp_headline(arms)
        comparison_matrix(arms)
    end
    normalization_compare(arms)   # economic vs unscaled — only prints when both are present

    # SECONDARY — at the largest n_rp, PROPOSED vs the k-means / k-medoids dirac baselines.
    grid = sort(unique(arms.n_rep_periods)); nfocus = maximum(grid)
    show_combos = [PROPOSED, ("k_means", "dirac"), ("k_medoids", "dirac")]
    # pin to the canonical config so ablation twins (unscaled, other tol) that share a
    # clustering×weight do not pool into these metrics; a no-op on comparison files (economic/0.01 only).
    canon = arms[.!occursin.("nocache", arms.name) .&
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
            @printf("  %-26s %11s %11s %8s %11s %11s\n", combo_label(cl, w),
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
            @printf("  %-26s %s\n", combo_label(cl, w), ismissing(d) ? "-" : @sprintf("%.1f%%", d))
        end
    end

    # RUNTIME (PROPOSED) vs n_rp — per stage + PGD iterations + cache hit-rate.
    prop = arms[(arms.clustering_type .== PROPOSED[1]) .& (arms.weight_type .== PROPOSED[2]) .&
                (arms.normalization .== "economic") .&
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
        "outputs/nrel/118bus.csv", "outputs/sienna/5bus.csv", "outputs/gridmod/rts.csv"] : ARGS
    found = false
    for f in files
        if isfile(f); found = true; summarize_file(f)
        else; println("(skipping $f — not found; run its case study first)"); end
    end
    found || println("\nNo result CSVs. Run e.g. run_case_studies(\"configs/5bus.toml\").")
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
