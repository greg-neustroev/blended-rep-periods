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

using CSV, DataFrames, Printf, Statistics
using Arrow
using HypothesisTests

const METHOD_ORDER = ["k_means", "k_medoids", "hierarchical", "chronological", "convex_hull", "convex_hull_with_null", "conical_hull"]
const WEIGHT_ORDER = ["dirac", "convex", "conical_bounded", "conical"]   # displayed: dirac, convex, subunit, conical
const NORM_ORDER = ["economic", "unscaled", "minmax"]
const PROPOSED = ("conical_hull", "convex")   # the method whose runtime/secondary metrics are reported
const HULL = ("convex_hull", "convex_hull_with_null", "conical_hull")          # the blended-hull family
const CONVENTIONAL = ("k_means", "k_medoids", "hierarchical", "chronological")  # the baselines it competes with

# Fixed field widths so wide cells never break alignment. Each n_rp splits into a regret and a
# time sub-column; the label holds the longest combo "convex_hull_0 / subunit" (23 chars).
const LABELW = 28
const REGW = 14    # regret sub-column; widest cell "441.1±138.5%"
const TIMEW = 13   # time sub-column;   widest cell "146.7±20.0s"
# stages that sum to the reported total time (matches the RUNTIME breakdown below).
const TIME_STAGES = (:time_to_cluster, :time_to_fit_weights, :time_to_solve)

# Shorter display names for two verbose category labels (the data keys are unchanged).
const DISPLAY = Dict("convex_hull_with_null" => "convex_hull_0", "conical_bounded" => "subunit")
disp(x) = get(DISPLAY, x, x)
combo_label(cl, w) = "$(disp(cl)) / $(disp(w))"

# normalizations present in `df`, in canonical order, with any unrecognized ones appended.
function ordered_norms(df)
    norms = [nrm for nrm in NORM_ORDER if any(df.normalization .== nrm)]
    append!(norms, sort(unique(df.normalization[.!in.(df.normalization, Ref(NORM_ORDER))])))
    return norms
end

getcol(r, c) = hasproperty(r, c) ? r[!, c] : fill(missing, nrow(r))

# mean of `col` over rows (nothing if empty), and its std as a "±x" suffix when >1 seed.
function meanstd(sub, col)
    vals = collect(skipmissing(getcol(sub, col)))
    isempty(vals) && return (nothing, "")
    (mean(vals), length(vals) > 1 && std(vals) > 1e-9 ? @sprintf("±%.1f", std(vals)) : "")
end

# regret / total-time cells, each "mean ± std over seeds"; "-" when there is no matching row.
regret_cell(sub) = ((m, s) = meanstd(sub, :regret_pct); m === nothing ? "-" : @sprintf("%.1f%s%%", m, s))
time_cell(sub)   = ((m, s) = meanstd(sub, :total_time); m === nothing ? "-" : @sprintf("%.1f%ss", m, s))

# Significance level for the paired t-test. Raised to 0.10 (from the usual 0.05) because there are
# only ~5 seeds per cell: at that sample size the test has little power, so a stricter α would fail
# to flag all but the most extreme differences. 0.10 trades a higher false-positive rate for the
# ability to actually distinguish methods given the few runs we have.
const ALPHA = 0.1
const MARGIN_ABS = 1.0    # regret practical-equivalence margin (percentage points), and
const MARGIN_REL = 0.25   # as a fraction of the better cell's mean; the larger of the two applies

# per-seed regret and time samples for a cell, keyed by seed so two cells sharing a seed set
# can be paired; nothing when the cell has no evaluated regret.
function samples(sub)
    isempty(sub) && return nothing
    reg = Dict(r.seed => r.regret_pct for r in eachrow(sub) if !ismissing(r.regret_pct))
    isempty(reg) && return nothing
    tim = Dict(r.seed => r.total_time for r in eachrow(sub) if !ismissing(r.total_time))
    (; reg, tim)
end

# the two value-vectors of `a` and `b` restricted to their shared seeds, in a common order.
function paired(a, b)
    ks = sort(collect(intersect(keys(a), keys(b))))
    ([a[k] for k in ks], [b[k] for k in ks])
end

# one-sided paired t-test p that sample `y` is worse (larger regret) than `x`: p that mean(y-x)>0.
# Unlike a rank test it weighs the magnitude of each per-seed difference, so it stays decisive at
# n=5. Returns 1.0 (no evidence) when the samples are tied or empty.
function worse_pvalue(x, y)
    d = y .- x
    (length(d) < 2 || all(iszero, d)) && return 1.0
    pvalue(OneSampleTTest(d); tail = :right)
end

# given one n_rp column of cells, the indices "reasonable to use": the lowest-mean-regret cell
# plus every cell that is not BOTH meaningfully worse (mean regret beyond the equivalence margin)
# AND statistically worse (paired t-test p<ALPHA). The margin handles deterministic near-ties that
# a test flags as significant only because they have zero seed variance; the test excludes noisy
# methods whose large mean gap survives their spread.
function tied_with_best(samps)
    valid = [i for i in eachindex(samps) if samps[i] !== nothing]
    isempty(valid) && return Set{Int}()
    best = argmin(i -> mean(values(samps[i].reg)), valid)
    bmean = mean(values(samps[best].reg))
    margin = max(MARGIN_ABS, MARGIN_REL * bmean)
    tied = Set{Int}()
    for i in valid
        br, ir = paired(samps[best].reg, samps[i].reg)
        gap = mean(ir) - mean(br)
        (gap <= margin || worse_pvalue(br, ir) >= ALPHA) && push!(tied, i)
    end
    return tied
end

# --- Pareto (non-domination) marking across clustering × weight × n_rp × normalization ---
# `dominates(a, b)`: cell `a` dominates cell `b` iff `a` is not meaningfully slower on average AND
#   (i)  `a` has meaningfully/significantly lower MEAN regret, or
#   (ii) `a` has COMPARABLE mean regret but significantly lower regret VARIANCE (an F-test; a
#        deterministic zero-variance method beats a seed-noisy one at equal mean).
# Cells are compared on their shared seeds (paired). A cell is coloured cyan when nothing dominates it.
const TIME_MARGIN_ABS = 0.1    # time practical-equivalence margin (seconds, ~one display unit), and
const TIME_MARGIN_REL = 0.15   # as a fraction of the faster cell's mean; the larger of the two applies

# `a` has meaningfully lower mean regret than `b`: strictly lower mean, the gap clears the
# practical-equivalence margin, AND it is statistically significant (or untestable at a single
# seed). Requiring BOTH — like the bold "reasonable" rule — is essential: deterministic methods
# repeat identically across seeds, so a paired t-test on them reads any nonzero gap as significant
# (zero variance ⇒ p≈0); without the margin gate a sub-margin gap would spuriously "dominate".
function reg_meanbetter(a, b)
    ra, rb = paired(a.reg, b.reg)
    isempty(ra) && return false
    ma, mb = mean(ra), mean(rb)
    ma >= mb && return false
    margin = max(MARGIN_ABS, MARGIN_REL * min(ma, mb))
    (mb - ma > margin) || return false                 # must clear the practical margin, and
    length(ra) < 2 || worse_pvalue(ra, rb) < ALPHA     # be significant (single seed: trust the mean)
end

# neither cell is meaningfully better in mean regret ⇒ their mean regrets are comparable.
reg_comparable(a, b) = !reg_meanbetter(a, b) && !reg_meanbetter(b, a)

# `a` has significantly lower regret variance than `b`: a deterministic (zero-variance) `a` beats a
# varying `b` outright; when both vary, require the variance-ratio F-test to reject equality.
function var_lower(a, b)
    ra, rb = paired(a.reg, b.reg)
    length(ra) < 2 && return false
    va, vb = var(ra), var(rb)
    (va >= vb || vb <= 0) && return false
    va <= 0 && return true
    pvalue(VarianceFTest(ra, rb)) < ALPHA
end

# `a` does not take meaningfully longer than `b` on average (mean total time, shared seeds):
# equal-or-faster, or slower only within the practical-equivalence margin (an absolute sub-second
# floor OR a relative fraction), or not significantly slower by the paired t-test. The floor keeps
# trivial sub-second jitter from blocking a regret win (e.g. a 0.08s hull vs a 0.05s chronological),
# while a genuinely-slower higher-n_rp cell still counts as slower.
function time_not_longer(a, b)
    ta, tb = paired(a.tim, b.tim)
    isempty(ta) && return true
    ma, mb = mean(ta), mean(tb)
    ma <= mb && return true
    margin = max(TIME_MARGIN_ABS, TIME_MARGIN_REL * mb)
    (ma - mb <= margin) || (worse_pvalue(tb, ta) >= ALPHA)
end

function dominates(a, b)
    time_not_longer(a, b) || return false
    reg_meanbetter(a, b) && return true
    return reg_comparable(a, b) && var_lower(a, b)
end

# keys (normalization, n_rp, clustering, weight) of the cells no other cell dominates.
function nondominated_keys(df)
    ks = NTuple{4,Any}[]; cells = Dict{NTuple{4,Any},Any}()
    for nrm in unique(df.normalization), n in unique(df.n_rep_periods),
        cl in unique(df.clustering_type), w in unique(df.weight_type)
        sub = df[(df.normalization .== nrm) .& (df.n_rep_periods .== n) .&
                 (df.clustering_type .== cl) .& (df.weight_type .== w), :]
        s = samples(sub); s === nothing && continue
        k = (nrm, n, cl, w); push!(ks, k); cells[k] = s
    end
    Set(k for k in ks if !any(k2 -> k2 != k && dominates(cells[k2], cells[k]), ks))
end

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

# Generic pairwise normalization comparison — is `normA` meaningfully different from `normB`?
# Pairs the two normalizations cell-by-cell (same clustering × weight × n_rp, shared seeds) at the
# canonical tol, classifies each with the same margin+t-test rule used elsewhere, and pools every
# per-seed regret difference into one two-sided paired t-test. Returns `nothing` when either
# normalization is absent from `arms`; returns `n_total=0` when present but never co-occurring on
# the same (clustering, weight, n_rp) cell. Factored out of `normalization_compare` so the same
# logic backs both the economic-vs-unscaled console report AND the synthetic system's 3-way
# (economic/peak/min-max) comparison used by the sensitivity export.
function normalization_compare_pair(arms, normA, normB)
    df = arms[.!occursin.("nocache", arms.name) .& isapprox.(arms.tol, 0.01; atol=1e-12), :]
    (any(df.normalization .== normA) && any(df.normalization .== normB)) || return nothing
    clusts = unique(df.clustering_type); weights = unique(df.weight_type); nrps = sort(unique(df.n_rep_periods))
    cellsamp(nrm, cl, w, n) = samples(df[(df.normalization .== nrm) .& (df.clustering_type .== cl) .&
                                         (df.weight_type .== w) .& (df.n_rep_periods .== n), :])
    ncomp = nAb = nBb = ntot = 0
    diffs = Float64[]; worst = (d = 0.0, name = "")
    for cl in clusts, w in weights, n in nrps
        aA = cellsamp(normA, cl, w, n); aB = cellsamp(normB, cl, w, n)
        (aA === nothing || aB === nothing) && continue
        ntot += 1
        if reg_meanbetter(aA, aB); nAb += 1
        elseif reg_meanbetter(aB, aA); nBb += 1
        else; ncomp += 1
        end
        dA, dB = paired(aA.reg, aB.reg)
        append!(diffs, dA .- dB)                       # normA − normB, per shared seed
        md = isempty(dA) ? 0.0 : mean(dA) - mean(dB)
        abs(md) > abs(worst.d) && (worst = (d = md, name = "$(combo_label(cl, w)) n_rp=$n"))
    end
    ntot == 0 && return (; normA, normB, n_total = 0, n_comparable = 0, n_a_better = 0, n_b_better = 0,
                          p_value = NaN, mean_diff = NaN, worst_diff = NaN, worst_combo = "")
    p = (length(diffs) < 2 || all(iszero, diffs)) ? 1.0 : pvalue(OneSampleTTest(diffs))
    return (; normA, normB, n_total = ntot, n_comparable = ncomp, n_a_better = nAb, n_b_better = nBb,
            p_value = p, mean_diff = (isempty(diffs) ? 0.0 : mean(diffs)), worst_diff = worst.d,
            worst_combo = worst.name)
end

# Margin-only classification: `a` has meaningfully lower mean regret than `b` by the
# practical-equivalence margin ALONE — no paired significance test. Needed because 5 of the paper's
# 7 clustering methods (conical hull, convex hull, convex-hull-with-null, hierarchical,
# chronological) are fully deterministic, so their 5 "seeds" are bit-identical: a paired t-test on
# such a cell has zero within-cell variance, so it reads any nonzero mean gap as significant
# (p≈0) — a degenerate, statistically invalid test. `reg_meanbetter` (above) papers over this by
# using the margin to gate the test, but still runs/reports a p-value; this variant drops the test
# entirely so the same rule applies uniformly to deterministic and stochastic cells alike.
function reg_meanbetter_margin(a, b)
    ra, rb = paired(a.reg, b.reg)
    isempty(ra) && return false
    ma, mb = mean(ra), mean(rb)
    ma >= mb && return false
    margin = max(MARGIN_ABS, MARGIN_REL * min(ma, mb))
    mb - ma > margin
end

# Margin-only sibling of `normalization_compare_pair`: identical pairing/looping over
# (clustering, weight, n_rp) cells, but classifies each cell as comparable / normA-better /
# normB-better using ONLY the practical-equivalence margin (the larger of 1pp or 25% of the
# better arm's mean regret) — no paired t-test fallback, no pooled p-value. Per reviewer guidance
# ("descriptive win/comparable/loss counts may be sufficient without fragile significance
# language"), and because a paired t-test is statistically invalid on the 5 deterministic
# clustering methods (zero seed variance). Kept as a separate function (rather than changing
# `normalization_compare_pair`'s signature/return shape) since that function is still called
# elsewhere (`normalization_compare`, `export_sensitivity` in export_summary_csvs.jl) expecting
# p_value/worst_diff/worst_combo in its return. Returns `nothing` when either normalization is
# absent from `arms`; `n_total=0` when present but never co-occurring on the same cell.
function normalization_compare_pair_margin_only(arms, normA, normB)
    df = arms[.!occursin.("nocache", arms.name) .& isapprox.(arms.tol, 0.01; atol=1e-12), :]
    (any(df.normalization .== normA) && any(df.normalization .== normB)) || return nothing
    clusts = unique(df.clustering_type); weights = unique(df.weight_type); nrps = sort(unique(df.n_rep_periods))
    cellsamp(nrm, cl, w, n) = samples(df[(df.normalization .== nrm) .& (df.clustering_type .== cl) .&
                                         (df.weight_type .== w) .& (df.n_rep_periods .== n), :])
    ncomp = nAb = nBb = ntot = 0
    diffs = Float64[]
    for cl in clusts, w in weights, n in nrps
        aA = cellsamp(normA, cl, w, n); aB = cellsamp(normB, cl, w, n)
        (aA === nothing || aB === nothing) && continue
        ntot += 1
        if reg_meanbetter_margin(aA, aB); nAb += 1
        elseif reg_meanbetter_margin(aB, aA); nBb += 1
        else; ncomp += 1
        end
        dA, dB = paired(aA.reg, aB.reg)
        append!(diffs, dA .- dB)                       # normA − normB, per shared seed
    end
    ntot == 0 && return (; normA, normB, n_total = 0, n_comparable = 0, n_a_better = 0, n_b_better = 0,
                          mean_diff = NaN)
    return (; normA, normB, n_total = ntot, n_comparable = ncomp, n_a_better = nAb, n_b_better = nBb,
            mean_diff = (isempty(diffs) ? 0.0 : mean(diffs)))
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

# Ablation arm labels (5bus): PROPOSED with one component knocked out.
function ablation_label(ct, wt, nrm)
    if ct == "conical_hull" && wt == "convex" && nrm == "economic"; return "PROPOSED"
    elseif ct == "conical_hull" && wt == "convex" && nrm == "unscaled"; return "ablate: -economic"
    elseif ct == "convex_hull" && wt == "convex" && nrm == "economic"; return "ablate: -conic-selection"
    elseif ct == "conical_hull" && wt == "dirac" && nrm == "economic"; return "ablate: -convex-weights"
    end
    return missing
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
