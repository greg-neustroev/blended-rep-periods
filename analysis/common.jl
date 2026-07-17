#!/usr/bin/env julia
#
# Shared regret/mismatch/statistics primitives, included by BOTH the console diagnostic tool
# (summarize.jl, summarize_knobs.jl) and the paper data-export scripts (analysis/paper/*.jl, via
# analysis/paper/common.jl), so the paper's tables/figures and the console diagnostic can never
# silently disagree. Deliberately has NO opinion on ARGS or output paths (those differ between the
# console tool and the export scripts) and does not read any CSV itself.
#
# Not meant to be run directly.

using CSV, DataFrames, Printf, Statistics
using Arrow
using HypothesisTests

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

const METHOD_ORDER = ["k_means", "k_medoids", "hierarchical", "chronological", "convex_hull", "convex_hull_with_null", "conical_hull"]
const WEIGHT_ORDER = ["dirac", "convex", "conical_bounded", "conical"]   # displayed: dirac, convex, subunit, conical
const NORM_ORDER = ["economic", "unscaled", "minmax"]
const PROPOSED = ("conical_hull", "convex")   # the method whose runtime/secondary metrics are reported
const HULL = ("convex_hull", "convex_hull_with_null", "conical_hull")          # the blended-hull family
const CONVENTIONAL = ("k_means", "k_medoids", "hierarchical", "chronological")  # the baselines it competes with

# stages that sum to the reported total time (matches the RUNTIME breakdown in summarize.jl).
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

# Generic pairwise normalization comparison — is `normA` meaningfully different from `normB`?
# Pairs the two normalizations cell-by-cell (same clustering × weight × n_rp, shared seeds) at the
# canonical tol, classifies each with the same margin+t-test rule used elsewhere, and pools every
# per-seed regret difference into one two-sided paired t-test. Returns `nothing` when either
# normalization is absent from `arms`; returns `n_total=0` when present but never co-occurring on
# the same (clustering, weight, n_rp) cell.
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
# clustering methods (zero seed variance). Returns `nothing` when either normalization is absent
# from `arms`; `n_total=0` when present but never co-occurring on the same cell.
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

# Ablation arm labels (5bus): PROPOSED with one component knocked out.
function ablation_label(ct, wt, nrm)
    if ct == "conical_hull" && wt == "convex" && nrm == "economic"; return "PROPOSED"
    elseif ct == "conical_hull" && wt == "convex" && nrm == "unscaled"; return "ablate: -economic"
    elseif ct == "convex_hull" && wt == "convex" && nrm == "economic"; return "ablate: -conic-selection"
    elseif ct == "conical_hull" && wt == "dirac" && nrm == "economic"; return "ablate: -convex-weights"
    end
    return missing
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

