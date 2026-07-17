#!/usr/bin/env julia
#
# appendix_normalization.tex §Min-Max Normalization: the dedicated full 112-combination min-max
# sweep on 5-bus (configs/sienna_5bus_minmax.toml, sharing outputs/sienna/5bus.csv's namespace per
# that config's own comment) — mean/median regret per scheme (naive mean is dominated by a few
# catastrophic k-means/k-medoids-at-low-n_rp outliers under economic normalization; median removes
# that distortion), the peak-vs-min-max and economic-vs-min-max cell-by-cell classifications, and
# the same classification restricted to conical hull only (16 combos: the method this paper
# recommends), all via `normalization_compare_pair_margin_only`.
#
# NOTE: this sweep may not have been run yet in the current campaign (it is a late addition, split
# out from the main 5-bus comparison specifically so min-max could be dropped from the main
# campaign) — this script warns and writes nothing until minmax rows actually appear in
# outputs/sienna/5bus.csv.
#
# Usage: julia --project=. analysis/paper/normalization_minmax_5bus.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

const FIVEBUS_PATH = "outputs/sienna/5bus.csv"

function export_normalization_minmax_5bus()
    loaded = load_case(FIVEBUS_PATH)
    loaded === nothing && (@warn "no reference optimum for $FIVEBUS_PATH"; return)
    arms = canonical(loaded.arms)
    if !any(arms.normalization .== "minmax")
        @warn "no minmax rows in $FIVEBUS_PATH yet — run configs/sienna_5bus_minmax.toml first"
        return
    end

    means = DataFrame()
    for nrm in ordered_norms(arms)
        r = collect(skipmissing(arms[arms.normalization .== nrm, :regret_pct]))
        isempty(r) && continue
        push!(means, (normalization = nrm, mean_regret_pct = mean(r), median_regret_pct = median(r),
                      n = length(r)); cols=:union)
    end
    write_csv("normalization_minmax_5bus_means.csv", means)

    pairs = DataFrame()
    for (a, b) in [("unscaled", "minmax"), ("economic", "minmax")]
        r = normalization_compare_pair_margin_only(arms, a, b)
        (r === nothing || r.n_total == 0) && continue
        push!(pairs, (scope = "all_112", normA = r.normA, normB = r.normB, n_total = r.n_total,
                      n_comparable = r.n_comparable, n_a_better = r.n_a_better,
                      n_b_better = r.n_b_better, mean_diff_pp = r.mean_diff); cols=:union)
    end

    # conical-hull-only subset (16 combos: 4 weight types x 4 n_rep) — the method this paper
    # recommends, per the appendix's own callout.
    conical = arms[arms.clustering_type .== "conical_hull", :]
    conical_means = DataFrame()
    for nrm in ordered_norms(conical)
        r = collect(skipmissing(conical[conical.normalization .== nrm, :regret_pct]))
        isempty(r) && continue
        push!(conical_means, (normalization = nrm, mean_regret_pct = mean(r),
                              median_regret_pct = median(r), n = length(r)); cols=:union)
    end
    write_csv("normalization_minmax_5bus_conical_means.csv", conical_means)

    for (a, b) in [("unscaled", "minmax"), ("economic", "minmax")]
        r = normalization_compare_pair_margin_only(conical, a, b)
        (r === nothing || r.n_total == 0) && continue
        push!(pairs, (scope = "conical_hull_only", normA = r.normA, normB = r.normB, n_total = r.n_total,
                      n_comparable = r.n_comparable, n_a_better = r.n_a_better,
                      n_b_better = r.n_b_better, mean_diff_pp = r.mean_diff); cols=:union)
    end
    write_csv("normalization_minmax_5bus_pairs.csv", pairs)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_normalization_minmax_5bus()
