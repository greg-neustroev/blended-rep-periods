#!/usr/bin/env julia
#
# appendix_normalization.tex §Comparison of the Three Schemes, first paragraph: the synthetic
# development-system 3-way comparison that motivated economic scaling as the default — mean
# regret per scheme, plus the pairwise practical-equivalence classification (no significance test:
# several clustering methods are deterministic, so a paired t-test would be statistically invalid
# here — see analysis/common.jl's `normalization_compare_pair_margin_only`) for economic-vs-peak
# and peak-vs-min-max ("peak" = normalization "unscaled" in the raw data; the paper calls the
# `unscaled` clustering-matrix option "peak-scaled").
#
# Usage: julia --project=. analysis/paper/normalization_synthetic.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

function export_normalization_synthetic()
    loaded = load_case(SYNTH_FILE)
    loaded === nothing && (@warn "no reference optimum for $SYNTH_FILE"; return)
    arms = canonical(loaded.arms)

    means = DataFrame()
    for nrm in ordered_norms(arms)
        r = collect(skipmissing(arms[arms.normalization .== nrm, :regret_pct]))
        isempty(r) && continue
        push!(means, (normalization = nrm, mean_regret_pct = mean(r), n = length(r)); cols=:union)
    end
    write_csv("normalization_synthetic_means.csv", means)

    pairs = DataFrame()
    norms = ordered_norms(arms)
    for i in eachindex(norms), j in (i+1):length(norms)
        r = normalization_compare_pair_margin_only(arms, norms[i], norms[j])
        (r === nothing || r.n_total == 0) && continue
        push!(pairs, (normA = r.normA, normB = r.normB, n_total = r.n_total,
                      n_comparable = r.n_comparable, n_a_better = r.n_a_better,
                      n_b_better = r.n_b_better, mean_diff_pp = r.mean_diff); cols=:union)
    end
    write_csv("normalization_synthetic_pairs.csv", pairs)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_normalization_synthetic()
