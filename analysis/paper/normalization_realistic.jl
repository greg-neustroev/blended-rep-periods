#!/usr/bin/env julia
#
# appendix_normalization.tex's sanity-check paragraph: peak-scaled vs. economic, paired cell by
# cell (same clustering, weight type, n_rep, matched over 5 seeds) across all 112 combinations per
# case study, on the four REAL case studies (not the synthetic development system —
# normalization_synthetic.jl). Same practical-equivalence-margin-only classification (no
# significance test) as normalization_synthetic.jl, via
# `normalization_compare_pair_margin_only`.
#
# Usage: julia --project=. analysis/paper/normalization_realistic.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

function export_normalization_realistic()
    out = DataFrame()
    for (path, meta) in CASE_META
        loaded = load_case(path)
        loaded === nothing && (@warn "no reference optimum for $path — skipping"; continue)
        arms = canonical(loaded.arms)
        r = normalization_compare_pair_margin_only(arms, "unscaled", "economic")
        (r === nothing || r.n_total == 0) &&
            (@warn "no economic/peak overlap for $(meta.case_study) yet"; continue)
        push!(out, (
            case_study = meta.case_study, n_total = r.n_total, n_comparable = r.n_comparable,
            n_peak_better = r.n_a_better, n_economic_better = r.n_b_better,
            mean_diff_pp_peak_minus_economic = r.mean_diff,
        ); cols=:union)
    end
    write_csv("normalization_realistic.csv", out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_normalization_realistic()
