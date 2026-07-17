#!/usr/bin/env julia
#
# Table tab:knobs (case_studies.tex): parameter sweep on the synthetic dataset, regret range and
# per-axis spread. One row per (knob group, swept value), columns range_min_pct/range_max_pct
# (the table's "Range, %" column) and weight_spread_pp/clustering_spread_pp/normalization_spread_pp
# (the table's "Weight"/"Clust."/"Norm." spread-in-percentage-points columns) — mirrors
# tab:knobs's row/column layout directly (min/max kept as two numeric columns rather than a
# formatted "0.7--5.7" range string, since this is a data file, not the typeset table itself).
#
# Re-implemented (not `include`d) because summarize_knobs.jl calls `main()` unconditionally at
# file scope — including it here would execute its CLI entry point with the wrong ARGS/cwd.
#
# Usage: julia --project=. analysis/paper/knob_sensitivity.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

const KNOB_GROUPS = [
    ("base",  "baseline",                       "β=60, renew=0.95, all knobs at default"),
    ("beta",  "buffer horizon (beta_days)",     "gates the inter-period storage reconstruction"),
    ("renew", "storage reliance (renew_share)", "gates storage throughput / whether error costs"),
    ("hull",  "residual-load hull geometry",    "-> weight class (convex vs conical vs bounded)"),
    ("reg",   "intrinsic dimensionality",       "-> selection and n_rp"),
    ("scale", "cross-block feature scale",      "-> normalization"),
    ("noise", "day-to-day irregularity",        "-> out-of-hull / reconstruction stress"),
    ("sd",    "seasonal:diurnal amplitude",     "-> inter- vs intra-period stress"),
]
function knob_suffix_num(name)
    m = match(r"(\d+)$", name); m === nothing ? -1 : parse(Int, m.captures[1])
end
function knob_group_of(name)
    best = nothing
    for (pre, _, _) in KNOB_GROUPS
        if startswith(name, pre) && (best === nothing || length(pre) > length(best)); best = pre; end
    end
    best
end
function knob_axes(path)
    d = CSV.read(path, DataFrame)
    for c in ("clustering_type", "weight_type", "normalization"); d[!, c] = string.(d[!, c]); end
    refrows = d[d.n_rep_periods .== 1, :]
    (isempty(refrows) || all(ismissing, refrows.objective_value)) && return nothing
    ref = mean(skipmissing(refrows.objective_value))
    a = d[d.n_rep_periods .!= 1, :]
    reg(sub) = [100 * (e / ref - 1) for e in sub.evaluated_objective_value if !ismissing(e) && e > 0]
    spr(v) = isempty(v) ? NaN : maximum(v) - minimum(v)
    W = reg(a[(a.clustering_type .== "conical_hull") .& (a.normalization .== "economic"), :])
    S = reg(a[(a.weight_type .== "convex") .& (a.normalization .== "economic"), :])
    N = reg(a[(a.clustering_type .== "conical_hull") .& (a.weight_type .== "convex"), :])
    allreg = reg(a)
    return (weight = spr(W), selection = spr(S), norm = spr(N),
            lo = isempty(allreg) ? NaN : minimum(allreg), hi = isempty(allreg) ? NaN : maximum(allreg))
end

function export_knob_sensitivity()
    dir = joinpath(REPO_ROOT, "outputs", "synth", "knobs")
    isdir(dir) || (@warn "no knob-study output dir: $dir"; return)
    out = DataFrame()
    for f in filter(f -> endswith(f, ".csv"), readdir(dir))
        variant = replace(f, ".csv" => "")
        grp = knob_group_of(variant); grp === nothing && continue
        r = try knob_axes(joinpath(dir, f)) catch err; @warn "skipping $f" exception=err; nothing end
        r === nothing && continue
        label, note = only(g[2:3] for g in KNOB_GROUPS if g[1] == grp)
        push!(out, (knob_group = grp, knob_label = label, knob_note = note, variant = variant,
                    knob_value = knob_suffix_num(variant),
                    range_min_pct = r.lo, range_max_pct = r.hi,
                    weight_spread_pp = r.weight, clustering_spread_pp = r.selection,
                    normalization_spread_pp = r.norm); cols=:union)
    end
    write_csv("knob_sensitivity.csv", out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_knob_sensitivity()
