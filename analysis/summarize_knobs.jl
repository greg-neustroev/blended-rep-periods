#!/usr/bin/env julia
#
# Summarize the DATA-PROPERTY sensitivity study (configs/synth_knobs.toml). Each result file is
# one synth/hydro variant that differs from the baseline in ONE generator knob, run through the
# per-axis diagnostic panel. For every variant we report the regret spread WITHIN each method
# axis — so reading down a knob's block shows which method axis that data property moves:
#
#   weight-axis spread  = max−min regret over conical_hull × {dirac,convex,conical,conical_bounded} @ economic
#   selection-axis      = over {k_means,k_medoids,hierarchical,chronological,convex_hull,conical_hull} × convex @ economic
#   normalization-axis  = over conical_hull/convex × {economic,unscaled,minmax}
#   overall             = min..max regret over the whole panel (the dataset's total method spread)
#
# A knob the method is SENSITIVE to widens the relevant axis as its value grows; a knob it is
# ROBUST to leaves the axes flat. β and reliance gate the whole panel (overall spread → 0 below
# their thresholds). Regret = evaluated_objective_value / reference_optimum − 1 (n_rep=1 row).
#
# Usage: julia --project=BlendedClustering analysis/summarize_knobs.jl [outputs/synth/knobs]

using CSV, DataFrames, Printf, Statistics

const DIR = isempty(ARGS) ? "outputs/synth/knobs" : ARGS[1]

# knob groups: name prefix => (label, note). Order sets the print order.
const GROUPS = [
    ("base",  "baseline",                         "β=60, renew=0.95, all knobs at default"),
    ("beta",  "β — buffer horizon",               "GATE: inter-period storage reconstruction"),
    ("renew", "storage reliance (renew_share)",   "GATE: storage throughput / whether error costs"),
    ("hull",  "residual-load hull geometry",      "→ WEIGHT class (convex vs conical vs bounded)"),
    ("reg",   "intrinsic dimensionality",         "→ SELECTION and n_rp"),
    ("scale", "cross-block feature scale",        "→ NORMALIZATION"),
    ("noise", "day-to-day noise",                 "→ ROBUSTNESS (expect axes flat)"),
    ("sd",    "seasonal:diurnal amplitude",       "→ inter- vs intra-period stress"),
]

# trailing integer of a variant name, for numeric ordering within a group (base ⇒ -1 ⇒ first).
function suffix_num(name)
    m = match(r"(\d+)$", name)
    m === nothing ? -1 : parse(Int, m.captures[1])
end

# which group prefix this variant name belongs to (longest matching prefix wins: base vs beta).
function group_of(name)
    best = nothing
    for (pre, _, _) in GROUPS
        if startswith(name, pre) && (best === nothing || length(pre) > length(best))
            best = pre
        end
    end
    return best
end

# per-axis regret spreads for one variant file.
function axes(path)
    d = CSV.read(path, DataFrame)
    for c in ("clustering_type", "weight_type", "normalization")
        d[!, c] = string.(d[!, c])
    end
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

function main()
    isdir(DIR) || (println("no study output dir: $DIR (run configs/synth_knobs.toml first)"); return)
    files = filter(f -> endswith(f, ".csv"), readdir(DIR))
    isempty(files) && (println("no result CSVs under $DIR"); return)
    # variant name (basename without .csv) => axes
    data = Dict{String,Any}()
    for f in files
        r = try
            axes(joinpath(DIR, f))
        catch err
            @warn "skipping $f" exception=err
            nothing
        end
        r === nothing || (data[replace(f, ".csv" => "")] = r)
    end
    println("\n", "="^88)
    println("DATA-PROPERTY SENSITIVITY — per-axis regret spread (pp) by knob   [", DIR, "]")
    println("="^88)
    println("axis spread = how much regret varies WITHIN that method axis; large ⇒ method sensitive there.")
    for (pre, label, note) in GROUPS
        vs = sort([n for n in keys(data) if group_of(n) == pre], by = suffix_num)
        isempty(vs) && continue
        println("\n", label, "   [", note, "]")
        @printf("  %-12s %10s %10s %10s %14s\n", "variant", "weight", "selection", "normaliz.", "overall[lo..hi]")
        for n in vs
            r = data[n]
            @printf("  %-12s %9.1f %10.1f %10.1f  %6.1f .. %-6.1f\n",
                n, r.weight, r.selection, r.norm, r.lo, r.hi)
        end
    end
    println()
end

main()
