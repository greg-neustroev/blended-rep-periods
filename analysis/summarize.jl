#!/usr/bin/env julia
#
# Summarize proposed-method results into the three deliverable tables:
#   1. Comparison  — traditional k-means / k-medoids / hierarchical vs the PROPOSED method,
#                    regret (%) across the representative-period grid.
#   2. Ablation    — PROPOSED vs each single-component knockout, regret (%) across the grid.
#   3. Sensitivity — PROPOSED selection+weights at n_rp = 10, regret (%) across the PGD
#                    tolerance and normalization.
#
# Regret = evaluated_objective_value / reference_optimum − 1, where the reference optimum is
# the objective of the n_rep = 1 (full-horizon) reference row. Results are averaged over
# seeds (only k-means differs across seeds; the hull / k-medoids / hierarchical arms are
# deterministic).
#
# Usage:
#   julia --project=BlendedClustering analysis/summarize.jl [result.csv ...]
# With no arguments it summarizes outputs/tyndp/gep.csv and outputs/sienna/5bus.csv if present.

using CSV, DataFrames, Printf, Statistics

# Classify one arm (a set of pipeline options) into a human label. `is_storage` tells us
# whether the PROPOSED method carries a chain matrix on this dataset (storage) or not
# (investment), so the same (selection, weight, normalization) triple reads correctly.
function arm_label(ct, wt, nrm, has_chain, is_storage)
    traditional = Dict("k_means" => "k-means", "k_medoids" => "k-medoids", "hierarchical" => "hierarchical")
    if wt == "dirac" && haskey(traditional, ct) && nrm == "unscaled" && !has_chain
        return "$(traditional[ct]) (traditional)"
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

# Print a label × n_rep table of regret percentages, rows in the given `order`.
function print_grid_table(title, sub, order)
    println("\n", title)
    grid = sort(unique(sub.n_rep_periods))
    @printf("  %-28s", "arm \\ n_rp")
    for n in grid
        @printf("%10d", n)
    end
    println()
    for label in order
        rows = sub[sub.label .== label, :]
        isempty(rows) && continue
        @printf("  %-28s", label)
        for n in grid
            r = rows[rows.n_rep_periods .== n, :]
            if isempty(r) || all(ismissing, r.regret_pct)
                @printf("%10s", "-")
            else
                @printf("%9.1f%%", mean(skipmissing(r.regret_pct)))
            end
        end
        println()
    end
end

function summarize_file(path)
    println("\n", "="^78)
    println("RESULTS: ", path)
    println("="^78)
    df = CSV.read(path, DataFrame)
    for c in ("clustering_type", "weight_type", "normalization")
        df[!, c] = string.(df[!, c])
    end

    # Reference optimum: the n_rep == 1 full-horizon solve's reduced objective.
    ref_rows = df[df.n_rep_periods .== 1, :]
    if isempty(ref_rows) || all(ismissing, ref_rows.objective_value)
        println("  !! no n_rep=1 reference row with an objective — cannot compute regret.")
        return
    end
    ref_opt = mean(skipmissing(ref_rows.objective_value))
    @printf("  reference optimum (n_rep=1): %.6g\n", ref_opt)

    arms = df[df.n_rep_periods .!= 1, :]
    is_storage = any(arms.fix_every .> 1)   # storage arms grade with fix_every=30; investment=1
    arms.has_chain = occursin.("chain", arms.name)
    arms.label = [arm_label(r.clustering_type, r.weight_type, r.normalization, r.has_chain, is_storage)
                  for r in eachrow(arms)]
    # Regret %, blank when the evaluation model was not solved to optimality.
    arms.regret_pct = map(eachrow(arms)) do r
        (ismissing(r.evaluated_objective_value) || r.evaluated_objective_value <= 0) && return missing
        100 * (r.evaluated_objective_value / ref_opt - 1)
    end

    # 1 + 2: comparison and ablation use the main-sweep tolerance (tol = 0.01) only.
    main = arms[isapprox.(arms.tol, 0.01; atol=1e-12), :]
    comparison_order = ["k-means (traditional)", "k-medoids (traditional)",
        "hierarchical (traditional)", "PROPOSED"]
    print_grid_table("1. COMPARISON — regret (%) vs n_rp", main, comparison_order)

    ablation_order = ["PROPOSED", "ablate: -economic", "ablate: -conic-selection",
        "ablate: -convex-weights", "ablate: -chain-split"]
    print_grid_table("2. ABLATION — regret (%) vs n_rp (each row knocks out one component)", main, ablation_order)

    # 3: sensitivity — the proposed selection+weights at n_rp = 10, over tol × normalization.
    # Pin the proposed chain setting (chain=convex on storage, none on investment) so the
    # tol=0.01 cell is the proposed arm alone, not averaged with the -chain ablation.
    sens = arms[(arms.clustering_type .== "conical_hull") .& (arms.weight_type .== "convex") .&
                (arms.n_rep_periods .== 10) .& (arms.has_chain .== is_storage), :]
    if !isempty(sens)
        println("\n3. SENSITIVITY — PROPOSED (conical_hull + convex) at n_rp=10, regret (%) vs tol")
        tols = sort(unique(sens.tol); rev=true)
        @printf("  %-16s", "normalization \\ tol")
        for t in tols
            @printf("%12s", t)
        end
        println()
        for nrm in sort(unique(sens.normalization))
            @printf("  %-16s", nrm)
            for t in tols
                r = sens[(sens.normalization .== nrm) .& isapprox.(sens.tol, t; atol=1e-12), :]
                if isempty(r) || all(ismissing, r.regret_pct)
                    @printf("%12s", "-")
                else
                    @printf("%11.1f%%", mean(skipmissing(r.regret_pct)))
                end
            end
            println()
        end
    end
    return
end

function main()
    files = isempty(ARGS) ? ["outputs/tyndp/gep.csv", "outputs/sienna/5bus.csv"] : ARGS
    any_found = false
    for f in files
        if isfile(f)
            any_found = true
            summarize_file(f)
        else
            println("(skipping $f — not found; run the case study first)")
        end
    end
    any_found || println("\nNo result CSVs found. Run e.g.:\n" *
        "  julia --project=BlendedClustering -e 'using BlendedClustering; run_case_studies(\"configs/gep.toml\")'")
    return
end

main()
