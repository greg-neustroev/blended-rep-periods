#!/usr/bin/env julia
#
# curtailment.csv (flagged low-confidence; superseded by curtailment_renewable.jl for the actual
# renewable-curtailment prose — kept for the console diagnostic's SECONDARY METRICS section) and
# capacity_decisions.csv, which backs case_studies.tex §Additional Evaluation Metrics' "capacity
# decision deviation" prose (11-59% on 5-bus, 30-73% on GEP).
#
# Usage: julia --project=. analysis/paper/secondary.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

function df_ref_name(path)
    df = CSV.read(joinpath(REPO_ROOT, path), DataFrame)
    df[df.n_rep_periods .== 1, :name]
end

function export_secondary()
    curt = DataFrame(); cap = DataFrame()
    for (path, meta) in CASE_META
        loaded = load_case(path)
        loaded === nothing && continue
        arms = canonical(loaded.arms)
        grid = sort(unique(arms.n_rep_periods))
        for (cl, w) in [PROPOSED, ("k_means", "dirac"), ("k_medoids", "dirac")], n in grid
            sub = arms[(arms.normalization .== "economic") .& (arms.clustering_type .== cl) .&
                       (arms.weight_type .== w) .& (arms.n_rep_periods .== n), :]
            isempty(sub) && continue
            f(c) = (v = collect(skipmissing(getcol(sub, c))); isempty(v) ? missing : mean(v))
            push!(curt, (case_study = meta.case_study, clustering_type = cl, weight_type = w,
                        n_rep_periods = n, total_spillage = f(:total_spillage),
                        total_borrow = f(:total_borrow), max_weight_sum = f(:max_weight_sum),
                        eval_cost_of_operations = f(:eval_cost_of_operations),
                        eval_cost_of_borrow = f(:eval_cost_of_borrow)); cols=:union)
            if meta.problem_class == "investment"
                d = capacity_l1_diff(joinpath(REPO_ROOT, path), first(sub.name),
                                      first(df_ref_name(path)))
                push!(cap, (case_study = meta.case_study, clustering_type = cl, weight_type = w,
                            n_rep_periods = n, capacity_l1_diff_pct = d); cols=:union)
            end
        end
    end
    write_csv("curtailment.csv", curt)
    write_csv("capacity_decisions.csv", cap)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_secondary()
