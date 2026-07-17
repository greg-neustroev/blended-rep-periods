#!/usr/bin/env julia
#
# Table tab:case-studies (case_studies.tex): the #var./#constr. columns (the |N|,|G|,... structural
# columns are case_study_sizes.jl). One row per (case_study, n_rep_periods), including the n_rep=1
# reference. n_variables/n_constraints are invariant to clustering/weight choice at fixed n_rp;
# asserted with a @warn if a run ever breaks that invariant.
#
# Usage: julia --project=. analysis/paper/model_size.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

function export_model_size()
    out = DataFrame()
    for (path, meta) in CASE_META
        df = CSV.read(joinpath(REPO_ROOT, path), DataFrame)
        for n in sort(unique(df.n_rep_periods))
            g = df[df.n_rep_periods .== n, :]
            vs = unique(skipmissing(g.n_variables)); cs = unique(skipmissing(g.n_constraints))
            evs = unique(skipmissing(g.eval_n_variables)); ecs = unique(skipmissing(g.eval_n_constraints))
            (length(vs) > 1 || length(cs) > 1) &&
                @warn "model size not invariant at $(meta.case_study) n_rp=$n: n_variables=$vs n_constraints=$cs"
            push!(out, (case_study = meta.case_study, problem_class = meta.problem_class,
                        n_rep_periods = n, n_variables = isempty(vs) ? missing : first(vs),
                        n_constraints = isempty(cs) ? missing : first(cs),
                        eval_n_variables = isempty(evs) ? missing : first(evs),
                        eval_n_constraints = isempty(ecs) ? missing : first(ecs)); cols=:union)
        end
    end
    write_csv("model_size.csv", out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_model_size()
