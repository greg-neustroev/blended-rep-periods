#!/usr/bin/env julia
#
# Proposition 1's remark and methodology.tex's acceleration-ablation prose ("orders-of-magnitude
# fewer iterations than plain PGD" / "Section~\ref{sec:sensitivity-ablation}'s acceleration
# ablation confirms the resulting speedup empirically") — plain-PGD vs Gram-only vs Gram+FISTA+
# restart, synth dev system. Isolates the Gram-matrix reformulation's contribution from FISTA's.
# The "proposed" arm already exists in outputs/synth/hydro.csv (from the main sweep); the
# "gram_only"/"plain" arms come from configs/synth_acceleration.toml, sharing the same output file
# (their pgd_variant is encoded as a name suffix, per Types.jl's ExperimentData).
#
# Prose-supporting data (not itself a printed table) — kept as a tidy per-(variant, n_rep) frame.
#
# Usage: julia --project=. analysis/paper/acceleration_ablation.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

function export_acceleration_ablation()
    isfile(joinpath(REPO_ROOT, SYNTH_FILE)) || return
    df = CSV.read(joinpath(REPO_ROOT, SYNTH_FILE), DataFrame)
    for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
    refrows = df[df.n_rep_periods .== 1, :]
    isempty(refrows) && return
    ref_opt = mean(skipmissing(refrows.objective_value))

    variant_of(name) = occursin("_gram_only", name) ? "gram_only" :
                        occursin("_plain", name) ? "plain" : "proposed"
    base = df[(df.clustering_type .== "conical_hull") .& (df.weight_type .== "convex") .&
              (df.normalization .== "economic") .& isapprox.(df.tol, 0.01; atol=1e-12) .&
              .!occursin.("nocache", df.name) .& (df.n_rep_periods .!= 1), :]
    isempty(base) && (@warn "no acceleration-ablation rows found in $SYNTH_FILE"; return)
    base = copy(base)
    base.variant = variant_of.(base.name)
    base.regret_pct = map(eachrow(base)) do r
        (ismissing(r.evaluated_objective_value) || r.evaluated_objective_value <= 0) && return missing
        100 * (r.evaluated_objective_value / ref_opt - 1)
    end

    out = DataFrame()
    for variant in ("proposed", "gram_only", "plain"), n in sort(unique(base.n_rep_periods))
        sub = base[(base.variant .== variant) .& (base.n_rep_periods .== n), :]
        isempty(sub) && continue
        iters = collect(skipmissing(getcol(sub, :pgd_total_iters)))
        times = collect(skipmissing(getcol(sub, :time_to_fit_weights)))
        regret = collect(skipmissing(sub.regret_pct))
        proj = collect(skipmissing(getcol(sub, :projection_error)))
        push!(out, (variant = variant, n_rep_periods = n, n_seeds = nrow(sub),
                    iters_mean = isempty(iters) ? missing : mean(iters),
                    iters_sd = length(iters) > 1 ? std(iters) : 0.0,
                    time_mean_s = isempty(times) ? missing : mean(times),
                    time_sd_s = length(times) > 1 ? std(times) : 0.0,
                    regret_mean_pct = isempty(regret) ? missing : mean(regret),
                    projection_error_mean = isempty(proj) ? missing : mean(proj)); cols=:union)
    end
    write_csv("acceleration_ablation.csv", out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_acceleration_ablation()
