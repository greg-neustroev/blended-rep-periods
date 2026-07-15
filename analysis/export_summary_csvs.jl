#!/usr/bin/env julia
#
# Stage 1 of the paper's two-stage data pipeline: reduce the raw experiment output CSVs/Arrow
# dumps under outputs/ into tidy, long-format summary CSVs written into the PAPER repo's data/
# directory. Stage 2 (TikZ/pgfplots scripts under clustering/plots/) reads ONLY those CSVs, so
# regenerating figures never requires re-running this script or having a solver license.
#
# Reuses summarize.jl's regret/Pareto/statistics logic via `include` so the paper's figures and
# the console diagnostic tool can never silently disagree. summarize.jl's own `main()` is guarded
# (only runs when it is the top-level script), so including it here is safe.
#
# Usage: julia --project=. analysis/export_summary_csvs.jl [outdir]
# (run from the experiments/ repo root; outdir defaults to ../clustering/data relative to repo root)

using CSV, DataFrames, Printf, Statistics
using Arrow

include(joinpath(@__DIR__, "summarize.jl"))

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTDIR = length(ARGS) >= 1 ? ARGS[1] : normpath(joinpath(REPO_ROOT, "..", "clustering", "data"))

# Case study -> (case_study id, problem_class, source). RTS deliberately excluded here — the
# campaign hasn't finished it yet; add a row once outputs/gridmod/rts.csv is complete, per the plan.
const CASE_META = [
    ("outputs/tyndp/gep.csv",   (case_study="gep",    problem_class="investment", source="tyndp")),
    ("outputs/sienna/5bus.csv", (case_study="5bus",   problem_class="investment", source="sienna")),
    ("outputs/tyndp/p2x.csv",   (case_study="p2x",    problem_class="dispatch",   source="tyndp")),
    ("outputs/nrel/118bus.csv", (case_study="118bus", problem_class="dispatch",   source="nrel")),
]
const SYNTH_FILE = "outputs/synth/hydro.csv"

mkpath(OUTDIR)

# ---------------------------------------------------------------------------------------------
# Helpers shared across exporters
# ---------------------------------------------------------------------------------------------

# reference optimum (mean objective_value over n_rep_periods==1 rows) and the arms (n_rep>1),
# with the derived columns every exporter needs: regret_pct, mismatch_pct, total_time (stage sum,
# matching summarize.jl's TIME_STAGES), cache_hit_rate, plus string-coerced category columns.
function load_case(path)
    df = CSV.read(joinpath(REPO_ROOT, path), DataFrame)
    for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
    refrows = df[df.n_rep_periods .== 1, :]
    (isempty(refrows) || all(ismissing, refrows.objective_value)) && return nothing
    ref_opt = mean(skipmissing(refrows.objective_value))
    arms = df[df.n_rep_periods .!= 1, :]
    arms.regret_pct = map(eachrow(arms)) do r
        (ismissing(r.evaluated_objective_value) || r.evaluated_objective_value <= 0) && return missing
        100 * (r.evaluated_objective_value / ref_opt - 1)
    end
    # Solution mismatch: |objective_value - evaluated_objective_value| / evaluated_objective_value.
    # A reliability diagnostic distinct from regret: how much the reduced model's OWN reported
    # objective differs from the true cost obtained by fixing its decisions and re-solving at full
    # resolution, i.e. how much to trust the reduced model's self-reported objective as a predictor
    # of real performance (regret instead grades decision quality against the full-resolution
    # optimum). Computable directly from existing columns -- no new experiments needed.
    arms.mismatch_pct = map(eachrow(arms)) do r
        (ismissing(r.objective_value) || ismissing(r.evaluated_objective_value) ||
         r.evaluated_objective_value == 0) && return missing
        100 * abs(r.objective_value - r.evaluated_objective_value) / abs(r.evaluated_objective_value)
    end
    arms.total_time = reduce(.+, (coalesce.(getcol(arms, s), 0.0) for s in TIME_STAGES))
    arms.cache_hit_rate = map(eachrow(arms)) do r
        tot = coalesce(r.cache_hits, 0) + coalesce(r.cache_misses, 0)
        tot > 0 ? 100 * r.cache_hits / tot : missing
    end
    return (; ref_opt, arms)
end

canonical(df) = df[.!occursin.("nocache", df.name) .& isapprox.(df.tol, 0.01; atol=1e-12), :]

write_csv(name, df) = (p = joinpath(OUTDIR, name); CSV.write(p, df); println("  wrote $p ($(nrow(df)) rows)"))

# ---------------------------------------------------------------------------------------------
# 1-2. regret_by_seed.csv / regret_summary.csv (the four real case studies)
# ---------------------------------------------------------------------------------------------

function export_regret()
    by_seed = DataFrame()
    summary = DataFrame()
    for (path, meta) in CASE_META
        loaded = load_case(path)
        loaded === nothing && (@warn "no reference optimum for $path — skipping"; continue)
        arms = canonical(loaded.arms)   # canonical tol=0.01, exclude nocache twins

        sub = select(arms, :n_rep_periods, :clustering_type, :weight_type, :normalization, :seed,
                     :regret_pct, :mismatch_pct, :total_time => :total_time_s,
                     :time_to_preprocess, :time_to_cluster, :time_to_fit_weights,
                     :time_to_formulate_model, :time_to_solve, :time_to_read,
                     :projection_error, :max_weight_sum, :cache_hit_rate, :pgd_total_iters,
                     :termination_status, :evaluation_termination_status,
                     :n_variables, :n_constraints, :eval_n_variables, :eval_n_constraints)
        sub.case_study .= meta.case_study; sub.problem_class .= meta.problem_class; sub.source .= meta.source
        by_seed = vcat(by_seed, sub; cols=:union)

        # Pareto-frontier and per-group "reasonable" flags, computed with summarize.jl's own
        # functions so the figures can never disagree with the console diagnostic.
        nondom = nondominated_keys(arms)
        for nrm in ordered_norms(arms), n in sort(unique(arms.n_rep_periods))
            g = arms[(arms.normalization .== nrm) .& (arms.n_rep_periods .== n), :]
            isempty(g) && continue
            clusts = unique(g.clustering_type); weights = unique(g.weight_type)
            keys_ = [(cl, w) for cl in clusts for w in weights]
            samps = [samples(g[(g.clustering_type .== cl) .& (g.weight_type .== w), :]) for (cl, w) in keys_]
            tied = Set(keys_[i] for i in tied_with_best(samps))
            for (cl, w) in keys_
                rows = g[(g.clustering_type .== cl) .& (g.weight_type .== w), :]
                isempty(rows) && continue
                (m_r, s_r) = meanstd(rows, :regret_pct)
                (m_t, s_t) = meanstd(rows, :total_time)
                (m_m, s_m) = meanstd(rows, :mismatch_pct)
                sdr = isempty(s_r) ? 0.0 : parse(Float64, replace(s_r, "±" => ""))
                push!(summary, (
                    case_study = meta.case_study, problem_class = meta.problem_class, source = meta.source,
                    n_rep_periods = n, clustering_type = cl, weight_type = w, normalization = nrm,
                    n_seeds = nrow(rows), regret_mean_pct = something(m_r, NaN), regret_sd_pct = sdr,
                    mismatch_mean_pct = something(m_m, NaN),
                    time_mean_s = something(m_t, NaN),
                    is_deterministic = sdr < 1e-9,
                    pareto_frontier = (nrm, n, cl, w) in nondom,
                    reasonable_in_group = (cl, w) in tied,
                ); cols=:union)
            end
        end
    end
    write_csv("regret_by_seed.csv", by_seed)
    write_csv("regret_summary.csv", summary)
end

# ---------------------------------------------------------------------------------------------
# 3. model_size.csv — one row per (case_study, n_rep_periods), including the n_rep=1 reference.
#    n_variables/n_constraints are invariant to clustering/weight choice at fixed n_rp; assert it.
# ---------------------------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------------------------
# 4. cache_hit_rate.csv — hit rate per real case study x hull type; cached-vs-uncached wall time
#    only exists for the synth dev system (the only dataset with matched _nocache twins).
# ---------------------------------------------------------------------------------------------

function export_cache_hit_rate()
    out = DataFrame()
    for (path, meta) in CASE_META
        loaded = load_case(path)
        loaded === nothing && continue
        arms = canonical(loaded.arms)
        for cl in HULL
            g = arms[arms.clustering_type .== cl, :]
            isempty(g) && continue
            hits = sum(skipmissing(getcol(g, :cache_hits))); misses = sum(skipmissing(getcol(g, :cache_misses)))
            tot = hits + misses
            push!(out, (case_study = meta.case_study, clustering_type = cl,
                        cache_hits = hits, cache_misses = misses,
                        hit_rate_pct = tot > 0 ? 100 * hits / tot : missing,
                        cluster_time_s_cached = missing, cluster_time_s_uncached = missing); cols=:union)
        end
    end
    # dev-system cached-vs-uncached wall time (explicitly labeled as such downstream)
    if isfile(joinpath(REPO_ROOT, SYNTH_FILE))
        df = CSV.read(joinpath(REPO_ROOT, SYNTH_FILE), DataFrame)
        for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
        nocache = df[occursin.("nocache", df.name), :]
        for r in eachrow(nocache)
            cached = df[df.name .== replace(r.name, "_nocache" => ""), :]
            isempty(cached) && continue
            cr = first(eachrow(cached))
            hits = coalesce(cr.cache_hits, 0); misses = coalesce(cr.cache_misses, 0); tot = hits + misses
            push!(out, (case_study = "synth_dev", clustering_type = r.clustering_type,
                        cache_hits = hits, cache_misses = misses,
                        hit_rate_pct = tot > 0 ? 100 * hits / tot : missing,
                        cluster_time_s_cached = cr.time_to_cluster, cluster_time_s_uncached = r.time_to_cluster);
                  cols=:union)
        end
    end
    write_csv("cache_hit_rate.csv", out)
end

# ---------------------------------------------------------------------------------------------
# 5-6. sensitivity_tol.csv / pgd_iters_by_tol.csv (synth dev system, n_rp=10 tol sweep)
# ---------------------------------------------------------------------------------------------

function export_sensitivity()
    isfile(joinpath(REPO_ROOT, SYNTH_FILE)) || return
    df = CSV.read(joinpath(REPO_ROOT, SYNTH_FILE), DataFrame)
    for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
    refrows = df[df.n_rep_periods .== 1, :]
    isempty(refrows) && return
    ref_opt = mean(skipmissing(refrows.objective_value))
    # Exclude both cache and pgd_variant knockout-ablation arms (configs/synth_knockout.toml shares
    # this same output file): this sweep is about tolerance/normalization sensitivity for the
    # PROPOSED method only, not acceleration variants, which need far more iterations and would
    # otherwise skew the pooled iteration-count statistics below.
    sens = df[(df.n_rep_periods .== 10) .& .!occursin.("nocache", df.name) .&
              .!occursin.("_gram_only", df.name) .& .!occursin.("_plain", df.name), :]
    (isempty(sens) || length(unique(sens.tol)) < 2) && return
    sens.regret_pct = map(eachrow(sens)) do r
        (ismissing(r.evaluated_objective_value) || r.evaluated_objective_value <= 0) && return missing
        100 * (r.evaluated_objective_value / ref_opt - 1)
    end
    sens.mismatch_pct = map(eachrow(sens)) do r
        (ismissing(r.objective_value) || ismissing(r.evaluated_objective_value) ||
         r.evaluated_objective_value == 0) && return missing
        100 * abs(r.objective_value - r.evaluated_objective_value) / abs(r.evaluated_objective_value)
    end
    sens.total_time = reduce(.+, (coalesce.(getcol(sens, s), 0.0) for s in TIME_STAGES))

    tol_out = DataFrame()
    for nrm in ordered_norms(sens), cl in METHOD_ORDER, w in WEIGHT_ORDER
        mask = (sens.normalization .== nrm) .& (sens.clustering_type .== cl) .& (sens.weight_type .== w)
        g = sens[mask, :]; length(unique(g.tol)) < 2 && continue
        for t in sort(unique(g.tol))
            gt = g[isapprox.(g.tol, t; atol=1e-12), :]
            v = collect(skipmissing(gt.regret_pct)); isempty(v) && continue
            tt = collect(skipmissing(gt.total_time))
            push!(tol_out, (normalization = nrm, clustering_type = cl, weight_type = w, tol = t,
                            regret_min_pct = minimum(v), regret_median_pct = median(v),
                            regret_max_pct = maximum(v), regret_range_pct = maximum(v) - minimum(v),
                            median_time_s = isempty(tt) ? missing : median(tt)); cols=:union)
        end
    end
    write_csv("sensitivity_tol.csv", tol_out)

    iters_out = DataFrame()
    for t in sort(unique(sens.tol); rev=true)
        it = collect(skipmissing(sens[isapprox.(sens.tol, t; atol=1e-12), :pgd_total_iters]))
        isempty(it) && continue
        push!(iters_out, (tol = t, median_pgd_iters = median(it)); cols=:union)
    end
    write_csv("pgd_iters_by_tol.csv", iters_out)

    # Per-seed tol × normalization grid for PROPOSED (conical_hull/convex), with mismatch --
    # keeps the seed column (unlike sensitivity_tol.csv's min/median/max) so a paired t-test across
    # tol or across normalization is possible downstream. Also computes the two paired tests
    # (tol: loosest vs tightest; normalization: every pair present) and writes them to
    # sensitivity_stats.csv so the numbers a paper table cites are never hand-copied from a console.
    prop = sens[(sens.clustering_type .== PROPOSED[1]) .& (sens.weight_type .== PROPOSED[2]), :]
    by_seed = select(prop, :normalization, :tol, :seed, :regret_pct, :mismatch_pct)
    write_csv("sensitivity_by_seed.csv", by_seed)

    stats = DataFrame()
    tols = sort(unique(prop.tol))
    if length(tols) >= 2
        tlo, thi = minimum(tols), maximum(tols)
        for nrm in ordered_norms(prop)
            g = prop[prop.normalization .== nrm, :]
            glo = Dict(r.seed => r.regret_pct for r in eachrow(g[isapprox.(g.tol, tlo; atol=1e-15), :]) if !ismissing(r.regret_pct))
            ghi = Dict(r.seed => r.regret_pct for r in eachrow(g[isapprox.(g.tol, thi; atol=1e-15), :]) if !ismissing(r.regret_pct))
            ks = sort(collect(intersect(keys(glo), keys(ghi))))
            length(ks) < 2 && continue
            d = [ghi[k] - glo[k] for k in ks]   # tightest − loosest (negative ⇒ tighter tol helps)
            p = all(iszero, d) ? 1.0 : pvalue(OneSampleTTest(d))
            push!(stats, (comparison = "tol", normalization = nrm, a = thi, b = tlo,
                          mean_diff = mean(d), p_value = p, n_paired = length(ks)); cols=:union)
        end
    end
    norms_present = ordered_norms(prop)
    for i in eachindex(norms_present), j in (i+1):length(norms_present)
        r = normalization_compare_pair(prop, norms_present[i], norms_present[j])
        (r === nothing || r.n_total == 0) && continue
        push!(stats, (comparison = "normalization", normalization = "$(r.normA) vs $(r.normB)",
                      a = missing, b = missing, mean_diff = r.mean_diff, p_value = r.p_value,
                      n_paired = r.n_total, n_a_better = r.n_a_better, n_b_better = r.n_b_better,
                      n_comparable = r.n_comparable); cols=:union)
    end
    write_csv("sensitivity_stats.csv", stats)
end

# ---------------------------------------------------------------------------------------------
# 7. ablation_synth.csv — PROPOSED vs. single-component knockouts, on the synth dev system.
# ---------------------------------------------------------------------------------------------

function export_ablation()
    isfile(joinpath(REPO_ROOT, SYNTH_FILE)) || return
    df = CSV.read(joinpath(REPO_ROOT, SYNTH_FILE), DataFrame)
    for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
    arms = df[df.n_rep_periods .!= 1, :]
    refrows = df[df.n_rep_periods .== 1, :]
    isempty(refrows) && return
    ref_opt = mean(skipmissing(refrows.objective_value))
    arms.regret_pct = map(eachrow(arms)) do r
        (ismissing(r.evaluated_objective_value) || r.evaluated_objective_value <= 0) && return missing
        100 * (r.evaluated_objective_value / ref_opt - 1)
    end
    arms.total_time = reduce(.+, (coalesce.(getcol(arms, s), 0.0) for s in TIME_STAGES))
    arms.ablabel = [ablation_label(r.clustering_type, r.weight_type, r.normalization) for r in eachrow(arms)]
    sub = arms[.!ismissing.(arms.ablabel), :]
    isempty(sub) && return
    out = DataFrame()
    order = ["PROPOSED", "ablate: -economic", "ablate: -conic-selection", "ablate: -convex-weights"]
    for label in order, n in sort(unique(sub.n_rep_periods))
        rows = sub[(sub.ablabel .=== label) .& (sub.n_rep_periods .== n), :]
        isempty(rows) && continue
        (m_r, s_r) = meanstd(rows, :regret_pct); (m_t, s_t) = meanstd(rows, :total_time)
        sdr = isempty(s_r) ? 0.0 : parse(Float64, replace(s_r, "±" => ""))
        push!(out, (arm_label = label, n_rep_periods = n, n_seeds = nrow(rows),
                    regret_mean_pct = something(m_r, NaN), regret_sd_pct = sdr,
                    time_mean_s = something(m_t, NaN)); cols=:union)
    end
    write_csv("ablation_synth.csv", out)
end

# ---------------------------------------------------------------------------------------------
# 8. acceleration_ablation.csv — plain-PGD vs Gram-only vs Gram+FISTA+restart, synth dev system.
#    Isolates the Gram-matrix reformulation's contribution from FISTA's (see plan §1.1). The
#    "proposed" arm already exists in outputs/synth/hydro.csv (from the main ablation sweep); the
#    "gram_only"/"plain" arms come from configs/synth_acceleration.toml, sharing the same file
#    (their pgd_variant is encoded as a name suffix, per Types.jl's ExperimentData).
# ---------------------------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------------------------
# 8b. knockout_ablation.csv — cache / FISTA / Gram-matrix knockouts, crossed over multiple
#     hull x weight combinations (configs/synth_knockout.toml), unlike acceleration_ablation.csv's
#     single fixed combination. Four arms per combination: proposed (cache=true, pgd=proposed),
#     no-cache (cache=false, pgd=proposed), gram_only (cache=true, pgd=gram_only), plain
#     (cache=true, pgd=plain). Regret + mismatch + total runtime, mean ± sd over 5 seeds.
# ---------------------------------------------------------------------------------------------

const KNOCKOUT_COMBOS = [("convex_hull", "convex"), ("convex_hull", "conical"),
                          ("conical_hull", "convex"), ("conical_hull", "conical")]
const KNOCKOUT_ARMS = ["proposed", "no_cache", "gram_only", "plain"]

# Name-suffix order is base_n_periodlen_clustering_weight_tol[_norm][_pgd_variant][_nocache]
# (see ExperimentData's constructor in Types.jl) — "nocache" is always the last suffix when
# present, so it cannot collide with a pgd_variant suffix (our sweep never combines cache=false
# with a non-proposed pgd_variant, but the check order below is robust even if it did).
function knockout_arm_of(name)
    occursin("_nocache", name) && return "no_cache"
    occursin("_gram_only", name) && return "gram_only"
    occursin("_plain", name) && return "plain"
    return "proposed"
end

function export_knockout_ablation()
    isfile(joinpath(REPO_ROOT, SYNTH_FILE)) || return
    df = CSV.read(joinpath(REPO_ROOT, SYNTH_FILE), DataFrame)
    for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
    refrows = df[df.n_rep_periods .== 1, :]
    isempty(refrows) && return
    ref_opt = mean(skipmissing(refrows.objective_value))

    is_knockout_combo = [(cl, w) in KNOCKOUT_COMBOS for (cl, w) in zip(df.clustering_type, df.weight_type)]
    arms = df[(df.n_rep_periods .!= 1) .& (df.normalization .== "economic") .&
              isapprox.(df.tol, 0.01; atol=1e-12) .& is_knockout_combo, :]
    arms = copy(arms)
    arms.arm = knockout_arm_of.(arms.name)
    isempty(arms) && (@warn "no knockout-ablation rows found in $SYNTH_FILE — run configs/synth_knockout.toml first"; return)
    arms.regret_pct = map(eachrow(arms)) do r
        (ismissing(r.evaluated_objective_value) || r.evaluated_objective_value <= 0) && return missing
        100 * (r.evaluated_objective_value / ref_opt - 1)
    end
    arms.mismatch_pct = map(eachrow(arms)) do r
        (ismissing(r.objective_value) || ismissing(r.evaluated_objective_value) ||
         r.evaluated_objective_value == 0) && return missing
        100 * abs(r.objective_value - r.evaluated_objective_value) / abs(r.evaluated_objective_value)
    end
    arms.total_time = reduce(.+, (coalesce.(getcol(arms, s), 0.0) for s in TIME_STAGES))

    out = DataFrame()
    for (cl, w) in KNOCKOUT_COMBOS, arm in KNOCKOUT_ARMS, n in sort(unique(arms.n_rep_periods))
        sub = arms[(arms.clustering_type .== cl) .& (arms.weight_type .== w) .&
                   (arms.arm .== arm) .& (arms.n_rep_periods .== n), :]
        isempty(sub) && continue
        regret = collect(skipmissing(sub.regret_pct))
        mismatch = collect(skipmissing(sub.mismatch_pct))
        times = collect(skipmissing(sub.total_time))
        iters = collect(skipmissing(getcol(sub, :pgd_total_iters)))
        push!(out, (clustering_type = cl, weight_type = w, arm = arm, n_rep_periods = n,
                    n_seeds = nrow(sub),
                    regret_mean_pct = isempty(regret) ? missing : mean(regret),
                    regret_sd_pct = length(regret) > 1 ? std(regret) : 0.0,
                    mismatch_mean_pct = isempty(mismatch) ? missing : mean(mismatch),
                    mismatch_sd_pct = length(mismatch) > 1 ? std(mismatch) : 0.0,
                    time_mean_s = isempty(times) ? missing : mean(times),
                    time_sd_s = length(times) > 1 ? std(times) : 0.0,
                    iters_mean = isempty(iters) ? missing : mean(iters),
                    iters_sd = length(iters) > 1 ? std(iters) : 0.0); cols=:union)
    end
    write_csv("knockout_ablation.csv", out)
end

# ---------------------------------------------------------------------------------------------
# 9. curtailment.csv / capacity_decisions.csv (flagged low-confidence) — secondary metrics
# ---------------------------------------------------------------------------------------------

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

function df_ref_name(path)
    df = CSV.read(joinpath(REPO_ROOT, path), DataFrame)
    df[df.n_rep_periods .== 1, :name]
end

# ---------------------------------------------------------------------------------------------
# 9b. curtailment_renewable.csv — TRUE renewable curtailment (availability minus dispatch in the
#     full-resolution EVALUATION model), as a % of total available renewable energy. Distinct from
#     curtailment.csv's total_spillage/total_borrow (which are seasonal-storage water-balance
#     slacks, not generic renewable curtailment) -- this is what R5.11 actually asks for. Scoped to
#     the two dispatch case studies (p2x, 118bus), where capacity is fixed (not investable), so a
#     single input-side availability/capacity read applies to every arm.
# ---------------------------------------------------------------------------------------------

const CURTAILMENT_DATASETS = [("tyndp/p2x", "p2x"), ("nrel/118bus", "118bus")]

function renewable_availability(ds)
    base = joinpath(REPO_ROOT, "inputs", ds)
    assets = CSV.read(joinpath(base, "assets.csv"), DataFrame)
    cap = Dict(string(r.id) => Float64(r.unit_capacity) * Float64(r.initial_units) for r in eachrow(assets))
    aprof = CSV.read(joinpath(base, "assets-profiles.csv"), DataFrame)
    aprof = aprof[string.(aprof.profile_type) .== "availability", :]
    profiles = CSV.read(joinpath(base, "profiles.csv"), DataFrame)
    avail = Dict{String,Float64}()
    for r in eachrow(aprof)
        asset = string(r.asset); col = Symbol(string(r.profile))
        (haskey(cap, asset) && cap[asset] > 0 && hasproperty(profiles, col)) || continue
        avail[asset] = sum(skipmissing(profiles[!, col])) * cap[asset]   # hourly data, timestep_duration=1h
    end
    return avail
end

function export_curtailment_renewable()
    out = DataFrame()
    for (ds, case_study) in CURTAILMENT_DATASETS
        avail = renewable_availability(ds)
        total_avail = sum(values(avail))
        total_avail <= 0 && (@warn "no renewable availability found for $ds"; continue)
        ids = Set(keys(avail))
        dir = joinpath(REPO_ROOT, "outputs", split(ds, "/")[1])
        isdir(dir) || continue
        entries = readdir(dir)
        dsname = split(ds, "/")[2]
        for (cl, w) in [PROPOSED, ("k_means", "dirac"), ("k_medoids", "dirac")], n in (10, 20, 40, 80)
            pat = Regex("^" * dsname * "_$(n)_24_$(cl)_$(w)_0\\.01_economic\$")
            idx = findfirst(e -> occursin(pat, e), entries)
            idx === nothing && continue
            f = joinpath(dir, entries[idx], "eval_power_out.arrow")
            isfile(f) || continue
            d = DataFrame(Arrow.Table(f))
            d = d[in.(d.id, Ref(ids)), :]
            isempty(d) && continue
            by_seed = combine(groupby(d, :seed), :value => (v -> sum(skipmissing(v))) => :dispatched)
            curtailed_pct = [100 * (total_avail - r.dispatched) / total_avail for r in eachrow(by_seed)]
            push!(out, (case_study = case_study, clustering_type = cl, weight_type = w, n_rep_periods = n,
                        curtailment_mean_pct = mean(curtailed_pct),
                        curtailment_sd_pct = length(curtailed_pct) > 1 ? std(curtailed_pct) : 0.0,
                        n_seeds = length(curtailed_pct)); cols=:union)
        end
    end
    write_csv("curtailment_renewable.csv", out)
end

# ---------------------------------------------------------------------------------------------
# 9c. storage_tracking.csv — decision-fidelity metric for the STORAGE case studies (p2x, 118bus),
#     the analogue of capacity_decisions.csv for investment case studies: how closely does the
#     reduced model's reconstructed daily inter-period storage trajectory track the full-horizon
#     (n_rep=1) model's actual trajectory? Reported as the mean absolute deviation across every
#     seasonal asset and every day of the year, each expressed as a % of THAT asset's own storage
#     capacity (so heterogeneous reservoir sizes don't swamp one another), then averaged.
# ---------------------------------------------------------------------------------------------

function asset_capacities(ds)
    path = joinpath(REPO_ROOT, "inputs", ds, "assets-storage.csv")
    isfile(path) || return Dict{String,Float64}()
    df = CSV.read(path, DataFrame)
    Dict(string(r.id) => Float64(r.capacity_storage_energy) for r in eachrow(df))
end

function export_storage_tracking()
    out = DataFrame()
    for (ds, case_study) in CURTAILMENT_DATASETS   # same two storage case studies
        cap = asset_capacities(ds)
        isempty(cap) && continue
        dir = joinpath(REPO_ROOT, "outputs", split(ds, "/")[1])
        isdir(dir) || continue
        entries = readdir(dir)
        dsname = split(ds, "/")[2]
        idx0 = findfirst(e -> occursin(Regex("^" * dsname * "_1_8760_"), e), entries)
        idx0 === nothing && (@warn "no n_rep=1 arm dir for $ds storage tracking"; continue)
        basefile = joinpath(dir, entries[idx0], "reduced_state_of_charge_intra.arrow")
        isfile(basefile) || continue
        base = DataFrame(Arrow.Table(basefile))
        base = base[base.timestep .% 24 .== 0, :]
        base.day = base.timestep .÷ 24
        # full-model daily trajectory per (asset, day), averaged over seeds first (decisions should
        # be seed-invariant for the n_rep=1 reference; average handles any residual solver noise).
        fullavg = combine(groupby(base, [:id, :day]), :value => (v -> mean(skipmissing(v))) => :full)
        full = Dict((string(r.id), Int(r.day)) => r.full for r in eachrow(fullavg))
        for (cl, w) in [PROPOSED, ("k_means", "dirac"), ("k_medoids", "dirac")], n in (10, 20, 40, 80)
            pat = Regex("^" * dsname * "_$(n)_24_$(cl)_$(w)_0\\.01_economic\$")
            idx = findfirst(e -> occursin(pat, e), entries)
            idx === nothing && continue
            f = joinpath(dir, entries[idx], "reduced_state_of_charge_inter.arrow")
            isfile(f) || continue
            d = DataFrame(Arrow.Table(f))
            davg = combine(groupby(d, [:id, :period]), :value => (v -> mean(skipmissing(v))) => :reduced)
            relsum = 0.0; n_obs = 0
            for r in eachrow(davg)
                key = (string(r.id), Int(r.period)); haskey(full, key) || continue
                c = get(cap, string(r.id), 0.0); c <= 0 && continue
                relsum += abs(r.reduced - full[key]) / c
                n_obs += 1
            end
            n_obs == 0 && continue
            push!(out, (case_study = case_study, clustering_type = cl, weight_type = w, n_rep_periods = n,
                        mean_tracking_error_pct_of_capacity = 100 * relsum / n_obs, n_obs = n_obs); cols=:union)
        end
    end
    write_csv("storage_tracking.csv", out)
end

# ---------------------------------------------------------------------------------------------
# 10. runtime_breakdown.csv — all 5 timing stages (keeps model formulation as its own segment),
#     PROPOSED vs a k_means/dirac baseline, per case study, vs n_rp; plus the n_rp=1 reference time.
# ---------------------------------------------------------------------------------------------

const STAGE_COLS = [:time_to_preprocess, :time_to_cluster, :time_to_fit_weights,
                    :time_to_formulate_model, :time_to_solve]
const STAGE_NAMES = ["read_preprocess", "cluster", "fit_weights", "formulate_model", "solve"]

function export_runtime_breakdown()
    out = DataFrame()
    for (path, meta) in CASE_META
        df = CSV.read(joinpath(REPO_ROOT, path), DataFrame)
        for c in ("clustering_type", "weight_type", "normalization"); df[!, c] = string.(df[!, c]); end
        refrows = df[df.n_rep_periods .== 1, :]
        methods = [("PROPOSED", PROPOSED[1], PROPOSED[2]), ("k_means", "k_means", "dirac"),
                   ("k_medoids", "k_medoids", "dirac"), ("hierarchical", "hierarchical", "dirac")]
        for (label, cl, w) in methods
            for n in sort(unique(df.n_rep_periods))
                n == 1 && continue
                sub = df[(df.normalization .== "economic") .& isapprox.(df.tol, 0.01; atol=1e-12) .&
                         .!occursin.("nocache", df.name) .&
                         (df.clustering_type .== cl) .& (df.weight_type .== w) .& (df.n_rep_periods .== n), :]
                isempty(sub) && continue
                for (col, stage) in zip(STAGE_COLS, STAGE_NAMES)
                    v = collect(skipmissing(getcol(sub, col))); isempty(v) && continue
                    push!(out, (case_study = meta.case_study, method_label = label, n_rep_periods = n,
                                stage = stage, time_mean_s = mean(v),
                                time_sd_s = length(v) > 1 ? std(v) : 0.0); cols=:union)
                end
            end
        end
        # n_rep=1 full-model reference solve time (read_preprocess + solve; no clustering/fitting)
        if !isempty(refrows)
            rt = collect(skipmissing(getcol(refrows, :time_to_solve)))
            isempty(rt) || push!(out, (case_study = meta.case_study, method_label = "full_reference",
                                       n_rep_periods = 1, stage = "solve", time_mean_s = mean(rt),
                                       time_sd_s = length(rt) > 1 ? std(rt) : 0.0); cols=:union)
        end
    end
    write_csv("runtime_breakdown.csv", out)
end

# ---------------------------------------------------------------------------------------------
# 11. reservoir_trajectory_p2x.csv
# ---------------------------------------------------------------------------------------------

function export_reservoir_trajectory(; top_n=1)
    dir = joinpath(REPO_ROOT, "outputs", "tyndp")
    isdir(dir) || return
    # find the n_rep=1 baseline arm directory to identify the largest reservoir(s)
    entries = readdir(dir)
    base_dir = findfirst(e -> occursin(r"^p2x_1_8760_", e), entries)
    base_dir === nothing && (@warn "no p2x n_rep=1 arm directory found for reservoir trajectory"; return)
    basepath = joinpath(dir, entries[base_dir], "reduced_state_of_charge_intra.arrow")
    isfile(basepath) || (@warn "missing $basepath"; return)
    base = DataFrame(Arrow.Table(basepath))
    base = base[base.timestep .% 24 .== 0, :]
    base.day = base.timestep .÷ 24
    peak = combine(groupby(base, :id), :value => (v -> mean(skipmissing(v))) => :m)
    top_ids = first(sort(peak, :m, rev=true), top_n).id

    out = DataFrame()
    for id in top_ids
        b = base[base.id .== id, :]
        agg = combine(groupby(b, :day), :value => (v -> mean(skipmissing(v))) => :m,
                                        :value => (v -> length(v) > 1 ? std(skipmissing(v)) : 0.0) => :s)
        for r in eachrow(agg)
            push!(out, (reservoir_id = id, day_of_year = r.day, source_label = "full_reference",
                        clustering_type = missing, weight_type = missing, n_rep_periods = 1,
                        storage_mwh_mean = r.m, storage_mwh_sd = r.s); cols=:union)
        end
    end

    for (cl, w) in [("k_means", "dirac"), ("k_medoids", "dirac"), PROPOSED]
        for n in (10, 20, 40, 80)
            pattern = Regex("^p2x_$(n)_24_$(cl)_$(w)_0\\.01_economic\$")
            idx = findfirst(e -> occursin(pattern, e), entries)
            idx === nothing && continue
            f = joinpath(dir, entries[idx], "reduced_state_of_charge_inter.arrow")
            isfile(f) || continue
            d = DataFrame(Arrow.Table(f))
            for id in top_ids
                di = d[d.id .== id, :]
                isempty(di) && continue
                agg = combine(groupby(di, :period), :value => (v -> mean(skipmissing(v))) => :m,
                                                     :value => (v -> length(v) > 1 ? std(skipmissing(v)) : 0.0) => :s)
                for r in eachrow(agg)
                    push!(out, (reservoir_id = id, day_of_year = r.period, source_label = "reduced",
                                clustering_type = cl, weight_type = w, n_rep_periods = n,
                                storage_mwh_mean = r.m, storage_mwh_sd = r.s); cols=:union)
                end
            end
        end
    end
    write_csv("reservoir_trajectory_p2x.csv", out)
end

# ---------------------------------------------------------------------------------------------
# 12. knob_sensitivity.csv — data-property study (per-axis regret spread per knob variant)
# ---------------------------------------------------------------------------------------------

# Re-implemented (not `include`d) because summarize_knobs.jl calls `main()` unconditionally at
# file scope — including it here would execute its CLI entry point with the wrong ARGS/cwd.
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
                    knob_value = knob_suffix_num(variant), weight_spread_pp = r.weight,
                    selection_spread_pp = r.selection, normalization_spread_pp = r.norm,
                    overall_lo_pct = r.lo, overall_hi_pct = r.hi); cols=:union)
    end
    write_csv("knob_sensitivity.csv", out)
end

# ---------------------------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------------------------

function export_provenance()
    envfile = joinpath(REPO_ROOT, "outputs", "environment.txt")
    content = isfile(envfile) ? read(envfile, String) : "(no outputs/environment.txt found)"
    commit = try strip(read(`git -C $REPO_ROOT rev-parse HEAD`, String)) catch; "unknown" end
    open(joinpath(OUTDIR, "PROVENANCE.md"), "w") do io
        println(io, "# Provenance")
        println(io, "\nGenerated by `analysis/export_summary_csvs.jl`.")
        println(io, "\nexperiments repo commit: `$commit`")
        println(io, "\n## outputs/environment.txt\n\n```\n$content\n```")
    end
    println("  wrote $(joinpath(OUTDIR, "PROVENANCE.md"))")
end

function main()
    println("Exporting summary CSVs to $OUTDIR")
    export_regret()
    export_model_size()
    export_cache_hit_rate()
    export_sensitivity()
    export_ablation()
    export_acceleration_ablation()
    export_knockout_ablation()
    export_secondary()
    export_curtailment_renewable()
    export_storage_tracking()
    export_runtime_breakdown()
    export_reservoir_trajectory()
    export_knob_sensitivity()
    export_provenance()
    println("Done.")
end

main()
