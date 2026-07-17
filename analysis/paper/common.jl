#!/usr/bin/env julia
#
# Shared helpers for the paper's data-export scripts: each script under analysis/paper/ is Stage 1
# of the paper's two-stage data pipeline, reducing the raw experiment output CSVs/Arrow dumps under
# outputs/ into tidy CSVs written into the PAPER repo's data/ directory, shaped to match the .tex
# table/prose paragraph they back (see the header comment of each script for its target). Stage 2
# (R scripts under analysis/R/, or hand-transcription into a table) reads ONLY those CSVs, so
# regenerating a figure never requires re-running this pipeline or having a solver license.
#
# `include`s analysis/common.jl for the regret/mismatch/statistics primitives shared with the
# console diagnostic tool (summarize.jl), and adds the export-pipeline-specific pieces: output
# path handling, case-study metadata, and the raw-CSV/Arrow loading helpers every export script
# needs.
#
# Not meant to be run directly — `include`d by every analysis/paper/*.jl script.

isdefined(Main, :METHOD_ORDER) || include(joinpath(@__DIR__, "..", "common.jl"))

const OUTDIR = length(ARGS) >= 1 ? ARGS[1] : normpath(joinpath(REPO_ROOT, "..", "clustering", "data"))
mkpath(OUTDIR)

# Case study -> (case_study id, problem_class, source). RTS deliberately excluded here — the
# campaign hasn't finished it yet; add a row once outputs/gridmod/rts.csv is complete, per the plan.
const CASE_META = [
    ("outputs/tyndp/gep.csv",   (case_study="gep",    problem_class="investment", source="tyndp")),
    ("outputs/sienna/5bus.csv", (case_study="5bus",   problem_class="investment", source="sienna")),
    ("outputs/tyndp/p2x.csv",   (case_study="p2x",    problem_class="dispatch",   source="tyndp")),
    ("outputs/nrel/118bus.csv", (case_study="118bus", problem_class="dispatch",   source="nrel")),
]
const SYNTH_FILE = "outputs/synth/hydro.csv"

# reference optimum (mean objective_value over n_rep_periods==1 rows) and the arms (n_rep>1),
# with the derived columns every exporter needs: regret_pct, mismatch_pct, total_time (stage sum,
# matching TIME_STAGES), cache_hit_rate, plus string-coerced category columns.
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

# arm directory lookup shared by every exporter that reads Arrow dumps directly (curtailment,
# storage tracking, feasibility slack, reservoir trajectory): the canonical (tol=0.01, economic,
# cached) run-directory name for a given (dataset, clustering, weight, n_rep).
function arm_dir(dataset_dir, dsname, cl, w, n)
    isdir(dataset_dir) || return nothing
    entries = readdir(dataset_dir)
    pat = Regex("^" * dsname * "_$(n)_24_$(cl)_$(w)_0\\.01_economic\$")
    idx = findfirst(e -> occursin(pat, e), entries)
    idx === nothing ? nothing : joinpath(dataset_dir, entries[idx])
end

# ids of assets whose `technology` column equals `tech` (e.g. the virtual "ENS" loss-of-load
# generators, excluded from Table tab:case-studies' |G| count and summed for loss-of-load).
function asset_ids_by_technology(ds, tech)
    path = joinpath(REPO_ROOT, "inputs", ds, "assets.csv")
    isfile(path) || return Set{String}()
    df = CSV.read(path, DataFrame)
    Set(string.(df.id[string.(df.technology) .== tech]))
end

# scalar value from a dataset's scalars.csv (e.g. "timestep_duration"); `default` when absent.
function scalar_value(ds, name, default)
    path = joinpath(REPO_ROOT, "inputs", ds, "scalars.csv")
    isfile(path) || return default
    df = CSV.read(path, DataFrame)
    r = df[string.(df.scalar) .== name, :value]
    isempty(r) ? default : Float64(first(r))
end

# total annual energy (MWh) dispatched against a set of demand ids: demand ids are columns of
# profiles.csv directly (no assets-profiles join, unlike generation availability), each scaled by
# its own peak value in demand.csv's `peak_demand` column and the dataset's timestep duration.
function total_demand_mwh(ds)
    base = joinpath(REPO_ROOT, "inputs", ds)
    dem = CSV.read(joinpath(base, "demand.csv"), DataFrame)
    profiles = CSV.read(joinpath(base, "profiles.csv"), DataFrame)
    tau = scalar_value(ds, "timestep_duration", 1.0)
    tot = 0.0
    for r in eachrow(dem)
        col = Symbol(string(r.id))
        hasproperty(profiles, col) || continue
        tot += sum(skipmissing(profiles[!, col])) * Float64(r.peak_demand) * tau
    end
    return tot
end
