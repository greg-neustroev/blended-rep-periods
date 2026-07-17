#!/usr/bin/env julia
#
# case_studies.tex §Additional Evaluation Metrics: storage-tracking fidelity (reservoir-trajectory
# gap, as a % of capacity) for the STORAGE case studies (p2x, 118bus) — the analogue of
# secondary.jl's capacity_decisions.csv for investment case studies. How closely does the reduced
# model's reconstructed daily inter-period storage trajectory track the full-horizon (n_rep=1)
# model's actual trajectory? Reported as the mean absolute deviation across every seasonal asset
# and every day of the year, each expressed as a % of THAT asset's own storage capacity (so
# heterogeneous reservoir sizes don't swamp one another), then averaged.
#
# Usage: julia --project=. analysis/paper/storage_tracking.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

const STORAGE_DATASETS = [("tyndp/p2x", "p2x"), ("nrel/118bus", "118bus")]

function asset_capacities(ds)
    path = joinpath(REPO_ROOT, "inputs", ds, "assets-storage.csv")
    isfile(path) || return Dict{String,Float64}()
    df = CSV.read(path, DataFrame)
    Dict(string(r.id) => Float64(r.capacity_storage_energy) for r in eachrow(df))
end

function export_storage_tracking()
    out = DataFrame()
    for (ds, case_study) in STORAGE_DATASETS
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
        # "regardless of clustering method, weight type, or n_rep" (case_studies.tex) needs the
        # FULL grid, not just PROPOSED + 2 baselines.
        for (cl, w) in [(cl, w) for cl in METHOD_ORDER for w in WEIGHT_ORDER], n in (10, 20, 40, 80)
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

abspath(PROGRAM_FILE) == (@__FILE__) && export_storage_tracking()
