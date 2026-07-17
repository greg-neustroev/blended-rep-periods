#!/usr/bin/env julia
#
# case_studies.tex §Additional Evaluation Metrics: TRUE renewable curtailment (availability minus
# dispatch in the full-resolution EVALUATION model), as a % of total available renewable energy.
# Distinct from secondary.jl's curtailment.csv (total_spillage/total_borrow — seasonal-storage
# water-balance slacks, not generic renewable curtailment). Scoped to the two dispatch case
# studies (p2x, 118bus), where capacity is fixed (not investable), so a single input-side
# availability/capacity read applies to every arm.
#
# Usage: julia --project=. analysis/paper/curtailment_renewable.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

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

abspath(PROGRAM_FILE) == (@__FILE__) && export_curtailment_renewable()
