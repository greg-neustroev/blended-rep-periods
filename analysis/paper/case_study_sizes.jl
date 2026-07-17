#!/usr/bin/env julia
#
# Table tab:case-studies (case_studies.tex): the structural columns |N|,|G|,|S^ST|,|S^seas|,|C|,
# |L|,|A^inv|,n_T, regret type, and the "max buff. (days)" column — everything except the
# #var./#constr. columns (model_size.jl). One row per case study, matching the table's own row
# order and column layout directly.
#
# |G| excludes the virtual "ENS" loss-of-load generators, per the table's own caption. Asset class
# (generation/storage/conversion) comes from joining assets.csv's `technology` column to
# technologies.csv's `type`; short-term vs. seasonal storage is assets-storage.csv's row count
# minus assets-storage-seasonal.csv's (the seasonal-only subset carries extra fields — peak
# inflow, spillage/borrow cost — so its id set is a subset of assets-storage.csv's).
#
# "max buff. (days)" — the reservoir buffer described in the table caption as "capacity in days of
# mean inflow" — is computed per seasonal-storage asset as capacity_storage_energy / (mean hourly
# inflow x 24h), reported as the min-max range across that case study's seasonal reservoirs (a
# single reservoir collapses min==max, printed as one number by whatever consumes this CSV). This
# is a best-effort reproduction of a descriptive characterization column, not a statistical claim
# — cross-check the synthetic-dataset row by hand (it should read as the 2-120 day range spanned by
# the reservoir-buffer knob sweep, tab:knobs, not just the 60-day default configuration).
#
# Usage: julia --project=. analysis/paper/case_study_sizes.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

const CASE_INPUT_DIRS = [
    ("synth/hydro", "synthetic", "storage"),
    ("sienna/5bus", "5bus", "investment"),
    ("tyndp/gep", "gep", "investment"),
    ("tyndp/p2x", "p2x", "storage"),
    ("nrel/118bus", "118bus", "storage"),
]

function asset_class_counts(ds)
    base = joinpath(REPO_ROOT, "inputs", ds)
    assets = CSV.read(joinpath(base, "assets.csv"), DataFrame)
    techs = CSV.read(joinpath(base, "technologies.csv"), DataFrame)
    type_of = Dict(string(r.id) => string(r.type) for r in eachrow(techs))
    non_ens = assets[string.(assets.technology) .!= "ENS", :]
    n_g = count(r -> get(type_of, string(r.technology), "") == "generation", eachrow(non_ens))
    n_c = count(r -> get(type_of, string(r.technology), "") == "conversion", eachrow(non_ens))
    return (g = n_g, c = n_c)
end

function storage_counts(ds)
    base = joinpath(REPO_ROOT, "inputs", ds)
    stpath = joinpath(base, "assets-storage.csv")
    isfile(stpath) || return (short = 0, seasonal = 0)
    st = CSV.read(stpath, DataFrame)
    seaspath = joinpath(base, "assets-storage-seasonal.csv")
    n_seasonal = isfile(seaspath) ? nrow(CSV.read(seaspath, DataFrame)) : 0
    (short = nrow(st) - n_seasonal, seasonal = n_seasonal)
end

function n_nodes(ds)
    base = joinpath(REPO_ROOT, "inputs", ds)
    locs = Set{String}()
    assets = CSV.read(joinpath(base, "assets.csv"), DataFrame)
    union!(locs, string.(assets.location))
    dpath = joinpath(base, "demand.csv")
    if isfile(dpath)
        dem = CSV.read(dpath, DataFrame)
        union!(locs, string.(dem.location))
    end
    length(locs)
end

function n_lines(ds)
    path = joinpath(REPO_ROOT, "inputs", ds, "transmission-lines.csv")
    isfile(path) ? nrow(CSV.read(path, DataFrame)) : 0
end

function n_investable(ds)
    path = joinpath(REPO_ROOT, "inputs", ds, "investments.csv")
    isfile(path) ? nrow(CSV.read(path, DataFrame)) : 0
end

# reservoir buffer, in days of mean inflow, for every seasonal-storage asset: capacity / (mean
# hourly inflow rate x 24h). Mean hourly inflow rate = mean(inflow profile fraction) x peak_inflow,
# the same peak-times-profile convention used elsewhere (e.g. curtailment_renewable.jl's
# renewable_availability).
function reservoir_buffer_days(ds)
    base = joinpath(REPO_ROOT, ds)
    seaspath = joinpath(base, "assets-storage-seasonal.csv")
    isfile(seaspath) || return Float64[]
    seas = CSV.read(seaspath, DataFrame)
    cap = Dict(string(r.id) => Float64(r.capacity_storage_energy)
               for r in eachrow(CSV.read(joinpath(base, "assets-storage.csv"), DataFrame)))
    aprof = CSV.read(joinpath(base, "assets-profiles.csv"), DataFrame)
    aprof = aprof[string.(aprof.profile_type) .== "inflows", :]
    profiles = CSV.read(joinpath(base, "profiles.csv"), DataFrame)
    days = Float64[]
    for r in eachrow(seas)
        id = string(r.id)
        haskey(cap, id) || continue
        prow = aprof[string.(aprof.asset) .== id, :]
        isempty(prow) && continue
        col = Symbol(string(first(prow).profile))
        hasproperty(profiles, col) || continue
        mean_frac = mean(skipmissing(profiles[!, col]))
        mean_hourly_inflow = mean_frac * Float64(r.peak_inflow)
        mean_hourly_inflow <= 0 && continue
        push!(days, cap[id] / (mean_hourly_inflow * 24))
    end
    days
end

function export_case_study_sizes()
    out = DataFrame()
    for (ds, case_study, regret_type) in CASE_INPUT_DIRS
        gc = asset_class_counts(ds)
        st = storage_counts(ds)
        buf = reservoir_buffer_days(joinpath("inputs", ds))
        # nT (period_length at n_rep=1), read off the matching outputs/ CSV when available (falls
        # back to `missing` if that case study hasn't produced output yet).
        meta_match = findfirst(m -> m.case_study == case_study, last.(CASE_META))
        out_csv = meta_match === nothing ? nothing : first(CASE_META[meta_match])
        out_csv === nothing && case_study == "synthetic" && (out_csv = SYNTH_FILE)
        n_t = missing
        if out_csv !== nothing && isfile(joinpath(REPO_ROOT, out_csv))
            df = CSV.read(joinpath(REPO_ROOT, out_csv), DataFrame)
            refrows = df[df.n_rep_periods .== 1, :]
            isempty(refrows) || (n_t = first(refrows.period_length))
        end
        push!(out, (
            case_study = case_study, n_nodes = n_nodes(ds), n_generation = gc.g,
            n_storage_short_term = st.short, n_storage_seasonal = st.seasonal,
            n_conversion = gc.c, n_lines = n_lines(ds), n_investable = n_investable(ds),
            n_timesteps = n_t, regret_type = regret_type,
            max_buffer_days_min = isempty(buf) ? missing : minimum(buf),
            max_buffer_days_max = isempty(buf) ? missing : maximum(buf),
        ); cols=:union)
    end
    write_csv("case_study_sizes.csv", out)
end

abspath(PROGRAM_FILE) == (@__FILE__) && export_case_study_sizes()
