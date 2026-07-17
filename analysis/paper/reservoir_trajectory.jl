#!/usr/bin/env julia
#
# reservoir_trajectory_p2x.csv — feeds only the currently-unused plots/ (candidate reservoir
# trajectory figure, not `\includegraphics`'d anywhere in the compiled paper). Migrated as-is and
# not otherwise re-verified in this pass; revisit if a future draft adds a figure that uses it.
#
# Usage: julia --project=. analysis/paper/reservoir_trajectory.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

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

abspath(PROGRAM_FILE) == (@__FILE__) && export_reservoir_trajectory()
