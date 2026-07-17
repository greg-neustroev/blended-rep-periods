#!/usr/bin/env julia
#
# case_studies.tex §Performance of the Projection Cache: hit rate per real case study x hull type;
# cached-vs-uncached wall time only exists for the synth dev system (the only dataset with matched
# _nocache twins).
#
# Usage: julia --project=. analysis/paper/cache_hit_rate.jl [outdir]

isdefined(Main, :write_csv) || include(joinpath(@__DIR__, "common.jl"))

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

abspath(PROGRAM_FILE) == (@__FILE__) && export_cache_hit_rate()
