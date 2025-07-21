export ExperimentData, RunData, ExperimentResult


"""
Data needed to run all of the experiments
"""
struct RunData
    # Sets
    times::AbstractVector{Int}
    locations::AbstractVector{Symbol}
    transmission_lines::AbstractVector{Tuple{Symbol,Symbol}}
    generators::AbstractVector{Tuple{Symbol,Symbol}}
    generation_technologies::AbstractVector{Symbol}

    # Dataframes
    demand::AbstractDataFrame
    generation_availability::AbstractDataFrame
    generation::AbstractDataFrame
    transmission_capacities::AbstractDataFrame

    # Scalars
    value_of_lost_load::Float64
    representative_period_duration::Int

    # Clustering
    clusters::AbstractVector{Int}
    distances::AbstractVector{SemiMetric}
    weight_types::AbstractVector{Symbol}
    scalings::AbstractVector{Bool}
    clustering_types::AbstractVector{Symbol}

    function RunData(config_dict::Dict{Symbol,Any})
        sets = config_dict[:sets]
        clustering = config_dict[:clustering]
        data = config_dict[:data]
        scalars = data[:scalars]

        return new(
            sets[:times],
            sets[:locations],
            sets[:transmission_lines],
            sets[:generators],
            sets[:generation_technologies],
            data[:demand],
            data[:generation_availability],
            data[:generation],
            data[:transmission_lines],
            scalars[:value_of_lost_load],
            scalars[:representative_period_duration],
            clustering[:clusters],
            clustering[:distances],
            clustering[:weights],
            clustering[:scalings],
            clustering[:types],
        )
    end
end


"""
Data needed to run a single experiment (i.e., a single optimization model)
"""
struct ExperimentData
    # Sets
    times::Vector{Int}
    locations::Vector{Symbol}
    transmission_lines::Vector{Tuple{Symbol,Symbol}}
    generators::Vector{Tuple{Symbol,Symbol}}
    generation_technologies::Vector{Symbol}

    # Dataframes
    demand::AbstractDataFrame
    generation_availability::AbstractDataFrame
    generation::AbstractDataFrame
    transmission_capacities::AbstractDataFrame

    # Scalars
    value_of_lost_load::Float64
    representative_period_duration::Int

    # Clustering
    n_clusters::Int
    distance::SemiMetric
    weight_type::Symbol
    scaling::Bool
    clustering_type::Symbol

    function ExperimentData(run_data::RunData, n_clusters::Int, distance::SemiMetric, weight_type::Symbol, scaling::Bool, clustering_type::Symbol)
        return new(
            run_data.times,
            run_data.locations,
            run_data.transmission_lines,
            run_data.generators,
            run_data.generation_technologies,
            run_data.demand,
            run_data.generation_availability,
            run_data.generation,
            run_data.transmission_capacities,
            run_data.value_of_lost_load,
            run_data.representative_period_duration,
            n_clusters,
            distance,
            weight_type,
            scaling,
            clustering_type,
        )
    end
end

"""
Structure to hold the time series used in clustering together with some
summary statistics on the data.
"""
mutable struct AuxiliaryClusteringData
    key_columns::AbstractVector{Symbol}
    period_duration::Int
    last_period_duration::Int
    n_periods::Int

    function AuxiliaryClusteringData(key_columns, period_duration, last_period_duration, n_periods)
        return new(key_columns, period_duration, last_period_duration, n_periods)
    end
end

"""
Structure to hold the clustering result.
"""
mutable struct ClusteringResult
    profiles::AbstractDataFrame
    weight_matrix::Union{SparseMatrixCSC{Float64,Int64},Matrix{Float64}}
    clustering_matrix::Union{Matrix{Float64},Nothing}
    rp_matrix::Union{Matrix{Float64},Nothing}

    function ClusteringResult(profiles, weight_matrix, clustering_matrix, rp_matrix)
        return new(profiles, weight_matrix, clustering_matrix, rp_matrix)
    end

    function ClusteringResult(profiles, weight_matrix)
        return new(profiles, weight_matrix, nothing, nothing)
    end
end

struct ExperimentResult
    objective::Float64
    eval_objective::Float64
    runtime::Float64

    function ExperimentResult(objective::Float64, eval_objective::Float64, runtime::Float64)
        return new(objective, eval_objective, runtime)
    end
end
