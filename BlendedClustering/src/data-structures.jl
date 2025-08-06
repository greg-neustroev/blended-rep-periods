export ExperimentData, RunData, ExperimentResult


"""
Data needed to run all of the experiments
"""
struct RunData <: AbstractDataFrame
    function RunData(path)
        df = CSV.read(path, DataFrame)
        # Ensure the DataFrame has the correct column names
        required_columns = [
            "name",
            "n_rep_periods",
            "period_length",
            "clustering_type",
            "distance",
            "weight_type",
            "niters",
            "learning_rate",
            "evaluation_type",
        ] |> Set
        df_columns = df |> names |> Set
        if !issubset(required_columns, df_columns)
            error("Input CSV file is missing required columns: $required_columns")
        end
        symbol_columns = [:clustering_type, :weight_type, :evaluation_type]
        transform!(df, symbol_columns .=> ByRow(Symbol) .=> symbol_columns)
        transform!(df, :distance => ByRow(string_to_semimetric) => :distance)
        return df
    end
end

function string_to_semimetric(s::AbstractString)
    if s == "Euclidean"
        return Euclidean()
    elseif s == "cityblock" || s == "manhattan"
        return Cityblock()
    elseif s == "cosine"
        return CosineDist()
    elseif s == "chebyshev"
        return Chebyshev()
    else
        error("Unknown distance metric: $s")
    end
end

"""
Data needed to run a single experiment (i.e., a single optimization model)
"""
struct ExperimentData
    name::String

    # Clustering
    n_rep_periods::Int
    period_length::Int
    clustering_type::Symbol
    distance::SemiMetric
    weight_type::Symbol
    niters::Int
    learning_rate::Float64
    evaluation_type::Symbol

    function ExperimentData(run_data_row::DataFrameRow{DataFrame,DataFrames.Index}; database::AbstractString=":memory:")
        return new(
            run_data_row.name,
            run_data_row.n_rep_periods,
            run_data_row.period_length,
            run_data_row.clustering_type,
            run_data_row.distance,
            run_data_row.weight_type,
            run_data_row.niters,
            run_data_row.learning_rate,
            run_data_row.evaluation_type
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
    name::String
    seed::Int
    termination_status::String
    objective_value::Float64
    time_to_preprocess::Float64
    time_to_cluster::Float64
    time_to_fit_weights::Float64
    time_to_formulate_model::Float64
    time_to_solve::Float64

    function ExperimentResult(
        name::String,
        seed::Int,
        solved_model::JuMP.AbstractModel,
        time_to_preprocess::Float64,
        time_to_cluster::Float64,
        time_to_fit_weights::Float64,
        time_to_formulate_model::Float64,
        time_to_solve::Float64,
    )
        if JuMP.is_solved_and_feasible(solved_model)
            return new(
                name,
                seed,
                solved_model |> termination_status |> string,
                solved_model |> objective_value,
                time_to_preprocess,
                time_to_cluster,
                time_to_fit_weights,
                time_to_formulate_model,
                time_to_solve,
            )
        else
            return new(
                name,
                seed,
                solved_model |> termination_status |> string,
                Float64.Inf,
                time_to_read,
                time_to_preprocess,
                time_to_cluster,
                time_to_fit_weights,
                time_to_formulate_model,
                time_to_solve,
            )
        end

    end
end
