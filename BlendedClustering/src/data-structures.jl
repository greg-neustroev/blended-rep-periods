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
            "input_dir",
            "n_rep_periods",
            "period_length",
            "clustering_type",
            "distance",
            "weight_type",
            "niters",
            "learning_rate"
        ] |> Set
        df_columns = df |> names |> Set
        if !issubset(required_columns, df_columns)
            error("Input CSV file is missing required columns: $required_columns")
        end
        base_dir = dirname(path)
        transform!(df, :input_dir => ByRow(s -> joinpath(base_dir, s)) => :input_dir)
        transform!(df, [:clustering_type, :weight_type] .=> ByRow(Symbol) .=> [:clustering_type, :weight_type])
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

    # Database connection
    db_connection::DuckDB.DB
    input_dir::String

    # Clustering
    n_rep_periods::Int
    period_length::Int
    clustering_type::Symbol
    distance::SemiMetric
    weight_type::Symbol
    niters::Int
    learning_rate::Float32

    function ExperimentData(run_data_row::DataFrameRow{DataFrame,DataFrames.Index}; database::AbstractString=":memory:")
        return new(
            run_data_row.name,
            DBInterface.connect(DuckDB.DB, database),
            run_data_row.input_dir,
            run_data_row.n_rep_periods,
            run_data_row.period_length,
            run_data_row.clustering_type,
            run_data_row.distance,
            run_data_row.weight_type,
            run_data_row.niters,
            run_data_row.learning_rate
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

"""
    ModelSets

Structure to hold all the model sets in one place.
"""
struct ModelSets
    N::Vector{Symbol}        # locations
    X::Vector{Symbol}        # carriers
    L::Vector{Symbol}        # transmission_lines
    T::Vector{Symbol}        # technologies
    A::Vector{Symbol}        # assets
    A_inv::Vector{Symbol}    # investable_assets
    A_not::Vector{Symbol}    # non_investable_assets
    G::Vector{Symbol}        # generation_assets
    S::Vector{Symbol}        # storage_assets
    S_ST::Vector{Symbol}     # short_term_storage_assets
    S_seas::Vector{Symbol}   # seasonal_storage_assets
    S_seas_in::Vector{Symbol} # seasonal_storage_assets_can_charge
    C::Vector{Symbol}        # conversion_assets
    R::Vector{Int}           # representative periods
    H::Vector{Int}           # timesteps
    D::Vector{Int}           # periods
end
