export ExperimentData, RunData, ExperimentResult


"""
Data needed to run all of the experiments
"""
struct RunData <: AbstractDataFrame
    function RunData(path)
        df = CSV.read(path, DataFrame)
        # Ensure the DataFrame has the correct column names
        required_columns = [
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
    if s == "euclidean"
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

function semimetric_to_string(distance::SemiMetric)
    if distance isa Euclidean
        return "euclidean"
    elseif distance isa Cityblock
        return "cityblock"
    elseif distance isa CosineDist
        return "cosine"
    elseif distance isa Chebyshev
        return "chebyshev"
    else
        return string(distance)
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

    function ExperimentData(run_data_row::DataFrameRow{DataFrame,DataFrames.Index}, base_name::String)
        name = join([
                base_name,
                run_data_row.n_rep_periods,
                run_data_row.period_length,
                string(run_data_row.clustering_type),
                semimetric_to_string(run_data_row.distance),
                string(run_data_row.weight_type),
                run_data_row.niters,
                run_data_row.learning_rate,
            ],
            "_"
        )
        return new(
            name,
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
    n_rep_periods::Int
    period_length::Int
    clustering_type::Symbol
    distance::SemiMetric
    weight_type::Symbol
    projection_error::Float64
    termination_status::String
    objective_value::Union{Float64,Missing}
    evaluation_termination_status::String
    evaluated_objective_value::Union{Float64,Missing}
    total_spillage::Float64
    total_borrow::Float64
    time_to_preprocess::Float64
    time_to_cluster::Float64
    time_to_fit_weights::Float64
    time_to_formulate_model::Float64
    time_to_solve::Float64

    function ExperimentResult(
        data::ExperimentData,
        seed::Int,
        solved_model::JuMP.AbstractModel,
        eval_model::Union{JuMP.AbstractModel,Nothing},
        projection_error::Float64,
        time_to_preprocess::Float64,
        time_to_cluster::Float64,
        time_to_fit_weights::Float64,
        time_to_formulate_model::Float64,
        time_to_solve::Float64,
    )
        name = data.name
        n_rep_periods = data.n_rep_periods
        period_length = data.period_length
        clustering_type = data.clustering_type
        distance = data.distance
        weight_type = data.weight_type
        termination_status = solved_model |> JuMP.termination_status |> string
        objective_value = if JuMP.is_solved_and_feasible(solved_model)
            solved_model |> JuMP.objective_value
        else
            missing
        end
        evaluation_termination_status = if !isnothing(eval_model)
            eval_model |> JuMP.termination_status |> string
        else
            "N/A"
        end
        evaluated_objective_value = if !isnothing(eval_model) && JuMP.is_solved_and_feasible(eval_model)
            eval_model |> JuMP.objective_value
        else
            missing
        end
        total_spillage = if !isnothing(eval_model) && haskey(eval_model, :spillage) && !isempty(eval_model[:spillage])
            value.(eval_model[:spillage]) |> sum
        else
            0.0
        end
        total_borrow = if !isnothing(eval_model) && haskey(eval_model, :borrow) && !isempty(eval_model[:borrow])
            value.(eval_model[:borrow]) |> sum
        else
            0.0
        end
        return new(
            name,
            seed,
            n_rep_periods,
            period_length,
            clustering_type,
            distance,
            weight_type,
            projection_error,
            termination_status,
            objective_value,
            evaluation_termination_status,
            evaluated_objective_value,
            total_spillage,
            total_borrow,
            time_to_preprocess,
            time_to_cluster,
            time_to_fit_weights,
            time_to_formulate_model,
            time_to_solve,
        )
    end
end

Tables.columns(res::ExperimentResult) = (;
    name=[res.name],
    seed=[res.seed],
    n_rep_periods=[res.n_rep_periods],
    period_length=[res.period_length],
    clustering_type=[string(res.clustering_type)],
    distance=[semimetric_to_string(res.distance)],
    weight_type=[string(res.weight_type)],
    projection_error=[res.projection_error],
    termination_status=[res.termination_status],
    objective_value=[res.objective_value],
    evaluation_termination_status=[res.evaluation_termination_status],
    evaluated_objective_value=[res.evaluated_objective_value],
    total_spillage=[res.total_spillage],
    total_borrow=[res.total_borrow],
    time_to_preprocess=[res.time_to_preprocess],
    time_to_cluster=[res.time_to_cluster],
    time_to_fit_weights=[res.time_to_fit_weights],
    time_to_formulate_model=[res.time_to_formulate_model],
    time_to_solve=[res.time_to_solve],
)
