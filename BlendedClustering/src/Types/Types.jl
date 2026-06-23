"""
Core data types shared across the package: experiment configuration, clustering
results, and the per-experiment result record.
"""
module Types

using CSV
using DataFrames
using JuMP
using SparseArrays
import Tables

export ExperimentData, ExperimentResult, read_run_data
export ClusteringResult, AuxiliaryClusteringData

# Default weight resolution tolerance used when a configuration CSV omits the
# optional `tol` column: weights are fit and reported only down to this precision.
const DEFAULT_PGD_TOL = 1e-2

# Default clustering-feature normalization used when a configuration CSV omits the
# optional `normalization` column. `:unscaled` reproduces the historical behavior:
# cluster directly on the dimensionless profiles as stored, where each value is a
# fraction of the asset's peak (0 means physically zero, 1 means the peak), with no
# scaling of its own. `:economic` instead rescales the features into common
# physical/economic units (peak × profile, per-block economic weights) before
# clustering. The features are never centered, so zero keeps meaning zero.
const DEFAULT_NORMALIZATION = :unscaled


"""
    read_run_data(path) -> DataFrame

Read the experiment configuration at `path` (a CSV with one experiment per row)
into a `DataFrame`, checking that the required columns are present and converting
the categorical columns (`clustering_type`, `weight_type`, `evaluation_type`, and
the optional `normalization`) to `Symbol`s. The projected gradient descent
tolerance column `tol` is optional and defaults to `DEFAULT_PGD_TOL` when absent;
the `normalization` column is optional and defaults to `DEFAULT_NORMALIZATION`.
"""
function read_run_data(path)
    df = CSV.read(path, DataFrame)
    required_columns = Set([
        "n_rep_periods",
        "period_length",
        "clustering_type",
        "weight_type",
        "evaluation_type",
    ])
    if !issubset(required_columns, Set(names(df)))
        error("Input CSV file is missing required columns: $required_columns")
    end
    symbol_columns = [:clustering_type, :weight_type, :evaluation_type]
    # `normalization` is optional; convert it to a Symbol only when present.
    if "normalization" in names(df)
        push!(symbol_columns, :normalization)
    end
    transform!(df, symbol_columns .=> ByRow(Symbol) .=> symbol_columns)
    return df
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
    weight_type::Symbol
    tol::Float64
    normalization::Symbol
    evaluation_type::Symbol

    function ExperimentData(run_data_row::DataFrameRow{DataFrame,DataFrames.Index}, base_name::String)
        # `tol` (the PGD tolerance ε) is optional; fall back to the default when
        # the configuration CSV does not provide the column.
        tol = hasproperty(run_data_row, :tol) ? Float64(run_data_row.tol) : DEFAULT_PGD_TOL
        # `normalization` is optional and defaults to the historical `:minmax`.
        normalization = hasproperty(run_data_row, :normalization) ?
                        Symbol(run_data_row.normalization) : DEFAULT_NORMALIZATION
        # Keep the experiment name (and thus output paths) unchanged for the
        # default normalization; only the non-default arm gets a suffix, so the
        # same RP grid can be run both ways without colliding.
        name_parts = [
            base_name,
            run_data_row.n_rep_periods,
            run_data_row.period_length,
            string(run_data_row.clustering_type),
            string(run_data_row.weight_type),
            tol,
        ]
        if normalization ≠ DEFAULT_NORMALIZATION
            push!(name_parts, string(normalization))
        end
        name = join(name_parts, "_")
        return new(
            name,
            run_data_row.n_rep_periods,
            run_data_row.period_length,
            run_data_row.clustering_type,
            run_data_row.weight_type,
            tol,
            normalization,
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
end

"""
Structure to hold the clustering result.
"""
mutable struct ClusteringResult
    profiles::AbstractDataFrame
    weight_matrix::Union{SparseMatrixCSC{Float64,Int64},Matrix{Float64}}
    clustering_matrix::Union{Matrix{Float64},Nothing}
    rp_matrix::Union{Matrix{Float64},Nothing}
end

# Convenience constructor: the clustering/representative-period matrices are not
# always available (e.g. the single-period fast path), so default them to nothing.
ClusteringResult(profiles, weight_matrix) =
    ClusteringResult(profiles, weight_matrix, nothing, nothing)

struct ExperimentResult
    name::String
    seed::Int
    n_rep_periods::Int
    period_length::Int
    clustering_type::Symbol
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

end # module Types
