export run_experiment

function run_experiment(
    experiment_data::ExperimentData,
    model::JuMP.GenericModel{Float64},
    eval_model::JuMP.GenericModel{Float64},
    connection,
    seed::Int
)
    @info styled"{info:Running experiment: $(experiment_data.name) with seed $seed}"
    Random.seed!(seed)

    # Extract parameters from the experiment data into local variables
    n_rep_periods = experiment_data.n_rep_periods
    clustering_type = experiment_data.clustering_type
    weight_type = experiment_data.weight_type
    period_length = experiment_data.period_length
    evaluation_type = experiment_data.evaluation_type

    @info "Preprocessing data"
    # Create the database views for quering constraints and objective data
    time_to_preprocess = @elapsed create_views(connection, period_length)

    @info styled"{bold:Finding $n_rep_periods $(string(clustering_type)) representative periods}"

    @info "Performing $(string(clustering_type)) clustering"
    time_to_cluster = @elapsed clustering_result = cluster_using_experiment_data(experiment_data, connection)

    projection_errors = clustering_result.rp_matrix * clustering_result.weight_matrix' - clustering_result.clustering_matrix
    initial_projection_error = norm(projection_errors)
    @info "Initial norm of the projection error: $initial_projection_error"

    @info "Fitting $(string(weight_type)) weights"
    # Fit the weights 
    time_to_fit_weights = @elapsed fit_rep_period_weights!(
        clustering_result;
        weight_type=experiment_data.weight_type,
        learning_rate=experiment_data.learning_rate,
        niters=experiment_data.niters,
    )

    projection_errors = clustering_result.rp_matrix * clustering_result.weight_matrix' - clustering_result.clustering_matrix
    final_projection_error = norm(projection_errors)
    @info "Norm of the projection error after fitting: $final_projection_error"

    # Register the representative period profiles in DuckDB
    # so that the views use the correct data
    DuckDB.register_data_frame(
        connection,
        clustering_result.profiles,
        "rp_profiles"
    )

    # Create model
    @info styled"{bold:Creating the model}"

    time_to_formulate_model = @elapsed create_optimization_model!(connection, model, clustering_result)

    # Solve
    @info styled"{bold:Solving the model}"
    time_to_solve = @elapsed optimize!(model)
    @info "Model termination status: $(JuMP.termination_status(model))"

    if JuMP.termination_status(model) == MOI.OPTIMAL && evaluation_type != :none
        # Evaluate the solution. Resolve the full problem with some of the
        # decision variables fixed to the values from the clustering model
        # to get the true objective value.
        @info styled"{bold:Evaluating the model}"
        evaluate_solution!(connection, model, eval_model, period_length, evaluation_type)
    end

    # Return the results
    @info styled"{bold:Returning the results}"
    result = ExperimentResult(
        experiment_data,
        seed,
        model,
        evaluation_type == :none ? nothing : eval_model,
        final_projection_error,
        time_to_preprocess,
        time_to_cluster,
        time_to_fit_weights,
        time_to_formulate_model,
        time_to_solve
    )
    return result
end
