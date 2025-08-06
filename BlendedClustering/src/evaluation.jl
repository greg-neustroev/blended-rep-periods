function evaluate_solution!(connection, model, eval_model, period_length, evaluation_type)
    n_rep_periods = 1
    raw_period_length = DBInterface.execute(
                            connection,
                            "SELECT MAX(timestep) FROM profiles_raw"
                        ) |> first |> first
    create_views(connection, raw_period_length)
    profiles = DBInterface.execute(
        connection,
        "SELECT period AS rep_period, * EXCLUDE period FROM profiles"
    ) |> DataFrame
    DuckDB.register_data_frame(connection, profiles, "rp_profiles")
    dummy_clustering_result = ClusteringResult(
        profiles,
        Matrix{Float64}(I, n_rep_periods, n_rep_periods),
    )
    create_optimization_model!(connection, eval_model, dummy_clustering_result)
    # Fix the variables
    @info "Fixing the variables in the evaluation model"
    if evaluation_type == :storage_regret
        fixed_eval_df = Containers.rowtable(
            eval_model[:state_of_charge_intra]; header=[:id, :rep_period, :timestep, :variable]
        ) |> DataFrame
        fixed_solution_df = Containers.rowtable(
            model[:state_of_charge_inter]; header=[:id, :period, :variable]
        ) |> DataFrame
        fixed_solution_df.value = value.(fixed_solution_df.variable)
        fixed_solution_df.timestep = fixed_solution_df.period .* period_length
        select!(fixed_solution_df, Not(:variable, :period))
        for row in innerjoin(fixed_eval_df, fixed_solution_df; on=[:id, :timestep]) |> eachrow
            JuMP.fix(row.variable, row.value; force=true)
        end
    elseif evaluation_type == :investment_regret
        fixed_eval_df = Containers.rowtable(
            eval_model[:invested_units]; header=[:id, :variable]
        ) |> DataFrame
        fixed_solution_df = Containers.rowtable(
            model[:invested_units]; header=[:id, :variable]
        ) |> DataFrame
        fixed_solution_df.value = value.(fixed_solution_df.variable)
        select!(fixed_solution_df, Not(:variable))
        for row in innerjoin(fixed_eval_df, fixed_solution_df; on=:id) |> eachrow
            JuMP.fix(row.variable, row.value; force=true)
        end
    end
    @info "Solving the evaluation model"
    optimize!(eval_model)
    # if termination_status(eval_model) == JuMP.INFEASIBLE
    #     compute_conflict!(eval_model)
    #     iis_model, reference_map = JuMP.copy_conflict(eval_model)
    #     print(iis_model)
    # end
end
