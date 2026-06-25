"""
    evaluate_solution!(connection, model, eval_model, period_length, evaluation_type; fix_every=1)

Measure the true cost of the clustered solution by solving the full-horizon model
`eval_model` with the relevant decisions fixed to the values from the clustered
`model`. `evaluation_type` selects which decisions are fixed:

  - `:investment_regret` fixes the investment decisions (`invested_units`);
  - `:storage_regret` fixes the inter-period storage levels (`state_of_charge_inter`).

`fix_every` controls the cadence of the storage-regret endpoint fixing: with the
default `1` every base-period boundary state of charge is pinned (the original,
strictest grading). With `fix_every = k > 1` only every `k`-th period boundary is
pinned, so the full model re-optimizes the storage trajectory between checkpoints.
Day-exact pinning forces the full horizon to honour each period's reconstructed end
state even when the real period's inflow/demand cannot supply it, closing the gap
through penalised borrow — an effect that grows with the number of representative
periods. A coarser cadence grades the *seasonal* trajectory the RP method is meant
to deliver instead. `fix_every` is ignored for `:investment_regret`.
"""
function evaluate_solution!(connection, model, eval_model, period_length, evaluation_type; fix_every::Int=1)
    fix_every ≥ 1 || error("fix_every must be ≥ 1, got $fix_every")
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
        # Sub-sample which period boundaries are pinned. `fix_every = 1` keeps every
        # boundary (the original behaviour); a coarser cadence pins only every k-th.
        if fix_every > 1
            filter!(:period => p -> p % fix_every == 0, fixed_solution_df)
        end
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
end
