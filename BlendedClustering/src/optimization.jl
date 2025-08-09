function create_optimization_model!(connection, model, clustering_result)
    # Create scalars
    operations_weight = get_scalar(connection, "operations_weight")
    timestep_duration = get_scalar(connection, "timestep_duration")

    # Create indexing sets; we use the same notation as in the paper here
    N = get_index_set(connection, "locations")
    X = get_index_set(connection, "carriers")
    L = get_index_set(connection, "transmission_lines")
    T = get_index_set(connection, "technologies")
    A = get_index_set(connection, "assets")
    A_inv = get_index_set(connection, "investable_assets")
    A_not = get_index_set(connection, "non_investable_assets")
    G = get_index_set(connection, "generation_assets")
    S = get_index_set(connection, "storage_assets")
    S_ST = get_index_set(connection, "short_term_storage_assets")
    S_seas = get_index_set(connection, "seasonal_storage_assets")
    S_seas_in = get_index_set(connection, "seasonal_storage_assets_can_charge")
    C = get_index_set(connection, "conversion_assets")
    H = get_index_set(connection, "timesteps")
    D = get_index_set(connection, "periods")
    R = collect(1:size(clustering_result.weight_matrix, 2))

    # Compute the representative period weights in the operations costs
    rp_weight = sum(clustering_result.weight_matrix, dims=1)
    rp_weight .*= operations_weight

    # Create variables
    @info "Creating variables"
    @variable(model, invested_units[A_inv] ≥ 0)
    @variable(model, power_out[A, R, H] ≥ 0)
    @variable(model, power_in[S_ST∪S_seas_in∪C, R, H] ≥ 0)
    @variable(model, state_of_charge_intra_0[S, R] ≥ 0)
    @variable(model, state_of_charge_intra[S, R, H] ≥ 0)
    @variable(model, state_of_charge_inter_0[S_seas] ≥ 0)
    @variable(model, state_of_charge_inter[S_seas, D] ≥ 0)
    @variable(model, spillage[S_seas, R, H] ≥ 0)
    @variable(model, borrow[S_seas, R, H] ≥ 0)
    @variable(model, flow[L, R, H])

    @info "Creating objective"
    # First build an expression for the investment cost;
    # start with a zero cost, query costs from `investment_cost_objective_view`
    # and add them if there are any
    @expression(model, cost_of_investment, AffExpr(0.0))
    investment_cost_data = DBInterface.execute(
        connection,
        "SELECT * FROM investment_cost_objective_view"
    )
    if !isempty(investment_cost_data)
        cost_of_investment += sum([
            row.cost * row.unit_capacity * invested_units[row.id]
            for row in rows(investment_cost_data)
        ])
    end

    # Next build an expression for the investment cost in the same way
    @expression(model, cost_of_operations, AffExpr(0.0))
    operations_cost_data = DBInterface.execute(
        connection,
        "SELECT * FROM operations_cost_objective_view"
    )
    if !isempty(operations_cost_data)
        cost_of_operations +=
            sum([
                rp_weight[r] *
                sum([
                    row.variable_cost * power_out[row.id, r, h]
                    for row in rows(operations_cost_data)
                ])
                for r in R, h in H
            ])
    end
    spillage_cost_data = DBInterface.execute(
        connection,
        "SELECT * FROM spillage_cost_objective_view"
    )
    if !isempty(spillage_cost_data)
        cost_of_operations +=
            sum([
                rp_weight[r] *
                sum([
                    row.spillage_cost * spillage[row.id, r, h]
                    for row in rows(spillage_cost_data)
                ])
                for r in R, h in H
            ])
    end
    borrow_cost_data = DBInterface.execute(
        connection,
        "SELECT * FROM borrow_cost_objective_view"
    )
    if !isempty(borrow_cost_data)
        cost_of_operations +=
            sum([
                rp_weight[r] *
                sum([
                    row.borrow_cost * borrow[row.id, r, h]
                    for row in rows(borrow_cost_data)
                ])
                for r in R, h in H
            ])
    end

    # Finally, formulate the objective function as the sum of the costs
    @objective(model, Min, cost_of_investment + cost_of_operations)

    @info "Creating constraints"

    @info "- Adding balance constraints"
    # First build the expressions for total power in/out and flows
    @expression(model, total_power_out[N, X, R, H], AffExpr(0.0))
    power_out_data = DBInterface.execute(
        connection,
        "SELECT * FROM power_out_expression_view"
    )
    for row in rows(power_out_data)
        for r in R, h in H
            total_power_out[row.location, row.carrier_out, r, h] += power_out[row.id, r, h]
        end
    end

    @expression(model, total_power_in[N, X, R, H], AffExpr(0.0))
    power_in_data = DBInterface.execute(
        connection,
        "SELECT * FROM power_in_expression_view"
    )
    for row in rows(power_in_data)
        for r in R, h in H
            total_power_in[row.location, row.carrier_in, r, h] += power_in[row.id, r, h]
        end
    end

    @expression(model, total_flow_in[N, X, R, H], AffExpr(0.0))
    @expression(model, total_flow_out[N, X, R, H], AffExpr(0.0))
    transmission_line_data = DBInterface.execute(
        connection,
        "SELECT * FROM transmission_line_expression_view"
    )
    for row in rows(transmission_line_data)
        for r in R, h in H
            total_flow_out[row.from, row.carrier, r, h] += flow[row.id, r, h]
            total_flow_in[row.to, row.carrier, r, h] += flow[row.id, r, h]
        end
    end

    # Now for each location-carrier combination, make the balance constraint
    balance_data = DBInterface.execute(
        connection,
        "SELECT * FROM balance_constraint_view"
    )
    for row in rows(balance_data)
        @constraint(model,
            total_power_out[row.location, row.carrier, row.rep_period, row.timestep]
            -
            total_power_in[row.location, row.carrier, row.rep_period, row.timestep]
            -
            total_flow_out[row.location, row.carrier, row.rep_period, row.timestep]
            +
            total_flow_in[row.location, row.carrier, row.rep_period, row.timestep]
            ==
            row.demand_profile * row.peak_demand
        )
    end

    @info "- Adding intra-period short-term storage constraints"
    intra_period_short_term_storage_data = DBInterface.execute(
        connection,
        "SELECT * FROM intra_period_short_term_storage_constraint_view"
    )
    for row in rows(intra_period_short_term_storage_data)
        if row.timestep == 1
            @constraint(model,
                state_of_charge_intra[row.id, row.rep_period, 1]
                -
                state_of_charge_intra_0[row.id, row.rep_period]
                ==
                (
                    row.efficiency_in * power_in[row.id, row.rep_period, 1]
                    -
                    power_out[row.id, row.rep_period, 1] / row.efficiency_out
                ) * timestep_duration
            )
        else
            @constraint(model,
                state_of_charge_intra[row.id, row.rep_period, row.timestep]
                -
                state_of_charge_intra[row.id, row.rep_period, row.timestep-1]
                ==
                (
                    row.efficiency_in * power_in[row.id, row.rep_period, row.timestep]
                    -
                    power_out[row.id, row.rep_period, row.timestep] / row.efficiency_out
                ) * timestep_duration
            )
        end
    end

    @info "- Adding intra-period seasonal storage constraints"
    intra_period_seasonal_storage_can_charge_data = DBInterface.execute(
        connection,
        """
        SELECT * FROM intra_period_seasonal_storage_can_charge_constraint_view
        """
    )
    for row in rows(intra_period_seasonal_storage_can_charge_data)
        if row.timestep == 1
            @constraint(model,
                state_of_charge_intra[row.id, row.rep_period, 1]
                -
                state_of_charge_intra_0[row.id, row.rep_period]
                ==
                (
                    row.efficiency_in * power_in[row.id, row.rep_period, 1]
                    -
                    power_out[row.id, row.rep_period, 1] / row.efficiency_out
                ) * timestep_duration
                -
                spillage[row.id, row.rep_period, 1]
                +
                borrow[row.id, row.rep_period, 1]
                +
                row.inflow_profile * row.peak_inflow
            )
        else
            @constraint(model,
                state_of_charge_intra[row.id, row.rep_period, row.timestep]
                -
                state_of_charge_intra[row.id, row.rep_period, row.timestep-1]
                ==
                (
                    row.efficiency_in * power_in[row.id, row.rep_period, row.timestep]
                    -
                    power_out[row.id, row.rep_period, row.timestep] / row.efficiency_out
                ) * timestep_duration
                -
                spillage[row.id, row.rep_period, row.timestep]
                +
                borrow[row.id, row.rep_period, row.timestep]
                +
                row.inflow_profile * row.peak_inflow
            )
        end
    end

    intra_period_seasonal_storage_cannot_charge_data = DBInterface.execute(
        connection,
        """
        SELECT * FROM intra_period_seasonal_storage_cannot_charge_constraint_view
        """
    )
    for row in rows(intra_period_seasonal_storage_cannot_charge_data)
        if row.timestep == 1
            @constraint(model,
                state_of_charge_intra[row.id, row.rep_period, 1]
                -
                state_of_charge_intra_0[row.id, row.rep_period]
                ==
                -power_out[row.id, row.rep_period, 1] / row.efficiency_out * timestep_duration
                -
                spillage[row.id, row.rep_period, 1]
                +
                borrow[row.id, row.rep_period, 1]
                +
                row.inflow_profile * row.peak_inflow
            )
        else
            @constraint(model,
                state_of_charge_intra[row.id, row.rep_period, row.timestep]
                -
                state_of_charge_intra[row.id, row.rep_period, row.timestep-1]
                ==
                -power_out[row.id, row.rep_period, row.timestep] / row.efficiency_out * timestep_duration
                -
                spillage[row.id, row.rep_period, row.timestep]
                +
                borrow[row.id, row.rep_period, row.timestep]
                +
                row.inflow_profile * row.peak_inflow
            )
        end
    end

    @info "- Adding inter-period storage constraints"
    # There is no additional data to query here, so we just iterate
    # through the seasonal storage assets and base periods, with the special
    # treatment of the first period
    @constraint(model, [s in S_seas],
        state_of_charge_inter[s, 1] - state_of_charge_inter_0[s]
        ==
        sum(
            clustering_result.weight_matrix[1, r]
            *
            (state_of_charge_intra[s, r, H[end]] - state_of_charge_intra_0[s, r])
            for r in R
        )
    )
    @constraint(model, [s in S_seas, d in D[2:end]],
        state_of_charge_inter[s, d] - state_of_charge_inter[s, d-1]
        ==
        sum(
            clustering_result.weight_matrix[d, r]
            *
            (state_of_charge_intra[s, r, H[end]] - state_of_charge_intra_0[s, r])
            for r in R
        )
    )

    @info "- Adding initial state of charge constraints"
    # The state of charge at the end of the last period is equal to the state of
    # charge at the beginning of the first period for seasonal storage assets so
    # that the model does not create extra state of charge in the first period
    # that can be used for free by discharging the storage through the periods.
    @constraint(model, [s in S_seas],
        state_of_charge_inter[s, D[end]] == state_of_charge_inter_0[s]
    )
    # Similarly, the state of charge at the end of each representative period
    # is equal to the state of charge at the beginning of that period
    # for short-term storage assets.
    @constraint(model, [s in S_ST, r in R],
        state_of_charge_intra[s, r, H[end]] == state_of_charge_intra_0[s, r]
    )

    initial_storage_data = DBInterface.execute(
        connection,
        "SELECT * FROM initial_storage_constraint_view"
    )
    for row in rows(initial_storage_data)
        @constraint(model, state_of_charge_inter_0[row.id] == row.initial_storage_level)
        @constraint(model,
            state_of_charge_inter[row.id, D[end]]
            ==
            sum(
                clustering_result.weight_matrix[D[end], r]
                *
                state_of_charge_intra[row.id, r, H[end]]
                for r in R
            )
        )
    end

    @info "- Adding conversion constraints"
    conversion_asset_data = DBInterface.execute(
        connection,
        """
        SELECT * FROM conversion_constraint_view
        """
    )
    for row in rows(conversion_asset_data)
        @constraint(model,
            power_in[row.id, row.rep_period, row.timestep] * row.efficiency_in
            ==
            power_out[row.id, row.rep_period, row.timestep] / row.efficiency_out
        )
    end

    @info "- Adding capacity constraints"
    @expression(model, accumulated_capacity[A], AffExpr(0.0))
    capacity_data = DBInterface.execute(
        connection,
        """
        SELECT * FROM capacity_constraint_view
        """
    )
    for row in rows(capacity_data)
        accumulated_capacity[row.id] = row.unit_capacity * (
            if row.investable
                row.initial_units + invested_units[row.id]
            else
                row.initial_units
            end
        )
        @constraint(model, power_out[row.id, row.rep_period, row.timestep] <= row.availability_profile * accumulated_capacity[row.id])
    end
    for s in S_ST, r in R, h in H
        @constraint(model, power_in[s, r, h] <= accumulated_capacity[s])
    end
    for s in S_seas_in, r in R, h in H
        @constraint(model, power_in[s, r, h] <= accumulated_capacity[s])
    end

    @info "- Adding intra-period ramping constraints"
    intra_period_ramping_data = DBInterface.execute(
        connection,
        "SELECT * FROM intra_period_ramping_constraint_view"
    )
    for row in rows(intra_period_ramping_data)
        @constraint(model,
            power_out[row.id, row.rep_period, row.timestep]
            -
            power_out[row.id, row.rep_period, row.timestep-1]
            <=
            row.ramping_rate * timestep_duration * accumulated_capacity[row.id]
        )
        @constraint(model,
            power_out[row.id, row.rep_period, row.timestep-1]
            -
            power_out[row.id, row.rep_period, row.timestep]
            <=
            row.ramping_rate * timestep_duration * accumulated_capacity[row.id]
        )
    end

    # @info "- Adding inter-period ramping constraints"
    # inter_period_ramping_data = DBInterface.execute(
    #     connection,
    #     "SELECT * FROM inter_period_ramping_constraint_view"
    # )
    # @expression(model, power_out_inter_start[g in G, d in D[2:end]], sum(
    #     clustering_result.weight_matrix[d, r] * power_out[g, r, 1]
    #     for r in R
    # ))
    # @expression(model, power_out_inter_end[g in G, d in D[1:(end-1)]], sum(
    #     clustering_result.weight_matrix[d, r] * power_out[g, r, H[end]]
    #     for r in R
    # ))
    # for row in rows(inter_period_ramping_data)
    #     @constraint(model,
    #         power_out_inter_start[row.id, row.period]
    #         -
    #         power_out_inter_end[row.id, row.period-1]
    #         <=
    #         row.ramping_rate * timestep_duration * accumulated_capacity[row.id]
    #     )
    #     @constraint(model,
    #         power_out_inter_end[row.id, row.period-1]
    #         -
    #         power_out_inter_start[row.id, row.period]
    #         <=
    #         row.ramping_rate * timestep_duration * accumulated_capacity[row.id]
    #     )
    # end

    @info "- Adding intra-period maximum state of charge constraints"
    intraperiod_storage_capacity_data = DBInterface.execute(
        connection,
        """
        SELECT * FROM intra_period_storage_capacity_constraint_view
        """
    )
    for row in rows(intraperiod_storage_capacity_data)
        JuMP.set_upper_bound(
            state_of_charge_intra[row.id, row.rep_period, row.timestep],
            row.capacity_storage_energy
        )
    end

    @info "- Adding inter-period maximum state of charge constraints"
    interperiod_storage_capacity_data = DBInterface.execute(
        connection,
        "SELECT * FROM inter_period_storage_capacity_constraint_view"
    )
    for row in rows(interperiod_storage_capacity_data)
        JuMP.set_lower_bound(
            state_of_charge_inter[row.id, row.period],
            row.min_storage_level * row.capacity_storage_energy
        )
        JuMP.set_upper_bound(
            state_of_charge_inter[row.id, row.period],
            row.max_storage_level * row.capacity_storage_energy
        )
    end

    @info "- Adding transmission capacity constraints"
    line_capacity_data = DBInterface.execute(
        connection,
        """
        SELECT * FROM transmission_line_capacity_constraint_view
        """
    )
    for row in rows(line_capacity_data)
        JuMP.set_lower_bound(flow[row.id, row.rep_period, row.timestep], -row.import_capacity)
        JuMP.set_upper_bound(flow[row.id, row.rep_period, row.timestep], row.export_capacity)
    end
end
