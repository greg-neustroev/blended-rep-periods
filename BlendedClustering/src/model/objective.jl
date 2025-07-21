"""
    create_objective!(model, con, sets, rp_weight)

Create the objective function for the energy system model.
"""
function create_objective!(model, con, sets, rp_weight)
    @info "Creating objective"

    # Investment costs
    cost_of_investment = create_investment_cost_expression(model, con, sets)

    # Operations costs  
    cost_of_operations = create_operations_cost_expression(model, con, sets, rp_weight)

    # Set objective
    @objective(model, Min, cost_of_investment + cost_of_operations)

    return nothing
end

function create_investment_cost_expression(model, con, sets)
    investment_data = DBInterface.execute(
        con,
        "SELECT id, unit_capacity, cost FROM investable_assets"
    )

    return @expression(
        model,
        isempty(investment_data) ?
        0.0 :
        sum([
            row.cost * row.unit_capacity * model[:invested_units][row.id]
            for row in rows(investment_data)
        ])
    )
end

function create_operations_cost_expression(model, con, sets, rp_weight)
    operations_data = DBInterface.execute(
        con,
        "SELECT id, variable_cost FROM generation_assets"
    )
    spillage_data = DBInterface.execute(
        con,
        "SELECT id, spillage_cost FROM seasonal_storage_assets"
    )

    return @expression(
        model,
        isempty(operations_data) ?
        0.0 :
        sum([
            rp_weight[r] * (
                sum([
                    row.variable_cost * model[:power_out][row.id, r, h]
                    for row in rows(operations_data)
                ])
                +
                sum([
                    row.spillage_cost * model[:spillage][row.id, r, h]
                    for row in rows(spillage_data)
                ])
            )
            for r in sets.R, h in sets.H
        ])
    )
end
