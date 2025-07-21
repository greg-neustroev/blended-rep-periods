"""
    create_variables!(model, sets)

Create all optimization variables for the energy system model.
"""
function create_variables!(model, sets)
    @info "Creating variables"

    # Extract sets for readability
    A_inv, A, S_ST, S_seas_in, C, S, S_seas, L, R, H, D =
        sets.A_inv, sets.A, sets.S_ST, sets.S_seas_in, sets.C,
        sets.S, sets.S_seas, sets.L, sets.R, sets.H, sets.D

    # Investment variables
    @variable(model, invested_units[A_inv] ≥ 0)

    # Power flow variables
    @variable(model, power_out[A, R, H] ≥ 0)
    @variable(model, power_in[S_ST∪S_seas_in∪C, R, H] ≥ 0)

    # Storage state variables
    @variable(model, state_of_charge_intra_0[S, R] ≥ 0)
    @variable(model, state_of_charge_intra[S, R, H] ≥ 0)
    @variable(model, state_of_charge_inter_0[S_seas] ≥ 0)
    @variable(model, state_of_charge_inter[S_seas, D] ≥ 0)

    # Spillage and transmission variables
    @variable(model, spillage[S_seas, R, H] ≥ 0)
    @variable(model, flow[L, R, H])

    return nothing
end
