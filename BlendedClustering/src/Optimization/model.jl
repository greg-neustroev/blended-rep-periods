"""
    create_optimization_model!(connection, model, clustering_result) -> Dict{String,Float64}

Build the investment-and-operations problem in `model` from the DuckDB views
reachable through `connection` and the representative-period weights in
`clustering_result`. Adds all variables, the objective, and every constraint
family in place, and returns a dictionary of per-block formulation timings (in
seconds), including the `"duckdb_queries"` total spent querying DuckDB.
"""
function create_optimization_model!(connection, model, clustering_result)
    # By default JuMP builds a String name for every variable and constraint
    # (e.g. "power_out[asset,1,4567]"); on the full-horizon models that is
    # millions of string interpolations and allocations, and they dominate the
    # build. We never read the MOI names (variables are accessed through the
    # `model[:symbol]` registry, not by name), so disable them — per the JuMP
    # performance tips this roughly halves build time and allocations.
    set_string_names_on_creation(model, false)

    # Per-block formulation timings (seconds). `duckdb_queries` accumulates the
    # time spent executing and *materializing* every DuckDB query, so it can be
    # compared against the JuMP-build time captured by the per-block timers.
    timings = Dict{String,Float64}()
    duckdb_time = Ref(0.0)
    function run_query(sql)
        local result
        duckdb_time[] += @elapsed begin
            result = DBInterface.execute(connection, sql)
            columns(result)  # force full materialization now, not lazily later
        end
        return result
    end

    @timed_step timings "setup" "Creating index sets" begin
        # Create scalars
        operations_weight = get_scalar(connection, "operations_weight")
        timestep_duration = get_scalar(connection, "timestep_duration")

        # Create indexing sets, named with the usual single-letter set notation
        # (L lines, A assets, G generators, S storage, C conversion, ...).
        L = get_index_set(connection, "transmission_lines")
        A = get_index_set(connection, "assets")
        A_inv = get_index_set(connection, "investable_assets")
        G = get_index_set(connection, "generation_assets")
        S = get_index_set(connection, "storage_assets")
        S_ST = get_index_set(connection, "short_term_storage_assets")
        S_seas = get_index_set(connection, "seasonal_storage_assets")
        S_seas_in = get_index_set(connection, "seasonal_storage_assets_can_charge")
        C = get_index_set(connection, "conversion_assets")
        # R, H and D index the representative periods, intra-period timesteps and
        # base periods — contiguous 1:n integer sets. Keep them as ranges so the
        # DenseAxisArray dimensions they index use O(1) arithmetic lookups rather
        # than Dict hashing (these dimensions dominate the indexing work).
        H = as_range(get_index_set(connection, "timesteps"))
        D = as_range(get_index_set(connection, "periods"))
        R = Base.OneTo(size(clustering_result.weight_matrix, 2))

        # Compute the representative period weights in the operations costs
        rp_weight = sum(clustering_result.weight_matrix, dims=1)
        rp_weight .*= operations_weight

        # The inter-period storage chain reconstructs each base period's seasonal
        # storage *increment* from the representatives — the "prolongation" role of the
        # weights. When a separate signed chain matrix W^ch was fit (chain split), it is
        # used for those increment dynamics; the objective/aggregation (rp_weight above),
        # the inter-period ramping, and the *absolute* end-of-horizon reconstruction (the
        # affine S^0 anchor in the storage_initial block) all keep the operational W^op
        # (`weight_matrix`). With no split the chain reuses `weight_matrix`, exactly the
        # historical single-matrix model. Both matrices share `rp_matrix`/`R`, so they
        # are indexed identically.
        chain_weight_matrix = clustering_result.chain_weight_matrix ≡ nothing ?
                              clustering_result.weight_matrix :
                              clustering_result.chain_weight_matrix
    end

    @timed_step timings "variables" "Creating variables" begin
        @variable(model, invested_units[A_inv] ≥ 0)
        @variable(model, power_out[A, R, H] ≥ 0)
        @variable(model, power_in[S_ST∪S_seas_in∪C, R, H] ≥ 0)
        @variable(model, state_of_charge_intra_0[S, R] ≥ 0)
        @variable(model, state_of_charge_intra[S, R, H] ≥ 0)
        @variable(model, state_of_charge_inter_0[S_seas] ≥ 0)
        @variable(model, state_of_charge_inter[S_seas, D] ≥ 0)
        @variable(model, spillage[S_seas, R, H] ≥ 0)
        @variable(model, borrow[S_seas, R, H] ≥ 0)
        # Elastic slacks on the inter-period state-of-charge band (see the
        # `soc_cap_inter` block): `soc_band_over` absorbs over-shoot above the upper
        # level, `soc_band_under` absorbs under-shoot below the lower level. They make
        # the seasonal band a penalised soft constraint rather than a hard one, so the
        # reduced model is feasible-by-construction for any weight class (a conical
        # prolongation that over-scales the band can no longer make it infeasible).
        @variable(model, soc_band_over[S_seas, D] ≥ 0)
        @variable(model, soc_band_under[S_seas, D] ≥ 0)
        @variable(model, flow[L, R, H])
    end

    @timed_step timings "objective" "Creating objective" begin
        # First build an expression for the investment cost;
        # start with a zero cost, query costs from `investment_cost_objective_view`
        # and add them if there are any
        @expression(model, cost_of_investment, AffExpr(0.0))
        investment_cost_data = run_query("SELECT * FROM investment_cost_objective_view")
        if !isempty(investment_cost_data)
            cost_of_investment += sum([
                row.cost * row.unit_capacity * invested_units[row.id]
                for row in rows(investment_cost_data)
            ])
        end

        # Operations, spillage and borrow costs are accumulated in a single pass
        # over each (small) cost table, fanning out across (r, h) with
        # `add_to_expression!` so we never materialise an R×H array of terms nor
        # re-iterate the DuckDB result once per timestep. Each cost component is a
        # separately named expression so the objective can be decomposed after the
        # solve via `value(model[:cost_of_operations])` etc. — kept distinct (rather
        # than lumped into one `cost_of_operations`) so capex/opex/spillage/ENS can
        # be reported without re-solving.
        @expression(model, cost_of_operations, AffExpr(0.0))
        operations_cost_data = run_query("SELECT * FROM operations_cost_objective_view")
        for row in rows(operations_cost_data)
            c = row.variable_cost
            for r in R
                w = rp_weight[r] * c
                for h in H
                    add_to_expression!(cost_of_operations, w, power_out[row.id, r, h])
                end
            end
        end

        @expression(model, cost_of_spillage, AffExpr(0.0))
        spillage_cost_data = run_query("SELECT * FROM spillage_cost_objective_view")
        for row in rows(spillage_cost_data)
            c = row.spillage_cost
            for r in R
                w = rp_weight[r] * c
                for h in H
                    add_to_expression!(cost_of_spillage, w, spillage[row.id, r, h])
                end
            end
        end

        @expression(model, cost_of_borrow, AffExpr(0.0))
        borrow_cost_data = run_query("SELECT * FROM borrow_cost_objective_view")
        for row in rows(borrow_cost_data)
            c = row.borrow_cost
            for r in R
                w = rp_weight[r] * c
                for h in H
                    add_to_expression!(cost_of_borrow, w, borrow[row.id, r, h])
                end
            end
        end

        # Penalty on the elastic inter-period state-of-charge band. The band slack is
        # indexed per base period d (not per representative period), so its per-unit
        # price is `operations_weight * borrow_cost` — exactly the per-base-period rate
        # at which `borrow` is charged (rp_weight already folds operations_weight and
        # the column sum of W into the intra-period borrow). At this rate the penalty
        # is an *exact* penalty: violating the band to free up a unit of stored energy
        # never beats the genuine alternatives (whose marginal value is bounded by the
        # borrow backstop at VOLL), so feasible models leave the slacks at zero and are
        # unchanged, while an otherwise-infeasible band reports finite regret instead.
        @expression(model, cost_of_soc_band, AffExpr(0.0))
        soc_band_cost_data = run_query("SELECT * FROM borrow_cost_objective_view")
        for row in rows(soc_band_cost_data)
            w = operations_weight * row.borrow_cost
            for d in D
                add_to_expression!(cost_of_soc_band, w, soc_band_over[row.id, d])
                add_to_expression!(cost_of_soc_band, w, soc_band_under[row.id, d])
            end
        end

        # Finally, formulate the objective function as the sum of the costs
        @objective(model, Min, cost_of_investment + cost_of_operations + cost_of_spillage + cost_of_borrow + cost_of_soc_band)
    end

    @info "Creating constraints"

    @timed_step timings "balance" "- Adding balance constraints" begin
        # Build the total power in/out and flow expressions as *sparse* maps from
        # (location, carrier, rep_period, timestep) to an AffExpr. Only the
        # combinations that actually carry generation/flow get an entry; absent
        # keys are treated as zero by `add_term!` below. This avoids allocating a
        # dense |N|×|X|×|R|×|H| array of (mostly zero) expressions.
        # Keyed by (location, carrier, rep_period, timestep). The key type is left
        # abstract on purpose: `location`/`carrier` come straight from the input
        # data and are dataset-dependent (e.g. integer bus ids in `sienna` vs
        # string country codes in `tyndp`), so no single concrete tuple type fits.
        total_power_out = Dict{Tuple,AffExpr}()
        total_power_in = Dict{Tuple,AffExpr}()
        total_flow_in = Dict{Tuple,AffExpr}()
        total_flow_out = Dict{Tuple,AffExpr}()

        power_out_data = run_query("SELECT * FROM power_out_expression_view")
        for row in rows(power_out_data), r in R, h in H
            e = get!(() -> AffExpr(0.0), total_power_out, (row.location, row.carrier_out, r, h))
            add_to_expression!(e, power_out[row.id, r, h])
        end

        power_in_data = run_query("SELECT * FROM power_in_expression_view")
        for row in rows(power_in_data), r in R, h in H
            e = get!(() -> AffExpr(0.0), total_power_in, (row.location, row.carrier_in, r, h))
            add_to_expression!(e, power_in[row.id, r, h])
        end

        transmission_line_data = run_query("SELECT * FROM transmission_line_expression_view")
        for row in rows(transmission_line_data), r in R, h in H
            e_out = get!(() -> AffExpr(0.0), total_flow_out, (row.from, row.carrier, r, h))
            add_to_expression!(e_out, flow[row.id, r, h])
            e_in = get!(() -> AffExpr(0.0), total_flow_in, (row.to, row.carrier, r, h))
            add_to_expression!(e_in, flow[row.id, r, h])
        end

        # Now for each location-carrier combination, make the balance constraint.
        # Retain the constraint references keyed by (location, carrier, rep_period,
        # timestep) in `model.ext` so the nodal/marginal prices (their duals) can be
        # exported after the solve without re-deriving which constraint priced which
        # node — these are anonymous constraints, so this is the only handle on them.
        balance_constraints = Dict{Tuple,Any}()
        balance_data = run_query("SELECT * FROM balance_constraint_view")
        for row in rows(balance_data)
            key = (row.location, row.carrier, row.rep_period, row.timestep)
            lhs = AffExpr(0.0)
            add_term!(lhs, total_power_out, key, 1.0)
            add_term!(lhs, total_power_in, key, -1.0)
            add_term!(lhs, total_flow_out, key, -1.0)
            add_term!(lhs, total_flow_in, key, 1.0)
            balance_constraints[key] = @constraint(model, lhs == row.demand_profile * row.peak_demand)
        end
        model.ext[:balance_constraints] = balance_constraints
    end

    @timed_step timings "storage_short_term" "- Adding intra-period short-term storage constraints" begin
        intra_period_short_term_storage_data = run_query(
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
    end

    @timed_step timings "storage_seasonal" "- Adding intra-period seasonal storage constraints" begin
        intra_period_seasonal_storage_can_charge_data = run_query(
            "SELECT * FROM intra_period_seasonal_storage_can_charge_constraint_view"
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

        intra_period_seasonal_storage_cannot_charge_data = run_query(
            "SELECT * FROM intra_period_seasonal_storage_cannot_charge_constraint_view"
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
    end

    @timed_step timings "storage_inter" "- Adding inter-period storage constraints" begin
        # There is no additional data to query here, so we just iterate
        # through the seasonal storage assets and base periods, with the special
        # treatment of the first period
        @constraint(model, [s in S_seas],
            state_of_charge_inter[s, 1] - state_of_charge_inter_0[s]
            ==
            sum(
                chain_weight_matrix[1, r]
                *
                (state_of_charge_intra[s, r, H[end]] - state_of_charge_intra_0[s, r])
                for r in R
            )
        )
        @constraint(model, [s in S_seas, d in D[2:end]],
            state_of_charge_inter[s, d] - state_of_charge_inter[s, d-1]
            ==
            sum(
                chain_weight_matrix[d, r]
                *
                (state_of_charge_intra[s, r, H[end]] - state_of_charge_intra_0[s, r])
                for r in R
            )
        )
    end

    @timed_step timings "storage_initial" "- Adding initial state of charge constraints" begin
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

        initial_storage_data = run_query("SELECT * FROM initial_storage_constraint_view")
        for row in rows(initial_storage_data)
            @constraint(model, state_of_charge_inter_0[row.id] == row.initial_storage_level)
            # This reconstructs an *absolute* end-of-horizon level that must hit the
            # nonzero storage tether S^0 (via the cyclic constraint above), so it needs
            # partition-of-unity weights and stays on the operational `weight_matrix`
            # (W^op). A signed chain matrix W^ch is fit for *increments* (its column
            # sums vanish), so it cannot reconstruct an absolute level inside the
            # bounded intra SoC — using it here is infeasible whenever S^0 ≠ 0. The
            # chain split therefore routes W^ch onto the increment dynamics above and
            # leaves this affine anchor with W^op.
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
    end

    @timed_step timings "conversion" "- Adding conversion constraints" begin
        conversion_asset_data = run_query("SELECT * FROM conversion_constraint_view")
        for row in rows(conversion_asset_data)
            @constraint(model,
                power_in[row.id, row.rep_period, row.timestep] * row.efficiency_in
                ==
                power_out[row.id, row.rep_period, row.timestep] / row.efficiency_out
            )
        end
    end

    @timed_step timings "capacity" "- Adding capacity constraints" begin
        @expression(model, accumulated_capacity[A], AffExpr(0.0))
        capacity_data = run_query("SELECT * FROM capacity_constraint_view")
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
    end

    @timed_step timings "ramping_intra" "- Adding intra-period ramping constraints" begin
        intra_period_ramping_data = run_query("SELECT * FROM intra_period_ramping_constraint_view")
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
    end

    @timed_step timings "ramping_inter" "- Adding inter-period ramping constraints" begin
        inter_period_ramping_data = run_query("SELECT * FROM inter_period_ramping_constraint_view")
        @expression(model, power_out_inter_start[g in G, d in D[2:end]], sum(
            clustering_result.weight_matrix[d, r] * power_out[g, r, 1]
            for r in R
        ))
        @expression(model, power_out_inter_end[g in G, d in D[1:(end-1)]], sum(
            clustering_result.weight_matrix[d, r] * power_out[g, r, H[end]]
            for r in R
        ))
        for row in rows(inter_period_ramping_data)
            @constraint(model,
                power_out_inter_start[row.id, row.period]
                -
                power_out_inter_end[row.id, row.period-1]
                <=
                row.ramping_rate * timestep_duration * accumulated_capacity[row.id]
            )
            @constraint(model,
                power_out_inter_end[row.id, row.period-1]
                -
                power_out_inter_start[row.id, row.period]
                <=
                row.ramping_rate * timestep_duration * accumulated_capacity[row.id]
            )
        end
    end

    @timed_step timings "soc_cap_intra" "- Adding intra-period maximum state of charge constraints" begin
        intraperiod_storage_capacity_data = run_query(
            "SELECT * FROM intra_period_storage_capacity_constraint_view"
        )
        for row in rows(intraperiod_storage_capacity_data)
            if row.investable && row.initial_units > 0
                # Investing in storage scales the energy reservoir together with the
                # power rating (a fixed-duration battery): the SoC cap grows in
                # proportion to the total units, so storage expansion adds MWh as well
                # as MW. Non-investable storage keeps a constant bound (the historical
                # behaviour, since `invested_units` is then absent).
                @constraint(model,
                    state_of_charge_intra[row.id, row.rep_period, row.timestep] <=
                    row.capacity_storage_energy *
                    (row.initial_units + invested_units[row.id]) / row.initial_units
                )
            else
                JuMP.set_upper_bound(
                    state_of_charge_intra[row.id, row.rep_period, row.timestep],
                    row.capacity_storage_energy
                )
            end
        end
    end

    @timed_step timings "soc_cap_inter" "- Adding inter-period maximum state of charge constraints" begin
        # The band on the inter-period state of charge is enforced *softly*: rather
        # than hard variable bounds, the level may breach the band by paying the
        # `cost_of_soc_band` penalty through the elastic slacks `soc_band_under`
        # (below the lower level) and `soc_band_over` (above the upper level). The
        # variable keeps its physical `≥ 0` floor (a deficit there is met by the
        # intra-period `borrow` slack, never infeasible); only the band's min/max
        # *levels* are relaxed. That is the single thing a conical prolongation can
        # over-scale into infeasibility, so softening it makes the model feasible for
        # any weight class while the exact penalty leaves feasible cells untouched.
        interperiod_storage_capacity_data = run_query(
            "SELECT * FROM inter_period_storage_capacity_constraint_view"
        )
        for row in rows(interperiod_storage_capacity_data)
            @constraint(model,
                state_of_charge_inter[row.id, row.period]
                ≥
                row.min_storage_level * row.capacity_storage_energy
                -
                soc_band_under[row.id, row.period]
            )
            @constraint(model,
                state_of_charge_inter[row.id, row.period]
                ≤
                row.max_storage_level * row.capacity_storage_energy
                +
                soc_band_over[row.id, row.period]
            )
        end
    end

    @timed_step timings "transmission" "- Adding transmission capacity constraints" begin
        line_capacity_data = run_query("SELECT * FROM transmission_line_capacity_constraint_view")
        for row in rows(line_capacity_data)
            JuMP.set_lower_bound(flow[row.id, row.rep_period, row.timestep], -row.import_capacity)
            JuMP.set_upper_bound(flow[row.id, row.rep_period, row.timestep], row.export_capacity)
        end
    end

    timings["duckdb_queries"] = duckdb_time[]
    @info "  (DuckDB query + materialize total: $(round(duckdb_time[]; digits=3))s)"
    return timings
end
