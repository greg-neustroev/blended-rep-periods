"""
    create_optimization_model!(connection, model, clustering_result) -> Dict{String,Float64}

Build the investment-and-operations problem in `model` from the DuckDB views
reachable through `connection` and the representative-period weights in
`clustering_result`. Adds all variables, the objective, and every constraint
family in place, and returns a dictionary of per-block formulation timings (in
seconds), including the `"duckdb_queries"` total spent querying DuckDB.

The inter-period storage chain uses exact-inflow injection (for any genuine multi-period
reduction): it reconstructs each base period's seasonal increment from only the *dispatch
response* of the representatives (the net increment minus each RP's exogenous inflow) and
adds the base period's own *real* inflow exactly, closing the per-base-day ledger with
conservation slacks and enforcing the seasonal band as a hard [min,max]·cap constraint.
This removes the integrated inflow-approximation drift and restores the true annual water
balance; the absolute-level reconstruction clause is unnecessary and dropped. A single base
period (the n_rep=1 reference / full-horizon evaluation) has no chain to reconstruct, so
injection is inert there and the model reduces to the plain formulation.
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

        # The inter-period storage chain reconstructs each base period's seasonal storage
        # *increment* from the representatives (the "prolongation" role of the weights),
        # using the operational weight matrix — the same weights that drive the
        # objective/aggregation (`rp_weight`) and the inter-period ramping.
        weight_matrix = clustering_result.weight_matrix
        # Exact-inflow injection is a multi-period *reconstruction* device: it only bites when
        # a chain of base periods is rebuilt from fewer representatives. A single base period
        # (the n_rp=1 reference optimum and the full-horizon evaluation model) has no chain to
        # reconstruct, so injection there would only let its conservation slacks perturb the
        # annual balance. It is therefore active only for a genuine reduction, so those models
        # stay the plain formulation and the regret denominator is unchanged.
        inject = length(D) > 1
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
        @variable(model, flow[L, R, H])
        # Per-base-day conservation slacks on the inter-period storage chain, added only
        # under exact-inflow injection. `soc_chain_spill` removes mass (spill) and
        # `soc_chain_borrow` adds it (energy-not-served), both at base-day resolution — the
        # resolution the injected forcing E_d already lives at. They are the recourse that
        # lets the seasonal band be *hard* (see `soc_cap_inter`) while staying feasible for
        # any dispatch: a run-of-river reservoir whose real daily inflow overtops the buffer
        # spills the excess here (correcting every downstream day at once) instead of letting
        # the reconstructed trajectory ratchet out of the reservoir. Priced at the asset's own
        # spillage/borrow cost (the full-model prices), so the regime is arbitrated by price,
        # not a classifier: a seasonal reservoir with a large buffer leaves them at zero.
        if inject
            @variable(model, soc_chain_spill[S_seas, D] ≥ 0)
            @variable(model, soc_chain_borrow[S_seas, D] ≥ 0)
        end
    end

    @timed_step timings "objective" "Creating objective" begin
        # First build an expression for the investment cost;
        # start with a zero cost, query costs from `investment_cost_objective_view`
        # and add them if there are any
        @expression(model, cost_of_investment, AffExpr(0.0))
        investment_cost_data = run_query("SELECT * FROM investment_cost_objective_view")
        # Mutate the *registered* expression in place (like every other cost term below).
        # A rebinding `cost_of_investment += …` would leave `model[:cost_of_investment]` at
        # zero, so the objective would be correct but the reported capex component (read via
        # `value(model[:cost_of_investment])`) would spuriously record 0.
        for row in rows(investment_cost_data)
            add_to_expression!(cost_of_investment, row.cost * row.unit_capacity, invested_units[row.id])
        end

        # Operations, spillage and borrow costs are accumulated in a single pass
        # over each (small) cost table, fanning out across (r, h) with
        # `add_to_expression!` so we never materialise an R×H array of terms nor
        # re-iterate the DuckDB result once per timestep. Each cost component is a
        # separately named expression so the objective can be decomposed after the
        # solve via `value(model[:cost_of_operations])` etc. — kept distinct (rather
        # than lumped into one `cost_of_operations`) so capex/opex/spillage/ENS can
        # be reported without re-solving.
        # Unit convention in the objective: all operational cost terms are put on a common
        # per-power ($/h) basis so MW and MWh variables agree. `power_out` is already a power
        # (MW), so its cost uses `variable_cost` directly; the energy-based slacks (spillage /
        # borrow below, which are per-step energy MWh — see the seasonal storage balance) are
        # divided by the timestep duration τ to convert MWh→MW. τ is applied to each unit
        # *price* (never to `rp_weight`, which is shared across terms and carries the period
        # count only). For hourly data (τ=1) this is a no-op; it matters only for sub-hourly
        # datasets such as RTS (τ = 5/60 h), where the objective becomes (true cost ÷ τ) — a
        # global scalar that leaves the optimum and the regret ratio unchanged.
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
            c = row.spillage_cost / timestep_duration   # MWh→MW: match the power-basis operations term
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
            c = row.borrow_cost / timestep_duration   # MWh→MW: match the power-basis operations term
            for r in R
                w = rp_weight[r] * c
                for h in H
                    add_to_expression!(cost_of_borrow, w, borrow[row.id, r, h])
                end
            end
        end

        # Cost of the per-base-day chain conservation slacks (injection only; zero otherwise).
        # Same per-base-day pricing convention as the band above (operations_weight × unit
        # cost): spill at the asset's spillage_cost, borrow at its borrow_cost — the prices the
        # full model pays. Kept as their own expressions (not folded into cost_of_spillage /
        # cost_of_borrow, which accumulate the intra-period slacks over (r, h)) so the chain
        # ledger can be reported separately from the intra-period physics.
        @expression(model, cost_of_chain_spill, AffExpr(0.0))
        @expression(model, cost_of_chain_borrow, AffExpr(0.0))
        if inject
            spill_price = Dict(string(r.id) => r.spillage_cost
                               for r in rows(run_query("SELECT id, spillage_cost FROM spillage_cost_objective_view")))
            borrow_price = Dict(string(r.id) => r.borrow_cost
                                for r in rows(run_query("SELECT id, borrow_cost FROM borrow_cost_objective_view")))
            for s in S_seas, d in D
                # `/ timestep_duration`: same per-power objective basis as the other slacks.
                add_to_expression!(cost_of_chain_spill, operations_weight * spill_price[string(s)] / timestep_duration, soc_chain_spill[s, d])
                add_to_expression!(cost_of_chain_borrow, operations_weight * borrow_price[string(s)] / timestep_duration, soc_chain_borrow[s, d])
            end
        end

        # Finally, formulate the objective function as the sum of the costs
        @objective(model, Min, cost_of_investment + cost_of_operations + cost_of_spillage +
                   cost_of_borrow + cost_of_chain_spill + cost_of_chain_borrow)
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
        # Exact-inflow injection. The net RP increment decomposes as Δσ_r = δ_r + Ē_{s,r},
        # where Ē_{s,r} = Σ_h inflow_profile·peak_inflow is representative period r's total
        # (exogenous) inflow and δ_r is the dispatch response. With injection on, the chain
        # blends only δ_r = Δσ_r − Ē_{s,r} and the base period's *real* total inflow
        # E_{s,d} = Σ_{h∈d} inflow·peak_inflow enters exactly on the RHS. Both are pure data
        # (no τ — inflow is already per-step energy, matching the intra balance). `inflow_r`
        # is keyed (id, rep_period) from the seasonal views; `inflow_d` is keyed (id, period)
        # from the raw `profiles` split, so it carries every base day's actual inflow. Both
        # default to 0 for a seasonal asset without an inflow profile (injection is then a
        # no-op for it). When injection is off, the maps are unused and the chain blends the
        # whole increment, reproducing the historical single-matrix / chain-split model.
        inflow_r = Dict{Tuple{String,Int},Float64}()
        inflow_d = Dict{Tuple{String,Int},Float64}()
        if inject && !isempty(S_seas)
            for row in rows(run_query("""
                SELECT id, rep_period, SUM(inflow_profile * peak_inflow) AS e
                FROM (
                    SELECT id, rep_period, inflow_profile, peak_inflow
                    FROM intra_period_seasonal_storage_can_charge_constraint_view
                    UNION ALL
                    SELECT id, rep_period, inflow_profile, peak_inflow
                    FROM intra_period_seasonal_storage_cannot_charge_constraint_view
                )
                GROUP BY id, rep_period
                """))
                inflow_r[(string(row.id), Int(row.rep_period))] = row.e
            end
            for row in rows(run_query("""
                SELECT p.id AS id, p.period AS period, SUM(p.value * s.peak_inflow) AS e
                FROM (SELECT id, period, value FROM profiles WHERE profile_type = 'inflows') p
                JOIN seasonal_storage_assets s ON p.id = s.id
                GROUP BY p.id, p.period
                """))
                inflow_d[(string(row.id), Int(row.period))] = row.e
            end
        end
        ebar(s, r) = get(inflow_r, (string(s), r), 0.0)
        ereal(s, d) = get(inflow_d, (string(s), d), 0.0)

        # Blended dispatch increment for base period d (whole increment when injection is off).
        # Under injection the base-day conservation slacks close the ledger at base-day
        # resolution: spill removes over-accumulated inflow, borrow supplies a shortfall.
        increment(s, d) = sum(
            weight_matrix[d, r]
            *
            (state_of_charge_intra[s, r, H[end]] - state_of_charge_intra_0[s, r]
             - (inject ? ebar(s, r) : 0.0))
            for r in R
        ) + (inject ? ereal(s, d) - soc_chain_spill[s, d] + soc_chain_borrow[s, d] : 0.0)

        # Special treatment of the first period (anchors to state_of_charge_inter_0).
        @constraint(model, [s in S_seas],
            state_of_charge_inter[s, 1] - state_of_charge_inter_0[s] == increment(s, 1)
        )
        @constraint(model, [s in S_seas, d in D[2:end]],
            state_of_charge_inter[s, d] - state_of_charge_inter[s, d-1] == increment(s, d)
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
            # The end-of-horizon tether σ^inter_D = S^0 does two jobs the single-matrix chain
            # fuses: cyclic closure, and gauge-fixing the per-RP additive freedom in σ^intra,0.
            # Its partition-of-unity weights both close the cycle and reconstruct S^0, so the
            # historical (non-injected) model keeps this absolute-reconstruction clause. Exact-
            # inflow injection makes the chain increment-only: column-sum closure telescopes the
            # increment dynamics to σ^inter_D = σ^inter_0 for *any* dispatch (with the cyclic
            # constraint above + the S^0 pin), and the σ^intra,0 gauge then couples to no
            # observable (σ^inter, slacks and cost depend only on the increments), so the clause
            # is dropped — not made feasible, made unnecessary.
            if !inject
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
        interperiod_storage_capacity_data = run_query(
            "SELECT * FROM inter_period_storage_capacity_constraint_view"
        )
        if inject
            # Under injection the band is enforced *hard* on every seasonal asset × base
            # period: the base-day conservation slacks (`soc_chain_spill`/`soc_chain_borrow`)
            # are the recourse, so the hard band is always feasible, and — unlike the soft
            # VOLL slack, which lets the level breach capacity pointwise without touching
            # conservation — spill/borrow correct the mass ledger and therefore every
            # downstream day at once, which is both cheaper and physically correct. The band
            # covers all base periods (not only those with a reservoir-level profile), so the
            # profile-less datasets get a real [0,1]·cap band instead of an unbounded level.
            cap = Dict(string(r.id) => r.capacity_storage_energy
                       for r in rows(run_query("SELECT id, capacity_storage_energy FROM seasonal_storage_assets")))
            bounds = Dict{Tuple{String,Int},Tuple{Float64,Float64}}()
            for row in rows(interperiod_storage_capacity_data)
                bounds[(string(row.id), Int(row.period))] = (row.min_storage_level, row.max_storage_level)
            end
            for s in S_seas, d in D
                mn, mx = get(bounds, (string(s), d), (0.0, 1.0))   # default [0,1] when no reservoir-level profile
                @constraint(model, state_of_charge_inter[s, d] ≥ mn * cap[string(s)])
                @constraint(model, state_of_charge_inter[s, d] ≤ mx * cap[string(s)])
            end
        else
            # Hard inter-period band, applied only when there is no genuine reduction
            # (length(D) == 1: the full-horizon reference / evaluation model). The physical
            # trajectory respects the reservoir bounds by construction, so no slack is needed.
            # For a real reduction the `inject` branch above enforces the band with its own
            # chain-slack recourse. (No rows ⇒ no band, for datasets without reservoir-level
            # profiles.)
            for row in rows(interperiod_storage_capacity_data)
                @constraint(model,
                    state_of_charge_inter[row.id, row.period]
                    ≥
                    row.min_storage_level * row.capacity_storage_energy
                )
                @constraint(model,
                    state_of_charge_inter[row.id, row.period]
                    ≤
                    row.max_storage_level * row.capacity_storage_energy
                )
            end
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
