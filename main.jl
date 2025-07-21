@info "Activating the environment"
cd(@__DIR__)
using Pkg: Pkg
Pkg.activate(".")
#Pkg.update()
Pkg.instantiate()

@info "Loading the packages"
using DuckDB
using Glob
using Gurobi
using DataFrames
using JuMP
using CSV
import Tables.columns, Tables.rows
using BlendedClustering
using Distances

include("io.jl")
include("helper.jl")

input_dir = "inputs"
n_rep_periods = 1
period_length = 8760
clustering_type = :k_means
distance = Euclidean()
weight_type = :dirac

con = DBInterface.connect(DuckDB.DB, ":memory:")
# read all the data from the CSV files
@info "Reading data from CSV files"
read_data_from_dir(con, input_dir, period_length)

# Creating scalars
operations_weight = DBInterface.execute(con,
                        "SELECT value FROM scalars WHERE scalar='operations_weight'"
                    ) |> first |> first
timestep_duration = DBInterface.execute(con,
                        "SELECT value FROM scalars WHERE scalar='timestep_duration'"
                    ) |> first |> first
# Create indexing sets
N = get_index_set(con, "locations")
X = get_index_set(con, "carriers")
L = get_index_set(con, "transmission_lines")
T = get_index_set(con, "technologies")
A = get_index_set(con, "assets")
A_inv = get_index_set(con, "investable_assets")
A_not = get_index_set(con, "non_investable_assets")
G = get_index_set(con, "generation_assets")
S = get_index_set(con, "storage_assets")
S_ST = get_index_set(con, "short_term_storage_assets")
S_seas = get_index_set(con, "seasonal_storage_assets")
S_seas_in = get_index_set(con, "seasonal_storage_assets_can_charge")
C = get_index_set(con, "conversion_assets")
R = collect(1:n_rep_periods)
H = get_index_set(con, "timesteps")
D = get_index_set(con, "periods")

# Clustering
clustering_df = DuckDB.query(con, "SELECT * FROM profiles ORDER BY period, timestep, id") |> DataFrame

@info "Finding $n_rep_periods $clustering_type representative periods"
clustering_method = if clustering_type ≡ :hull
    if weight_type ≡ :conical
        :conical_hull
    elseif weight_type ≡ :conical_bounded
        :convex_hull_with_null
    else
        :convex_hull
    end
else
    clustering_type
end

clustering_result = find_representative_periods(
    clustering_df, n_rep_periods;
    drop_incomplete_last_period=true,
    method=clustering_method,
    distance=distance,
    init=:kmcen
)

@info "Fitting $(string(weight_type)) weights"
@time begin
    fit_rep_period_weights!(clustering_result; weight_type=weight_type)
end

@info "Reinterpreting the clustering results"
weight = clustering_result.weight_matrix
rp_weight = sum(weight, dims=1)
rp_weight .*= operations_weight
DuckDB.register_data_frame(con, clustering_result.profiles, "rp_profiles")

# Create model
model = Model(Gurobi.Optimizer)

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
@variable(model, flow[L, R, H])

@info "Creating objective"
# Build expressions for costs
investment_data = DBInterface.execute(
    con,
    "SELECT id, unit_capacity, cost FROM investable_assets"
)
cost_of_investment = @expression(
    model,
    isempty(investment_data) ?
    0.0 :
    sum([
        row.cost * row.unit_capacity * invested_units[row.id]
        for row in rows(investment_data)
    ])
)

operations_data = DBInterface.execute(
    con,
    "SELECT id, variable_cost FROM generation_assets"
)
spillage_data = DBInterface.execute(
    con,
    "SELECT id, spillage_cost FROM seasonal_storage_assets"
)

cost_of_operations = @expression(
    model,
    isempty(operations_data) ?
    0.0 :
    sum([
        rp_weight[r] * (
            sum([
                row.variable_cost * power_out[row.id, r, h]
                for row in rows(operations_data)
            ])
            +
            sum([
                row.spillage_cost * spillage[row.id, r, h]
                for row in rows(spillage_data)
            ])
        )
        for r in R, h in H
    ])
)

# Formulate the objective function
@objective(model, Min, cost_of_investment + cost_of_operations)

# Add the constraints
@info "Adding constraints"
## balance constraint
@info "Adding balance constraints"
### first build expressions
power_out_data = DBInterface.execute(
    con,
    """
    SELECT a.id, a.location, t.carrier_out
    FROM assets AS a
    JOIN technologies AS t
    ON a.technology = t.id
    """
)
@expression(model, total_power_out[N, X, R, H], AffExpr(0.0))
for row in rows(power_out_data)
    for r in R, h in H
        total_power_out[row.location, row.carrier_out, r, h] += power_out[row.id, r, h]
    end
end
power_in_data = DBInterface.execute(
    con,
    """
    SELECT a.id, a.location, t.carrier_out as carrier_in
    FROM seasonal_storage_assets_can_charge AS a
    JOIN technologies AS t
    ON a.technology = t.id
    UNION ALL
    SELECT a.id, a.location, t.carrier_out as carrier_in
    FROM short_term_storage_assets AS a
    JOIN technologies AS t
    ON a.technology = t.id
    UNION ALL
    SELECT a.id, a.location, t.carrier_in
    FROM conversion_assets AS a
    JOIN
    (SELECT * FROM technologies NATURAL JOIN technologies_conversion) as t
    ON a.technology = t.id
    """
)
@expression(model, total_power_in[N, X, R, H], AffExpr(0.0))
for row in rows(power_in_data)
    for r in R, h in H
        total_power_in[row.location, row.carrier_in, r, h] += power_in[row.id, r, h]
    end
end
@expression(model, total_flow_in[N, X, R, H], AffExpr(0.0))
@expression(model, total_flow_out[N, X, R, H], AffExpr(0.0))
transmission_line_data = DBInterface.execute(
    con,
    """
    SELECT id, "from", "to", carrier
    FROM transmission_lines
    """
)
for row in rows(transmission_line_data)
    for r in R, h in H
        total_flow_out[row.from, row.carrier, r, h] += flow[row.id, r, h]
        total_flow_in[row.to, row.carrier, r, h] += flow[row.id, r, h]
    end
end

demand_data = DBInterface.execute(
    con,
    """
    SELECT
    d.location, d.carrier, t.rep_period, t.timestep,
    COALESCE(rp.value, 1.0) AS demand_profile,
    d.peak_demand
    FROM (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
    CROSS JOIN
    (SELECT demand.id, lc.location, lc.carrier,
    COALESCE(demand.peak_demand, 0.0) AS peak_demand
    FROM all_locations_and_carriers lc
    LEFT JOIN demand
    ON lc.carrier=demand.carrier AND lc.location = demand.location) AS d
    LEFT JOIN (SELECT * FROM rp_profiles WHERE profile_type = 'demand') AS rp
    ON rp.rep_period = t.rep_period AND rp.timestep = t.timestep AND rp.id = d.id
    ORDER BY
    d.location, d.carrier, t.rep_period, t.timestep
    """
)
for row in rows(demand_data)
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

@info "Adding storage constraints"
short_term_storage_asset_data = DBInterface.execute(
    con,
    """
    SELECT
    s.id, t.rep_period, t.timestep, s.efficiency_in, s.efficiency_out,
    FROM (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
    CROSS JOIN short_term_storage_assets AS s
    ORDER BY
    s.id, t.rep_period, t.timestep
    """
)
@info "Adding intra-period short-term storage constraints"
for row in rows(short_term_storage_asset_data)
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
seasonal_storage_can_charge_asset_data = DBInterface.execute(
    con,
    """
    SELECT
    s.id, t.rep_period, t.timestep, s.efficiency_in, s.efficiency_out, 
    COALESCE(rp.value, 0.0) AS inflow_profile,
    s.peak_inflow,
    FROM (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
    CROSS JOIN seasonal_storage_assets_can_charge AS s
    LEFT JOIN (SELECT * FROM rp_profiles WHERE profile_type = 'inflows') AS rp
    ON rp.rep_period = t.rep_period AND rp.timestep = t.timestep AND rp.id = s.id
    ORDER BY
    s.id, t.rep_period, t.timestep
    """
)
@info "Adding intra-period seasonal storage constraints"
for row in rows(seasonal_storage_can_charge_asset_data)
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
            row.inflow_profile * row.peak_inflow
        )
    end
end
seasonal_storage_cannot_charge_asset_data = DBInterface.execute(
    con,
    """
    SELECT
    s.id, t.rep_period, t.timestep, s.efficiency_out, 
    COALESCE(rp.value, 0.0) AS inflow_profile,
    s.peak_inflow,
    FROM (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
    CROSS JOIN seasonal_storage_assets_cannot_charge AS s
    LEFT JOIN (SELECT * FROM rp_profiles WHERE profile_type = 'inflows') AS rp
    ON rp.rep_period = t.rep_period AND rp.timestep = t.timestep AND rp.id = s.id
    ORDER BY
    s.id, t.rep_period, t.timestep
    """
)
@info "Adding intra-period seasonal storage constraints"
for row in rows(seasonal_storage_cannot_charge_asset_data)
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
            row.inflow_profile * row.peak_inflow
        )
    end
end
@info "Adding inter-period storage constraints"
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
@info "Adding cyclic state of charge constraints"
@constraint(model, [s in S_seas],
    state_of_charge_inter[s, D[end]] == state_of_charge_inter_0[s]
)
@constraint(model, [s in S_ST, r in R],
    state_of_charge_intra[s, r, H[end]] == state_of_charge_intra_0[s, r]
)
@info "Adding initial state of charge constraints"
initial_storage_data = DBInterface.execute(
    con,
    """
    SELECT id, initial_storage_level
    FROM seasonal_storage_assets
    """
)
for row in rows(initial_storage_data)
    @constraint(model, state_of_charge_inter_0[row.id] == row.initial_storage_level)
    @constraint(model, sum(clustering_result.weight_matrix[D[end], r] * state_of_charge_intra_0[row.id, r] for r in R) == row.initial_storage_level)
end
@info "Adding conversion constraints"
conversion_asset_data = DBInterface.execute(
    con,
    """
    SELECT
    id, rep_period, timestep, efficiency_in, efficiency_out
    FROM (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
    CROSS JOIN conversion_assets AS s
    ORDER BY
    id, rep_period, timestep
    """
)
for row in rows(conversion_asset_data)
    @constraint(model,
        power_in[row.id, row.rep_period, row.timestep] * row.efficiency_in
        ==
        power_out[row.id, row.rep_period, row.timestep] / row.efficiency_out
    )
end
@info "Adding maximum power output constraints"
all_assets_data = DBInterface.execute(
    con,
    """
    SELECT a.id, t.rep_period, t.timestep, a.investable,
        COALESCE(rp.value, 1.0) AS availability_profile,
        a.unit_capacity, a.initial_units
    FROM
    (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
    CROSS JOIN
    (SELECT *, i.id IS NOT NULL AS investable FROM assets LEFT JOIN investments i ON assets.id = i.id) AS a
    LEFT JOIN
    (SELECT * FROM rp_profiles WHERE profile_type = 'availability') AS rp
    ON
    rp.rep_period = t.rep_period AND rp.timestep = t.timestep AND rp.id = a.id
    ORDER BY
    a.id, t.rep_period, t.timestep
    """
)
@expression(model, accumulated_capacity[A], AffExpr(0.0))
for row in rows(all_assets_data)
    accumulated_capacity[row.id] = row.unit_capacity * (
        if row.investable
            row.initial_units + invested_units[row.id]
        else
            row.initial_units
        end
    )
    @constraint(model, power_out[row.id, row.rep_period, row.timestep] <= row.availability_profile * accumulated_capacity[row.id])
end
for s in S_ST ∪ S_seas_in, r in R, h in H
    @constraint(model, power_in[s, r, h] <= accumulated_capacity[s])
end
@info "Adding intraperiod maximum state of charge constraints"
intraperiod_storage_capacity_data = DBInterface.execute(
    con,
    """
    SELECT id, t.rep_period, t.timestep, capacity_storage_energy
    FROM
    (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
    CROSS JOIN
    storage_assets
    ORDER BY
    id, t.rep_period, t.timestep
    """
)
for row in rows(intraperiod_storage_capacity_data)
    JuMP.set_upper_bound(
        state_of_charge_intra[row.id, row.rep_period, row.timestep],
        row.capacity_storage_energy
    )
end
@info "Adding interperiod maximum state of charge constraints"
interperiod_storage_capacity_data = DBInterface.execute(
    con,
    """
    WITH
    -- 1. Join profiles to profile names and include period
    raw_profiles AS (
        SELECT
            arl.asset,
            arl.profile_type,
            prl.period,
            prl.value
        FROM
            assets_storage_seasonal_reservoir_levels arl
        JOIN
            profiles_reservoir_levels prl
        ON
            arl.profile_name = prl.id
    ),

    -- 2. Compute average value per asset, profile_type, and period
    avg_profiles AS (
        SELECT
            asset,
            profile_type,
            period,
            AVG(value) AS avg_value
        FROM
            raw_profiles
        GROUP BY
            asset, profile_type, period
    ),

    -- 3. Pivot: get min and max levels side by side per (asset, period)
    pivoted_profiles AS (
        SELECT
            asset AS id,
            period,
            MAX(CASE WHEN profile_type = 'min_storage_level' THEN avg_value END) AS min_storage_level,
            MAX(CASE WHEN profile_type = 'max_storage_level' THEN avg_value END) AS max_storage_level
        FROM
            avg_profiles
        GROUP BY id, period
    ),

    -- 4. Get all combinations of seasonal_storage_assets × all periods (from profiles table)
    all_asset_periods AS (
        SELECT
            sa.id, prl.period, sa.capacity_storage_energy,
        FROM
            seasonal_storage_assets sa
        CROSS JOIN (
            SELECT DISTINCT period FROM profiles_reservoir_levels
        ) prl
    )

    -- 5. Final: join and fill missing min/max with defaults
    SELECT
        ap.id,
        ap.period,
        COALESCE(pp.min_storage_level, 0.0) AS min_storage_level,
        COALESCE(pp.max_storage_level, 1.0) AS max_storage_level,
        ap.capacity_storage_energy
    FROM
        all_asset_periods ap
    LEFT JOIN
        pivoted_profiles pp
    ON
        ap.id = pp.id AND ap.period = pp.period
    ORDER BY
    ap.id, ap.period
    """
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
@info "Adding transmission capacity constraints"
line_capacity_data = DBInterface.execute(
    con,
    """
    SELECT id, t.rep_period, t.timestep, export_capacity, import_capacity
    FROM (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
    CROSS JOIN
    transmission_lines
    ORDER BY
    id, t.rep_period, t.timestep
    """
)
for row in rows(line_capacity_data)
    JuMP.set_lower_bound(flow[row.id, row.rep_period, row.timestep], -row.import_capacity)
    JuMP.set_upper_bound(flow[row.id, row.rep_period, row.timestep], row.export_capacity)
end

# Write the model to a file for debugging
# @info "Saving the model"
# write_to_file(model, "model.lp")

# Solve
@info "Solving the model"
optimize!(model)
