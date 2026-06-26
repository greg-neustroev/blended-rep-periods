"""
    get_index_set(connection, table_name) -> Vector

Return the sorted, distinct `id` column of `table_name` — the index set used to
build the optimization model's variables and constraints over that table.
"""
function get_index_set(connection, table_name)
    query = "SELECT DISTINCT id FROM $table_name ORDER BY id"
    return columns(DBInterface.execute(connection, query)).id
end

"""
    get_scalar(connection, scalar_name)

Return the single value stored under `scalar_name` in the `scalars` table.
"""
function get_scalar(connection, scalar_name)
    query = "SELECT value FROM scalars WHERE scalar = '$scalar_name'"
    return DBInterface.execute(connection, query) |> first |> first
end

"""
    get_clustering_profiles(connection) -> DataFrame

Return the full long-format profile table (`period`, `timestep`, `id`,
`profile_type`, `value`) used as the clustering input.
"""
function get_clustering_profiles(connection)
    return DBInterface.execute(connection, "SELECT * FROM profiles") |> DataFrame
end

"""
    count_profile_periods(connection) -> Int

Return the number of base periods present in the profile data.
"""
function count_profile_periods(connection)
    return DBInterface.execute(connection, "SELECT max(period) AS n FROM profiles") |> first |> first
end

"""
    get_inflow_peaks(connection) -> Dict{String,Float64}

Return the peak inflow `E^max_s` of every seasonal-storage asset, keyed by asset id.
Used to scale the per-period inflow profile into physical net-inflow energy when
building the storage-increment data for the signed chain-weight fit.
"""
function get_inflow_peaks(connection)
    result = DBInterface.execute(connection, "SELECT id, peak_inflow FROM seasonal_storage_assets")
    return Dict(string(row.id) => Float64(row.peak_inflow) for row in result)
end

"""
    get_single_period_profiles(connection) -> DataFrame

Return the profile data relabelled as a single representative period (`rep_period`
= 1), ordered by `timestep, id`. Used for the single-representative-period fast path.
"""
function get_single_period_profiles(connection)
    return DBInterface.execute(
        connection,
        "SELECT 1 AS rep_period, timestep, id, profile_type, value FROM profiles ORDER BY timestep, id",
    ) |> DataFrame
end

"""
    get_economic_scaling_data(connection) -> NamedTuple

Fetch the raw parameter inputs needed to build the economic clustering normalization
(see `build_economic_feature_scale`). This is pure data access — it applies no
economics (no slack-technology filtering, no weights). Returns a `NamedTuple` with:

  - `is_investment::Bool`: whether the dataset has investable assets;
  - `timestep_duration::DataFrame`: the `timestep_duration` scalar row(s);
  - `demand`, `inflow`: `(id, peak)` for the demand / inflow profiles;
  - `availability`: `(id, technology, unit_capacity, initial_units, cost)` per
    availability profile (`cost` is the investment cost `I_g`, missing when the asset
    is not investable);
  - `generation_costs`: distinct `(technology, variable_cost)` over generation assets.
"""
function get_economic_scaling_data(connection)
    q(sql) = DBInterface.execute(connection, sql) |> DataFrame
    is_investment = (q("SELECT count(*) AS n FROM investable_assets").n[1]) > 0
    return (
        is_investment = is_investment,
        timestep_duration = q("SELECT value FROM scalars WHERE scalar = 'timestep_duration'"),
        demand = q("""
            SELECT ap.asset AS id, d.peak_demand AS peak
            FROM assets_profiles ap
            LEFT JOIN demand d ON ap.asset = d.id
            WHERE ap.profile_type = 'demand'
        """),
        inflow = q("""
            SELECT ap.asset AS id, s.peak_inflow AS peak
            FROM assets_profiles ap
            LEFT JOIN assets_storage_seasonal s ON ap.asset = s.id
            WHERE ap.profile_type = 'inflows'
        """),
        availability = q("""
            SELECT ap.asset AS id, a.technology AS technology,
                   a.unit_capacity AS unit_capacity, a.initial_units AS initial_units,
                   inv.cost AS cost
            FROM assets_profiles ap
            JOIN assets a ON ap.asset = a.id
            LEFT JOIN investments inv ON a.id = inv.id
            WHERE ap.profile_type = 'availability'
        """),
        generation_costs = q("SELECT DISTINCT technology, variable_cost FROM generation_assets"),
    )
end
