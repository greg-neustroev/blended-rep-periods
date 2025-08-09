export read_data_from_dir, save_result_to_csv, save_variables_to_csv

function read_data_from_dir(connection, input_dir)
    files = glob("*.csv", input_dir)
    for file_name in files
        table_name = replace(replace(basename(file_name), r"\.csv$" => ""), "-" => "_")
        if startswith(table_name, "profiles")
            # if table name starts with 'profiles', read it as profile data
            read_profile_data(connection, table_name, file_name)
        else
            # otherwise, read it as regular data
            read_data(connection, table_name, file_name)
        end
    end

    create_common_views(connection)
end

function read_data(connection, table_name, file_pattern)
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE TABLE $table_name
        AS SELECT *
        FROM read_csv('$file_pattern', null_padding = true, header = true, union_by_name = true)
        """
    )
end

function read_profile_data(connection, table_name, file_pattern)
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE TABLE $(table_name)_raw
        AS SELECT *
        FROM read_csv('$file_pattern', null_padding = true, header = true, union_by_name = true)
        """
    )
end

function create_dummy_rp_profiles_view(connection)
    rp_profiles = DataFrame(
        rep_period=Int[],
        timestep=Int[],
        id=String[],
        profile_type=String[],
        value=Float64[]
    )
    DuckDB.register_data_frame(connection, rp_profiles, "rp_profiles")
end

function create_index_views(connection)
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW locations AS
        SELECT DISTINCT location AS id
        FROM assets
        """
    )

    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW carriers AS
        SELECT DISTINCT carrier_out AS id FROM technologies
        UNION BY NAME
        SELECT DISTINCT carrier_in AS id FROM technologies_conversion
        """
    )

    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW all_locations_and_carriers AS
        SELECT locations.id AS location, carriers.id AS carrier
        FROM locations, carriers
        """
    )
end

function create_non_storage_views(connection)
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW investable_assets AS
        SELECT *
        FROM assets
        NATURAL JOIN investments
        """
    )
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW non_investable_assets AS
        SELECT assets.*
        FROM assets
        LEFT JOIN investments USING (id)
        WHERE investments.id IS NULL
        """
    )
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW conversion_assets AS
        SELECT *
        FROM assets
        NATURAL JOIN
            (SELECT id AS technology, * EXCLUDE id
            FROM technologies
            NATURAL JOIN technologies_conversion
            WHERE type = 'conversion')
        """
    )
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW generation_assets AS
        SELECT *
        FROM assets
        NATURAL JOIN
            (SELECT id AS technology, * EXCLUDE id
            FROM technologies
            NATURAL JOIN technologies_generation
            WHERE type = 'generation')
        """
    )
end

function create_storage_views(connection)
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW storage_assets AS
        SELECT *
        FROM assets
        NATURAL JOIN
            (SELECT id AS technology, * EXCLUDE id
            FROM technologies
            NATURAL JOIN technologies_storage
            WHERE type = 'storage')
        NATURAL JOIN
        assets_storage
        """
    )
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW seasonal_storage_assets AS
        SELECT *
        FROM
        (SELECT * FROM storage_assets WHERE is_seasonal)
        NATURAL JOIN
        assets_storage_seasonal
        """
    )
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW seasonal_storage_assets_can_charge AS
        SELECT *
        FROM
        (SELECT * FROM storage_assets WHERE is_seasonal)
        NATURAL JOIN
        (SELECT * FROM assets_storage_seasonal WHERE can_charge)
        """
    )
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW seasonal_storage_assets_cannot_charge AS
        SELECT *
        FROM
        (SELECT * FROM storage_assets WHERE is_seasonal)
        NATURAL JOIN
        (SELECT * FROM assets_storage_seasonal WHERE NOT can_charge)
        """
    )
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW short_term_storage_assets AS
        SELECT *
        FROM storage_assets
        WHERE NOT is_seasonal
        """
    )
end

function get_index_set(con, table_name)
    query = "SELECT DISTINCT id FROM $table_name ORDER BY id"
    return columns(DBInterface.execute(con, query)).id
end

function get_scalar(con, scalar_name)
    query = "SELECT value FROM scalars WHERE scalar = '$scalar_name'"
    return DBInterface.execute(con, query) |> first |> first
end

function create_common_views(connection)
    create_index_views(connection)
    create_non_storage_views(connection)
    create_storage_views(connection)
    create_period_independent_views(connection)
end

function create_views(connection, period_length::Int=8760)
    create_profile_data_views(connection, period_length)
    create_reservoir_profile_data_views(connection, period_length)

    create_dummy_rp_profiles_view(connection)
    create_rp_dependent_views(connection)
end

function create_profile_data_views(connection, period_length)
    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW profiles AS
        SELECT
        p.period, p.timestep, a.asset AS id, a.profile_type, p.value
        FROM
        (SELECT 
        ((timestep - 1) // $period_length) + 1 AS period,
        ((timestep - 1) % $period_length) + 1 AS timestep,
        id, value
        FROM
        (UNPIVOT
            (SELECT * FROM profiles_raw)
        ON COLUMNS
            (* EXCLUDE(timestep))
        INTO
            NAME id
            VALUE value
        )) p
        JOIN
        assets_profiles a
        ON
        p.id = a.profile
        ORDER BY
        p.period, p.timestep, a.asset
        """
    )

    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW timesteps AS
        SELECT DISTINCT timestep AS id
        FROM profiles
        """
    )

    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW periods AS
        SELECT DISTINCT period AS id
        FROM profiles
        """
    )
end

function create_reservoir_profile_data_views(connection, period_length)
    has_reservoir_levels =
        DBInterface.execute(
            connection,
            "SELECT EXISTS (SELECT 1 FROM profiles_reservoir_levels_raw) AS has_rows"
        ) |> first |> first
    if has_reservoir_levels
        DBInterface.execute(
            connection,
            """
            CREATE OR REPLACE VIEW profiles_reservoir_levels
            AS 
            (SELECT 
            ((timestep - 1) // $period_length) + 1 AS period,
            ((timestep - 1) % $period_length) + 1 AS timestep,
            id, value
            FROM
            (UNPIVOT
                (SELECT * FROM profiles_reservoir_levels_raw)
            ON COLUMNS
                (* EXCLUDE(timestep))
            INTO
                NAME id
                VALUE value
            ))
            ORDER BY
                period, timestep, id
            """
        )
    else
        DBInterface.execute(
            connection,
            """
            CREATE OR REPLACE VIEW profiles_reservoir_levels AS
            SELECT
                CAST(NULL AS INTEGER) AS period,
                CAST(NULL AS INTEGER) AS timestep,
                CAST(NULL AS TEXT) AS id,
                CAST(NULL AS DOUBLE) AS value
            WHERE FALSE
            """
        )
    end

    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW avg_profiles AS 
        SELECT asset, profile_type, period, AVG(value) AS avg_value
        FROM
        (SELECT arl.asset, arl.profile_type, prl.period, prl.value
        FROM assets_storage_seasonal_reservoir_levels arl
        JOIN profiles_reservoir_levels prl
        ON arl.profile_name = prl.id)
        GROUP BY asset, profile_type, period
        """
    )

    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW pivoted_profiles AS
        SELECT
            asset AS id,
            period,
            MAX(CASE WHEN profile_type = 'min_storage_level' THEN avg_value END)
                AS min_storage_level,
            MAX(CASE WHEN profile_type = 'max_storage_level' THEN avg_value END)
                AS max_storage_level
        FROM avg_profiles
        GROUP BY id, period
        """
    )

    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW all_asset_periods AS
        SELECT sa.id, prl.period, sa.capacity_storage_energy,
        FROM seasonal_storage_assets sa
        CROSS JOIN
        (SELECT DISTINCT period FROM profiles_reservoir_levels) prl
        """
    )

    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW inter_period_storage_capacity_constraint_view AS
        SELECT
            ap.id,
            ap.period,
            COALESCE(pp.min_storage_level, 0.0) AS min_storage_level,
            COALESCE(pp.max_storage_level, 1.0) AS max_storage_level,
            ap.capacity_storage_energy
        FROM all_asset_periods ap
        LEFT JOIN pivoted_profiles pp
        ON ap.id = pp.id AND ap.period = pp.period
        ORDER BY ap.id, ap.period
        """
    )
end

function create_period_independent_views(connection)
    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW investment_cost_objective_view AS
        SELECT id, unit_capacity, cost
        FROM investable_assets
        """
    )

    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW operations_cost_objective_view AS
        SELECT id, variable_cost
        FROM generation_assets
        """
    )

    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW spillage_cost_objective_view AS
        SELECT id, spillage_cost
        FROM seasonal_storage_assets
        """
    )

    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW borrow_cost_objective_view AS
        SELECT id, borrow_cost
        FROM seasonal_storage_assets
        """
    )

    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW power_out_expression_view AS
        SELECT a.id, a.location, t.carrier_out
        FROM assets AS a
        JOIN technologies AS t
        ON a.technology = t.id
        """
    )

    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW power_in_expression_view AS
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

    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW transmission_line_expression_view AS
        SELECT id, "from", "to", carrier
        FROM transmission_lines
        """
    )

    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW initial_storage_constraint_view AS
        SELECT id, initial_storage_level
        FROM seasonal_storage_assets
        """
    )
end

function create_rp_dependent_views(connection)
    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW balance_constraint_view AS
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
    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW intra_period_short_term_storage_constraint_view AS
        SELECT
        s.id, t.rep_period, t.timestep, s.efficiency_in, s.efficiency_out,
        FROM (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
        CROSS JOIN short_term_storage_assets AS s
        ORDER BY
        s.id, t.rep_period, t.timestep
        """
    )
    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW intra_period_seasonal_storage_can_charge_constraint_view AS
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
    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW intra_period_seasonal_storage_cannot_charge_constraint_view AS
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
    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW conversion_constraint_view AS
        SELECT
        id, rep_period, timestep, efficiency_in, efficiency_out
        FROM (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
        CROSS JOIN conversion_assets AS s
        ORDER BY
        id, rep_period, timestep
        """
    )
    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW capacity_constraint_view AS
        SELECT a.id, t.rep_period, t.timestep, a.investable,
            COALESCE(rp.value, 1.0) AS availability_profile,
            a.unit_capacity, a.initial_units
        FROM
        (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
        CROSS JOIN
            (SELECT *, i.id IS NOT NULL AS investable
            FROM assets
            LEFT JOIN investments i
            ON assets.id = i.id)
            AS a
        LEFT JOIN
        (SELECT * FROM rp_profiles WHERE profile_type = 'availability') AS rp
        ON
        rp.rep_period = t.rep_period AND rp.timestep = t.timestep AND rp.id = a.id
        ORDER BY
        a.id, t.rep_period, t.timestep
        """
    )
    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW intra_period_storage_capacity_constraint_view AS
        SELECT id, t.rep_period, t.timestep, capacity_storage_energy
        FROM
        (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
        CROSS JOIN
        storage_assets
        ORDER BY
        id, t.rep_period, t.timestep
        """
    )
    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW intra_period_ramping_constraint_view AS
        SELECT id, t.rep_period, t.timestep, ramping_rate
        FROM (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
        CROSS JOIN
        generation_assets
        WHERE ramping_rate < 1.0 AND t.timestep > 1
        ORDER BY
        id, t.rep_period, t.timestep
        """
    )
    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW inter_period_ramping_constraint_view AS
        SELECT g.id, p.id AS period, ramping_rate
        FROM
        generation_assets g
        CROSS JOIN
        periods p
        WHERE ramping_rate < 1.0 AND p.id > 1
        ORDER BY
        g.id, p.id
        """
    )
    DBInterface.execute(
        connection,
        """
        CREATE OR REPLACE VIEW transmission_line_capacity_constraint_view AS
        SELECT id, t.rep_period, t.timestep, export_capacity, import_capacity
        FROM (SELECT DISTINCT rep_period, timestep FROM rp_profiles) AS t
        CROSS JOIN
        transmission_lines
        ORDER BY
        id, t.rep_period, t.timestep
        """
    )
end

function clustering_type_to_method(clustering_type, weight_type)
    if clustering_type ≡ :hull
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
end

function save_result_to_csv(path::String, result::ExperimentResult, time_to_read::Float64)
    row = result |> DataFrame
    row.time_to_read .= time_to_read
    write_header = !isfile(path)
    CSV.write(path, row; append=true, writeheader=write_header)
end

function save_variable_to_csv(
    model,
    varname::Symbol,
    index_names::Vector{Symbol},
    filename::String,
    outputs_dir::AbstractString,
    result_name::AbstractString,
    seed::Int
)
    subdir = joinpath(outputs_dir, result_name)
    mkpath(subdir)  # ensure directory exists
    header = [:id, index_names..., :variable]
    df = Containers.rowtable(model[varname]; header=header) |> DataFrame
    if isempty(df)
        return
    end

    df.value = value.(df.variable)
    select!(df, Not(:variable))
    df.seed .= seed
    select!(df, Cols(:seed, Not(:seed))) # move seed to the first column

    path = joinpath(subdir, filename)
    write_header = !isfile(path)

    CSV.write(path, df; append=true, writeheader=write_header)
end

function save_variables_to_csv(model, outputs_dir::AbstractString, result_name::AbstractString, seed::Int)
    save_variable_to_csv(
        model,
        :state_of_charge_inter,
        [:period],
        "inter_period_storage_values.csv",
        outputs_dir,
        result_name,
        seed
    )
    save_variable_to_csv(
        model,
        :state_of_charge_intra,
        [:rep_period, :timestep],
        "intra_period_storage_values.csv",
        outputs_dir,
        result_name,
        seed
    )
    save_variable_to_csv(
        model,
        :invested_units,
        Symbol[],
        "invested_units.csv",
        outputs_dir,
        result_name,
        seed
    )
end
