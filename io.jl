function read_data(connection, table_name, file_pattern)
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE TABLE $table_name
        AS SELECT *
        FROM read_csv('$file_pattern', null_padding = true, header = true, union_by_name = true)
        """
    )
end

function read_profile_data(connection, table_name, file_pattern, period_length)
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE TABLE $table_name
        AS 
        (SELECT 
        ((timestep - 1) // $period_length) + 1 AS period,
        ((timestep - 1) % $period_length) + 1 AS timestep,
        id, value
        FROM
        (UNPIVOT
            (SELECT * FROM read_csv_auto('$file_pattern'))
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
end

function read_data_from_dir(connection, input_dir, period_length::Int=8760)
    files = glob("*.csv", input_dir)
    for file_name in files
        table_name = replace(replace(basename(file_name), r"\.csv$" => ""), "-" => "_")
        if startswith(table_name, "profiles")
            # if table name starts with 'profiles', read it as profile data
            read_profile_data(connection, table_name, file_name, period_length)
        else
            # otherwise, read it as regular data
            read_data(connection, table_name, file_name)
        end
    end
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE TABLE profiles AS
        SELECT
        p.period, p.timestep, a.asset AS id, a.profile_type, p.value
        FROM
        profiles p
        JOIN
        assets_profiles a
        ON
        p.id = a.profile
        """
    )
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
        (SELECT id AS technology, * EXCLUDE id FROM technologies NATURAL JOIN technologies_conversion WHERE type = 'conversion')
        """
    )
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW generation_assets AS
        SELECT *
        FROM assets
        NATURAL JOIN
        (SELECT id AS technology, * EXCLUDE id FROM technologies NATURAL JOIN technologies_generation WHERE type = 'generation')
        """
    )
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW storage_assets AS
        SELECT *
        FROM assets
        NATURAL JOIN
        (SELECT id AS technology, * EXCLUDE id FROM technologies NATURAL JOIN technologies_storage WHERE type = 'storage')
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
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE VIEW all_locations_and_carriers AS
        SELECT locations.id AS location, carriers.id AS carrier FROM locations, carriers
        """
    )
end

function get_index_set(con, table_name)
    query = "SELECT DISTINCT id FROM $table_name ORDER BY id"
    return columns(DBInterface.execute(con, query)).id
end
