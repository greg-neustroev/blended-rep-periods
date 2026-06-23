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
