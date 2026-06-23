"""
    read_data_from_dir(connection, input_dir)

Load every `*.csv` file in `input_dir` into a DuckDB table named after the file
(dashes become underscores). Files whose name starts with `profiles` are loaded
as raw profile data (into a `<name>_raw` staging table); the rest are loaded
directly. Finally builds the period-independent ("common") views.
"""
function read_data_from_dir(connection, input_dir)
    files = glob("*.csv", input_dir)
    for file_name in files
        table_name = replace(replace(basename(file_name), r"\.csv$" => ""), "-" => "_")
        # Profile data is staged into `<name>_raw` so the period-aware views can
        # reshape it later; all other data is loaded directly.
        is_profile = startswith(table_name, "profiles")
        read_table(connection, table_name, file_name; raw=is_profile)
    end

    create_common_views(connection)
end

"""
    read_table(connection, table_name, file_pattern; raw=false)

Create (or replace) the DuckDB table for the CSV file(s) matching `file_pattern`.
The table is named `<table_name>_raw` when `raw=true` and `<table_name>`
otherwise.
"""
function read_table(connection, table_name, file_pattern; raw::Bool=false)
    table = raw ? "$(table_name)_raw" : table_name
    DBInterface.execute(connection,
        """
        CREATE OR REPLACE TABLE $table
        AS SELECT *
        FROM read_csv('$file_pattern', null_padding = true, header = true, union_by_name = true)
        """
    )
end
