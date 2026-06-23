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
