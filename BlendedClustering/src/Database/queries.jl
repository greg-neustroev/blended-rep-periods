function get_index_set(con, table_name)
    query = "SELECT DISTINCT id FROM $table_name ORDER BY id"
    return columns(DBInterface.execute(con, query)).id
end

function get_scalar(con, scalar_name)
    query = "SELECT value FROM scalars WHERE scalar = '$scalar_name'"
    return DBInterface.execute(con, query) |> first |> first
end
