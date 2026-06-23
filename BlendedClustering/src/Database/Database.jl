"""
Database layer: ingest the CSV inputs into DuckDB, build the SQL views that
express the optimization data model, and provide small query accessors.
"""
module Database

using CSV
using DataFrames
using DuckDB
using Glob
import Tables: columns

export read_data_from_dir
export create_views, create_common_views
export get_index_set, get_scalar

include("ingestion.jl")   # CSV files -> DuckDB tables
include("schema.jl")      # CREATE VIEW definitions (the data model)
include("queries.jl")     # get_index_set, get_scalar

end # module Database
