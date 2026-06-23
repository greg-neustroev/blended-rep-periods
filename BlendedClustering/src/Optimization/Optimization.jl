"""
Optimization layer: build the JuMP investment-and-operations model from the
DuckDB views produced by the `Database` module.
"""
module Optimization

using ..Types
using ..Database
using ..Utils

using DuckDB
using JuMP
import Tables: columns, rows

export create_optimization_model!

include("jump_helpers.jl")   # add_term!, as_range
include("model.jl")          # create_optimization_model!

end # module Optimization
