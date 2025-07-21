module BlendedClustering

using CSV
using DataFrames
using Distances
using DuckDB
using Clustering
using LinearAlgebra
using SparseArrays
using Statistics
using JuMP
using Glob
import Tables.columns, Tables.rows

include("data-structures.jl")
include("weight_fitting.jl")
include("cluster.jl")
include("io.jl")
include("optimization.jl")

end # module BlendedClustering
