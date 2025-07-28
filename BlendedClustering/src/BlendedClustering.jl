module BlendedClustering

using Clustering
using CSV
using DataFrames
using Distances
using DuckDB
using Glob
using JuMP
using LinearAlgebra
using Random
using SparseArrays
using Statistics
using StyledStrings
import Tables.columns, Tables.rows

include("data-structures.jl")
include("weight_fitting.jl")
include("cluster.jl")
include("io.jl")
include("optimization.jl")
include("experiment.jl")

end # module BlendedClustering
