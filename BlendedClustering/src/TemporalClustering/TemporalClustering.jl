"""
Temporal clustering: reduce the full time horizon to a small set of
representative periods and fit their weights. Named `TemporalClustering` rather
than `Clustering` to avoid clashing with the external `Clustering.jl` dependency.
"""
module TemporalClustering

using ..Types
using ..Database  # data-access layer: all DuckDB queries live here, not in this module

using Clustering: kmeans, kmedoids, hclust, cutree   # restrict import: Clustering also exports `ClusteringResult`, which would clash with `Types.ClusteringResult`
using DataFrames
using Distances
using LinearAlgebra
using SparseArrays
using Statistics

export find_representative_periods, split_into_periods!
export fit_rep_period_weights!, projected_gradient_descent!, project_onto_simplex
export cluster_using_experiment_data
export build_economic_feature_scale

include("weight_fitting.jl")          # projected gradient descent + weight fitting
include("economic_normalization.jl")  # economic feature scaling for clustering
include("clustering.jl")              # representative-period selection

end # module TemporalClustering
