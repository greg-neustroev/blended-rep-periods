"""
Experiment orchestration: run a single clustering-and-optimization experiment,
evaluate its solution against the full horizon, and persist the results.
"""
module Experiments

using ..Types
using ..Database
using ..TemporalClustering
using ..Optimization

using Arrow
using CSV
using DataFrames
using DuckDB
using Gurobi
using JuMP
using LinearAlgebra
using Pkg
using Random
using StyledStrings
using TOML

export run_experiment, run_experiments, run_case_studies
export save_result_to_csv, save_solution_to_arrow, save_clustering_artifacts

include("results.jl")      # save_result_to_csv, save_solution_to_arrow, save_clustering_artifacts
include("evaluation.jl")   # evaluate_solution!
include("experiment.jl")   # run_experiment
include("runner.jl")       # run_experiments, run_case_studies

end # module Experiments
