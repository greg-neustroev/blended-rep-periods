"""
Experiment orchestration: run a single clustering-and-optimization experiment,
evaluate its solution against the full horizon, and persist the results.
"""
module Experiments

using ..Types
using ..Database
using ..TemporalClustering
using ..Optimization

using CSV
using DataFrames
using DuckDB
using JuMP
using LinearAlgebra
using Random
using StyledStrings

export run_experiment
export save_result_to_csv, save_variables_to_csv

include("results.jl")      # save_result_to_csv, save_variable(s)_to_csv
include("evaluation.jl")   # evaluate_solution!
include("experiment.jl")   # run_experiment

end # module Experiments
