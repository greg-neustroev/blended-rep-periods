@info "Activating the environment"
cd(@__DIR__)
using Pkg: Pkg
Pkg.activate(".")
Pkg.instantiate()

@info "Loading packages"
using BlendedClustering
using DuckDB
using Gurobi
using JuMP
using Random

n_random_seeds = 5
Random.seed!(123)
seeds = rand(1:1000, n_random_seeds)
input_dir = "tyndp2024"
experiments_file = "experiments.csv"

@info "Reading experiment configuration"
base_dir = joinpath(dirname(@__FILE__), "inputs")
run_data = joinpath(base_dir, experiments_file) |> RunData

@info "Reading data shared across all experiments"
connection = DBInterface.connect(DuckDB.DB, ":memory:")
read_data_from_dir(connection, joinpath(base_dir, input_dir))

for seed in seeds
    for run_data_row in run_data |> eachrow
        experiment_data = run_data_row |> ExperimentData
        model = Gurobi.Optimizer |> Model
        run_experiment(experiment_data, model, connection, seed)
    end
end
