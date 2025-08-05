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

@info "Reading experiment configuration"
base_dir = joinpath(dirname(@__FILE__), "inputs")
run_data = joinpath(base_dir, "experiments.csv") |> RunData
input_dir = joinpath(base_dir, "tyndp2024")

@info "Reading data from CSV files"
connection = DBInterface.connect(DuckDB.DB, ":memory:")
read_data_from_dir(connection, input_dir)

for seed in seeds
    for run_data_row in eachrow(run_data)
        experiment_data = run_data_row |> ExperimentData
        model = Model(Gurobi.Optimizer)
        run_experiment(experiment_data, model, connection, seed)
    end
end
