@info "Activating the environment"
cd(@__DIR__)
using Pkg: Pkg
Pkg.activate(".")
Pkg.instantiate()

@info "Loading packages"
using BlendedClustering
using Gurobi
using JuMP
using Random

n_random_seeds = 5
Random.seed!(123)
seeds = rand(1:1000, n_random_seeds)

@info "Reading experiment configuration"
run_data = RunData("inputs/experiments.csv")

for seed in seeds
    for run_data_row in eachrow(run_data)
        experiment_data = run_data_row |> ExperimentData
        model = Model(Gurobi.Optimizer)
        run_experiment(experiment_data, model, seed)
    end
end
