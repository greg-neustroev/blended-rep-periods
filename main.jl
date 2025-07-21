@info "Activating the environment"
cd(@__DIR__)
using Pkg: Pkg
Pkg.activate(".")
Pkg.instantiate()

@info "Loading the packages"
using BlendedClustering
using Gurobi
using JuMP
using StyledStrings

@info "Reading the experiment configuration"
run_data = RunData("inputs/experiments.csv")

for run_data_row in eachrow(run_data)
    experiment_data = run_data_row |> ExperimentData
    @info styled"{info:Running experiment: $(experiment_data.name)}"
    model = Model(Gurobi.Optimizer)
    run_experiment(experiment_data, model)
end
