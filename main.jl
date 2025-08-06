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
inputs = ["tyndp2024",]
inputs_dir = joinpath(dirname(@__FILE__), "inputs")
outputs_dir = joinpath(dirname(@__FILE__), "outputs")
if !isdir(outputs_dir)
    mkpath(outputs_dir)
end

for input in inputs
    @info "Reading experiment configuration"
    run_data = joinpath(inputs_dir, "$(input).csv") |> RunData

    @info "Reading data shared across all experiments"
    connection = DBInterface.connect(DuckDB.DB, ":memory:")
    read_data_from_dir(connection, joinpath(inputs_dir, input))

    output_file = joinpath(outputs_dir, "$(input).csv")

    for seed in seeds
        for run_data_row in run_data |> eachrow
            experiment_data = ExperimentData(run_data_row, input)
            model = Gurobi.Optimizer |> Model
            eval_model = Gurobi.Optimizer |> Model
            result = run_experiment(experiment_data, model, eval_model, connection, seed)
            save_result_to_csv(output_file, result)
        end
    end
end
