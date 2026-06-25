@info "Activating the environment"
cd(@__DIR__)
using Pkg: Pkg
Pkg.activate(".")
Pkg.instantiate()

@info "Loading packages"
using BlendedClustering

# Run the case studies described by a TOML config under `configs/`. Pass the config
# path as the first argument; with no argument the whole suite runs. Each config
# (one per case study, plus `sensitivity.toml` and `full_suite.toml`) selects which
# datasets, seeds, and solver are used — no code changes required. Examples:
#   julia main.jl configs/sensitivity.toml   # settle tol/normalization/init first
#   julia main.jl configs/gep.toml           # a single case study
#   julia main.jl                            # the full suite
# Runs resume automatically: experiments already recorded in the output CSV are skipped.
config = isempty(ARGS) ? joinpath(@__DIR__, "configs", "full_suite.toml") : ARGS[1]
@info "Running case studies from $config"
run_case_studies(config)
