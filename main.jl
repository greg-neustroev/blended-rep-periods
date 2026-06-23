@info "Activating the environment"
cd(@__DIR__)
using Pkg: Pkg
Pkg.activate(".")
Pkg.instantiate()

@info "Loading packages"
using BlendedClustering

# Run the case studies described in `case_studies.toml`. Edit that file to change
# which datasets, seeds, or solver are used; no code changes required. To run a
# single experiment interactively instead, call e.g.
#   run_experiments(["sienna/5bus"]; seeds=[123])
run_case_studies(joinpath(@__DIR__, "case_studies.toml"))
