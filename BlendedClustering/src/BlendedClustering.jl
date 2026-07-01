module BlendedClustering

# Submodules are included in dependency order; each `using` brings the
# submodule's exported names into this parent scope so they can be re-exported
# as the package's public API.
include("Utils/Utils.jl")
using .Utils

include("Types/Types.jl")
using .Types

include("Database/Database.jl")
using .Database

include("TemporalClustering/TemporalClustering.jl")
using .TemporalClustering

include("Optimization/Optimization.jl")
using .Optimization

include("Experiments/Experiments.jl")
using .Experiments

# Public API (re-exported from the submodules above).
export ExperimentData, ExperimentResult, read_run_data
export find_representative_periods, split_into_periods!
export fit_rep_period_weights!, projected_gradient_descent!, project_onto_simplex
export read_data_from_dir, save_result_to_csv
export run_experiment, run_experiments, run_case_studies

# Internal names exercised by the test suite as `BlendedClustering.<name>`.
using .TemporalClustering: greedy_convex_hull
using .Types: DEFAULT_PGD_TOL, DEFAULT_NORMALIZATION, DEFAULT_INFLOW_INTEGRAL_WEIGHT, DEFAULT_FIX_EVERY, DEFAULT_CHAIN_WEIGHT_TYPE

end # module BlendedClustering
