# BlendedClustering

Representative-period clustering for energy-system **investment + operations**
optimization, with *blended* (convex / conical) period weights.

Long-horizon capacity-expansion models are expensive to solve at full temporal
resolution, so the year is reduced to a handful of **representative periods**
(e.g. typical days). This project goes beyond assigning each period to a single
representative: it represents every period as a **blended combination** of
representatives — a convex or conical mix whose weights are fitted by projected
gradient descent — and it can pick representatives by k-means/k-medoids **or** by
a greedy convex/conical-hull selection. It then builds the reduced linear program
in [JuMP](https://jump.dev), solves it, and optionally measures the **regret** of
the reduced solution by re-solving the full-horizon model with the reduced
model's investment or storage decisions fixed.

The repository is both a reusable Julia package (`BlendedClustering/`) and an
experiment harness that runs it across several public power-system datasets.

## How it works

For each experiment the pipeline is:

1. **Ingest** the dataset CSVs into an in-memory [DuckDB](https://duckdb.org)
   database and build the SQL views that express the optimization data model.
2. **Cluster** the time series into `n_rep_periods` representative periods
   (`k_means`, `k_medoids`, or hull-based selection).
3. **Fit weights** mapping every base period to a blend of representatives
   (`dirac` = single assignment, `convex`, `conical`, or `conical_bounded`).
4. **Formulate** the JuMP investment-and-operations LP over the representative
   periods and **solve** it.
5. **Evaluate** (optional) the true cost by fixing the reduced model's
   investment (`investment_regret`) or storage (`storage_regret`) decisions in a
   full-horizon model and re-solving.

## Repository layout

```
.
├── main.jl                 # entry point: runs the case studies in case_studies.toml
├── case_studies.toml       # which datasets / seeds / solver to run
├── Project.toml            # top-level environment (depends on the package below)
├── inputs/                 # datasets: <source>/<name>.csv config + <source>/<name>/ data
│                           #   (see inputs/README.md for dataset details)
├── outputs/                # results written here (git-ignored)
└── BlendedClustering/      # the Julia package
    ├── Project.toml
    ├── src/
    │   ├── BlendedClustering.jl     # parent module: includes & re-exports submodules
    │   ├── Utils/                   # generic helpers (@timed_step)
    │   ├── Types/                   # core data types + read_run_data
    │   ├── Database/                # DuckDB ingestion, SQL views, query accessors
    │   ├── TemporalClustering/      # representative-period selection + weight fitting
    │   ├── Optimization/            # JuMP model construction
    │   └── Experiments/             # orchestration, evaluation, result export, runner
    └── test/                        # unit tests
```

## Requirements

- **[Julia](https://julialang.org) ≥ 1.11** (developed and tested on 1.12).
- **[Gurobi](https://www.gurobi.com)** and a valid license. Gurobi is commercial
  but offers a free [academic license](https://www.gurobi.com/academia/). The
  solver factory is configurable (see `solver` below), but Gurobi is the only one
  wired up out of the box.

All Julia package dependencies are pinned in the checked-in `Manifest.toml`
files and installed by the instantiate step below.

## Installation

```bash
git clone <repository-url>
cd experiments
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Ensure Gurobi is installed and licensed (`Gurobi.jl` picks up your license via
the standard `GUROBI_HOME` / license-file environment). `main.jl` also activates
and instantiates the environment on startup, so a first `julia main.jl` will
install dependencies if you skip the step above.

## Running

### Run the configured case studies

```bash
julia main.jl
```

This runs every case study listed in [`case_studies.toml`](case_studies.toml).
Each entry in `inputs` is a dataset; its `<input>.csv` file defines the parameter
sweep (one experiment per row), and the experiments run for every random seed.

```toml
# case_studies.toml
inputs = ["tyndp/gep", "tyndp/p2x", "sienna/5bus", "sienna/118bus", "gridmod/rts"]

# Draw n_seeds seeds from 1:seed_max, reproducibly:
#   Random.seed!(master_seed); rand(1:seed_max, n_seeds)
# Or set `seeds = [...]` to use fixed seeds instead.
master_seed = 123
n_seeds     = 5
seed_max    = 1000

solver      = "gurobi"      # see BlendedClustering.Experiments.SOLVERS
inputs_dir  = "inputs"      # resolved relative to this file when not absolute
outputs_dir = "outputs"
```

To run a different subset, edit `inputs` (or point `main.jl` at another TOML
file). No code changes are required.

### Run interactively

```julia
using BlendedClustering

# A whole campaign from a config file:
run_case_studies("case_studies.toml")

# Or a single dataset / seed directly, with full control:
run_experiments(["sienna/5bus"]; seeds = [123])
run_experiments(["sienna/5bus"]; seeds = [123], outputs_dir = "/tmp/out")
```

### Experiment configuration (`inputs/<name>.csv`)

Each row is one experiment:

| Column            | Meaning | Values |
|-------------------|---------|--------|
| `n_rep_periods`   | number of representative periods | integer (`1` = unclustered reference) |
| `period_length`   | length of a period, in timesteps | integer (e.g. `24` for daily) |
| `clustering_type` | how representatives are selected | `k_means`, `k_medoids`, `hull` |
| `weight_type`     | how periods blend representatives | `dirac`, `convex`, `conical`, `conical_bounded` |
| `niters`          | max weight-fitting iterations | integer |
| `evaluation_type` | full-horizon regret evaluation | `none`, `investment_regret`, `storage_regret` |

See the files under [`inputs/`](inputs/) for complete examples and
[`inputs/README.md`](inputs/README.md) for the dataset descriptions and modelling
conventions.

## Outputs

For each case study, results are written under `outputs/`:

- `outputs/<input>.csv` — one row per experiment with the objective value,
  projection error, termination status, spillage/borrow totals, and per-stage
  timings (preprocess, cluster, fit weights, formulate, solve, read).
- `outputs/<experiment name>/` — the solved decision variables
  (`invested_units.csv`, `inter_period_storage_values.csv`,
  `intra_period_storage_values.csv`).

## Datasets

Five datasets across three independent sources (European TYNDP, US Sienna,
US RTS-GMLC), ranging from a 5-bus illustrative system to a 73-bus case at
5-minute resolution (105 408 timesteps). See [`inputs/README.md`](inputs/README.md)
for sizes, provenance, and per-dataset modelling notes.

## Tests

```bash
cd BlendedClustering
julia --project=. -e 'using Pkg; Pkg.test()'
```

The unit tests cover the clustering, weight-fitting, and configuration code and
do not require a solver license.

## License

[Apache License 2.0](LICENSE).
