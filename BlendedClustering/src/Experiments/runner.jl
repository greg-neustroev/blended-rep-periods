# Solver name (as written in a case-studies file) -> JuMP optimizer factory.
const SOLVERS = Dict("gurobi" => Gurobi.Optimizer)

"""
    run_experiments(inputs::AbstractVector; seeds, optimizer=Gurobi.Optimizer,
                    inputs_dir="inputs", outputs_dir="outputs")

Run every experiment for each dataset in `inputs`. For each input, the matching
`<input>.csv` configuration (under `inputs_dir`) defines the parameter sweep:
one experiment is run per (seed, configuration row) pair. Results are appended to
`<outputs_dir>/<input>.csv` and the solved variables are written under
`<outputs_dir>/<experiment name>/`.

  - `inputs`: dataset identifiers, e.g. `["tyndp/gep", "sienna/5bus"]`.
  - `seeds`: random seeds to run each configuration with.
  - `optimizer`: a JuMP-compatible optimizer factory.
  - `inputs_dir`, `outputs_dir`: input/output root directories.
"""
function run_experiments(
    inputs::AbstractVector{<:AbstractString};
    seeds::AbstractVector{<:Integer},
    optimizer=Gurobi.Optimizer,
    inputs_dir::AbstractString="inputs",
    outputs_dir::AbstractString="outputs",
)
    isdir(outputs_dir) || mkpath(outputs_dir)
    record_environment(outputs_dir)

    for input in inputs
        @info "Reading experiment configuration"
        run_data = read_run_data(joinpath(inputs_dir, "$(input).csv"))

        @info "Reading data shared across all experiments"
        connection = DBInterface.connect(DuckDB.DB, ":memory:")
        output_file = joinpath(outputs_dir, "$(input).csv")
        mkpath(dirname(output_file))

        # Resume support: skip any (experiment name, seed) already recorded in the
        # output CSV so an interrupted chain picks up where it left off without
        # re-solving finished runs. This is most valuable for the single-period
        # reference, which is expensive and identical across the whole sweep.
        completed = load_completed_runs(output_file)

        # Every configuration is run for every seed. Hull selection and the
        # single-period reference are deterministic, so their regret/selection is
        # identical across seeds — but wall-clock timing is not, so the repeats
        # supply the runtime variance needed for the runtime/Pareto error bars.
        for seed in seeds
            experiments = [ExperimentData(row, input) for row in eachrow(run_data)]
            pending = [ed for ed in experiments if (ed.name, seed) ∉ completed]
            if isempty(pending)
                @info "All experiments for seed $seed already completed; skipping"
                continue
            end
            # Read the dataset once per seed, but only when something still needs running.
            time_to_read = @elapsed read_data_from_dir(connection, joinpath(inputs_dir, input))
            for experiment_data in pending
                model = Model(optimizer)
                eval_model = Model(optimizer)
                result, clustering_result = run_experiment(experiment_data, model, eval_model, connection, seed)
                save_result_to_csv(output_file, result, time_to_read)
                save_variables_to_csv(model, outputs_dir, result.name, seed)
                # Full primal solution + nodal-price duals for both the reduced and
                # the full-horizon evaluation model (the latter is a no-op when the
                # experiment skips evaluation, since it is then never solved).
                save_solution_to_arrow(model, outputs_dir, result.name, seed; kind="reduced")
                save_solution_to_arrow(eval_model, outputs_dir, result.name, seed; kind="eval")
                # Clustering-side artifacts (weights, assignment, selected RPs,
                # spectrum, residuals, PGD/cache diagnostics) for offline analysis.
                save_clustering_artifacts(clustering_result, outputs_dir, result.name, seed)
                push!(completed, (result.name, seed))
            end
        end
    end
    return nothing
end

"""
    run_case_studies(config_file::AbstractString)

Run the case studies described by the TOML file at `config_file`. Each entry in
`inputs` is a case study (a dataset, e.g. `"tyndp/gep"`); the experiments for
each are driven by its `<input>.csv` parameter sweep.

Recognised keys: `inputs` (required), `solver`, `inputs_dir`, `outputs_dir`, and
either an explicit `seeds` array or the trio `master_seed` / `n_seeds` /
`seed_max` used to draw `n_seeds` reproducible seeds from `1:seed_max`. Relative
`inputs_dir`/`outputs_dir` are resolved against the config file's directory.
"""
function run_case_studies(config_file::AbstractString)
    config = TOML.parsefile(config_file)
    base_dir = dirname(abspath(config_file))

    inputs = config["inputs"]
    seeds = haskey(config, "seeds") ? config["seeds"] : draw_seeds(config)

    solver = lowercase(get(config, "solver", "gurobi"))
    haskey(SOLVERS, solver) ||
        error("Unsupported solver \"$solver\"; known solvers: $(join(sort!(collect(keys(SOLVERS))), ", "))")
    optimizer = SOLVERS[solver]

    inputs_dir = resolve_path(base_dir, get(config, "inputs_dir", "inputs"))
    outputs_dir = resolve_path(base_dir, get(config, "outputs_dir", "outputs"))

    return run_experiments(inputs; seeds, optimizer, inputs_dir, outputs_dir)
end

# Draw `n_seeds` reproducible seeds from `1:seed_max`, seeded with `master_seed`.
function draw_seeds(config::AbstractDict)
    master_seed = config["master_seed"]
    n_seeds = config["n_seeds"]
    seed_max = get(config, "seed_max", 1000)
    Random.seed!(master_seed)
    return rand(1:seed_max, n_seeds)
end

# Resolve `path` against `base` unless it is already absolute.
resolve_path(base::AbstractString, path::AbstractString) =
    isabspath(path) ? path : normpath(joinpath(base, path))

# Read the (experiment name, seed) pairs already recorded in `output_file`, so a
# resumed run can skip them. Each finished experiment appends one such row, so a
# present (name, seed) means that run completed. Returns an empty set when the
# file does not exist yet.
function load_completed_runs(output_file::AbstractString)
    completed = Set{Tuple{String,Int}}()
    isfile(output_file) || return completed
    df = CSV.read(output_file, DataFrame)
    ("name" in names(df) && "seed" in names(df)) || return completed
    for row in eachrow(df)
        (ismissing(row.name) || ismissing(row.seed)) && continue
        push!(completed, (String(row.name), Int(row.seed)))
    end
    return completed
end

# Record the run environment — Julia/solver/package versions, machine, and the
# experiments-repo git commit — to `<outputs_dir>/environment.txt` for
# reproducibility (so reported runtimes are tied to a known machine and code).
function record_environment(outputs_dir::AbstractString)
    isdir(outputs_dir) || mkpath(outputs_dir)
    cpu = isempty(Sys.cpu_info()) ? "unknown" : Sys.cpu_info()[1].model
    commit = try
        strip(read(`git -C $(@__DIR__) rev-parse HEAD`, String))
    catch
        "unknown"
    end
    open(joinpath(outputs_dir, "environment.txt"), "w") do io
        println(io, "[julia]")
        println(io, "version = \"", VERSION, "\"")
        println(io, "\n[system]")
        println(io, "machine = \"", Sys.MACHINE, "\"")
        println(io, "cpu = \"", strip(cpu), "\"")
        println(io, "threads = ", Sys.CPU_THREADS)
        println(io, "total_memory_gb = ", round(Sys.total_memory() / 2^30; digits=1))
        println(io, "hostname = \"", gethostname(), "\"")
        println(io, "\n[git]")
        println(io, "experiments_commit = \"", commit, "\"")
        println(io, "\n[packages]")
        for (_, dep) in Pkg.dependencies()
            dep.name in ("JuMP", "Gurobi", "DuckDB", "Arrow", "Clustering", "Distances") || continue
            println(io, dep.name, " = \"", dep.version, "\"")
        end
    end
    return
end
