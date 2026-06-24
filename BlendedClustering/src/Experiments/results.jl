# Variables exported per experiment: (model variable, index column names, output file).
const EXPORTED_VARIABLES = [
    (:state_of_charge_inter, [:period], "inter_period_storage_values.csv"),
    (:state_of_charge_intra, [:rep_period, :timestep], "intra_period_storage_values.csv"),
    (:invested_units, Symbol[], "invested_units.csv"),
]

"""
    save_result_to_csv(path, result, time_to_read)

Append `result` (with `time_to_read` filled in) as a row to the CSV at `path`,
writing the header only when the file does not yet exist.
"""
function save_result_to_csv(path::String, result::ExperimentResult, time_to_read::Float64)
    row = result |> DataFrame
    row.time_to_read .= time_to_read
    write_header = !isfile(path)
    CSV.write(path, row; append=true, writeheader=write_header)
end

"""
    save_variable_to_csv(model, varname, index_names, filename, outputs_dir, result_name, seed)

Append the solved values of `model[varname]` to `<outputs_dir>/<result_name>/<filename>`,
labelling the index columns with `index_names` and tagging each row with `seed`.
Does nothing if the model is not solved or the variable is empty.
"""
function save_variable_to_csv(
    model,
    varname::Symbol,
    index_names::Vector{Symbol},
    filename::String,
    outputs_dir::AbstractString,
    result_name::AbstractString,
    seed::Int
)
    if JuMP.is_solved_and_feasible(model)
        subdir = joinpath(outputs_dir, result_name)
        mkpath(subdir)  # ensure directory exists
        header = [:id, index_names..., :variable]
        df = Containers.rowtable(model[varname]; header=header) |> DataFrame
        if isempty(df)
            return
        end

        df.value = value.(df.variable)
        select!(df, Not(:variable))
        df.seed .= seed
        select!(df, Cols(:seed, Not(:seed))) # move seed to the first column

        path = joinpath(subdir, filename)
        write_header = !isfile(path)

        CSV.write(path, df; append=true, writeheader=write_header)
    end
end

"""
    save_variables_to_csv(model, outputs_dir, result_name, seed)

Write each variable listed in [`EXPORTED_VARIABLES`](@ref) to its own CSV under
`<outputs_dir>/<result_name>/`.
"""
function save_variables_to_csv(model, outputs_dir::AbstractString, result_name::AbstractString, seed::Int)
    for (varname, index_names, filename) in EXPORTED_VARIABLES
        save_variable_to_csv(model, varname, index_names, filename, outputs_dir, result_name, seed)
    end
end

# Every primal variable container worth dumping in full, with the names of its
# index dimensions (the first axis is always the asset/line id). Saving these
# (plus the duals below) means downstream metrics — curtailment, generation mix,
# storage dispatch, ENS, congestion — are derivable offline without re-solving.
const SOLUTION_VARIABLES = [
    (:invested_units, Symbol[]),
    (:power_out, [:rep_period, :timestep]),
    (:power_in, [:rep_period, :timestep]),
    (:state_of_charge_intra, [:rep_period, :timestep]),
    (:state_of_charge_intra_0, [:rep_period]),
    (:state_of_charge_inter, [:period]),
    (:state_of_charge_inter_0, Symbol[]),
    (:spillage, [:rep_period, :timestep]),
    (:borrow, [:rep_period, :timestep]),
    (:flow, [:rep_period, :timestep]),
]

# Tidy DataFrame of a solved variable container: one row per index tuple, columns
# `[:seed, :id, index_names..., :value]`. Returns `nothing` when the model lacks
# the variable or the container is empty (e.g. `invested_units` on a dispatch-only
# dataset).
function _variable_dataframe(model, varname::Symbol, index_names::Vector{Symbol}, seed::Int)
    haskey(model, varname) || return nothing
    header = [:id, index_names..., :variable]
    df = Containers.rowtable(model[varname]; header=header) |> DataFrame
    isempty(df) && return nothing
    df.value = value.(df.variable)
    select!(df, Not(:variable))
    df.seed .= seed
    select!(df, Cols(:seed, Not(:seed)))
    return df
end

"""
    save_solution_to_arrow(model, outputs_dir, result_name, seed; kind)

Dump the full primal solution of `model` — every variable in
[`SOLUTION_VARIABLES`](@ref) — and, when available, the nodal balance duals
(marginal/clearing prices) to Arrow files under `<outputs_dir>/<result_name>/`,
each prefixed with `kind` (`"reduced"` or `"eval"`). Arrow keeps the full-horizon
arrays compact. Does nothing if `model` is not solved to a feasible point.
"""
function save_solution_to_arrow(
    model,
    outputs_dir::AbstractString,
    result_name::AbstractString,
    seed::Int;
    kind::AbstractString,
)
    JuMP.is_solved_and_feasible(model) || return
    subdir = joinpath(outputs_dir, result_name)
    mkpath(subdir)
    for (varname, index_names) in SOLUTION_VARIABLES
        df = _variable_dataframe(model, varname, index_names, seed)
        df === nothing && continue
        Arrow.write(joinpath(subdir, "$(kind)_$(varname).arrow"), df)
    end
    # Nodal/marginal prices: duals of the retained balance constraints (LMPs).
    # Collect keys and constraints into aligned vectors before reading duals.
    if JuMP.has_duals(model) && haskey(model.ext, :balance_constraints)
        bc = model.ext[:balance_constraints]
        if !isempty(bc)
            ks = collect(keys(bc))
            price_df = DataFrame(
                location=[string(k[1]) for k in ks],
                carrier=[string(k[2]) for k in ks],
                rep_period=[k[3] for k in ks],
                timestep=[k[4] for k in ks],
                price=[JuMP.dual(bc[k]) for k in ks],
            )
            price_df.seed .= seed
            Arrow.write(joinpath(subdir, "$(kind)_nodal_prices.arrow"), price_df)
        end
    end
    return
end
