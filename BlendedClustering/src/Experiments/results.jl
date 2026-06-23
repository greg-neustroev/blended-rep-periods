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
