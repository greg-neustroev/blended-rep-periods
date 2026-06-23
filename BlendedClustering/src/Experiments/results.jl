function save_result_to_csv(path::String, result::ExperimentResult, time_to_read::Float64)
    row = result |> DataFrame
    row.time_to_read .= time_to_read
    write_header = !isfile(path)
    CSV.write(path, row; append=true, writeheader=write_header)
end

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

function save_variables_to_csv(model, outputs_dir::AbstractString, result_name::AbstractString, seed::Int)
    save_variable_to_csv(
        model,
        :state_of_charge_inter,
        [:period],
        "inter_period_storage_values.csv",
        outputs_dir,
        result_name,
        seed
    )
    save_variable_to_csv(
        model,
        :state_of_charge_intra,
        [:rep_period, :timestep],
        "intra_period_storage_values.csv",
        outputs_dir,
        result_name,
        seed
    )
    save_variable_to_csv(
        model,
        :invested_units,
        Symbol[],
        "invested_units.csv",
        outputs_dir,
        result_name,
        seed
    )
end
