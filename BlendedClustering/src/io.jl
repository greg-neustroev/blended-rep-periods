export read_config, dataframe_to_dict

"""
    keys_to_symbols(dict::AbstractDict{String,Any}; recursive::Bool=true)::Dict{Symbol,Any}

Create a new dictionary that is identical to `dict`, except all of the keys 
are converted from strings to
[symbols](https://docs.julialang.org/en/v1/manual/metaprogramming/#Symbols).
Symbols are [interned](https://en.wikipedia.org/wiki/String_interning) and are
faster to work with when they serve as unique identifiers.
"""
function keys_to_symbols(
    dict::AbstractDict{String,Any};
    recursive::Bool=true
)::Dict{Symbol,Any}
    return Dict(Symbol(k) =>
        if recursive && typeof(v) <: Dict
            keys_to_symbols(v)
        else
            v
        end
                for (k, v) in dict
    )
end

"""
    read_config(config_path::AbstractString)::Dict{Symbol,Any}

Parse the contents of the config file `config_path` into a dictionary. The file
must be in TOML format.
"""
function read_config(config_path::AbstractString)::Dict{Symbol,Any}
    current_dir = pwd()  # current working directory
    full_path = (current_dir, config_path) |> joinpath |> abspath  # full path to the config file

    # read config to a dictionary and change keys to symbols
    config = full_path |> TOML.parsefile |> keys_to_symbols

    # aliases for input config dictionaries 
    data_config = config[:inputs][:data]
    sets_config = config[:inputs][:sets]
    clustering_config = config[:inputs][:clustering]

    # find the input directory
    config_dir = full_path |> dirname  # directory where the config is located
    input_dir = (config_dir, data_config[:dir]) |> joinpath |> abspath  # input data directory

    # read the dataframes from files

    function read_file!(path::AbstractString, key::Symbol, format::Symbol)
        if format ≡ :CSV
            data_config[key] = (path, data_config[key]) |> joinpath |> CSV.File |> DataFrame
            string_columns = findall(col -> eltype(col) <: AbstractString, eachcol(data_config[key]))
            data_config[key][!, string_columns] = Symbol.(data_config[key][!, string_columns])
        elseif format ≡ :TOML
            data_config[key] = (path, data_config[key]) |> joinpath |> TOML.parsefile |> keys_to_symbols
        end
    end

    read_file!(input_dir, :demand, :CSV)
    read_file!(input_dir, :generation_availability, :CSV)
    read_file!(input_dir, :generation, :CSV)
    read_file!(input_dir, :transmission_lines, :CSV)
    read_file!(input_dir, :scalars, :TOML)

    # remove the directory entry as it has been added to the file paths
    # and therefore it is no longer needed
    delete!(data_config, :dir)

    # resolve the sets

    if sets_config[:times] == "auto"
        t_min = min(minimum(data_config[:demand].timestep), minimum(data_config[:generation_availability].timestep))
        t_max = max(maximum(data_config[:demand].timestep), maximum(data_config[:generation_availability].timestep))
        sets_config[:times] = t_min:t_max
    end

    if sets_config[:locations] == "auto"
        sets_config[:locations] =
            data_config[:demand].location ∪
            data_config[:generation_availability].location ∪
            data_config[:generation].location ∪
            data_config[:transmission_lines].from ∪
            data_config[:transmission_lines].to
    end

    if sets_config[:generators] == "auto"
        sets_config[:generators] =
            Tuple.(map(collect, zip(data_config[:generation].location, data_config[:generation].technology)))
    end
    sets_config[:generation_technologies] = unique([g[2] for g ∈ sets_config[:generators]])

    if sets_config[:transmission_lines] == "auto"
        sets_config[:transmission_lines] =
            Tuple.(map(collect, zip(data_config[:transmission_lines].from, data_config[:transmission_lines].to)))
    end

    period_duration = data_config[:scalars][:representative_period_duration]
    split_into_periods!(data_config[:demand]; period_duration)
    split_into_periods!(data_config[:generation_availability]; period_duration)

    clustering_config[:types] = Symbol.(clustering_config[:types])
    clustering_config[:weights] = Symbol.(clustering_config[:weights])
    clustering_config[:distances] = Symbol.(clustering_config[:distances])
    clustering_config[:distances] = map(
        x -> if x ≡ :Euclidean
            Distances.Euclidean()
        elseif x ≡ :cosine
            Distances.CosineDist()
        elseif x ≡ :corr
            Distances.CorrDist()
        elseif x ≡ :cityblock
            Distances.Cityblock()
        elseif x ≡ :Chebyshev
            Distances.Chebyshev()
        else
            throw(DomainError(x, "Unknown distance semimetric"))
        end,
        clustering_config[:distances]
    )

    return config
end

"""
    dataframe_to_dict(
        df::AbstractDataFrame,
        keys::Union{Symbol, Vector{Symbol}},
        value::Symbol
    ) -> Dict

Convert the dataframe `df` to a dictionary using columns `keys` as keys and
`value` as values. `keys` can contain more than one column symbol, in which case
a tuple key is constructed.
"""
function dataframe_to_dict(
    df::AbstractDataFrame,
    keys::Union{Symbol,Vector{Symbol}},
    value::Symbol
)::Dict
    return if typeof(keys) <: AbstractVector
        Dict(Tuple.(eachrow(df[!, keys])) .=> Vector(df[!, value]))
    else
        Dict(Vector(df[!, keys]) .=> Vector(df[!, value]))
    end
end

"""
Returns a `DataFrame` with the values of the variables from the JuMP container
`container`. The column names of the `DataFrame` can be specified for the
indexing columns in `dim_names`, and the name of the data value column by a
Symbol `value_col` e.g. `:Value`.
"""
function jump_container_to_df(container::JuMP.Containers.DenseAxisArray;
    dim_names=Vector{Symbol}(),
    value_col::Symbol=:value)

    if isempty(container)
        return DataFrame()
    end

    if length(dim_names) == 0
        dim_names = [Symbol("dim$i") for i in 1:length(container.axes)]
    end

    if length(dim_names) != length(container.axes)
        throw(ArgumentError("Length of given name list does not fit the number of variable dimensions"))
    end

    dim_names = [dim_names..., value_col]
    df = DataFrame(Containers.rowtable(container))

    for (a, b) in zip(dim_names, names(df))
        if typeof(a) <: Tuple
            select!(df, Not(b), b => [a...]) # expand tuple columns
        else
            select!(df, Not(b), b => a) # rename other columns
        end
    end

    # move the value column to the end
    select!(df, Not(value_col), value_col)

    return df
end
