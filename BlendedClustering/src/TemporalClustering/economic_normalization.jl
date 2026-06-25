# Economic normalization of the clustering features.
#
# Under `normalization = :economic`, the dimensionless [0,1] profiles are rescaled
# into common physical/economic units before clustering, instead of being clustered
# as-is (the `:unscaled` default). Each profile value is a fraction of its asset's
# peak (0 = physically zero, 1 = peak). Features are never centered, so a zero
# profile value stays physically zero — the correct semantics for the origin-anchored
# conic / sub-unit-conic hulls.
#
# THE FEATURE SPACE IS ENERGY. Every feature is the energy delivered/available over
# one time step (MWh), so the space is dimensionally homogeneous and there is no
# power-vs-energy reconciliation between blocks. Demand and availability are powers,
# so they are converted to per-step energy by the step duration τ; inflow is already
# per-step energy (it enters the storage balance without ×τ), so it carries no τ:
#
#   demand_energy        = τ · peak_demand · profile        (MWh)
#   availability_energy  = τ · cap_g · profile              (MWh)   [operations]
#   inflow_energy        = peak_inflow · profile            (MWh)   [already energy]
#
# Energy is then weighted by marginal operating value V̄ (€/MWh) — a single scalar
# across all three blocks, with NO τ in the weighting (τ lives inside the features).
#
#   - Operations regime (no investable assets, e.g. P2X / 5bus): V̄ is common to all
#     blocks and cancels, leaving a price-free, parameter-free normalization. τ is a
#     global scalar over the whole (energy) feature space, so it cancels too.
#   - Investment regime (capacity is a decision variable, e.g. GEP): demand/inflow are
#     valued by operating cost V̄ (× their energy); availability is instead valued by
#     annualized capital I_g (€/MW·yr), which prices *capacity* (power), so it is NOT
#     converted to energy and carries no τ:
#
#       demand:        scale = V̄ · (τ · peak_demand)
#       inflow:        scale = V̄ · peak_inflow
#       availability:  scale = I_g · unit_capacity        (capital × capacity)
#
# Pulling out the common V̄ leaves a single cross-block ratio governing availability
# vs. demand: κ_g = I_g / (V̄·τ). Here τ is NOT a unit reconciliation — it is the
# annualization reconciling capital (€/MW·yr) against per-step operating value
# (€/MWh·step), the same annualization the objective uses for operations weighting.
# So τ appears in exactly one place with a clear economic meaning, never as a free-
# floating power/energy fudge. This makes the normalization robust to τ ≠ 1 (e.g.
# RTS) with no guard, because the feature space is energy by construction.

# Generation "technologies" that are feasibility/penalty slacks (energy-not-served /
# value-of-lost-load), not real supply. They are excluded from the reference
# operational cost V̄: their cost is the price of *failing* to serve (a penalty that
# caps the supply curve), not a dispatch cost on it, so it must not set the economic
# scale of the availability features. Kept as an explicit allow-list — never an
# outlier heuristic — so no genuinely expensive real generator is dropped and no
# VOLL-type asset silently re-enters the reference cost.
const SLACK_GENERATION_TECHNOLOGIES = ("ENS",)

"""
    build_economic_feature_scale(connection) -> (scale, info)

Build the per-feature economic scale applied to the clustering matrix under
`normalization = :economic`.

Returns `(scale, info)` where `scale::Dict{Tuple{String,String},Float64}` maps each
`(asset id, profile_type)` to the absolute-unit factor multiplying its dimensionless
profile, and `info` is a `NamedTuple` of diagnostics (regime, the reference
operational cost V̄ and its [avg, max] real-generator range, τ, and the
capital-to-operational ratio κ_g per technology in the investment regime).

Every value is read from the existing parameter inputs; nothing is invented. If a
profile that needs scaling lacks its parameter (a demand asset without `peak_demand`,
an inflow asset without `peak_inflow`, or — in the investment regime — an
availability asset without an investment cost `I_g`), the function throws rather than
guessing.
"""
function build_economic_feature_scale(connection)
    # All DuckDB access is delegated to the Database layer; this function applies
    # only the economics (slack filtering, weights, τ, κ_g) to the raw parameters.
    data = get_economic_scaling_data(connection)
    is_investment = data.is_investment

    # τ must be a uniform positive scalar. A per-timestep τ_h would turn the
    # power-block τ into a per-row weight (a different, real generalization), so we
    # assert the scalar assumption rather than silently averaging it away.
    if nrow(data.timestep_duration) ≠ 1
        error("Economic normalization: expected exactly one scalar timestep_duration " *
              "(uniform τ), found $(nrow(data.timestep_duration)).")
    end
    τ = Float64(data.timestep_duration.value[1])
    (isfinite(τ) && τ > 0) ||
        error("Economic normalization: timestep_duration τ = $τ is not a positive scalar.")

    scale = Dict{Tuple{String,String},Float64}()
    missing_params = String[]

    # Reference operational cost V̄ = max variable_cost over real (non-slack)
    # generators, with the [avg, max] range over distinct real technologies kept for
    # the sensitivity note. Only needed — and only meaningful — in the investment
    # regime; in operations it is the common factor that cancels.
    reference_operational_cost = NaN
    reference_operational_cost_avg = NaN
    if is_investment
        real_costs = filter(:technology => t -> !(t in SLACK_GENERATION_TECHNOLOGIES), data.generation_costs)
        if nrow(real_costs) == 0
            error("Economic normalization (investment regime): no real (non-slack) " *
                  "generators found to define the reference cost V̄.")
        end
        reference_operational_cost = maximum(real_costs.variable_cost)
        reference_operational_cost_avg = sum(real_costs.variable_cost) / nrow(real_costs)
        reference_operational_cost > 0 ||
            error("Economic normalization (investment regime): reference cost " *
                  "V̄ = $reference_operational_cost ≤ 0.")
    end

    # Marginal operating value V̄ (€/MWh): the single weight on every energy block.
    # In the investment regime it is the real-generator reference cost; in operations
    # it is common to all blocks and cancels, so we use 1 (a parameter-free
    # normalization). τ is NOT a weight — it lives inside the energy features below.
    operating_value = is_investment ? reference_operational_cost : 1.0

    # ---- demand block: energy = τ · peak_demand · profile, valued at V̄ ----------
    for r in eachrow(data.demand)
        if ismissing(r.peak)
            push!(missing_params, "demand asset $(r.id) has no peak_demand")
        else
            scale[(String(r.id), "demand")] = operating_value * τ * Float64(r.peak)
        end
    end

    # ---- inflow block: already per-step energy (no τ), valued at V̄ --------------
    for r in eachrow(data.inflow)
        if ismissing(r.peak)
            push!(missing_params, "inflow asset $(r.id) has no peak_inflow")
        else
            scale[(String(r.id), "inflows")] = operating_value * Float64(r.peak)
        end
    end

    # ---- availability block -----------------------------------------------------
    # operations: operational energy = τ · (unit_capacity · initial_units) · profile,
    #             valued at V̄ (cancels) — same energy basis as demand/inflow.
    # investment: valued by annualized capital I_g (€/MW·yr), which prices *capacity*
    #             (power), so it is NOT converted to energy and carries no τ.
    capital_to_operational_ratio_by_tech = Dict{String,Float64}()
    for r in eachrow(data.availability)
        if is_investment
            if ismissing(r.cost)
                push!(missing_params,
                      "availability asset $(r.id) ($(r.technology)) has no investment " *
                      "cost I_g, required in the investment regime")
                continue
            end
            scale[(String(r.id), "availability")] = Float64(r.cost) * Float64(r.unit_capacity)
            # κ_g = I_g / (V̄·τ): the lone cross-block ratio. τ here is the annualization
            # reconciling capital (€/MW·yr) with per-step operating value (€/MWh·step) —
            # not a power/energy reconciliation (the feature space is energy throughout).
            capital_to_operational_ratio_by_tech[String(r.technology)] =
                Float64(r.cost) / (reference_operational_cost * τ)
        else
            scale[(String(r.id), "availability")] =
                operating_value * τ * Float64(r.unit_capacity) * Float64(r.initial_units)
        end
    end

    if !isempty(missing_params)
        error("Economic normalization cannot proceed; missing parameters:\n  - " *
              join(missing_params, "\n  - "))
    end

    info = (
        regime = is_investment ? :investment : :operations,
        timestep_duration = τ,
        reference_operational_cost = reference_operational_cost,
        reference_operational_cost_avg = reference_operational_cost_avg,
        slack_technologies = SLACK_GENERATION_TECHNOLOGIES,
        capital_to_operational_ratio_by_tech = capital_to_operational_ratio_by_tech,
        has_inflow = nrow(data.inflow) > 0,
    )
    return scale, info
end

"""
    economic_scale_vector(feature_scale, keys) -> Vector{Float64}

Expand the `(asset id, profile_type) => scale` dictionary into the per-row scale
vector aligned with the clustering matrix, using the `keys` table (one row per
feature, with `:id` and `:profile_type` columns). Throws — rather than guessing —
if any feature has no scale entry, listing the offending `(id, profile_type)` pairs.
"""
function economic_scale_vector(feature_scale::AbstractDict, keys::AbstractDataFrame)
    scale_vector = Vector{Float64}(undef, nrow(keys))
    missing_keys = Set{Tuple{String,String}}()
    for (i, row) in enumerate(eachrow(keys))
        k = (String(row.id), String(row.profile_type))
        s = get(feature_scale, k, nothing)
        if s ≡ nothing
            push!(missing_keys, k)
        else
            scale_vector[i] = s
        end
    end
    isempty(missing_keys) ||
        error("Economic normalization: no scale for features $(sort(collect(missing_keys)))")
    return scale_vector
end

"""
    log_economic_scale_info(info)

Log the dataset-level economic-normalization parameters (regime, the reference
operational cost V̄ and its [avg, max] real-generator range, τ, and the
capital-to-operational ratio κ_g per technology). Computed by
[`build_economic_feature_scale`](@ref); nothing here is fabricated.
"""
function log_economic_scale_info(info::NamedTuple)
    @info "Economic normalization" regime = info.regime τ = info.timestep_duration V̄ = info.reference_operational_cost V̄_avg_to_max = (info.reference_operational_cost_avg, info.reference_operational_cost)
    if !isempty(info.capital_to_operational_ratio_by_tech)
        ranked = sort(collect(info.capital_to_operational_ratio_by_tech); by = x -> x[2])
        @info "  κ_g = I_g/(V̄·τ) per technology" pairs = ranked
    end
    return nothing
end

"""
    log_scaled_selection_diagnostics(profile_types, keep, selection_matrix)

Emit (never fabricate) the geometry diagnostics for the scaled selection space: how
many constant feature rows were dropped, and — per block and overall — the Frobenius
share (which block dominates the hull) and the distance to the origin (minimum entry
and minimum column norm, the conic-vs-convex tell).

`profile_types` is the block label of every original feature row; `keep` marks the
rows that vary across periods (constant rows are dropped); `selection_matrix` is the
scaled, pruned matrix actually clustered (its rows correspond to `profile_types[keep]`).
"""
function log_scaled_selection_diagnostics(
    profile_types::AbstractVector,
    keep::AbstractVector{Bool},
    selection_matrix::AbstractMatrix{Float64},
)
    n_total = length(keep)
    n_dropped = n_total - count(keep)
    @info "  constant feature rows dropped (zero range)" dropped = n_dropped of = n_total

    kept_types = collect(profile_types)[keep]
    total_sq = sum(abs2, selection_matrix)
    for pt in sort(unique(kept_types))
        rows = findall(==(pt), kept_types)
        block = @view selection_matrix[rows, :]
        share = total_sq > 0 ? sum(abs2, block) / total_sq : 0.0
        col_norms = [norm(@view selection_matrix[rows, c]) for c in axes(selection_matrix, 2)]
        @info "  block $pt" rows = length(rows) frobenius_share = round(share; digits = 3) min_entry = round(minimum(block); digits = 3) min_col_norm = round(minimum(col_norms); digits = 2)
    end
    overall = [norm(@view selection_matrix[:, c]) for c in axes(selection_matrix, 2)]
    @info "  distance to origin (column norm)" min = round(minimum(overall); digits = 2) median = round(sort(overall)[cld(end, 2)]; digits = 2)
    return nothing
end
