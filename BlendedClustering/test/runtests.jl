using BlendedClustering
using Test
using Random
using LinearAlgebra
using DataFrames
using DuckDB
using DuckDB: DBInterface

const BC = BlendedClustering
const TC = BlendedClustering.TemporalClustering

# A tall synthetic clustering frame: (n_timesteps * length(assets)) features per
# period, kept larger than the RP count so the weight-fitting matrices stay
# well-conditioned (the experiments are always in this tall regime).
function synthetic_clustering_df(; n_periods=12, n_timesteps=4, assets=["a", "b"], seed=0)
  Random.seed!(seed)
  DataFrame([
    (period=p, timestep=t, asset=a, value=rand())
    for p in 1:n_periods for t in 1:n_timesteps for a in assets
  ])
end

@testset "BlendedClustering.jl" begin

  @testset "project_onto_simplex" begin
    # A point already on the simplex is unchanged.
    v = [0.2, 0.3, 0.5]
    @test project_onto_simplex(v) ≈ v atol = 1e-10
    # Arbitrary vectors land on the simplex (nonnegative, sum to one).
    for seed in 1:20
      Random.seed!(seed)
      q = project_onto_simplex(randn(rand(2:10)))
      @test sum(q) ≈ 1
      @test all(q .>= -1e-12)
    end
    # Singleton.
    @test project_onto_simplex([5.0]) == [1.0]
  end

  @testset "projected_gradient_descent! minimizes the projection objective" begin
    Random.seed!(1)
    R = rand(8, 3)
    c = rand(8)
    g = x -> R' * (R * x - c)
    step = 1 / opnorm(R, 2)^2  # the principled 1/L step used in the package
    w, n_iter = projected_gradient_descent!(
      rand(3); gradient=g, projection=project_onto_simplex,
      learning_rate=step, tol=1e-8,
    )
    @test n_iter >= 1                    # reports the realized iteration count
    @test sum(w) ≈ 1
    @test all(w .>= -1e-9)
    # The fitted point is no worse than an arbitrary feasible point.
    feasible = project_onto_simplex(rand(3))
    @test norm(R * w - c) <= norm(R * feasible - c) + 1e-6
  end

  @testset "greedy hull cache is sound (cached == uncached)" begin
    # The projection-cached run must select exactly the same representatives,
    # in the same order, as recomputing every projection from scratch.
    for seed in 1:30
      Random.seed!(seed)
      M = rand(rand(5:18), rand(20:50))
      for k in (3, 5, 8)
        @test BC.greedy_convex_hull(M; n_points=k, cache=true) ==
              BC.greedy_convex_hull(M; n_points=k, cache=false)
      end
    end
  end

  @testset "greedy hull selection is well-formed" begin
    Random.seed!(2)
    M = rand(10, 30)
    idx = BC.greedy_convex_hull(M; n_points=5)
    @test length(idx) == 5
    @test length(unique(idx)) == 5          # no repeats
    @test all(1 .<= idx .<= 30)
    # The mean-furthest point is always picked first.
    mean_vec = vec(sum(M, dims=2) ./ size(M, 2))
    @test idx[1] == argmax([norm(mean_vec - M[:, j]) for j in axes(M, 2)])
  end

  @testset "find_representative_periods: every method" begin
    df = synthetic_clustering_df(seed=3)
    for m in (:k_means, :k_medoids, :hierarchical, :convex_hull, :convex_hull_with_null, :conical_hull)
      res = find_representative_periods(df, 3; method=m)
      @test size(res.rp_matrix, 2) == 3
      @test size(res.weight_matrix, 2) == 3
      @test size(res.clustering_matrix, 2) == 12   # one column per base period
    end
  end

  @testset "fit_rep_period_weights!: weight-class structure" begin
    df = synthetic_clustering_df(seed=4)
    # Convex: nonnegative rows summing to one.
    W = Matrix(fit_rep_period_weights!(
      find_representative_periods(df, 3; method=:k_means); weight_type=:convex))
    @test all(W .>= -1e-8)
    @test all(abs.(sum(W, dims=2) .- 1) .< 1e-5)
    # Conic: nonnegative.
    W = Matrix(fit_rep_period_weights!(
      find_representative_periods(df, 3; method=:k_means); weight_type=:conical))
    @test all(W .>= -1e-8)
    # Sub-unit conic: nonnegative rows summing to at most one.
    W = Matrix(fit_rep_period_weights!(
      find_representative_periods(df, 3; method=:k_means); weight_type=:conical_bounded))
    @test all(W .>= -1e-8)
    @test all(sum(W, dims=2) .<= 1 + 1e-6)
  end

  @testset "sub-unit conic fitting is robust when features == RPs + 1" begin
    # Here rp_matrix is square and, once a zero column is appended for the
    # null-point trick, singular; the pseudoinverse initial guess keeps fitting
    # well-defined where `rp_matrix \\ clustering_matrix` would throw.
    df = synthetic_clustering_df(n_timesteps=2, assets=["a", "b"], seed=7)  # 4 features
    res = find_representative_periods(df, 3; method=:k_means)               # 4x3 RP matrix
    W = Matrix(fit_rep_period_weights!(res; weight_type=:conical_bounded))
    @test all(W .>= -1e-8)
    @test all(sum(W, dims=2) .<= 1 + 1e-6)
  end

  @testset "experiment config schema: optional tol defaults" begin
    mktempdir() do dir
      path = joinpath(dir, "run.csv")
      open(path, "w") do io
        println(io, "n_rep_periods,period_length,clustering_type,weight_type,evaluation_type")
        println(io, "5,24,hull,convex,investment_regret")
      end
      rd = read_run_data(path)
      @test "distance" ∉ names(rd)
      @test "learning_rate" ∉ names(rd)
      @test "niters" ∉ names(rd)
      ed = ExperimentData(first(eachrow(rd)), "gep")
      # The optional `tol` column is absent, so it falls back to the default.
      @test ed.tol == BC.DEFAULT_PGD_TOL
      @test ed.name == "gep_5_24_hull_convex_$(BC.DEFAULT_PGD_TOL)"
    end
  end

  @testset "experiment config schema: explicit tol" begin
    mktempdir() do dir
      path = joinpath(dir, "run.csv")
      open(path, "w") do io
        println(io, "n_rep_periods,period_length,clustering_type,weight_type,tol,evaluation_type")
        println(io, "5,24,hull,convex,0.001,investment_regret")
      end
      ed = ExperimentData(first(eachrow(read_run_data(path))), "gep")
      @test ed.tol == 0.001
      @test ed.name == "gep_5_24_hull_convex_0.001"
    end
  end

  @testset "experiment config schema: inflow_integral_weight default off / suffixes name" begin
    mktempdir() do dir
      # Absent column -> default off, name unchanged.
      path = joinpath(dir, "run.csv")
      open(path, "w") do io
        println(io, "n_rep_periods,period_length,clustering_type,weight_type,evaluation_type")
        println(io, "5,24,hull,convex,investment_regret")
      end
      ed = ExperimentData(first(eachrow(read_run_data(path))), "gep")
      @test ed.inflow_integral_weight == BC.DEFAULT_INFLOW_INTEGRAL_WEIGHT
      @test ed.name == "gep_5_24_hull_convex_$(BC.DEFAULT_PGD_TOL)"
      # Explicit non-zero λ -> field set and name carries the suffix.
      path2 = joinpath(dir, "run2.csv")
      open(path2, "w") do io
        println(io, "n_rep_periods,period_length,clustering_type,weight_type,tol,evaluation_type,inflow_integral_weight")
        println(io, "5,24,hull,convex,0.01,investment_regret,0.5")
      end
      ed2 = ExperimentData(first(eachrow(read_run_data(path2))), "gep")
      @test ed2.inflow_integral_weight == 0.5
      @test ed2.name == "gep_5_24_hull_convex_0.01_inflowint0.5"
    end
  end

  @testset "experiment config schema: fix_every default 1 / suffixes name" begin
    mktempdir() do dir
      # Absent column -> default cadence, name unchanged.
      path = joinpath(dir, "run.csv")
      open(path, "w") do io
        println(io, "n_rep_periods,period_length,clustering_type,weight_type,evaluation_type")
        println(io, "5,24,hull,convex,storage_regret")
      end
      ed = ExperimentData(first(eachrow(read_run_data(path))), "gep")
      @test ed.fix_every == BC.DEFAULT_FIX_EVERY
      @test ed.name == "gep_5_24_hull_convex_$(BC.DEFAULT_PGD_TOL)"
      # Explicit coarser cadence -> field set and name carries the suffix.
      path2 = joinpath(dir, "run2.csv")
      open(path2, "w") do io
        println(io, "n_rep_periods,period_length,clustering_type,weight_type,tol,evaluation_type,fix_every")
        println(io, "5,24,hull,convex,0.01,storage_regret,7")
      end
      ed2 = ExperimentData(first(eachrow(read_run_data(path2))), "gep")
      @test ed2.fix_every == 7
      @test ed2.name == "gep_5_24_hull_convex_0.01_fixevery7"
    end
  end

  @testset "experiment config schema: chain_weights default off / suffixes name" begin
    mktempdir() do dir
      # Absent column -> chain split off, name unchanged.
      path = joinpath(dir, "run.csv")
      open(path, "w") do io
        println(io, "n_rep_periods,period_length,clustering_type,weight_type,evaluation_type")
        println(io, "5,24,hull,convex,storage_regret")
      end
      ed = ExperimentData(first(eachrow(read_run_data(path))), "gep")
      @test ed.chain_weights == BC.DEFAULT_CHAIN_WEIGHTS
      @test ed.name == "gep_5_24_hull_convex_$(BC.DEFAULT_PGD_TOL)"
      # Explicit chain split -> field set and name carries the `_chain` suffix.
      path2 = joinpath(dir, "run2.csv")
      open(path2, "w") do io
        println(io, "n_rep_periods,period_length,clustering_type,weight_type,tol,evaluation_type,chain_weights")
        println(io, "5,24,hull,convex,0.01,storage_regret,true")
      end
      ed2 = ExperimentData(first(eachrow(read_run_data(path2))), "gep")
      @test ed2.chain_weights == true
      @test ed2.name == "gep_5_24_hull_convex_0.01_chain"
    end
  end

  @testset "fit_chain_weights: signed closed-form chain fit" begin
    fcw = BC.TemporalClustering.fit_chain_weights
    # 2 assets x 6 base periods, 3 representatives.
    G = [1.0 3.0 2.0 5.0 0.0 4.0;
         2.0 1.0 0.0 3.0 1.0 2.0]
    sel = [1, 4, 6]
    W = fcw(G, sel)
    @test size(W) == (6, 3)                       # D x R, same shape as W^op
    # Cyclic closure: column sums vanish (1' W^ch = 0), so the reconstructed
    # inter-period trajectory returns to its start for any intra solution.
    @test all(abs.(vec(sum(W, dims=1))) .< 1e-10)
    # Reconstruction matches the least-squares projection of the centered increments.
    Gc = G .- sum(G; dims=2) ./ size(G, 2)
    @test maximum(abs.(Gc[:, sel] * transpose(W) - Gc)) < 1e-9
    # Signed: the fit is free to use negative weights (unlike convex/conical W^op).
    @test any(W .< -1e-6)
  end

  @testset "resolve_input + cross-sweep experiment identity" begin
    ri = BC.Experiments.resolve_input
    @test ri("tyndp/gep") == ("tyndp/gep", "tyndp/gep")
    @test ri(Dict("sweep" => "tyndp/sens_gep", "data" => "tyndp/gep")) == ("tyndp/sens_gep", "tyndp/gep")
    @test ri(Dict("sweep" => "x")) == ("x", "x")          # data defaults to sweep
    # A sensitivity-sweep row and the matching main-sweep row produce the SAME
    # experiment name when both key off the same data, so resume dedups them.
    main = DataFrame(n_rep_periods=40, period_length=24, clustering_type=:hull,
      weight_type=:convex, tol=0.01, evaluation_type=:investment_regret, normalization=:unscaled)
    sens = DataFrame(n_rep_periods=40, period_length=24, clustering_type=:hull,
      weight_type=:convex, tol=0.01, evaluation_type=:investment_regret, normalization=:unscaled,
      inflow_integral_weight=0.0)
    @test ExperimentData(first(eachrow(main)), "tyndp/gep").name ==
          ExperimentData(first(eachrow(sens)), "tyndp/gep").name
  end

  @testset "load_completed_runs (resume support)" begin
    lcr = BC.Experiments.load_completed_runs
    mktempdir() do dir
      f = joinpath(dir, "out.csv")
      @test isempty(lcr(f))                      # missing file -> nothing completed
      open(f, "w") do io
        println(io, "name,seed,objective_value")
        println(io, "ds_5_24_hull_convex_0.01,123,1.0")
        println(io, "ds_5_24_hull_convex_0.01,456,2.0")
      end
      done = lcr(f)
      @test ("ds_5_24_hull_convex_0.01", 123) in done
      @test ("ds_5_24_hull_convex_0.01", 456) in done
      @test ("ds_5_24_hull_convex_0.01", 999) ∉ done   # un-run seed not skipped
    end
  end

  @testset "split_into_periods!" begin
    df = DataFrame([:timestep => 1:4, :value => 5:8])
    split_into_periods!(df; period_duration=2)
    @test df.period == [1, 1, 2, 2]
    @test df.timestep == [1, 2, 1, 2]
  end

  @testset "original_space_centroids" begin
    M = [1.0 2.0 3.0 4.0; 10.0 20.0 30.0 40.0]   # 2 features, 4 periods
    C = TC.original_space_centroids(M, [1, 1, 2, 2], 2)
    @test C ≈ [1.5 3.5; 15.0 35.0]               # per-cluster means of the original columns
    # An empty cluster throws rather than emitting a fabricated all-zero representative.
    @test_throws ArgumentError TC.original_space_centroids(M, [1, 1, 1, 1], 2)
  end

  @testset "economic_scale_vector" begin
    keys = DataFrame(
      timestep=[1, 2, 1], id=["A", "A", "B"], profile_type=["demand", "demand", "availability"])
    fs = Dict(("A", "demand") => 2.0, ("B", "availability") => 3.0)
    @test TC.economic_scale_vector(fs, keys) == [2.0, 2.0, 3.0]
    # A feature with no scale entry is an error, not a silent guess.
    @test_throws ErrorException TC.economic_scale_vector(Dict(("A", "demand") => 2.0), keys)
  end

  # A clustering frame in the (id, profile_type) layout the economic path expects.
  function econ_clustering_df(; n_periods=12, n_timesteps=3, seed=0, constant_block=false)
    Random.seed!(seed)
    rows = NamedTuple[]
    for p in 1:n_periods, t in 1:n_timesteps
      push!(rows, (period=p, timestep=t, id="A", profile_type="demand", value=rand()))
      push!(rows, (period=p, timestep=t, id="B", profile_type="availability", value=rand()))
      if constant_block
        # "C" is identical in every period (constant across columns) -> droppable.
        push!(rows, (period=p, timestep=t, id="C", profile_type="availability", value=0.5))
      end
    end
    DataFrame(rows)
  end

  @testset "economic path: uniform scaling leaves selection invariant, profiles in original units" begin
    df = econ_clustering_df(seed=11)
    fs = Dict(("A", "demand") => 3.0, ("B", "availability") => 3.0)   # one global scalar
    for m in (:convex_hull, :k_medoids, :conical_hull)
      # Reseed so the randomized methods (k_medoids) start from the same state;
      # the invariance is about the scaling, not the RNG.
      Random.seed!(5); base = find_representative_periods(df, 4; method=m)
      Random.seed!(5); econ = find_representative_periods(df, 4; method=m, feature_scale=fs)
      # Same representatives chosen, and profiles stay in the original units (a global
      # scalar cannot change the geometry, so the selected period profiles are identical).
      @test isequal(base.profiles, econ.profiles)
      # The space weights are fitted in is the scaled feature space (×3 on every entry).
      @test econ.clustering_matrix ≈ 3.0 .* base.clustering_matrix
    end
  end

  @testset "economic path: constant feature rows are dropped from selection but kept in profiles" begin
    df = econ_clustering_df(seed=12, constant_block=true)
    fs = Dict(("A", "demand") => 2.0, ("B", "availability") => 2.0, ("C", "availability") => 2.0)
    base = find_representative_periods(df, 3; method=:convex_hull)
    econ = find_representative_periods(df, 3; method=:convex_hull, feature_scale=fs)
    # The constant "C" rows (one per timestep) are pruned from the scaled selection space.
    n_timesteps = 3
    @test size(econ.clustering_matrix, 1) == size(base.clustering_matrix, 1) - n_timesteps
    # But the representative profiles still carry every asset, including "C".
    @test "C" in econ.profiles.id
    @test Set(econ.profiles.id) == Set(base.profiles.id)
  end

  # Register a DataFrame as a DuckDB view, so build_economic_feature_scale can query it.
  reg!(conn, df, name) = DuckDB.register_data_frame(conn, df, name)

  function economic_scale_fixture(; is_investment, τ=1.0, initial_units=2, investment_ids=["g1", "g2"])
    conn = DBInterface.connect(DuckDB.DB, ":memory:")
    reg!(conn, DataFrame(scalar=["timestep_duration"], value=[τ]), "scalars")
    reg!(conn, DataFrame(
        asset=["d1", "g1", "g2", "s1"],
        profile_type=["demand", "availability", "availability", "inflows"]), "assets_profiles")
    reg!(conn, DataFrame(id=["d1"], peak_demand=[100.0]), "demand")
    reg!(conn, DataFrame(id=["s1"], peak_inflow=[50.0]), "assets_storage_seasonal")
    reg!(conn, DataFrame(
        id=["g1", "g2"], technology=["Wind", "Sun"],
        unit_capacity=[10.0, 5.0], initial_units=[initial_units, initial_units]), "assets")
    # ENS is a high-cost slack and must not set V̄.
    reg!(conn, DataFrame(
        technology=["ENS", "Wind", "Sun", "Gas"], variable_cost=[1000.0, 2.0, 0.0, 5.0]),
      "generation_assets")
    if is_investment
        costs = Dict("g1" => 20.0, "g2" => 8.0)
        reg!(conn, DataFrame(id=investment_ids, cost=[costs[i] for i in investment_ids]), "investments")
        reg!(conn, DataFrame(id=investment_ids), "investable_assets")
    else
        reg!(conn, DataFrame(id=String[], cost=Float64[]), "investments")
        reg!(conn, DataFrame(id=String[]), "investable_assets")
    end
    return conn
  end

  @testset "build_economic_feature_scale: investment regime (V̄ excludes ENS, κ_g)" begin
    conn = economic_scale_fixture(is_investment=true, τ=1.0)
    scale, info = TC.build_economic_feature_scale(conn)
    @test info.regime == :investment
    @test info.reference_operational_cost == 5.0                  # max(Wind 2, Sun 0, Gas 5); ENS 1000 excluded
    # demand carries V̄·τ; availability is capital I_g·P_g; inflow carries V̄.
    @test scale[("d1", "demand")] ≈ 5.0 * 1.0 * 100.0
    @test scale[("g1", "availability")] ≈ 20.0 * 10.0            # I_g · unit_capacity
    @test scale[("g2", "availability")] ≈ 8.0 * 5.0
    @test scale[("s1", "inflows")] ≈ 5.0 * 50.0                  # w_E = V̄
    @test info.capital_to_operational_ratio_by_tech["Wind"] ≈ 20.0 / (5.0 * 1.0)
  end

  @testset "build_economic_feature_scale: operations regime (price-free, τ on power blocks)" begin
    conn = economic_scale_fixture(is_investment=false, τ=2.0, initial_units=2)
    scale, info = TC.build_economic_feature_scale(conn)
    @test info.regime == :operations
    @test isnan(info.reference_operational_cost)                 # V̄ not used in operations
    @test scale[("d1", "demand")] ≈ 2.0 * 100.0                  # τ · peak_demand
    @test scale[("g1", "availability")] ≈ 2.0 * 10.0 * 2.0       # τ · unit_capacity · initial_units
    @test scale[("s1", "inflows")] ≈ 50.0                        # inflow is already energy: no τ
  end

  @testset "build_economic_feature_scale: halts on a missing parameter" begin
    # g2 carries an availability profile but has no investment cost I_g in the
    # investment regime, so the scale cannot be built and the function must throw.
    conn = economic_scale_fixture(is_investment=true, τ=1.0, investment_ids=["g1"])
    @test_throws ErrorException TC.build_economic_feature_scale(conn)
  end

  @testset "minmax_rescale" begin
    M = [1.0 3.0 5.0; 10.0 10.0 10.0]            # row 2 is constant across periods
    scaled, keep = TC.minmax_rescale(M)
    @test keep == [true, false]                  # constant row dropped (no range / avoids 0/0)
    @test scaled ≈ [0.0 0.5 1.0]                 # row 1 rescaled: min→0, max→1
  end

  @testset "minmax normalization: selection in [0,1], profiles in original units" begin
    # Original values well outside [0,1] so original units are distinguishable from
    # the min-max-rescaled selection space.
    Random.seed!(21)
    rows = NamedTuple[]
    for p in 1:12, t in 1:3
      push!(rows, (period=p, timestep=t, id="A", profile_type="demand", value=100 + 100rand()))
      push!(rows, (period=p, timestep=t, id="B", profile_type="availability", value=50rand()))
    end
    df = DataFrame(rows)
    res = find_representative_periods(df, 4; method=:convex_hull, minmax=true)
    # Selection space: every feature row rescaled to [0,1].
    @test all(isapprox.(minimum(res.clustering_matrix, dims=2), 0; atol=1e-9))
    @test all(isapprox.(maximum(res.clustering_matrix, dims=2), 1; atol=1e-9))
    # Representative profiles stay in the original units (values far outside [0,1]).
    @test maximum(res.profiles.value) > 1
    # minmax and feature_scale are mutually exclusive.
    @test_throws ArgumentError find_representative_periods(
      df, 4; method=:convex_hull, minmax=true, feature_scale=Dict(("A", "demand") => 1.0))
  end

  # A clustering frame carrying a demand, an availability, and an inflow asset, in the
  # (id, profile_type) layout the inflow-integral augmentation reads.
  function inflow_clustering_df(; n_periods=12, n_timesteps=3, seed=0)
    Random.seed!(seed)
    rows = NamedTuple[]
    for p in 1:n_periods, t in 1:n_timesteps
      push!(rows, (period=p, timestep=t, id="A", profile_type="demand", value=rand()))
      push!(rows, (period=p, timestep=t, id="W", profile_type="availability", value=rand()))
      push!(rows, (period=p, timestep=t, id="H", profile_type="inflows", value=rand()))
    end
    DataFrame(rows)
  end

  @testset "build_inflow_integral_rows" begin
    C = [1.0 2.0; 0.5 0.5; 0.2 0.4]   # rows: one demand, two inflow timesteps of "h"
    keys = DataFrame(timestep=[1, 1, 2], id=["d", "h", "h"],
      profile_type=["demand", "inflows", "inflows"])
    # Unscaled / min-max: dimensionless period mean (Σ_h profile)/H, with τ, E^max
    # cancelling. H = 2 timesteps: (0.5+0.2)/2 = 0.35, (0.5+0.4)/2 = 0.45.
    @test TC.build_inflow_integral_rows(C, keys) ≈ [0.35 0.45]
    # Economic: physical cost-weighted energy τ·(V̄·E^max)·Σ_h. τ=2, V̄·E^max=10:
    # 2·10·(0.5+0.2)=14, 2·10·(0.5+0.4)=18.
    @test TC.build_inflow_integral_rows(C, keys; energy_scale=Dict("h" => 10.0),
      timestep_duration=2.0) ≈ [14.0 18.0]
    # No inflow rows -> an empty (0 × n_periods) block, never an error.
    keys_no_inflow = DataFrame(timestep=[1, 1, 2], id=["d", "d", "d"],
      profile_type=["demand", "demand", "demand"])
    @test size(TC.build_inflow_integral_rows(C, keys_no_inflow)) == (0, 2)
    # In the economic mode a missing inflow scale is an error, not a silent guess.
    @test_throws ErrorException TC.build_inflow_integral_rows(C, keys;
      energy_scale=Dict{String,Float64}())
  end

  @testset "inflow-integral augmentation (unscaled): period-mean row, profiles unchanged" begin
    df = inflow_clustering_df(seed=31)
    λ, H = 0.5, 3   # n_timesteps = 3
    base = find_representative_periods(df, 4; method=:convex_hull)
    aug = find_representative_periods(df, 4; method=:convex_hull, inflow_integral_weight=λ)
    # One integrated-inflow row (a single inflow asset) is appended to the selection/
    # weight-fitting space; everything else keeps its shape.
    @test size(aug.clustering_matrix, 1) == size(base.clustering_matrix, 1) + 1
    @test size(aug.clustering_matrix, 2) == size(base.clustering_matrix, 2)
    # That last row is exactly λ · (Σ_h profile_H)/H — the dimensionless period mean.
    for p in 1:12
      expected = λ * sum(df[(df.id.=="H").&(df.period.==p), :value]) / H
      @test aug.clustering_matrix[end, p] ≈ expected
    end
    # The augmented selection matrix is consistent with rp_matrix (same row count),
    # so weight fitting reconstructs the integrated rows too.
    @test size(aug.rp_matrix, 1) == size(aug.clustering_matrix, 1)
    # Representative profiles carry only the original features, in original units —
    # no fabricated integral rows leak into the model profiles.
    @test Set(zip(aug.profiles.id, aug.profiles.profile_type)) ==
          Set(zip(base.profiles.id, base.profiles.profile_type))
    @test nrow(aug.profiles) == nrow(base.profiles)
  end

  @testset "inflow-integral augmentation (economic): τ·V̄·E^max·Σ row from feature_scale" begin
    df = inflow_clustering_df(seed=34)
    λ, τ = 0.5, 2.0
    fs = Dict(("A", "demand") => 3.0, ("W", "availability") => 2.0, ("H", "inflows") => 5.0)
    aug = find_representative_periods(df, 4; method=:convex_hull,
      feature_scale=fs, inflow_integral_weight=λ, timestep_duration=τ)
    # The economic integral row carries the physical energy λ·τ·(V̄·E^max)·Σ_h, with
    # V̄·E^max read back from feature_scale[("H","inflows")] = 5.0.
    for p in 1:12
      expected = λ * τ * 5.0 * sum(df[(df.id.=="H").&(df.period.==p), :value])
      @test aug.clustering_matrix[end, p] ≈ expected
    end
  end

  @testset "inflow-integral augmentation: λ=0 and no inflows are no-ops" begin
    df = inflow_clustering_df(seed=32)
    base = find_representative_periods(df, 4; method=:convex_hull)
    # λ = 0 leaves the selection space (and everything else) untouched.
    zero_λ = find_representative_periods(df, 4; method=:convex_hull, inflow_integral_weight=0.0)
    @test zero_λ.clustering_matrix == base.clustering_matrix
    # A non-zero λ on data without inflow profiles is a no-op too.
    no_inflow = synthetic_clustering_df(seed=33)   # (timestep, asset) layout, no inflows
    a = find_representative_periods(no_inflow, 3; method=:convex_hull)
    b = find_representative_periods(no_inflow, 3; method=:convex_hull, inflow_integral_weight=0.7)
    @test a.clustering_matrix == b.clustering_matrix
  end

end
