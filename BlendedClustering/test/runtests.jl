using BlendedClustering
using Test
using Random
using LinearAlgebra
using DataFrames

const BC = BlendedClustering

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
    w = projected_gradient_descent!(
      rand(3); gradient=g, projection=project_onto_simplex,
      learning_rate=step, niters=5000,
    )
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
    for m in (:k_means, :k_medoids, :convex_hull, :convex_hull_with_null, :conical_hull)
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
      find_representative_periods(df, 3; method=:k_means); weight_type=:convex, niters=500))
    @test all(W .>= -1e-8)
    @test all(abs.(sum(W, dims=2) .- 1) .< 1e-5)
    # Conic: nonnegative.
    W = Matrix(fit_rep_period_weights!(
      find_representative_periods(df, 3; method=:k_means); weight_type=:conical, niters=500))
    @test all(W .>= -1e-8)
    # Sub-unit conic: nonnegative rows summing to at most one.
    W = Matrix(fit_rep_period_weights!(
      find_representative_periods(df, 3; method=:k_means); weight_type=:conical_bounded, niters=500))
    @test all(W .>= -1e-8)
    @test all(sum(W, dims=2) .<= 1 + 1e-6)
  end

  @testset "sub-unit conic fitting is robust when features == RPs + 1" begin
    # Here rp_matrix is square and, once a zero column is appended for the
    # null-point trick, singular; the pseudoinverse initial guess keeps fitting
    # well-defined where `rp_matrix \\ clustering_matrix` would throw.
    df = synthetic_clustering_df(n_timesteps=2, assets=["a", "b"], seed=7)  # 4 features
    res = find_representative_periods(df, 3; method=:k_means)               # 4x3 RP matrix
    W = Matrix(fit_rep_period_weights!(res; weight_type=:conical_bounded, niters=500))
    @test all(W .>= -1e-8)
    @test all(sum(W, dims=2) .<= 1 + 1e-6)
  end

  @testset "experiment config schema has no distance or learning_rate" begin
    mktempdir() do dir
      path = joinpath(dir, "run.csv")
      open(path, "w") do io
        println(io, "n_rep_periods,period_length,clustering_type,weight_type,niters,evaluation_type")
        println(io, "5,24,hull,convex,1000,investment_regret")
      end
      rd = read_run_data(path)
      @test "distance" ∉ names(rd)
      @test "learning_rate" ∉ names(rd)
      ed = ExperimentData(first(eachrow(rd)), "gep")
      @test ed.name == "gep_5_24_hull_convex_1000"   # no distance/learning_rate fields
    end
  end

  @testset "split_into_periods!" begin
    df = DataFrame([:timestep => 1:4, :value => 5:8])
    split_into_periods!(df; period_duration=2)
    @test df.period == [1, 1, 2, 2]
    @test df.timestep == [1, 2, 1, 2]
  end

end
