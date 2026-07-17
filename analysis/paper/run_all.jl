#!/usr/bin/env julia
#
# Runs every analysis/paper/*.jl export script in one process (each script also runs standalone,
# e.g. `julia --project=. analysis/paper/regret.jl`, for regenerating a single artifact). Replaces
# the old monolithic export_summary_csvs.jl's `main()` — same "one command regenerates everything"
# convenience, now split so each table/figure/prose paragraph has its own script.
#
# Order matters only where a later script reads an earlier script's OUTPUT csv rather than raw
# outputs/ (prose_callouts_regret.jl needs regret_summary.csv; prose_callouts_runtime.jl needs
# runtime_breakdown.csv) — both are placed after their dependency below.
#
# Usage: julia --project=. analysis/paper/run_all.jl [outdir]

include(joinpath(@__DIR__, "common.jl"))

println("Exporting paper data CSVs to $OUTDIR")

include(joinpath(@__DIR__, "regret.jl"));                    export_regret()
include(joinpath(@__DIR__, "model_size.jl"));                 export_model_size()
include(joinpath(@__DIR__, "case_study_sizes.jl"));           export_case_study_sizes()
include(joinpath(@__DIR__, "cache_hit_rate.jl"));             export_cache_hit_rate()
include(joinpath(@__DIR__, "sensitivity.jl"));                export_sensitivity()
include(joinpath(@__DIR__, "knockout_ablation.jl"));          export_knockout_ablation()
include(joinpath(@__DIR__, "acceleration_ablation.jl"));      export_acceleration_ablation()
include(joinpath(@__DIR__, "ablation_synth.jl"));             export_ablation()
include(joinpath(@__DIR__, "knob_sensitivity.jl"));           export_knob_sensitivity()
include(joinpath(@__DIR__, "secondary.jl"));                  export_secondary()
include(joinpath(@__DIR__, "curtailment_renewable.jl"));      export_curtailment_renewable()
include(joinpath(@__DIR__, "storage_tracking.jl"));           export_storage_tracking()
include(joinpath(@__DIR__, "runtime_breakdown.jl"));          export_runtime_breakdown()
include(joinpath(@__DIR__, "reservoir_trajectory.jl"));       export_reservoir_trajectory()
include(joinpath(@__DIR__, "feasibility_slack.jl"));          export_feasibility_slack()
include(joinpath(@__DIR__, "weight_mass_excess.jl"));         export_weight_mass_excess()
include(joinpath(@__DIR__, "loss_of_load_by_weight.jl"));     export_loss_of_load_by_weight()
include(joinpath(@__DIR__, "normalization_synthetic.jl"));    export_normalization_synthetic()
include(joinpath(@__DIR__, "normalization_realistic.jl"));    export_normalization_realistic()
include(joinpath(@__DIR__, "normalization_minmax_5bus.jl"));  export_normalization_minmax_5bus()
include(joinpath(@__DIR__, "prose_callouts_regret.jl"));      export_prose_callouts_regret()
include(joinpath(@__DIR__, "prose_callouts_runtime.jl"));     export_prose_callouts_runtime()
include(joinpath(@__DIR__, "provenance.jl"));                 export_provenance()

println("Done.")
