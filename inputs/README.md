# Experiment input datasets

Each dataset is a directory of CSVs in the format consumed by
`read_data_from_dir` (see `BlendedClustering/src/Database/ingestion.jl`), plus a sibling
`<name>.csv` experiment-configuration file. Datasets are grouped by **data
source** so the experiments draw on more than one independent origin:

```
inputs/
  tyndp/    gep/  p2x/                 European TYNDP 2024 (generation expansion; power-to-X storage)
  sienna/   5bus/                      US, Sienna/PowerSystems test data (PJM 5-bus greenfield GEP)
  nrel/     118bus/                    US, primary NREL-118 PLEXOS data (large-scale dispatch)
  gridmod/  rts/                       US, RTS-GMLC from github.com/gridmod/rts-gmlc
```

`case_studies.toml` lists the datasets to run in `inputs`. An input id is a
relative path such as `nrel/118bus`; the data lives in `inputs/nrel/118bus/`
and its config in `inputs/nrel/118bus.csv`. Every dataset needed to reproduce
the results is committed here — the repository is self-contained.

The NREL, Sienna and gridmod datasets were derived from their respective upstream
sources (see the per-dataset notes below).

## Dataset files

Each dataset directory contains the standard set: `scalars`, `technologies`
(+ `-generation`, `-storage`, `-conversion`), `assets` (+ `-storage`,
`-storage-seasonal`, `-storage-seasonal-reservoir-levels`), `demand`,
`assets-profiles`, `profiles`, `profiles-reservoir-levels`, `transmission-lines`,
and `investments`. Profiles are normalized to `[0,1]`; demand/availability/inflow
are scaled at use by `peak_demand` / `unit_capacity` / `peak_inflow`. All numeric
values are rounded to three decimals (nothing here is measured more finely).

The config CSV has required columns `n_rep_periods, period_length, clustering_type,
weight_type, evaluation_type` and optional columns: `tol` (the PGD tolerance ε,
default `0.01`), `normalization` (`unscaled` (default), `minmax`, or `economic`),
and `cache` (the greedy-hull projection cache, default `true`; set `false` only for
the cached-vs-uncached benchmark). `clustering_type` is one of `k_means`, `k_medoids`,
`hierarchical`, `chronological`, `convex_hull`, `convex_hull_with_null`, or
`conical_hull`; `weight_type` is `dirac` (single assignment), `convex`, `conical`, or
`conical_bounded`; `evaluation_type` is `none`, `investment_regret`, or `storage_regret`.
The first row (`1, <horizon>, …, none`) is the unclustered reference; the remaining rows
sweep clustering/weight types under the dataset's regret metric. Each run also records the environment (Julia/solver
versions, machine, git commit) to `<outputs_dir>/environment.txt`.

## The datasets

Counts are exact. **Gen.** is the number of dispatchable/renewable generators and
**excludes** the per-node energy-not-served (ENS) slack that every dataset carries
(one pseudo-generator per node, priced at VOLL). **Batt.** is short-term
(non-seasonal) storage; **Seas.** is seasonal storage (reservoirs / pumped hydro).
**Conv.** is conversion assets (e.g. electrolyzers, SMR). **Inv.** is the number of
investable assets.

| Source | Dataset | Nodes | Gen. | Batt. | Seas. | Conv. | Lines | Inv. | Time series | Regret |
|---|---|---|---|---|---|---|---|---|---|---|
| TYNDP | gep | 20 | 106 | 0 | 0 | 0 | 44 | 106 | 8760 hourly | investment |
| TYNDP | p2x | 32 | 165 | 25 | 47 | 31 | 92 | 0 | 8760 hourly | storage |
| Sienna | 5bus | 5 | 8 | 0 | 2 | 0 | 7 | 8 | 8760 hourly | investment |
| NREL | 118bus | 118 | 312 | 0 | 15 | 0 | 186 | 0 | 8760 hourly | storage |
| gridmod | rts | 73 | 134 | 1 | 20 | 0 | 120 | 0 | 105 408 (5-min) | storage |

The suite spans two problem classes — **generation expansion** (gep, 5bus) and
**operational storage dispatch** (p2x, 118bus, rts) — over three independent data
sources and a wide size range (5 to 118 nodes, 2 880 to 105 408 timesteps).

- **rts** is the long, high-resolution case: the full leap year at **5-minute**
  resolution (105 408 steps). Representative periods are days, but each day is
  288 five-minute steps, so unlike the other (hourly, 8760-step) datasets it
  exercises the sub-hourly regime. The series length is read from the data, so the
  model handles hourly and sub-hourly resolutions without any code changes.

## Modelling conventions and per-dataset notes

The model is a **pure linear economic dispatch** — no unit commitment, no
reserves, ramping off. It therefore leans on the cheapest/free generation more
than unit-commitment production-cost models do; expect the thermal share a little
higher and renewables/hydro a little lower than published UC results.

**Feasibility penalties** are calibrated consistently with the TYNDP studies:
`ENS variable_cost = borrow_cost` (a last-resort value of lost load) and
`spillage_cost = 0`. Spillage represents inflow passing to downstream unmodelled
reservoirs and is physically free. TYNDP uses ENS = borrow = 1 on its ~0.08
thermal-cost scale; the US datasets use **10 000 $/MWh** (≈ the 10 k EUR European
VOLL, also a standard US value), well above their 35–127 $/MWh thermal so the
optimum never sheds load or borrows.

- **5bus** — the PJM 5-bus system posed as a **greenfield generation-expansion**
  problem. As pure dispatch the system is uninformative: the loads (~71 MW peak) sit
  far below its ~570 MW fleet (free hydro inflow alone exceeds annual demand), so it is
  grossly over-supplied, storage never binds, and every representative-period choice is
  near-optimal (all methods within ~0.1% regret). Demand is scaled ×5 for capacity
  adequacy, and the 8 generators are made **greenfield-investable** (`initial_units = 0`,
  technology-calibrated annualized capex borrowed from the TYNDP gep costs: gas 23.3,
  solar 24.0, wind 48.4, hydro 48.4 k$/MW), while the 2 hydro reservoirs and the ENS
  slacks stay fixed. The build decision is then genuinely temporally sensitive — which
  days represent the year drives the invested mix — so investment regret spreads over
  three orders of magnitude across clustering/weight/normalization (≈1% for hull
  selection to >1000% for k-means at n_rp=10). It is the small, illustrative expansion
  case, and mirrors the held-out TYNDP gep (hull selection wins; unscaled beats economic).
- **118bus** — the NREL-118 system (Peña–Martinez-Anido–Hodge 2018), rebased on the
  **primary NREL PLEXOS input files** and posed as a **high-renewable decarbonization storage
  scenario** (see `nrel/118bus/README.md` for the full derivation). The merit order is first
  restored directly from the source — per-fuel prices and capacity-weighted multi-band
  incremental heat rates — replacing the earlier Sienna conversion's flattened single-price /
  band-1-only costs. Source-faithful, the system is temporally non-discriminating (firm
  thermal ≈ peak, seasonal hydro ~3% of demand, cheap gas backstop → every method within
  ~0.05% regret — genuine physics). So the fleet is placed in a coherent decarbonization
  scenario that makes seasonal storage pivotal: a **$100/t CO₂ price** (from the source
  emission rates), **renewables ×6** (27 GW; the dataset's stated high-penetration purpose),
  the **fossil fleet halved** (firm capacity below peak → scarcity), and the seasonal
  reservoirs turned into **enlarged pumped storage** (7.7 TWh, chargeable, 81% round-trip).
  With storage now binding, discrimination appears on the **weight axis**: `dirac` weights
  give ~3.5% storage regret, blended (`convex`) weights ~1%, and the clustering method is a
  wash — mirroring `tyndp/p2x` and complementing the generation-expansion cases (`tyndp/gep`,
  `sienna/5bus`), where discrimination instead lives on the *clustering* axis. It is the
  large-scale demonstration that blended representative-period weights are essential for
  storage-coupled operation.
- **rts** — converted from gridmod/rts-gmlc. Hydro reservoirs and CSP are seasonal
  storage (inflow-only); the battery is short-term storage. `timestep_duration` is
  5/60 h.
