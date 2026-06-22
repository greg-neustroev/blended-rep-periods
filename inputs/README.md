# Experiment input datasets

Each dataset is a directory of CSVs in the format consumed by
`read_data_from_dir` (see `BlendedClustering/src/io.jl`), plus a sibling
`<name>.csv` experiment-configuration file. Datasets are grouped by **data
source** so the experiments draw on more than one independent origin:

```
inputs/
  tyndp/    gep/  p2x/        European TYNDP 2024 (generation expansion; power-to-X storage)
  sienna/   5bus/ 118bus/     US, from the Sienna/PowerSystems test data (PJM 5-bus; NREL-118)
  gridmod/  rts/              US, RTS-GMLC from github.com/gridmod/rts-gmlc
```

`main.jl` lists the datasets to run in `inputs`. An input id is a relative path
such as `sienna/118bus`; the data lives in `inputs/sienna/118bus/` and its config
in `inputs/sienna/118bus.csv`.

The Sienna/gridmod datasets were produced from their upstream sources by a one-off
converter kept in `aux_data/` (not part of the package and not committed). The
upstream sources also live under `aux_data/` (gitignored).

## Dataset files

Each dataset directory contains the standard set: `scalars`, `technologies`
(+ `-generation`, `-storage`, `-conversion`), `assets` (+ `-storage`,
`-storage-seasonal`, `-storage-seasonal-reservoir-levels`), `demand`,
`assets-profiles`, `profiles`, `profiles-reservoir-levels`, `transmission-lines`,
and `investments`. Profiles are normalized to `[0,1]`; demand/availability/inflow
are scaled at use by `peak_demand` / `unit_capacity` / `peak_inflow`. All numeric
values are rounded to three decimals (nothing here is measured more finely).

The config CSV has columns `n_rep_periods, period_length, clustering_type,
weight_type, niters, evaluation_type`. The first row (`1, <horizon>, …, none`) is
the unclustered reference; the remaining rows sweep clustering/weight types under
the dataset's regret metric.

## The datasets

| Source | Dataset | Locations | Generators | Lines | Seasonal storage | Time series | Regret |
|---|---|---|---|---|---|---|---|
| TYNDP | gep | 20 countries | 125 | 43 | 0 | 8760 hourly | investment |
| TYNDP | p2x | 32 countries | 297 | 91 | 47 | 8760 hourly | storage |
| Sienna | 5bus | 5 buses | 10 | 7 | 2 | 8760 hourly | storage |
| Sienna | 118bus | 118 buses | ~310 | 186 | 15 | 8760 hourly | storage |
| gridmod | rts | 73 buses | ~205 | 120 | 20 | **105 408 (5-min)** | storage |

- **rts** is the long, high-resolution case: the full leap year at **5-minute**
  resolution (105 408 steps). Representative periods are days, but each day is
  288 five-minute steps — directly answering reviews that both prior case studies
  were exactly 8760 steps. The series length is read from the data, so pointing
  the converter at the hourly DAY_AHEAD files instead would just work.

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

- **5bus** — the PJM 5-bus loads (~71 MW peak) sit far below its ~570 MW fleet, so
  raw it is grossly over-supplied (free hydro inflow alone exceeds annual demand)
  and storage never binds. **Demand is scaled ×5** so the system is
  capacity-adequate; the dispatch is then ~57% thermal / 23% hydro / 19% wind+solar
  with <1% ENS, and the hydro reservoirs are actively cycled. It is the small,
  illustrative case.
- **118bus** — converted from the NREL-118 data. Assumptions, since the data lacks
  them: per-fuel prices use standard $/MMBTU values; each region's load is split
  equally across its buses (no participation factors in the source); the 15
  dispatchable hydro units (monthly energy budgets) are modelled as seasonal
  storage with monthly inflow, the 28 non-dispatchable units as generation. The
  resulting dispatch is gas-dominated (~80%), consistent with the NREL-118 paper.
- **rts** — converted from gridmod/rts-gmlc. Hydro reservoirs and CSP are seasonal
  storage (inflow-only); the battery is short-term storage. `timestep_duration` is
  5/60 h.
