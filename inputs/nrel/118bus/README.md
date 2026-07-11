# NREL-118 — high-renewable decarbonization scenario

**Base data.** The NREL-118 system — the IEEE 118-bus network with the high-renewable
generation fleet of Peña, Martinez-Anido and Hodge, *"An Extended IEEE 118-Bus Test
System With High Renewable Penetration"*, IEEE TPWRS 33(1):281–289, 2018
([doi:10.1109/TPWRS.2017.2695963](https://doi.org/10.1109/TPWRS.2017.2695963)). The
primary NREL PLEXOS input files (generators, per-fuel prices, multi-band heat rates,
hourly load / solar / wind / hydro time series, CO₂ emission rates) were taken from the
[`evgenytsydenov/ieee118_power_flow_data`](https://github.com/evgenytsydenov/ieee118_power_flow_data/tree/main/data/raw/nrel118)
mirror, which redistributes the original `input-files.zip` and `additional-files-mti-118.zip`
downloaded from NREL.

**Role.** The large-scale (118-node) **operational storage-dispatch** case — the second
storage problem in the base set alongside `tyndp/p2x` (the two generation-expansion cases
are `tyndp/gep` and `sienna/5bus`). Evaluated with `storage_regret`.

**Why a scenario, not the base dispatch.** As-converted, source-faithful NREL-118 is
*temporally non-discriminating*: firm thermal (≈17.3 GW) equals the coincident peak, the
seasonal hydro carries only ~3 % of demand and never binds, and a cheap gas fleet backstops
everything, so every representative-period method lands within ~0.05 % regret (the reduced
solutions differ but the full-year re-dispatch collapses them to the same cost — genuine
physics, not a conversion artifact). To exercise the method on a large system where seasonal
storage is genuinely pivotal, the base fleet is placed in a **deeply decarbonized,
high-renewable, storage-reliant scenario**. Every lever points the same way (decarbonization
+ renewables + storage), so the scenario is a single coherent story rather than a set of
unrelated knobs.

## Deviations from the source (data-only, all applied on top of the primary-source merit order)

| # | Change | From | To | Rationale |
|---|--------|------|----|-----------|
| 1 | **Carbon price $100/t CO₂** folded into thermal `variable_cost` using the source CO₂ emission rates (`Fuels and emission rates.csv`: coal 203.5, gas 118, biomass 130 lb/MMBTU) | — | +`heat_rate × CO₂-rate × $100/t` | Prices the cheap-gas backstop out of the money so stored clean energy has value. |
| 2 | **Renewables ×6** (`assets.csv` `unit_capacity`, solar + wind) | 4.5 GW | 27.1 GW | High renewable penetration — the dataset's explicit purpose — creates large seasonal/diurnal net-load swings that only storage can smooth. |
| 3 | **Fossil fleet halved** (`assets.csv` `unit_capacity`, CC/CT/ST/ICE ×0.5) | 17.3 GW firm | 8.7 GW firm | Retiring half the thermal fleet makes firm capacity fall well below peak, so the system must lean on renewables + stored energy (scarcity). |
| 4 | **Seasonal reservoirs → pumped storage, enlarged ×10** (`assets-storage.csv` `capacity_storage_energy` ×10; `assets-storage-seasonal.csv` `can_charge = true`; `technologies-storage.csv` efficiency 0.9/0.9) | 0.77 TWh, inflow-only | 7.7 TWh, chargeable (81 % round-trip) | A seasonal storage buildout that can absorb renewable surplus and shift it across seasons — the inter-period coupling representative-period selection must reproduce. |

**Restored from the primary source (independent of the scenario):** the merit order —
per-fuel prices and capacity-weighted multi-band incremental heat rates — was restored
directly from the primary NREL PLEXOS data, replacing the earlier Sienna conversion's
flattened single-price / band-1-only costs. See the git history of
`technologies-generation.csv`.

**Unchanged / inherited (faithful to NREL-118):** the 118-node / 186-line topology and line
limits, the hydro / solar / wind / demand profile *shapes* and the seasonal inflow pattern,
and generator nameplate *ratios* within each fuel. Feasibility penalties follow the suite
convention (`ENS variable_cost = borrow_cost = 10 000 $/MWh`, `spillage_cost = 0`).

## Effect

With storage genuinely pivotal, representative-period selection now matters, and the result
is consistent with the suite's storage cases: **the weight type dominates.** At n_rp = 10,
single-assignment (`dirac`) weights give ~3.5 % storage regret while blended (`convex`)
weights cut that to ~1 %; the clustering method (k-means vs the hull selectors) is a wash
within that. This mirrors `tyndp/p2x` and the synthetic development system, and complements
the generation-expansion cases (`tyndp/gep`, `sienna/5bus`) where discrimination instead
lives on the *clustering* axis (hull selection wins). NREL-118 is thus the **large-scale
demonstration that blended representative-period weights are essential for storage-coupled
operation.**

To recover the source-faithful non-discriminating dispatch case, revert deviations 1–4.
