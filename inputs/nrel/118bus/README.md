# NREL-118 — provenance and deviations

**Base data.** The NREL-118 system — the IEEE 118-bus network with the high-renewable
generation fleet of Peña, Martinez-Anido and Hodge, *"An Extended IEEE 118-Bus Test
System With High Renewable Penetration"*, IEEE TPWRS 33(1):281–289, 2018
([doi:10.1109/TPWRS.2017.2695963](https://doi.org/10.1109/TPWRS.2017.2695963)). The
primary NREL PLEXOS input files (generators, per-fuel prices, multi-band heat rates,
hourly load / solar / wind / hydro time series) were taken from the
[`evgenytsydenov/ieee118_power_flow_data`](https://github.com/evgenytsydenov/ieee118_power_flow_data/tree/main/data/raw/nrel118)
mirror, which redistributes the original `input-files.zip` and `additional-files-mti-118.zip`
downloaded from NREL.

**Role.** The large-scale (118-node) **operational storage-dispatch** case — the second
dispatch problem in the base set alongside `tyndp/p2x` (the two generation-expansion cases
are `tyndp/gep` and `sienna/5bus`). Evaluated with `storage_regret`.

## Fleet (as modelled, excluding the per-node ENS slack)

| | Nodes | Gen. | Seas. hydro | Lines | Peak demand | Annual demand |
|---|---|---|---|---|---|---|
| NREL-118 | 118 | 312 | 15 reservoirs | 186 | 17.3 GW (coincident) | 96.0 TWh (LF 63%) |

- **Firm thermal ≈ 17.3 GW** (CC/CT natural gas 14.6 GW, coal/other steam 2.5 GW, plus
  biomass, geothermal, ICE) — essentially equal to the coincident peak.
- **Hydro 18.65 GW nameplate** across 43 units; the 15 without a fixed dispatch profile are
  modelled as **seasonal storage** (monthly inflow, `is_seasonal = true`), the 28 with a
  PLEXOS dispatch profile as non-dispatchable generation.
- **VRE ≈ 4.5 GW** (solar 3.4 GW + wind 1.1 GW), ≈ 11 % of generation capacity.

## What was restored from the primary source (vs. the earlier Sienna conversion)

This dataset previously lived under `sienna/118bus`, converted from the
Sienna/PowerSystems repackaging of NREL-118. That conversion was faithful in topology,
capacities, transmission limits and profile *shapes* (verified: **all 312 generators match
the primary source by name**, and the per-bus peak demands sum to 18.5 GW, consistent with
the 17.3 GW coincident peak), **but it flattened the merit order**: it used a single fuel
price and only the first heat-rate band, which mis-ordered the supply curve (oil/ICE priced
~5× too low, coal too high, gas understated).

The **merit order is restored directly from the primary NREL PLEXOS data**, which is the
only economically material change:

| | Primary NREL source | Earlier Sienna conversion |
|---|---|---|
| Fuel prices | Coal 1.8, Biomass 2.4, **Natural Gas 5.4**, Oil 21 $/MMBTU (`Fuels and emission rates.csv`) | single ≈ 3.5 $/MMBTU |
| Heat rates | 5 rising incremental bands per unit (`Generators.csv`) | band 1 only |

Each generator's single-block LP `variable_cost` is set to its **capacity-weighted-average
incremental heat rate × its fuel price + VO&M charge**. The result is a correctly ordered,
steeper supply curve: geothermal ≈ 3 $/MWh, coal/steam ≈ 21–60, combined-cycle gas ≈ 29–60,
peaking gas turbines up to ≈ 233 $/MWh.

**Unchanged / inherited (faithful to NREL-118):** the 118-node / 186-line topology and line
limits, generator nameplate capacities, the hydro / solar / wind / demand profile shapes and
the seasonal inflow pattern, and the seasonal-storage formulation. Feasibility penalties
follow the suite convention (`ENS variable_cost = borrow_cost = 10 000 $/MWh`,
`spillage_cost = 0`).

## Why it is a non-discriminating control

Even with the corrected merit order, every representative-period method (k-means, k-medoids,
hierarchical, chronological, and the hull selectors) lands within ≈ 0.04 % `storage_regret`
at n_rp = 10. This is a **genuine physical property of NREL-118, not a conversion artifact**:
the seasonal hydro carries only **2.56 TWh/yr ≈ 3 % of the 96 TWh demand** (2–5 % capacity
factor), so there is no strong inter-period storage coupling for period selection to get
wrong; the reduced solutions differ, but the reconstructed full-year dispatch collapses to
the same cost. NREL-118 is therefore kept as the **large-scale, source-faithful, negative
control** — evidence that the method costs nothing when temporal aggregation is easy — while
storage discrimination is carried by `tyndp/p2x` and the synthetic development system, and
investment discrimination by `tyndp/gep` and `sienna/5bus`.

To recover the source-faithful *flattened* costs, revert `technologies-generation.csv`
`variable_cost` to the Sienna values.
