#!/usr/bin/env python3
# Synthetic SENSITIVITY dataset for BlendedClustering (replaces 5bus, which post-injection
# is run-of-river and sensitive to nothing). Design goals:
#   * storage BINDS: large seasonal reservoirs (β ≫ horizon) so the reduced solution's
#     reconstructed seasonal SoC trajectory drives deployed cost — the channel that must
#     be sensitive to the clustering/weight/normalization/tol choices.
#   * storage-RELIANT: cheap hydro + solar carry most of demand; gas is a costly backstop,
#     so mis-timed discharge genuinely costs (gas or reservoir borrow at VOLL).
#   * DIVERSE day-shapes: seasonal × diurnal demand vs seasonal(anti-phase) × diurnal solar,
#     plus staggered reservoir inflows, so representative *selection* and weight *interpolation*
#     change the reconstruction (post-injection the inflow channel is exact, so the signal
#     must come from the dispatch/day-shape channel — hence the rich intra/inter structure).
#   * FAST: one node, ~10 assets, D=120 days × 24 h = 2880 steps.
#
# `mismatch` (0..1) is the difficulty knob: it scales the seasonal supply/demand phase offset
# and amplitude — larger ⇒ more seasonal arbitrage the reduction must get right ⇒ more spread.
import os, math

D, H = 120, 24                      # base periods (days) × intra-period steps
NP = D                              # season index runs over days
K = 3                               # seasonal reservoirs, staggered peaks
KAPPA = 3.0                         # seasonal concentration of hydro inflow

def clamp(x, lo=0.0, hi=None): return max(lo, x if hi is None else min(hi, x))

def demand(day, h, mismatch):
    # winter-peaked season (peak at day 0/D) × evening-peaked diurnal, plus a mild
    # day-to-day ripple so consecutive days are not identical (helps selection bite).
    season = 1.0 + (0.25 + 0.35 * mismatch) * math.cos(2 * math.pi * day / D)
    diurnal = 1.0 + 0.35 * math.cos(2 * math.pi * (h - 19) / H)      # evening peak ~19h
    ripple = 1.0 + 0.08 * math.cos(2 * math.pi * day / 7.0)          # weekly ripple
    return season * diurnal * ripple

def solar(day, h, mismatch):
    # summer-peaked season (anti-phase to demand) × midday diurnal; zero at night.
    season = clamp(1.0 + (0.3 + 0.5 * mismatch) * math.cos(2 * math.pi * day / D - math.pi))
    diurnal = clamp(math.cos(2 * math.pi * (h - 12) / H))           # 0 at night, peak noon
    return season * diurnal

def hydro_season(day, phi):
    # von Mises bump over the year peaked at day-fraction phi (spring melt etc.).
    return math.exp(KAPPA * math.cos(2 * math.pi * (day / D - phi)))

def gen(mismatch, d):
    os.makedirs(d, exist_ok=True)
    W = lambda n, t: open(os.path.join(d, n), "w").write(t)

    # ---- profiles (normalized to peak = 1 over the horizon, i.e. fractions of peak_*) ----
    dem = [[demand(day, h, mismatch) for h in range(H)] for day in range(D)]
    sol = [[solar(day, h, mismatch) for h in range(H)] for day in range(D)]
    phis = [0.30 + 0.20 * (k / max(K - 1, 1)) for k in range(K)]     # staggered spring peaks
    inf = [[[hydro_season(day, phis[k]) for h in range(H)] for day in range(D)] for k in range(K)]
    dpk = max(max(r) for r in dem); spk = max(max(r) for r in sol)
    ipk = [max(max(r) for r in inf[k]) for k in range(K)]
    # peak_inflow_k in MW: total annual inflow energy budget / (annual demand) ~ storage-reliant.
    # size so hydro+solar ≈ annual demand; gas fills the rest (costly).
    ann_dem = sum(sum(r) for r in dem) / dpk                          # in "peak-hours" units
    # hydro delivers ~65% of demand energy, solar ~30%, gas the ~5% residual (peaks/gaps).
    hydro_share, solar_share = 0.65, 0.30
    cols = ["Load_demand", "Solar_avail"] + [f"Hydro{k+1}_inflow" for k in range(K)]
    lines = ["timestep," + ",".join(cols)]
    for day in range(D):
        for h in range(H):
            row = [dem[day][h] / dpk, sol[day][h] / spk] + [inf[k][day][h] / ipk[k] for k in range(K)]
            lines.append(f"{day*H+h+1}," + ",".join(f"{v:.5f}" for v in row))
    W("profiles.csv", "\n".join(lines) + "\n")
    W("scalars.csv", "scalar,value\noperations_weight,1.0\ntimestep_duration,1.0\n")

    PEAK_DEMAND = 100.0
    ann_dem_MWh = PEAK_DEMAND * ann_dem                               # annual demand energy
    # Solar capacity so its annual energy ≈ solar_share·demand.
    ann_sol_frac = sum(sum(r) for r in sol) / spk                     # peak-hours of solar
    SOLAR_CAP = solar_share * ann_dem_MWh / max(ann_sol_frac, 1e-9)
    # Each reservoir's annual inflow so total hydro energy ≈ hydro_share·demand.
    ann_inf_frac = [sum(sum(r) for r in inf[k]) / ipk[k] for k in range(K)]
    E_res = hydro_share * ann_dem_MWh / K                             # MWh/yr per reservoir
    peak_inflow = [E_res / max(ann_inf_frac[k], 1e-9) for k in range(K)]
    # Reservoir energy capacity = ~half the annual inflow ⇒ β ≫ intra-day, seasonal binding.
    RES_CAP = [0.5 * E_res for _ in range(K)]

    assets = "id,location,technology,unit_capacity,initial_units\n"
    techs = "id,type,carrier_out\n"
    tstor = "id,is_seasonal,efficiency_in,efficiency_out\n"
    sstor = "id,initial_storage_level,peak_inflow,spillage_cost,borrow_cost,can_charge\n"
    astor = "id,capacity_storage_energy,initial_storage_units\n"
    aprof = "asset,profile_type,profile\nLoad,demand,Load_demand\nSolar,availability,Solar_avail\n"
    for k in range(K):
        hid = f"Hydro{k+1}"
        # discharge power rating: enough to cover a good share of peak demand when needed.
        assets += f"{hid},1,{hid},40.0,1\n"; techs += f"{hid},storage,electricity\n"
        tstor += f"{hid},true,1.0,1.0\n"
        sstor += f"{hid},{0.5*RES_CAP[k]:.5f},{peak_inflow[k]:.5f},0.0,10000.0,false\n"
        astor += f"{hid},{RES_CAP[k]:.5f},1\n"; aprof += f"{hid},inflows,{hid}_inflow\n"
    assets += f"Solar,1,Solar,{SOLAR_CAP:.5f},1\n"; techs += "Solar,generation,electricity\n"
    assets += "Gas,1,Gas,150.0,1\n"; techs += "Gas,generation,electricity\n"

    W("assets.csv", assets); W("technologies.csv", techs); W("technologies-storage.csv", tstor)
    W("assets-storage-seasonal.csv", sstor); W("assets-storage.csv", astor)
    W("assets-profiles.csv", aprof)
    # Solar is a must-take renewable (variable_cost 0); gas is the costly backstop.
    W("technologies-generation.csv",
      "id,variable_cost,ramping_rate\nSolar,0.0,1.0\nGas,30.0,1.0\n")
    W("technologies-conversion.csv", "id,carrier_in,efficiency_in,efficiency_out\n")
    W("assets-storage-seasonal-reservoir-levels.csv", "asset,profile_type,profile_name\n")
    W("demand.csv", f"id,location,carrier,peak_demand\nLoad,1,electricity,{PEAK_DEMAND}\n")
    W("transmission-lines.csv", "id,from,to,carrier,export_capacity,import_capacity\n")
    W("investments.csv", "id,cost\n")
    W("profiles-reservoir-levels.csv", "timestep\n")   # empty ⇒ model uses hard [0,1]·cap band
    return dict(SOLAR_CAP=SOLAR_CAP, peak_inflow=peak_inflow, RES_CAP=RES_CAP)


def write_sweeps(data_name, sweep_dir):
    # Sensitivity: clustering × weight × normalization × tol at n_rp=10 (dirac has no PGD tol).
    clus = ["k_means", "k_medoids", "hierarchical", "chronological",
            "convex_hull", "convex_hull_with_null", "conical_hull"]
    norms = ["economic", "unscaled", "minmax"]
    tols = {"dirac": ["0.01"], "convex": ["0.01", "0.001", "0.0001", "0.00001"],
            "conical": ["0.01", "0.001", "0.0001", "0.00001"],
            "conical_bounded": ["0.01", "0.001", "0.0001", "0.00001"]}
    hdr = "n_rep_periods,period_length,clustering_type,weight_type,tol,evaluation_type,normalization,fix_every\n"
    sens = [hdr, f"1,{D*H},k_means,dirac,0.01,none,economic,1\n"]
    for cl in clus:
        for w in ["dirac", "convex", "conical", "conical_bounded"]:
            for nrm in norms:
                for t in tols[w]:
                    sens.append(f"10,{H},{cl},{w},{t},storage_regret,{nrm},1\n")
    open(os.path.join(sweep_dir, f"sens_{os.path.basename(data_name)}.csv"), "w").write("".join(sens))

    # Ablation: PROPOSED vs single-component knockouts, over an n_rp grid (D=120 ⇒ up to 40).
    abl = [hdr, f"1,{D*H},k_means,dirac,0.01,none,economic,1\n"]
    for n in [5, 10, 20, 40]:
        abl.append(f"{n},{H},conical_hull,convex,0.01,storage_regret,economic,1\n")     # PROPOSED
        abl.append(f"{n},{H},conical_hull,convex,0.01,storage_regret,unscaled,1\n")     # -economic
        abl.append(f"{n},{H},convex_hull,convex,0.01,storage_regret,economic,1\n")      # -conic-selection
        abl.append(f"{n},{H},conical_hull,dirac,0.01,storage_regret,economic,1\n")      # -convex-weights
    open(os.path.join(sweep_dir, f"{os.path.basename(data_name)}.csv"), "w").write("".join(abl))


if __name__ == "__main__":
    data_name = "synth/hydro"
    ddir = os.path.join("inputs", data_name)
    info = gen(0.6, ddir)                                # mismatch=0.6 (moderate difficulty)
    write_sweeps(data_name, os.path.dirname(ddir))
    print(f"generated {ddir}: solar_cap={info['SOLAR_CAP']:.1f} "
          f"peak_inflow={[round(x,2) for x in info['peak_inflow']]} res_cap={[round(x) for x in info['RES_CAP']]}")
    print(f"wrote sweeps inputs/{os.path.dirname(data_name)}/sens_hydro.csv (sensitivity) and {data_name}.csv (ablation)")
