#!/usr/bin/env python3
# Computes beta (reservoir buffer horizon, in days of mean daily inflow) for every seasonal
# storage asset in the real case studies, for the paper's case-study summary table.
#
# beta = capacity_storage_energy / mean_daily_inflow_energy, where mean_daily_inflow_energy is
# read directly off the input profiles (peak_inflow x the asset's inflow profile column, averaged
# over all timesteps and scaled to a day). This mirrors the synthetic generator's own definition
# (scripts/gen_synth_sensitivity.py: RES_CAP = beta_days * (E_res / D)), computed here from real
# input data instead of a knob. GEP has no storage (N/A); RTS is excluded from the paper for now.
#
# Usage: python3 scripts/compute_beta.py   (run from the experiments/ repo root)
# Writes ../clustering/data/beta_values.csv (per-asset) with summary stats printed to stdout.

import csv
import os
import statistics

DATASETS = ["tyndp/p2x", "nrel/118bus", "sienna/5bus"]
OUT = os.path.join(os.path.dirname(__file__), "..", "..", "clustering", "data", "beta_values.csv")


def compute_beta(ds):
    base = os.path.join(os.path.dirname(__file__), "..", "inputs", ds)
    cap = {}
    with open(os.path.join(base, "assets-storage.csv")) as f:
        for r in csv.DictReader(f):
            cap[r["id"]] = float(r["capacity_storage_energy"])
    peak_inflow = {}
    with open(os.path.join(base, "assets-storage-seasonal.csv")) as f:
        for r in csv.DictReader(f):
            peak_inflow[r["id"]] = float(r["peak_inflow"])
    prof_col = {}
    with open(os.path.join(base, "assets-profiles.csv")) as f:
        for r in csv.DictReader(f):
            if r["profile_type"] == "inflows":
                prof_col[r["asset"]] = r["profile"]
    with open(os.path.join(base, "profiles.csv")) as f:
        rdr = csv.DictReader(f)
        cols = [c for c in rdr.fieldnames if c != "timestep"]
        sums = {c: 0.0 for c in cols}
        counts = {c: 0 for c in cols}
        for row in rdr:
            for c in cols:
                v = row.get(c, "")
                if v not in ("", None):
                    sums[c] += float(v)
                    counts[c] += 1
    mean_profile = {c: (sums[c] / counts[c] if counts[c] else 0.0) for c in cols}

    betas = {}
    for asset, col in prof_col.items():
        if asset not in cap or asset not in peak_inflow or col not in mean_profile:
            continue
        mean_daily_inflow = mean_profile[col] * peak_inflow[asset] * 24  # hourly data -> per day
        if mean_daily_inflow <= 0:
            continue
        betas[asset] = cap[asset] / mean_daily_inflow
    return betas


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case_study", "asset_id", "beta_days"])
        for ds in DATASETS:
            case_study = ds.split("/")[-1]
            betas = compute_beta(ds)
            vals = list(betas.values())
            print(f"=== {case_study} ===  n_assets={len(vals)}")
            if vals:
                print(f"  min={min(vals):.1f}  median={statistics.median(vals):.1f}  "
                      f"mean={statistics.mean(vals):.1f}  max={max(vals):.1f}")
            else:
                print("  (no seasonal-storage assets with inflow)")
            for asset, beta in sorted(betas.items(), key=lambda kv: -kv[1]):
                w.writerow([case_study, asset, f"{beta:.2f}"])
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
