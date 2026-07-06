#!/usr/bin/env python3
# Synthetic SENSITIVITY dataset(s) for BlendedClustering — the development system that replaces
# 5bus (post exact-inflow injection 5bus is run-of-river and sensitive to nothing).
#
# The point of a *synthetic* dev system is orthogonal control: each generator knob varies ONE
# property of the data so we can map which properties the method is sensitive to and which it is
# robust to. Run from the repo root: `python3 scripts/gen_synth_sensitivity.py`.
#
# Fixed skeleton: one node, 3 seasonal hydro reservoirs (can_charge=false) + solar + a gas
# backstop, D=120 days × H=24 h = 2880 steps (small ⇒ the full sweep runs fast). Cheap hydro+
# solar carry the demand; gas is the costly backstop; reservoir borrow is VOLL.
#
# KNOBS (baseline in gen()'s defaults sits in the binding, storage-reliant, sensitive corner):
#   beta_days      buffer horizon: reservoir energy = beta_days · mean-daily-inflow. GATES the
#                  inter-period reconstruction — β≈1–3 d is run-of-river (storage carries nothing
#                  ⇒ every method ties), β≈months is seasonal (reconstruction drives regret).
#   gas_cost       €/MWh of the backstop. GATES storage-*reliance*: cheap gas ⇒ seasonal storage
#                  is skipped (spill in spring, burn gas in winter) ⇒ reconstruction irrelevant;
#                  dear gas ⇒ storage essential ⇒ reconstruction errors cost.
#   renew_share    annual share of demand energy from hydro+solar (rest from gas). Second facet
#                  of reliance (supply side).
#   hull_spread    per-day magnitude dispersion of the residual-load cloud. 0 ⇒ in-hull continuum
#                  (weight class irrelevant); large ⇒ out-of-(convex-)hull days ⇒ conical / sub-
#                  unit-conic weights help, convex/dirac struggle. (Anchored on residual load, not
#                  inflow: exact-inflow injection cancels the inflow-hull channel.)
#   n_regimes      0 ⇒ smooth seasonal continuum; k>0 ⇒ k discrete day archetypes. Sets intrinsic
#                  dimensionality ⇒ how much SELECTION and n_rp matter.
#   scale_ratio    raw magnitude of the inflow feature block (its profile peak = scale_ratio, while
#                  demand/solar are peak 1); peak_inflow rescales the energy back so only the
#                  cross-block clustering-feature scale changes ⇒ `unscaled` is dominated by inflow
#                  while `economic`/`minmax` are not ⇒ NORMALIZATION matters; 1.0 ⇒ homogeneous.
#   noise          per-(day,h) multiplicative irregularity on demand/solar. Empirically a SECOND
#                  out-of-hull stressor (erratic days push the residual-load cloud out of hull ⇒
#                  weight/reconstruction sensitivity), NOT a clean robustness dial — the robustness
#                  side of the story is read from the axes that stay flat + the per-seed variance.
#   season_diurnal seasonal:diurnal amplitude ratio. Shifts stress between the inter-period storage
#                  chain (seasonal) and intra-period dispatch (diurnal).
import os, math, random

D, H, K = 120, 24, 3
PEAK_DEMAND = 100.0


def _clamp(x, lo=0.0, hi=None):
    return max(lo, x if hi is None else min(hi, x))


def _build_profiles(beta_days, gas_cost, renew_share, hull_spread, n_regimes,
                    scale_ratio, noise, season_diurnal, seed):
    rng = random.Random(seed)
    a_diur = 0.35
    a_seas = 0.40 * season_diurnal

    # day-phase: continuum, or snapped to n_regimes discrete archetypes.
    def dayphase(day):
        f = day / D
        if n_regimes and n_regimes > 0:
            f = (math.floor(f * n_regimes) + 0.5) / n_regimes
        return f

    # per-day residual-load magnitude multipliers: lognormal-ish dispersion ⇒ out-of-hull cloud.
    day_mag = [math.exp(hull_spread * rng.gauss(0.0, 1.0)) for _ in range(D)]
    mean_mag = sum(day_mag) / D
    day_mag = [m / mean_mag for m in day_mag]                     # normalize mean to 1

    def demand(day, h):
        ph = dayphase(day)
        season = 1.0 + a_seas * math.cos(2 * math.pi * ph)         # winter (ph=0) peak
        diurnal = 1.0 + a_diur * math.cos(2 * math.pi * (h - 19) / H)
        eps = 1.0 + noise * rng.gauss(0.0, 1.0)
        return _clamp(day_mag[day] * season * diurnal * eps, 0.0)

    def solar(day, h):
        ph = dayphase(day)
        season = _clamp(1.0 + (0.3 + 0.5 * season_diurnal) * math.cos(2 * math.pi * ph - math.pi))
        diurnal = _clamp(math.cos(2 * math.pi * (h - 12) / H))
        eps = 1.0 + noise * rng.gauss(0.0, 1.0)
        return _clamp(season * diurnal * eps, 0.0)

    def hydro(day, k):
        phi = 0.30 + 0.20 * (k / max(K - 1, 1))                    # staggered spring peaks
        return math.exp(3.0 * math.cos(2 * math.pi * (dayphase(day) - phi)))

    dem = [[demand(day, h) for h in range(H)] for day in range(D)]
    sol = [[solar(day, h) for h in range(H)] for day in range(D)]
    inf = [[hydro(day, k) for day in range(D)] for k in range(K)]  # flat within day
    return dem, sol, inf


def gen(outdir, *, beta_days=60.0, gas_cost=30.0, renew_share=0.95, hull_spread=0.0,
        n_regimes=0, scale_ratio=1.0, noise=0.0, season_diurnal=1.0, seed=0):
    os.makedirs(outdir, exist_ok=True)
    W = lambda n, t: open(os.path.join(outdir, n), "w").write(t)
    dem, sol, inf = _build_profiles(beta_days, gas_cost, renew_share, hull_spread,
                                    n_regimes, scale_ratio, noise, season_diurnal, seed)

    dpk = max(max(r) for r in dem)
    spk = max(max(r) for r in sol) or 1.0
    ipk = [max(inf[k]) for k in range(K)]
    shape = [[inf[k][day] / ipk[k] for day in range(D)] for k in range(K)]   # inflow shape, peak 1
    # scale_ratio sets the inflow feature's raw magnitude (peak = scale_ratio) as the CLUSTERING
    # sees it: demand/solar features are peak 1, inflow is peak scale_ratio. peak_inflow below
    # rescales the energy back, so scale_ratio changes ONLY the cross-block feature scale — which
    # `unscaled` clustering is dominated by while `economic`/`minmax` are not — and never the
    # physics. scale_ratio=1 ⇒ homogeneous (normalization moot); large ⇒ normalization bites.
    cols = ["Load_demand", "Solar_avail"] + [f"Hydro{k+1}_inflow" for k in range(K)]
    lines = ["timestep," + ",".join(cols)]
    for day in range(D):
        for h in range(H):
            row = [dem[day][h] / dpk, sol[day][h] / spk] + [scale_ratio * shape[k][day] for k in range(K)]
            lines.append(f"{day*H+h+1}," + ",".join(f"{v:.5f}" for v in row))
    W("profiles.csv", "\n".join(lines) + "\n")
    W("scalars.csv", "scalar,value\noperations_weight,1.0\ntimestep_duration,1.0\n")

    # ---- sizing (energy balance) ----
    ann_dem = PEAK_DEMAND * sum(sum(r) for r in dem) / dpk          # annual demand energy (MWh)
    hydro_share, solar_share = renew_share * 0.68, renew_share * 0.32
    ann_sol = sum(sum(r) for r in sol) / spk
    SOLAR_CAP = solar_share * ann_dem / max(ann_sol, 1e-9)
    E_res = hydro_share * ann_dem / K                              # annual inflow per reservoir
    # peak_inflow_k rescales the inflow *shape* (× scale_ratio in the profile) back to energy E_res,
    # so the physical inflow is E_res independent of scale_ratio.
    ann_shape = [sum(shape[k]) * H for k in range(K)]              # ×H: flat within day
    peak_inflow = [E_res / (scale_ratio * max(a, 1e-9)) for a in ann_shape]
    RES_CAP = [beta_days * (E_res / D) for _ in range(K)]           # β days of mean daily inflow

    assets = "id,location,technology,unit_capacity,initial_units\n"
    techs = "id,type,carrier_out\n"
    tstor = "id,is_seasonal,efficiency_in,efficiency_out\n"
    sstor = "id,initial_storage_level,peak_inflow,spillage_cost,borrow_cost,can_charge\n"
    astor = "id,capacity_storage_energy,initial_storage_units\n"
    aprof = "asset,profile_type,profile\nLoad,demand,Load_demand\nSolar,availability,Solar_avail\n"
    for k in range(K):
        hid = f"Hydro{k+1}"
        assets += f"{hid},1,{hid},60.0,1\n"; techs += f"{hid},storage,electricity\n"
        tstor += f"{hid},true,1.0,1.0\n"
        sstor += f"{hid},{0.5*RES_CAP[k]:.5f},{peak_inflow[k]:.5f},0.0,10000.0,false\n"
        astor += f"{hid},{RES_CAP[k]:.5f},1\n"; aprof += f"{hid},inflows,{hid}_inflow\n"
    assets += f"Solar,1,Solar,{SOLAR_CAP:.5f},1\n"; techs += "Solar,generation,electricity\n"
    assets += f"Gas,1,Gas,{1.3*PEAK_DEMAND:.5f},1\n"; techs += "Gas,generation,electricity\n"

    W("assets.csv", assets); W("technologies.csv", techs); W("technologies-storage.csv", tstor)
    W("assets-storage-seasonal.csv", sstor); W("assets-storage.csv", astor)
    W("assets-profiles.csv", aprof)
    W("technologies-generation.csv",
      f"id,variable_cost,ramping_rate\nSolar,0.0,1.0\nGas,{gas_cost},1.0\n")
    W("technologies-conversion.csv", "id,carrier_in,efficiency_in,efficiency_out\n")
    W("assets-storage-seasonal-reservoir-levels.csv", "asset,profile_type,profile_name\n")
    W("demand.csv", f"id,location,carrier,peak_demand\nLoad,1,electricity,{PEAK_DEMAND}\n")
    W("transmission-lines.csv", "id,from,to,carrier,export_capacity,import_capacity\n")
    W("investments.csv", "id,cost\n")
    W("profiles-reservoir-levels.csv", "timestep\n")              # empty ⇒ hard [0,1]·cap band
    return dict(SOLAR_CAP=SOLAR_CAP, peak_inflow=peak_inflow, RES_CAP=RES_CAP, beta_days=beta_days)


# ---- sweep writers ------------------------------------------------------------------------
_HDR = "n_rep_periods,period_length,clustering_type,weight_type,tol,evaluation_type,normalization\n"

def write_full_sweeps(base):
    # Sensitivity: clustering × weight × normalization × tol at n_rp=10 (dirac has no PGD tol).
    clus = ["k_means", "k_medoids", "hierarchical", "chronological",
            "convex_hull", "convex_hull_with_null", "conical_hull"]
    norms = ["economic", "unscaled", "minmax"]
    tols = {"dirac": ["0.01"], "convex": ["0.01", "0.001", "0.0001", "0.00001"],
            "conical": ["0.01", "0.001", "0.0001", "0.00001"],
            "conical_bounded": ["0.01", "0.001", "0.0001", "0.00001"]}
    sens = [_HDR, f"1,{D*H},k_means,dirac,0.01,none,economic\n"]
    for cl in clus:
        for w in ["dirac", "convex", "conical", "conical_bounded"]:
            for nrm in norms:
                for t in tols[w]:
                    sens.append(f"10,{H},{cl},{w},{t},storage_regret,{nrm}\n")
    open(f"inputs/synth/sens_{base}.csv", "w").write("".join(sens))
    # Ablation: PROPOSED vs single-component knockouts over the n_rp grid {10,20,40,80}.
    abl = [_HDR, f"1,{D*H},k_means,dirac,0.01,none,economic\n"]
    for n in [10, 20, 40, 80]:
        abl.append(f"{n},{H},conical_hull,convex,0.01,storage_regret,economic\n")   # PROPOSED
        abl.append(f"{n},{H},conical_hull,convex,0.01,storage_regret,unscaled\n")   # -economic
        abl.append(f"{n},{H},convex_hull,convex,0.01,storage_regret,economic\n")    # -conic-selection
        abl.append(f"{n},{H},conical_hull,dirac,0.01,storage_regret,economic\n")    # -convex-weights
    open(f"inputs/synth/{base}.csv", "w").write("".join(abl))


def write_panel_sweep(path):
    # Diagnostic panel (n_rp=10): arms chosen so PER-AXIS regret spread is computable for each
    # dataset — the weight axis (conical_hull × 4 weights @ economic), the selection axis (6
    # selections × convex @ economic), and the normalization axis (conical_hull/convex × 3 norms).
    # Running it on each knob variant tells which METHOD axis a given DATA property moves.
    rows = [_HDR, f"1,{D*H},k_means,dirac,0.01,none,economic\n"]
    add = lambda cl, w, nrm: rows.append(f"10,{H},{cl},{w},0.01,storage_regret,{nrm}\n")
    for w in ["dirac", "convex", "conical", "conical_bounded"]:      # weight axis
        add("conical_hull", w, "economic")
    for cl in ["k_means", "k_medoids", "hierarchical", "chronological", "convex_hull"]:  # selection axis
        add(cl, "convex", "economic")                               # (conical_hull/convex already above)
    for nrm in ["unscaled", "minmax"]:                              # normalization axis
        add("conical_hull", "convex", nrm)                          # (economic already above)
    open(path, "w").write("".join(rows))


# Data-property study grid: one entry per (variant name, knobs), each varying ONE knob off the
# baseline (β=60, renew=0.95, gas=30, everything else default). Grouped by the method axis the
# knob is expected to move — β/reliance GATE all sensitivity; the rest each target one axis.
STUDY = [
    ("base",     {}),
    # β — buffer horizon: gates the inter-period storage reconstruction.
    ("beta2", dict(beta_days=2)), ("beta5", dict(beta_days=5)), ("beta15", dict(beta_days=15)),
    ("beta30", dict(beta_days=30)), ("beta120", dict(beta_days=120)),
    # reliance — storage throughput (gates whether reconstruction error costs).
    ("renew40", dict(renew_share=0.40)), ("renew60", dict(renew_share=0.60)),
    ("renew80", dict(renew_share=0.80)), ("renew100", dict(renew_share=1.00)),
    # hull geometry of the residual-load cloud → WEIGHT class (convex vs conical vs bounded).
    ("hull03", dict(hull_spread=0.3)), ("hull06", dict(hull_spread=0.6)), ("hull10", dict(hull_spread=1.0)),
    # intrinsic dimensionality (# day archetypes) → SELECTION and n_rp.
    ("reg3", dict(n_regimes=3)), ("reg6", dict(n_regimes=6)), ("reg12", dict(n_regimes=12)), ("reg30", dict(n_regimes=30)),
    # cross-block feature scale → NORMALIZATION (economic vs unscaled vs minmax).
    ("scale05", dict(scale_ratio=0.5)), ("scale2", dict(scale_ratio=2.0)), ("scale4", dict(scale_ratio=4.0)),
    # noise → day-shape irregularity: empirically a second out-of-hull / reconstruction stressor.
    ("noise10", dict(noise=0.1)), ("noise30", dict(noise=0.3)), ("noise60", dict(noise=0.6)),
    # seasonal:diurnal amplitude → inter- vs intra-period stress.
    ("sd03", dict(season_diurnal=0.3)), ("sd30", dict(season_diurnal=3.0)),
]


def write_study():
    # Regenerate every study variant into the ignorable inputs/synth/knobs/ subdir, the shared
    # diagnostic panel, and the consolidated config that runs the panel on all of them.
    os.makedirs("inputs/synth/knobs", exist_ok=True)
    for name, knobs in STUDY:
        gen(f"inputs/synth/knobs/{name}", **knobs)
    write_panel_sweep("inputs/synth/knobs/panel.csv")
    hdr = (
        "# Data-property sensitivity study. The diagnostic panel (per-axis method arms at\n"
        "# n_rp=10) is run on synth/hydro variants that each vary ONE generator knob, so the\n"
        "# per-axis regret spread of each variant says which METHOD axis that DATA property\n"
        "# moves (β & reliance GATE everything; hull→weight, dimensionality→selection,\n"
        "# scale→normalization, noise→robustness, season_diurnal→inter/intra). Regenerate the\n"
        "# variant datasets first:  python3 scripts/gen_synth_sensitivity.py study\n"
        "# (they live under the gitignored inputs/synth/knobs/, so they are not committed.)\n"
        "inputs = [\n"
    )
    body = "".join(f'  {{ sweep = "synth/knobs/panel", data = "synth/knobs/{name}" }},\n'
                   for name, _ in STUDY)
    tail = (']\n'
            'master_seed = 382507\n'
            'n_seeds     = 5          # deterministic methods repeat; k-means/medoids get variance bars\n'
            'solver      = "gurobi"\n'
            'inputs_dir  = "../inputs"\n'
            'outputs_dir = "../outputs"\n')
    open("configs/synth_knobs.toml", "w").write(hdr + body + tail)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "study":
        write_study()
        print(f"study: generated {len(STUDY)} variants under inputs/synth/knobs/ + panel + configs/synth_knobs.toml")
    else:
        info = gen("inputs/synth/hydro")                          # baseline dev dataset
        write_full_sweeps("hydro")
        print(f"baseline synth/hydro: solar_cap={info['SOLAR_CAP']:.1f} "
              f"peak_inflow={[round(x,2) for x in info['peak_inflow']]} res_cap={[round(x) for x in info['RES_CAP']]}")
        print("wrote inputs/synth/sens_hydro.csv (sensitivity) + inputs/synth/hydro.csv (ablation)")
