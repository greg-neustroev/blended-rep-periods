# Blended representative-period clustering — the proposed method

This note records the *theoretical* basis of the representative-period (RP) clustering
method implemented on the `proposed-method` branch, and the design of the comparison,
ablation, and sensitivity studies that test it. It is the deployable distillate of a long
investigation; only the components that earned their keep are kept (see
[§7 What was pruned](#7-what-was-pruned-and-why)).

The companion working notes with the full empirical record live in
`notes/usage-informed-clustering-findings.md`.

---

## 1. The problem regret grades

An energy-system planning LP over a full year (8760 h) is expensive. RP clustering replaces
the year's `D` base periods (days) by a small set of `R` representatives, solves the reduced
model, and reconstructs a full-horizon solution. We grade the reduced solution by its
**regret**: the true full-horizon cost incurred when the *committed block* of decisions is
fixed to the reduced solution's values and everything else is re-optimised as recourse.

The committed block depends on the problem:

| Problem | Committed block | Regret | Reconstruction |
|---|---|---|---|
| Investment (GEP) | capacity vector `invested_units` (period-independent) | `investment_regret` | copy capacity into the full model, re-optimise all operations |
| Storage / dispatch (5bus) | seasonal inter-period state-of-charge trajectory | `storage_regret` | prolong the reduced R-period trajectory to 365 periods |

Writing the year's value as `F(u) = c_I'u + Σ_d φ_d(u)`, where `φ_d` is period `d`'s
recourse value function (a convex, piecewise-linear function of the committed block `u`), RP
clustering is a **quadrature rule**: approximate `Σ_d φ_d` by `Σ_r w_r φ_r` over the selected
representatives `r` with weights `w_r`. Two facts drive every design choice below:

* **Regret is a *committed-block* quantity, not a representation norm.** Operational
  reconstruction error is recourse — regret re-optimises it away. So the pipeline must serve
  the committed block, not the Euclidean fit of the full operational vector.
* **The regret surface has a smooth part and a cliff.** Where the reduced model is broadly
  correct, `regret ≈ ½ δu' H δu` — a *quadratic form in the constraint metric*, not a plain
  norm of the feature residual. Where the reduced model under-represents scarce periods it
  *under-builds*, and regret jumps by `VOLL · Σ_t [peak·demand − Σ_g availability·capacity]_+`
  — an adequacy **hinge** over the worst days. Representing the *binding* periods matters far
  more than minimising average error.

---

## 2. The pipeline: five operators

```
 profiles ──0──▶ normalise ──1──▶ select reps ──2──▶ fit weights ──3──▶ solve ──4──▶ reconstruct
             (feature scale)     (which days)      (blend W, chain Wᶜʰ)  (Gurobi)   (→ full horizon)
```

Each step exposes several options in code; the **proposed method** picks one per step, and
the **ablation** knocks each back to a baseline one axis at a time. The steps and options are
implemented in `BlendedClustering/src/TemporalClustering/` (`economic_normalization.jl`,
`clustering.jl`, `weight_fitting.jl`) and driven by a per-experiment CSV row
(`ExperimentData`, columns `normalization`, `clustering_type`, `weight_type`,
`chain_weight_type`, `tol`, …).

---

## 3. Step 0 — Normalisation (the feature metric)

**Options:** `unscaled` · `minmax` · **`economic`** (proposed).

Clustering happens in a feature space; the metric of that space decides which days look
"extreme". The default `unscaled` clusters the dimensionless [0,1] profiles as-is (an
isotropic metric over ~2000 feature rows). `minmax` rescales each row to [0,1] — the classic
baseline; it *centres* (subtracts the row minimum), so it is used only for selection, never
for the profiles handed to the model.

**Proposed: `economic`.** Multiply each feature row by its constraint-incidence scale so the
space is in constraint-RHS (energy/cost) units:

* demand `× τ·peak_demand·V̄`,
* availability `× I_g·unit_capacity` (investment) or `× τ·capacity·V̄` (operations),
* inflow `× peak_inflow·V̄`.

**Why.** First-order regret is `⟨S'ν, ε⟩`, where `ε` is the feature residual, `ν` are prices,
and `S` is the *a-priori* incidence map sending a feature residual to a constraint-RHS
residual (peak on demand rows, capacity/`κ_g` on availability rows). The natural selection
metric is `SᵀΛS`. Its structural half `S` is free (no solve); only the prices `Λ` need a
solve. The economic scaling **is the diagonal of `S`** — it puts features in the units in
which they load constraints. This *replaces* the ad-hoc minmax normalisation with a
constraint-informed one, rather than adding a knob.

Diagonal scaling cannot, on its own, form the off-diagonal residual-load direction
(demand − renewables) that the adequacy hinge lives on; a subtraction between two different
rows is not a per-row rescale. That off-diagonal direction is recovered *geometrically* by the
conic hull in Step 1, which is why economic + conic works as a pair.

---

## 4. Step 1 — Selection (which days are representatives)

**Options:** `k_means` · `k_medoids` · `hierarchical` · `chronological` (traditional
baselines) · `convex_hull` · `convex_hull_with_null` · **`conical_hull`** (proposed).

Selection is now specified **independently of the weight class** — the historical `:hull`
alias, whose hull geometry was silently inferred from the weight type, has been removed. You
name `:convex_hull` or `:conical_hull` directly, so any selection geometry can be paired with
any weight class (this is what lets the proposed method use conic selection *and* convex
weights).

* **k-means** — synthetic centroids; the classic RP method. Not real days.
* **k-medoids / hierarchical (Ward)** — cluster on Euclidean distance; each cluster's
  representative is its medoid (a real day). The standard "traditional clustering" baselines.
* **chronological** — partition the year into `n_rp` contiguous, equal-length time blocks
  and take each block's medoid. An order-preserving "sequential linked periods" baseline
  that groups by time rather than feature similarity.
* **convex / conical hull** — greedy: repeatedly add the day furthest (in projection residual)
  from the current hull. The convex hull captures interior-spanning extremes; the **conical
  hull** anchors at the origin and captures extreme *rays* — the days that load a constraint
  hardest.

**Why conical hull.** Adequacy is a hinge: the reduced model must *see* the binding
(cold-calm-dark, high-residual-load) days or it under-builds and pays the VOLL cliff.
Traditional clustering spans the cloud and dilutes those few days; the conical hull, in the
economic (constraint-incidence) metric, sits on the extreme rays where capacity binds and so
grabs them. Empirically, on GEP the raw convex hull leaves 64 % regret at 5 RP while
conic+economic sits at ~13 %, and an explicit bespoke residual-load/price-ladder selection
metric was at best a *wash* and often worse — so the deployable realisation of "cluster in
the constraint metric" is exactly **economic normalisation + conic hull**, not a hand-built
metric. The bespoke ladder is also `κ`-unstable on greenfield investment, another reason it
is not part of the proposed method.

---

## 5. Step 2 — Weight fitting (the quadrature rule)

**Ops options:** `dirac` (nearest-rep, no blending) · **`convex`** (proposed) · `conical` ·
`conical_bounded` (sub-unit conic).
**Chain options (storage only):** `none` · **`convex`** (proposed) · `conical`.

Weights are fit by projected gradient descent (Condat's simplex projection for convex, the
non-negative-orthant projection for conical), one blend per base period, in the **original
physical feature space** — *not* the selection metric.

* **`dirac`** — each period keeps its single nearest representative; no blending. The "no
  weight fitting" baseline.
* **`convex`** — non-negative weights summing to one (partition of unity).
* **`conical`** — non-negative, unbounded above.

**Why convex.** As a quadrature rule for a convex PWL integrand, the weight class sets the
*bias direction* of the cost error via Jensen's inequality:

* **convex** (and any sub-unit non-negative blend, since `φ(0)=0`) over-charges,
  `Σ_r w_r φ(b_r) ≥ φ(b̄)` — it errs **safe-side** (conservative, over-build pressure);
* **dirac** with a missed extreme errs **unsafe** (under-build → the ENS cliff);
* general signed/unbounded blends certify nothing in either direction.

Safe-side is what a planning model wants. Non-negativity is the property that keeps the
reconstruction well-behaved; it is the proven weight-class mechanism.

**Why fit in physical, not selection, space.** The model consumes demand and availability
*separately* (balance vs capacity constraints), so the weights must reconstruct the physical
profiles faithfully. Selecting in the constraint metric but fitting the blend in the physical
metric is a deliberate restriction/prolongation split of *metrics*.

**Storage: a separate chain matrix `Wᶜʰ`.** For a per-period committed block (seasonal SoC),
the prolongation operator is fit *independently* of the operational weights (a
Petrov–Galerkin split: restriction `W` ≠ prolongation `Wᶜʰ`), on the seasonal net-inflow
increments, with the same non-negative class. This owns the seasonal-borrow channel; the
clustering metric does not fix it. On 5bus the chain split is the dominant, robust effect for
storage regret; a single matrix playing both roles is far worse.

---

## 6. Step 3 — Solve · Step 4 — Reconstruct

**Step 3.** Solve the reduced model with Gurobi. No method content here.

**Step 4 — reconstruction, by committed block:**

* **Investment (GEP): no reconstruction matrix.** The committed quantity is
  period-independent: fix `invested_units` in the full 8760-h model and re-optimise all
  operations as free recourse. Regret is the true operating cost at that fixed build minus
  the optimum. What matters upstream is therefore selection *adequacy*.
* **Storage (5bus): the chain matrix `Wᶜʰ`.** Prolong the reduced R-period seasonal SoC
  trajectory to 365 periods via `Wᶜʰ` (fit in Step 2). The seasonal-borrow channel lives
  here, not in the clustering metric.

---

## 7. What was pruned (and why)

"Implement only what the method needs." Removed from the code on this branch:

* **Chain classes `signed`, `l1_ball`, `clipped_affine`, `affine`** — all dead ends. The
  proven chain mechanism is non-negativity (interpolation); the chain class is `convex` /
  `conical`. The three projectors that served only the removed classes (`project_box_sum`,
  `project_subunit_conic`, `project_l1_ball`) are gone.
* **The `:hull` alias** — it coupled hull geometry to the weight class. Selection is now named
  explicitly (`:convex_hull`, `:convex_hull_with_null`, `:conical_hull`), independent of the
  weight class.
* **The bespoke usage-informed residual-load *seeding* and price-*ladder* selection metrics**
  are *not* part of the proposed method: they were a wash against economic + conic hull, and
  unstable/wrong-channel where distinct. Their lasting value is the *explanation* (§3–§4) for
  why economic + conic works, not a separate arm.

Retained as available (non-proposed) options for comparison: the sub-unit-conic weight class
`conical_bounded` and its null-augmented hull selection `convex_hull_with_null`.

---

## 8. The proposed method, concretely

| | 5bus (storage) | GEP (investment) |
|---|---|---|
| 0. normalise | `economic` | `economic` |
| 1. select | `conical_hull` | `conical_hull` |
| 2. weights | `convex` | `convex` |
| 2b. chain | `convex` (`Wᶜʰ`) | — (n/a) |
| 3. solve | Gurobi | Gurobi |
| 4. reconstruct | prolong SoC via `Wᶜʰ` | fix capacity, re-optimise ops |
| grade | `storage_regret` | `investment_regret` |

---

## 9. Experiment design (held-out methodology)

Datasets play distinct roles so the method is never tuned on the systems it is judged on. We
presuppose **neither the selection method nor the weight type** — both are free axes. Regret is
`evaluated_objective_value / reference_optimum − 1` (reference = the `n_rep=1` full-horizon
solve). Storage is graded **day-exact** (`fix_every=1`): the convex chain reconstructs every
seasonal boundary, so day-exact is the honest "reconstruct as best as possible" test — the
chain needs base-period anchors, so every real-period method carries it (chain = convex) while
k-means is the no-chain baseline. RP grids cap at 80 for the ~365-period systems; RTS scales
higher.

* **5bus — development system** (`configs/5bus.toml`):
  - **Sensitivity** — *every* clustering type × *weight type* × normalisation × PGD tolerance
    (`{1e-2,1e-3,1e-4,1e-5}`) at n_rp=10. The selection method, weight type, normalisation and
    tolerance are chosen here, and per-(clustering×weight) hyperparameter robustness is measured
    (we do not presuppose conic+convex wins — a method that is more hyperparameter-robust could
    be preferable). α is fixed at 1/L (Lipschitz-optimal) and the iteration count is not imposed
    — the realised N(ε) is reported; a cache-on/off pair validates the greedy-hull cache.
  - **Ablation** — proposed vs each single-component knockout (−economic, −conic-selection,
    −convex-weights, −chain-split) across n_rp ∈ {5,10,20,40,80}.
* **GEP, P2X, 118-bus — held-out test systems** (`configs/comparison.toml`): the full
  clustering × weight matrix (6 × 4 = 24 combinations) at the economic normalisation and the
  tolerance fixed on 5bus, across the RP grid — the generalisation test on systems the method
  was never tuned on. No sensitivity or ablation here.
* **RTS-GMLC — large-scale / scaling** (`configs/rts.toml`, run last): 5-minute data with
  6-hour periods (1464 base periods). Too large for the full weight matrix, so a reduced
  comparison (traditional dirac baselines vs proposed) + a proposed-only extension into the
  hundreds of RPs ({160,320,500}) for the runtime / weight-fit-bottleneck curve.

Reported metrics (all from recorded columns, or the invested-units dumps): regret with seed
variance, cost breakdown, curtailment & feasibility (borrow, weight-sum > 1), capacity-decision
differences (GEP), and per-stage runtime with N(ε) and cache hit-rate.

---

## 10. Scope and honest caveats

* The "lower representation error ⇒ lower regret" story **holds as a selection criterion**
  (represent the binding periods), **not as a norm** (the isotropic feature norm the plain
  hull minimises is the wrong quadratic).
* The regret win is **gated** by two conditions: (a) a **scarcity cliff** — a steep regret
  surface (real VOLL / tight adequacy); comfortable systems with a flat surface show no win;
  and (b) **system-level** scarcity that a period-aggregate signal can target — under
  transmission-*local* scarcity even a true-dual oracle mis-targets, and a nodal/PTDF-aware
  object would be required.
* The explicit constraint metric **explains** why economic + conic hull works; it does not
  *beat* it. The deployable contribution is the coherent, pruned pipeline: economic
  normalisation, conic-hull selection, convex weights, and — for storage — a separate
  non-negative chain prolongation.
