#!/usr/bin/env python3
from ablation import analyze, mean, std, grp
import statistics as st

rows, arms, ref = analyze("outputs/synth/hydro.csv")
print(f"synth/hydro  reference optimum (n_rp=1) = {ref:.6g}   arms={len(arms)}\n")

def cell(rs, key="_regret"):
    v=[r[key] for r in rs if r[key] is not None]
    if not v: return "   -   "
    m=st.mean(v); s=st.pstdev(v) if len(v)>1 else 0
    return f"{m:6.2f}±{s:4.1f}" if s>0.05 else f"{m:6.2f}    "

CLUS=["k_means","k_medoids","hierarchical","chronological","convex_hull","convex_hull_with_null","conical_hull"]
WGT=["dirac","convex","conical","conical_bounded"]
NRM=["economic","unscaled","minmax"]
NRP=[10,20,40,80]
TOLS=sorted({r["_tol"] for r in arms if r["_tol"]})

def sel(cl=None,w=None,nrm=None,nrp=None,tol=None):
    out=[]
    for r in arms:
        if cl and r["clustering_type"]!=cl: continue
        if w and r["weight_type"]!=w: continue
        if nrm and r["normalization"]!=nrm: continue
        if nrp and r["_nrp"]!=nrp: continue
        if tol is not None and abs((r["_tol"] or -9)-tol)>1e-12: continue
        out.append(r)
    return out

# ============ 1. NORMALIZATION: regret by normalization, at tol=0.01, per n_rp ============
# aggregate over ALL clustering x weight combos (the whole method family) and also for the
# canonical proposed combo conical_hull/convex.
print("="*78)
print("NORMALIZATION — regret (%) vs n_rp, tol=0.01")
print("="*78)
print("\n(a) pooled over ALL clustering×weight combos (family-wide effect)")
print(f"  {'norm':<10}"+"".join(f"{n:>14}" for n in NRP))
for nrm in NRM:
    print(f"  {nrm:<10}"+"".join(f"{cell(sel(nrm=nrm,nrp=n,tol=0.01)):>14}" for n in NRP))

print("\n(b) proposed combo  conical_hull / convex")
print(f"  {'norm':<10}"+"".join(f"{n:>14}" for n in NRP))
for nrm in NRM:
    print(f"  {nrm:<10}"+"".join(f"{cell(sel('conical_hull','convex',nrm,n,0.01)):>14}" for n in NRP))

print("\n(c) best clustering×weight PER normalization at n_rp=10 (does the winner change?)")
for nrm in NRM:
    best=None
    for cl in CLUS:
        for w in WGT:
            v=[r["_regret"] for r in sel(cl,w,nrm,10,0.01) if r["_regret"] is not None]
            if not v: continue
            m=st.mean(v)
            if best is None or m<best[0]: best=(m,cl,w)
    if best: print(f"  {nrm:<10} best = {best[1]}/{best[2]:<16} regret={best[0]:6.2f}%")

# ============ 2. NORMALIZATION time cost ============
print("\n"+"="*78)
print("NORMALIZATION — time (s), proposed combo conical_hull/convex, tol=0.01")
print("="*78)
print(f"  {'norm':<10}{'n_rp':>6}{'cluster':>10}{'weightfit':>11}{'solve':>10}{'total':>10}")
for nrm in NRM:
    for n in NRP:
        rs=sel('conical_hull','convex',nrm,n,0.01)
        print(f"  {nrm:<10}{n:>6}{mean([r['_tcluster'] for r in rs]) or 0:>10.3f}"
              f"{mean([r['_tweight'] for r in rs]) or 0:>11.3f}"
              f"{mean([r['_tsolve'] for r in rs]) or 0:>10.3f}"
              f"{mean([r['_ttotal'] for r in rs]) or 0:>10.3f}")

# ============ 3. CLUSTERING: hull vs medoids (best weight per, economic, tol .01) ============
print("\n"+"="*78)
print("CLUSTERING TYPE — regret (%), economic, tol=0.01, weight=convex")
print("="*78)
print(f"  {'clustering':<24}"+"".join(f"{n:>14}" for n in NRP))
for cl in CLUS:
    print(f"  {cl:<24}"+"".join(f"{cell(sel(cl,'convex','economic',n,0.01)):>14}" for n in NRP))

# ============ 4. WEIGHT TYPE (vs dirac), conical_hull, economic ============
print("\n"+"="*78)
print("WEIGHT TYPE — regret (%), conical_hull, economic, tol=0.01")
print("="*78)
print(f"  {'weight':<18}"+"".join(f"{n:>14}" for n in NRP))
for w in WGT:
    print(f"  {w:<18}"+"".join(f"{cell(sel('conical_hull',w,'economic',n,0.01)):>14}" for n in NRP))

# ============ 5. TOL sensitivity (conical_hull/convex/economic) ============
print("\n"+"="*78)
print("PGD TOLERANCE — regret(%) & weight-fit time(s), conical_hull/convex/economic, n_rp=10")
print("="*78)
print(f"  {'tol':<10}{'regret':>12}{'weightfit_s':>14}{'pgd_iters':>12}")
for t in sorted(TOLS,reverse=True):
    rs=sel('conical_hull','convex','economic',10,t)
    print(f"  {t:<10.0e}{cell(rs):>12}{mean([r['_tweight'] for r in rs]) or 0:>14.3f}"
          f"{mean([r['_pgd_total_iters'] for r in rs]) or 0:>12.0f}")

# ============ 6. CACHING (data availability) ============
print("\n"+"="*78)
print("CACHING — cache_hits / cache_misses observed (proposed, economic, tol=0.01)")
print("="*78)
for n in NRP:
    rs=sel('conical_hull','convex','economic',n,0.01)
    ch=mean([r['_cache_hits'] for r in rs]); cm=mean([r['_cache_misses'] for r in rs])
    print(f"  n_rp={n:<4} hits={ch} misses={cm}")
