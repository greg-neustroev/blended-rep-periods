#!/usr/bin/env python3
from ablation import analyze, mean
import statistics as st
rows, arms, ref = analyze("outputs/synth/hydro.csv")

CLUS=["k_means","k_medoids","hierarchical","chronological","convex_hull","convex_hull_with_null","conical_hull"]
WGT=["dirac","convex","conical","conical_bounded"]
NRM=["economic","unscaled","minmax"]

def sel(cl,w,nrm,nrp=10,tol=0.01):
    return [r for r in arms if r["clustering_type"]==cl and r["weight_type"]==w
            and r["normalization"]==nrm and r["_nrp"]==nrp
            and abs((r["_tol"] or -9)-tol)<1e-12]
def m(rs,k):
    v=[r[k] for r in rs if r[k] is not None]; return st.mean(v) if v else None

# Head-to-head at n_rp=10, tol=0.01: regret per combo, all three norms side by side
print("HEAD-TO-HEAD at n_rp=10, tol=0.01 — regret (%) [and weight-fit time (s)]")
print(f"{'clustering / weight':<34}{'economic':>18}{'unscaled':>18}{'minmax':>18}")
wins={n:0 for n in NRM}
for cl in CLUS:
    for w in WGT:
        line=f"{cl+' / '+w:<34}"
        vals={}
        for nrm in NRM:
            rs=sel(cl,w,nrm); r=m(rs,"_regret"); t=m(rs,"_tweight")
            vals[nrm]=r
            line+=f"{(f'{r:6.2f} [{t:5.2f}s]' if r is not None else '   -'):>18}"
        print(line)
        present={k:v for k,v in vals.items() if v is not None}
        if present:
            best=min(present,key=present.get); wins[best]+=1
print(f"\nper-combo winner tally (lowest regret): "+
      ", ".join(f"{n}={wins[n]}" for n in NRM))

# aggregate regret + total time over the full 7x4 grid, per norm (apples-to-apples)
print("\nGRID-WIDE mean over all 28 combos at n_rp=10, tol=0.01:")
print(f"{'norm':<10}{'mean regret%':>14}{'median regret%':>16}{'mean weightfit_s':>18}{'mean total_s':>14}")
for nrm in NRM:
    rr=[m(sel(cl,w,nrm),"_regret") for cl in CLUS for w in WGT]
    rr=[x for x in rr if x is not None]
    tt=[m(sel(cl,w,nrm),"_tweight") for cl in CLUS for w in WGT]; tt=[x for x in tt if x is not None]
    tot=[m(sel(cl,w,nrm),"_ttotal") for cl in CLUS for w in WGT]; tot=[x for x in tot if x is not None]
    print(f"{nrm:<10}{st.mean(rr):>14.2f}{st.median(rr):>16.2f}{st.mean(tt):>18.3f}{st.mean(tot):>14.3f}")

# PGD iterations per norm (why economic is slow) — conical_hull/convex across tols
print("\nPGD total iterations, conical_hull/convex, n_rp=10, by tol × norm:")
print(f"{'tol':<10}{'economic':>16}{'unscaled':>16}{'minmax':>16}")
for tol in sorted({r['_tol'] for r in arms if r['_tol']},reverse=True):
    line=f"{tol:<10.0e}"
    for nrm in NRM:
        rs=sel('conical_hull','convex',nrm,10,tol); v=m(rs,"_pgd_total_iters")
        line+=f"{(f'{v:,.0f}' if v else '-'):>16}"
    print(line)
