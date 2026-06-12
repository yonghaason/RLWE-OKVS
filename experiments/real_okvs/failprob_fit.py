#!/usr/bin/env python3
"""Fit lambda = -log2(Pr[fail]) ~ a*w + g per n (Appendix-F methodology) and extrapolate
to lambda=40. Visualization/fitting only; the matching experiment is in failprob.cpp."""
import csv, math, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
with open(os.path.join(HERE, "failprob.csv")) as f:
    for r in csv.DictReader(f):
        rows.append({k:(float(v) if k!="n" else int(float(v))) for k,v in r.items()})

ns = sorted(set(r["n"] for r in rows))
MIN_FAILS = 15          # need enough failures for a reliable point
# usable points: enough failures AND not saturated (failprob<0.7)
def usable(r): return r["fails"] >= MIN_FAILS and r["failprob"] < 0.7 and r["failprob"] > 0

# per-n independent fits
fits = {}
print(f"{'n':>10} {'slope a':>9} {'intercept g':>12} {'R^2':>6} {'w@2^-40':>9} {'#pts':>5}")
for n in ns:
    pts = [(r["w"], -math.log2(r["failprob"])) for r in rows if r["n"]==n and usable(r)]
    if len(pts) < 2:
        print(f"{n:>10}  (too few measurable points)"); continue
    w = np.array([p[0] for p in pts]); lam = np.array([p[1] for p in pts])
    a, g = np.polyfit(w, lam, 1)
    pred = a*w+g; ss=1-np.sum((lam-pred)**2)/np.sum((lam-lam.mean())**2)
    wstar = (40-g)/a
    fits[n] = (a,g,ss,wstar,w,lam)
    print(f"2^{int(math.log2(n)):<8} {a:9.4f} {g:12.4f} {ss:6.3f} {wstar:9.1f} {len(pts):5d}")

# shared-slope fit (eq. 4: slope a depends only on eps, not n) -> per-n intercepts
allpts = [(r["n"], r["w"], -math.log2(r["failprob"])) for r in rows if usable(r)]
if allpts:
    nidx = {n:i for i,n in enumerate(ns)}
    A = np.zeros((len(allpts), 1+len(ns))); y = np.zeros(len(allpts))
    for i,(n,w,lam) in enumerate(allpts):
        A[i,0]=w; A[i,1+nidx[n]]=1.0; y[i]=lam
    sol,_,_,_ = np.linalg.lstsq(A,y,rcond=None)
    a_sh = sol[0]
    print(f"\nShared-slope fit (a depends only on eps): a={a_sh:.4f}")
    print(f"{'n':>10} {'intercept':>10} {'w@2^-40':>9}")
    shared = {}
    for n in ns:
        g = sol[1+nidx[n]]; wstar=(40-g)/a_sh; shared[n]=(g,wstar)
        print(f"2^{int(math.log2(n)):<8} {g:10.4f} {wstar:9.1f}")

# plot
fig, ax = plt.subplots(figsize=(7,4.6))
colors = plt.cm.viridis(np.linspace(0,0.85,len(ns)))
for n,col in zip(ns,colors):
    pts = [(r["w"], -math.log2(r["failprob"])) for r in rows if r["n"]==n and usable(r)]
    if pts:
        w=[p[0] for p in pts]; lam=[p[1] for p in pts]
        ax.plot(w,lam,"o",color=col,label=f"n=2^{int(math.log2(n))}")
    if n in fits:
        a,g,ss,wstar,_,_=fits[n]
        ww=np.linspace(min(w),wstar,50); ax.plot(ww,a*ww+g,"--",color=col,lw=1)
ax.axhline(40,color="gray",ls=":",lw=1); ax.text(ax.get_xlim()[0],41,r"$\lambda=40$",fontsize=9,color="gray")
ax.set_xlabel("band width $w$"); ax.set_ylabel(r"$\lambda=-\log_2 \Pr[\mathrm{fail}]$")
ax.set_title("Real-field (Gaussian) encoding failure: $\\lambda$ linear in $w$ (Appendix-F style)\n"
             "$\\varepsilon=1$; dashed = linear fit extrapolated to $\\lambda=40$")
ax.legend(); fig.tight_layout()
fig.savefig(os.path.join(HERE,"figures","fig_failprob.png"),dpi=130)
print("\nwrote figures/fig_failprob.png")
