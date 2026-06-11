#!/usr/bin/env python3
"""Visualization-only: reads results.csv (produced by okvs_real C++) and renders figures.
No OKVS math here; the encoding/solving is all done in okvs_real.cpp."""
import csv, math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

rows = []
with open(os.path.join(HERE, "results.csv")) as f:
    for r in csv.DictReader(f):
        for k in r:
            if k not in ("sweep", "band"):
                r[k] = float(r[k])
        rows.append(r)

def sel(sweep, band=None, ok_min=0.5):
    out = [r for r in rows if r["sweep"] == sweep and r["ok"] >= ok_min]
    if band: out = [r for r in out if r["band"] == band]
    return out

def log10safe(x):
    return math.log10(x) if (x and math.isfinite(x) and x > 0) else 320.0  # overflow -> ~1e320

# ---- Fig 1: vs band width w (gauss) ----
d = sorted(sel("w", "gauss"), key=lambda r: r["w"])
w = [r["w"] for r in d]
fig, ax1 = plt.subplots(figsize=(6.4, 4.2))
ax1.plot(w, [r["ratio_min"] for r in d], "o-", color="tab:blue", label=r"$\|p^*\|/\|v\|$ (min-norm)")
ax1.set_xlabel("band width $w$"); ax1.set_ylabel(r"$\|p^*\|/\|v\|$", color="tab:blue")
ax1.tick_params(axis="y", labelcolor="tab:blue"); ax1.set_ylim(0, 1.0)
ax2 = ax1.twinx()
ax2.plot(w, [r["sig_min"] for r in d], "s--", color="tab:red", label=r"$\sigma_{\min}(M)$")
ax2.set_ylabel(r"$\sigma_{\min}(M)$", color="tab:red"); ax2.tick_params(axis="y", labelcolor="tab:red")
ax1.set_title("Min-norm conditioning vs band width $w$\n(Gauss band, $n=1024$, $\\varepsilon=1$)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_w.png"), dpi=130); plt.close(fig)

# ---- Fig 2: vs epsilon (gauss) ----
d = sorted(sel("eps", "gauss"), key=lambda r: r["eps"])
eps = [r["eps"] for r in d]
fig, ax1 = plt.subplots(figsize=(6.4, 4.2))
ax1.plot(eps, [r["ratio_min"] for r in d], "o-", color="tab:blue", label=r"$\|p^*\|/\|v\|$")
ax1.plot(eps, [r["g_min_rms"] for r in d], "^-", color="tab:green", label="garbage rms (min-norm)")
ax1.set_xlabel(r"expansion $\varepsilon$  ($m=(1+\varepsilon)n$)")
ax1.set_ylabel(r"$\|p^*\|/\|v\|$  and  garbage rms")
ax2 = ax1.twinx()
ax2.plot(eps, [r["sig_min"] for r in d], "s--", color="tab:red", label=r"$\sigma_{\min}(M)$")
ax2.set_ylabel(r"$\sigma_{\min}(M)$", color="tab:red"); ax2.tick_params(axis="y", labelcolor="tab:red")
ax1.legend(loc="upper right"); ax1.set_yscale("log")
ax1.set_title("Redundancy $\\varepsilon$ is the strong conditioning knob\n(Gauss band, $n=1024$, $w=24$)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_eps.png"), dpi=130); plt.close(fig)

# ---- Fig 3: vs n (scale invariance) ----
d = sorted(sel("n", "gauss"), key=lambda r: r["n"])
n = [r["n"] for r in d]
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(n, [r["ratio_min"] for r in d], "o-", label=r"$\|p^*\|/\|v\|$")
ax.plot(n, [r["sig_min"] for r in d], "s-", label=r"$\sigma_{\min}(M)$")
ax.plot(n, [r["cond"]/100 for r in d], "^-", label=r"cond$(M)/100$")
ax.set_xscale("log", base=2); ax.set_xlabel("$n$ (#keys)"); ax.set_ylim(0, 1.0)
ax.set_title("Min-norm conditioning is scale-invariant in $n$\n(Gauss band, $\\varepsilon=1$, $w=24$)")
ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_n.png"), dpi=130); plt.close(fig)

# ---- Fig 4: gauss vs binary bands ----
db = sorted(sel("band", "binary"), key=lambda r: r["w"])
dg = sorted(sel("band", "gauss"), key=lambda r: r["w"])
fig, (axA, axB) = plt.subplots(1, 2, figsize=(10, 4.2))
axA.plot([r["w"] for r in db], [r["ratio_min"] for r in db], "o-", label="binary")
axA.plot([r["w"] for r in dg], [r["ratio_min"] for r in dg], "s-", label="gauss")
axA.set_xlabel("band width $w$"); axA.set_ylabel(r"$\|p^*\|/\|v\|$"); axA.legend()
axA.set_title("min-norm size")
axB.plot([r["w"] for r in db], [r["sig_min"] for r in db], "o-", label="binary")
axB.plot([r["w"] for r in dg], [r["sig_min"] for r in dg], "s-", label="gauss")
axB.set_xlabel("band width $w$"); axB.set_ylabel(r"$\sigma_{\min}(M)$"); axB.legend()
axB.set_title(r"$\sigma_{\min}$")
fig.suptitle("Gauss bands condition better than binary (same $w$); binary is the addition-only choice\n($n=1024$, $\\varepsilon=1$)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_band.png"), dpi=130); plt.close(fig)

# ---- Fig 5: the chasm -- naive triangular blow-up vs min-norm (n sweep) ----
d = sorted(sel("n", "gauss"), key=lambda r: r["n"])
n = [int(r["n"]) for r in d]
xs = list(range(len(n)))
fig, ax1 = plt.subplots(figsize=(7.2, 4.4))
bars = ax1.bar([x-0.0 for x in xs], [log10safe(r["ratio_naive"]) for r in d],
               width=0.55, color="tab:red", alpha=0.75, label=r"naive triangular: $\log_{10}(\|p\|/\|v\|)$")
ax1.set_ylabel(r"$\log_{10}\,(\|p\|/\|v\|)$  -- naive (red bars)")
ax1.set_xticks(xs); ax1.set_xticklabels([str(x) for x in n]); ax1.set_xlabel("$n$ (#keys)")
for x, r in zip(xs, d):
    val = r["ratio_naive"]
    txt = "overflow" if not math.isfinite(val) else f"$10^{{{int(round(log10safe(val)))}}}$"
    ax1.text(x, min(log10safe(val), 320)+3, txt, ha="center", fontsize=8, color="darkred")
ax2 = ax1.twinx()
ax2.plot(xs, [r["ratio_min"] for r in d], "o-", color="tab:blue", lw=2,
         label=r"min-norm: $\|p^*\|/\|v\|$")
ax2.set_ylabel(r"$\|p^*\|/\|v\|$ -- min-norm (blue)", color="tab:blue")
ax2.tick_params(axis="y", labelcolor="tab:blue"); ax2.set_ylim(0, 1.0)
ax1.set_title("The chasm: naive finite-field triangulation ported to $\\mathbb{R}$ explodes;\n"
              "banded min-norm stays $\\|p^*\\| < \\|v\\|$  (Gauss, $\\varepsilon=1$, $w=24$)")
h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2, loc="upper left", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_chasm.png"), dpi=130); plt.close(fig)

print("wrote figures:", sorted(os.listdir(FIG)))
