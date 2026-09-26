import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"font.size": 9, "axes.linewidth": 0.6})
fig, axes = plt.subplots(1, 2, figsize=(3.4, 2.3))
B, R, G = "#1f4e79", "#c9563c", "#3a8f5c"
# (a) percolation thresholds of emergent hyperedges vs support radius
pc = json.load(open("../emergent/perc_thresholds.json"))
ax = axes[0]
for K, c, ls in ((2, B, "-"), (3, R, "-")):
    r = [0, 1, 2, 3]; p = [pc[f"{K},{x}"] for x in r]
    ax.plot(r, p, "o" + ls, color=c, ms=3.5, lw=1.2, label=f"$K={K}$")
# range of Borromean fractions found in the (U,W) scans, and of support radii
ax.axhspan(0.02, 0.5, color="#e8e8e8", lw=0)
ax.text(1.55, 0.13, "$f_B$ found", fontsize=7, color="#555555")
ax.set_yscale("log"); ax.set_xlabel("support radius $r$"); ax.set_ylabel("percolation threshold $p_c$")
ax.set_xticks([0, 1, 2, 3]); ax.set_ylim(1e-4, 1)
ax.legend(frameon=False, fontsize=7, loc="lower left")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.text(0.01, 0.95, "(a)", fontsize=9, fontweight="bold")
# (b) critical disorder of the composites vs U
loc = json.load(open("../emergent/localization.json"))
ax = axes[1]
for K, c, ls in ((2, B, "-"), (3, R, "-")):
    rows = loc[str(K)]; U = np.array([r["U"] for r in rows]); Wc3 = np.array([r["Wc3"] for r in rows]); Wc2 = np.array([r["Wc2"] for r in rows])
    U2u = 2 * (K - 1) / np.sqrt(K)
    ax.plot(U, Wc3, "-", color=c, lw=1.4)
    m = U > U2u * 1.02
    ax.plot(U[m], Wc2[m], "--", color=c, lw=1.0)
    ax.axhline({2: 17.39, 3: 32.26}[K], color=c, lw=0.6, ls=":")
ax.set_yscale("log"); ax.set_xlabel("$U/t$"); ax.set_ylabel("critical disorder $W_c/t$")
ax.set_ylim(0.3, 50)
ax.text(2.6, 21, "one particle", fontsize=6.5, color=B); ax.text(2.6, 38, "one particle", fontsize=6.5, color=R)
ax.text(3.2, 6.5, "pair", fontsize=7, color="#555555"); ax.text(3.6, 1.35, "trimer", fontsize=7, color="#555555")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.text(0.52, 0.95, "(b)", fontsize=9, fontweight="bold")
fig.tight_layout(w_pad=1.2); fig.savefig("fig4.pdf"); fig.savefig("fig4.png", dpi=150)
