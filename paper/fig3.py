import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"font.size": 9, "axes.linewidth": 0.6})
fig, axes = plt.subplots(1, 2, figsize=(3.4, 2.3))
thr = {2: (1.369, 1.414, 1.966), 3: (2.159, 2.309, 3.251)}  # U3 (L2), U2^u, U2 (L2)
xr = {2: (1.33, 1.85), 3: (2.12, 2.55)}
for ax, K in zip(axes, (2, 3)):
    r = sorted(json.load(open(f"../trimer/phase_onset_K{K}.json")))
    U = np.array([q[0] for q in r]); n1 = np.array([q[3] for q in r])
    U3, U2u, U2 = thr[K]
    ok = n1 < 0.35
    Uc = np.concatenate([[U3], U[ok]]); nc = np.concatenate([[0.0], n1[ok]])
    ax.fill_between(Uc, 0, nc, color="#e8e8e8", lw=0)
    ax.plot(Uc, nc, "-", color="#1f4e79", lw=1.4)
    ax.plot(U[ok], n1[ok], "o", color="#1f4e79", ms=3)
    for u, c in ((U3, "#1f4e79"), (U2u, "#c9563c"), (U2, "#c9563c")):
        if xr[K][0] < u < xr[K][1]: ax.axvline(u, color=c, lw=0.6, ls=":")
    ax.set_title(f"$K={K}$", fontsize=9)
    ax.set_xlabel("$U/t$"); ax.set_xlim(*xr[K]); ax.set_ylim(0, 0.3)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[0].set_ylabel("density per colour $n_1$")
axes[0].text(0.72, 0.18, "trimer\nliquid", transform=axes[0].transAxes, fontsize=8, color="#1f4e79", ha="center")
axes[0].text(0.32, 0.8, "with\nsuperfluid", transform=axes[0].transAxes, fontsize=8, color="#c9563c", ha="center")
axes[1].text(0.72, 0.18, "trimer\nliquid", transform=axes[1].transAxes, fontsize=8, color="#1f4e79", ha="center")
fig.tight_layout(w_pad=1.0); fig.savefig("fig3.pdf"); fig.savefig("fig3.png", dpi=150)
