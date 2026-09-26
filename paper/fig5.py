"""Phase diagram of the trimer liquid in the (U, W) plane: glass line W_glass(U) (every state of the
liquid, n < n_1(U), localised) and W_c^(3)(U) (every trimer state localised), K = 2 and 3."""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"font.size": 9, "axes.linewidth": 0.6})
fig, axes = plt.subplots(1, 2, figsize=(3.4, 2.3))
B, R = "#1f4e79", "#c9563c"
g = json.load(open("../emergent/glass_line.json"))
thr = {2: (1.369, 1.414), 3: (2.159, 2.309)}
for ax, K in zip(axes, (2, 3)):
    rows = [r for r in g[str(K)] if r["n1"] < 1]
    U = np.array([r["U"] for r in rows]); Wg = np.array([r["W_glass"] for r in rows]); Wc = np.array([r["Wc3"] for r in rows])
    ax.fill_between(U, 0, Wg, color="#dfe7ef", lw=0)
    ax.fill_between(U, Wg, Wc, color="#f3d9d2", lw=0)
    ax.fill_between(U, Wc, 3.2, color="#e8e8e8", lw=0)
    ax.plot(U, Wg, "o-", color=R, ms=3, lw=1.3); ax.plot(U, Wc, "-", color="black", lw=1.0)
    for u in thr[K]: ax.axvline(u, color="#888888", lw=0.6, ls=":")
    ax.set_title(f"$K={K}$", fontsize=9); ax.set_xlabel("$U/t$"); ax.set_ylim(0, 3.2)
    ax.set_xlim(U.min() - 0.02, U.max() + 0.02)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[0].set_ylabel("disorder $W/t$")
axes[0].text(0.5, 0.12, "liquid", transform=axes[0].transAxes, fontsize=7.5, color=B, ha="center")
axes[0].text(0.5, 0.5, "trimer\nglass", transform=axes[0].transAxes, fontsize=7.5, color=R, ha="center")
axes[0].text(0.5, 0.86, "all localised", transform=axes[0].transAxes, fontsize=7, color="#444444", ha="center")
axes[1].text(0.5, 0.12, "liquid", transform=axes[1].transAxes, fontsize=7.5, color=B, ha="center")
axes[1].text(0.5, 0.55, "trimer\nglass", transform=axes[1].transAxes, fontsize=7.5, color=R, ha="center")
fig.tight_layout(w_pad=1.0); fig.savefig("fig5.pdf"); fig.savefig("fig5.png", dpi=150)
