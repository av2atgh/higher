"""Phase diagram of the trimer liquid in the (U, W) plane: glass line W_glass(U) (every state of the
liquid, n < n_1(U), localised) and W_c^(3)(U) (every trimer state localised), K = 2 and 3."""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"font.size": 10, "axes.linewidth": 0.7})
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
B, R = "#1f4e79", "#c9563c"
g = json.load(open("../emergent/glass_line.json")); g1 = json.load(open("../emergent/glass_line_single.json")); Ucompact = {2: 1.7, 3: 2.5}
thr = {2: (1.369, 1.414), 3: (2.159, 2.309)}
for ax, K in zip(axes, (2, 3)):
    rows = [r for r in g[str(K)] if r["n1"] < 1 and r["U"] >= Ucompact[K] - 1e-9]
    r1 = [r for r in g1[str(K)] if r["U"] >= Ucompact[K] - 1e-9]; U1 = np.array([r["U"] for r in r1]); W1 = np.array([r["W1_glass"] for r in r1])
    U = np.array([r["U"] for r in rows]); Wg = np.array([r["W_glass"] for r in rows]); Wc = np.array([r["Wc3"] for r in rows])
    ax.fill_between(U, 0, Wg, color="#dfe7ef", lw=0)
    ax.fill_between(U, Wg, Wc, color="#f3d9d2", lw=0)
    ax.fill_between(U, Wc, 60, color="#e8e8e8", lw=0)
    ax.plot(U1, W1, "s-", color="#666666", ms=3.5, lw=1.2)
    ax.plot(U, Wg, "o-", color=R, ms=4, lw=1.6); ax.plot(U, Wc, "-", color="black", lw=1.3)
    for u in thr[K]: ax.axvline(u, color="#888888", lw=0.6, ls=":")
    ax.set_title(f"$K={K}$", fontsize=10); ax.set_xlabel("$U/t$"); ax.set_yscale("log"); ax.set_ylim(0.3, 60)
    ax.set_xlim(U.min() - 0.02, U.max() + 0.02)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[0].set_ylabel("disorder $W/t$")
for ax, K in zip(axes, (2, 3)):
    rows = [r for r in g[str(K)] if r["n1"] < 1 and r["U"] >= Ucompact[K] - 1e-9]
    U = np.array([r["U"] for r in rows]); Wg = np.array([r["W_glass"] for r in rows]); Wc = np.array([r["Wc3"] for r in rows])
    r1 = [r for r in g1[str(K)] if r["U"] >= Ucompact[K] - 1e-9]; W1 = np.array([r["W1_glass"] for r in r1])
    um = 0.5 * (U.min() + U.max())
    ax.text(um, W1.mean() * 1.45, "single particles at the same density", fontsize=9, color="#666666", ha="center", va="bottom")
    ax.text(um, 0.5, "liquid", fontsize=10, color=B, ha="center", va="center")
    ul = U.min() + 0.22 * (U.max() - U.min())
    ax.annotate("trimer glass", xy=(ul, np.sqrt(np.interp(ul, U, Wg) * np.interp(ul, U, Wc))), xytext=(ul, 5.5), fontsize=10, color=R, ha="center", va="bottom", arrowprops=dict(arrowstyle="-|>", color=R, lw=0.9, shrinkB=2))
    ax.text(U.max(), Wc[-1] * 1.22, "all trimer states localised", fontsize=8, color="#333333", ha="right", va="bottom")
fig.tight_layout(w_pad=2.0); fig.text(0.01, 0.95, "(a)", fontsize=10, fontweight="bold"); fig.text(0.51, 0.95, "(b)", fontsize=10, fontweight="bold"); fig.savefig("fig5.pdf"); fig.savefig("fig5.png", dpi=150)
