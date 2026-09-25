import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"font.size": 9, "axes.linewidth": 0.6})
fig, axes = plt.subplots(1, 2, figsize=(3.4, 2.3))
for ax, K in zip(axes, (2, 3)):
    ec = (K - 1) / np.sqrt(K)
    r = np.array(json.load(open(f"../trimer/scan_bs_K{K}_L40.json")))  # eps0, frac, L, u2, u3
    x, u2, u3 = r[:, 1], r[:, 3], r[:, 4]
    ax.fill_between(x, u3, u2, color="#e8e8e8", lw=0)
    ax.plot(x, u2, "-", color="#c9563c", lw=1.4)
    ax.plot(x, u3, "-", color="#1f4e79", lw=1.4)
    # uniform-sector thresholds of the clean tree, for reference
    U2u = 2 * (K - 1) / np.sqrt(K); U3u = {2: 1.310, 3: 2.073}[K]
    ax.plot([0], [U2u], "o", color="#c9563c", ms=4, mfc="white", mew=1.2)
    ax.plot([0], [U3u], "o", color="#1f4e79", ms=4, mfc="white", mew=1.2)
    ax.set_title(f"$K={K}$", fontsize=9)
    ax.set_xlabel("$|\\epsilon_0|/\\epsilon_c$")
    ax.set_xlim(-0.03, 1); ax.set_ylim(0, None)
    ax.set_xticks([0, 0.5, 1])
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[0].set_ylabel("threshold / $t$")
axes[0].text(0.35, 1.72, "pair $U_2$", color="#c9563c", fontsize=8)
axes[0].text(0.04, 0.95, "trimer $U_3$", color="#1f4e79", fontsize=8)
fig.tight_layout(w_pad=1.0); fig.savefig("fig2.pdf")
