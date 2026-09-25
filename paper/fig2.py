import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"font.size": 9, "axes.linewidth": 0.6})
d = json.load(open("../trimer/scan_imp.json"))
fig, axes = plt.subplots(1, 2, figsize=(3.4, 2.2), sharey=False)
cols = {"2": "#1f4e79", "3": "#c9563c"}
for ax, K in zip(axes, ("2", "3")):
    Kf = float(K); ec = (Kf - 1) / np.sqrt(Kf)
    rows = np.array(d[K])            # eps0, L1, u2, u3, L2, u2, u3
    x = -rows[:, 0] / ec
    u2, u3 = rows[:, -2], rows[:, -1]
    ax.plot(x, u2, "-", color="#c9563c", lw=1.4)
    ax.plot(x, u3, "-", color="#1f4e79", lw=1.4)
    ax.fill_between(x, u3, u2, where=u3 < u2, color="#e8e8e8", lw=0)
    ax.set_title(f"$K={K}$", fontsize=9)
    ax.set_xlabel("$|\\epsilon_0|/\\epsilon_c$")
    ax.set_xlim(0, 1); ax.set_ylim(0, None)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[0].set_ylabel("threshold / $t$")
axes[0].text(0.05, axes[0].get_ylim()[1] * 0.92, "pair $U_2$", color="#c9563c", fontsize=8)
axes[0].text(0.05, axes[0].get_ylim()[1] * 0.80, "trimer $U_3$", color="#1f4e79", fontsize=8)
fig.tight_layout(); fig.savefig("fig2.pdf")
