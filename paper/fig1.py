import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"font.size": 9, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
d = json.load(open("../trimer/scan.json"))["3"]
K = 3; U2 = d["U2"]; U3 = d["thr"]["180"]
g = np.array(d["grid"])  # U, E3, E2+edge1, E3anti
U, E3, E2p = g[:, 0], g[:, 1], g[:, 2]
b3 = np.minimum(E3 + 6 * np.sqrt(K), 0)
b2 = np.minimum(E2p + 2 * np.sqrt(K) + 4 * np.sqrt(K), 0)  # E2 + 4 sqrt K
fig, ax = plt.subplots(figsize=(3.4, 2.5))
ax.axvspan(U3, U2, color="#e8e8e8", lw=0, zorder=0)
ax.plot(U, b3, "-", color="#1f4e79", lw=1.6, label="trimer, $E_3+6\\sqrt{K}t$")
ax.plot(U, b2, "-", color="#c9563c", lw=1.6, label="pair, $E_2+4\\sqrt{K}t$")
ax.axhline(0, color="0.5", lw=0.6, ls="--")
ax.set_xlabel("$U/t$"); ax.set_ylabel("binding energy / $t$")
ax.set_xlim(0.5, 6); ax.set_ylim(-9, 0.4)
ax.text(U3 - 0.05, -8.6, "$U_3$", ha="right", va="bottom", fontsize=8)
ax.text(U2 + 0.05, -8.6, "$U_2$", ha="left", va="bottom", fontsize=8)
ax.text(4.6, -6.2, "trimer", color="#1f4e79", fontsize=8)
ax.text(4.9, -0.9, "pair", color="#c9563c", fontsize=8)
ax.set_title("$z=4$ ($K=3$)", fontsize=9)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("fig1.pdf")
