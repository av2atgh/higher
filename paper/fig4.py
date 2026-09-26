import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"font.size": 9, "axes.linewidth": 0.6})
fig, ax = plt.subplots(figsize=(3.4, 2.4))
B, R = "#1f4e79", "#c9563c"
loc = json.load(open("../emergent/localization.json"))
U2L2 = {2: 1.966, 3: 3.251}; Ucompact = {2: 1.7, 3: 2.5}
for K, c in ((2, B), (3, R)):
    rows = loc[str(K)]; U = np.array([r["U"] for r in rows]); Wc3 = np.array([r["Wc3"] for r in rows]); Wc2 = np.array([r["Wc2"] for r in rows])
    m3 = U >= Ucompact[K] - 1e-9; ax.plot(U[m3], Wc3[m3], "-", color=c, lw=1.4)
    m2 = U >= U2L2[K]; ax.plot(U[m2], Wc2[m2], "--", color=c, lw=1.0)
    ax.axhline({2: 17.39, 3: 32.26}[K], color=c, lw=0.6, ls=":")
ax.set_yscale("log"); ax.set_xlabel("$U/t$"); ax.set_ylabel("critical disorder $W_c/t$"); ax.set_ylim(0.3, 50); ax.set_xlim(1.5, 5.2)
ax.text(1.6, 20, "one particle, $K=2$", fontsize=7, color=B); ax.text(1.6, 37, "one particle, $K=3$", fontsize=7, color=R)
ax.text(4.2, 6.8, "pair", fontsize=7.5, color="#555555"); ax.text(4.2, 0.9, "trimer", fontsize=7.5, color="#555555")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("fig4.pdf"); fig.savefig("fig4.png", dpi=150)
