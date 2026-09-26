"""Schematic of the problem: three distinguishable particles on the Bethe lattice
(a) and the regimes separated by the two thresholds (b)."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch
from matplotlib import rcParams
rcParams.update({"font.size": 9, "axes.linewidth": 0.6})
fig = plt.figure(figsize=(3.4, 3.6))
ax = fig.add_axes([0.0, 0.36, 1.0, 0.64]); ax.set_aspect("equal"); ax.axis("off")
# --- (a) Cayley tree K=2, depth 3, radial layout
K, D = 2, 3
pos = {(): (0.0, 0.0)}; edges = []
def grow(path, ang0, ang1, depth):
    if depth == D: return
    n = K + 1 if depth == 0 else K
    for i in range(n):
        a0 = ang0 + (ang1 - ang0) * i / n; a1 = ang0 + (ang1 - ang0) * (i + 1) / n
        am = 0.5 * (a0 + a1); r = 1.0 + 1.0 * depth
        child = path + (i,); pos[child] = (r * np.cos(am), r * np.sin(am)); edges.append((path, child))
        grow(child, a0, a1, depth + 1)
grow((), np.pi / 2, np.pi / 2 + 2 * np.pi, 0)
for u, v in edges:
    ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="#b0b0b0", lw=0.8, zorder=1)
for v, (x, y) in pos.items():
    ax.plot(x, y, "o", color="white", mec="#808080", ms=3.2, mew=0.6, zorder=2)
# leaves fade: dotted stubs beyond depth 3
for v, (x, y) in pos.items():
    if len(v) == D:
        r = np.hypot(x, y); ax.plot([x, x * (r + 0.45) / r], [y, y * (r + 0.45) / r], color="#d0d0d0", lw=0.8, ls=":", zorder=0)
# particles: A blue at (0,1) depth 2 ; B red at (1,) depth 1 ; C green at (2,0,1) depth 3 ; median = root
cols = {"1": "#1f4e79", "2": "#c9563c", "3": "#3a8f5c"}
parts = {"1": (0, 1), "2": (1,), "3": (2, 0, 1)}
def path_to_root(v): return [v[:i] for i in range(len(v), -1, -1)]
for p, v in parts.items():
    pth = path_to_root(v)
    for u, w in zip(pth[:-1], pth[1:]):
        ax.plot([pos[u][0], pos[w][0]], [pos[u][1], pos[w][1]], color=cols[p], lw=1.6, zorder=3, alpha=0.6)
    x, y = pos[v]; ax.plot(x, y, "o", color=cols[p], ms=7, zorder=5)
    ax.text(x, y, p, color="white", fontsize=6, ha="center", va="center", zorder=6)
x0, y0 = pos[()]; ax.plot(x0, y0, "o", color="black", ms=3.8, zorder=5)
ax.text(x0 + 0.12, y0 - 0.32, "median", fontsize=7, ha="left", va="top")
# distance labels
lab = {"1": ("$a=2$", (0.55, 0.35)), "2": ("$b=1$", (-0.55, -0.15)), "3": ("$c=3$", (0.2, -0.2))}
for p, v in parts.items():
    pth = path_to_root(v); mid = pth[len(pth) // 2]; xm, ym = pos[mid]
    dx, dy = lab[p][1]; ax.text(xm + dx, ym + dy, lab[p][0], color=cols[p], fontsize=7, ha="center", va="center")
# hop arrow from particle 1 to its unused child, label t
v = parts["1"]; w = v + (0,); ax.annotate("", xy=pos[w], xytext=pos[v], arrowprops=dict(arrowstyle="-|>", color="#1f4e79", lw=1.0, shrinkA=5, shrinkB=2), zorder=6)
ax.text(0.5 * (pos[v][0] + pos[w][0]) + 0.28, 0.5 * (pos[v][1] + pos[w][1]), "$-t$", fontsize=7, color="#1f4e79", ha="left", va="center")
# on-site attraction inset: two particles on one site
xi, yi = 2.9, -2.6
ax.plot(xi, yi, "o", color="white", mec="#808080", ms=3.2, mew=0.6, zorder=2)
ax.plot(xi - 0.11, yi, "o", color="#1f4e79", ms=7, zorder=5); ax.plot(xi + 0.11, yi, "o", color="#c9563c", ms=7, zorder=5)
ax.text(xi, yi - 0.4, "$-U$ per pair\non a site", fontsize=7, ha="center", va="top")
ax.text(-3.6, 3.3, "(a)", fontsize=9, fontweight="bold")
ax.set_xlim(-3.8, 3.8); ax.set_ylim(-3.9, 3.7)
# --- (b) attraction axis with regimes
bx = fig.add_axes([0.0, 0.0, 1.0, 0.36]); bx.axis("off"); bx.set_xlim(0, 10.6); bx.set_ylim(0, 3.2)
bx.annotate("", xy=(9.8, 1.0), xytext=(0.5, 1.0), arrowprops=dict(arrowstyle="-|>", color="black", lw=0.8))
bx.text(9.85, 0.9, "$U$", fontsize=9, ha="left", va="top")
for x, name, c in ((4.0, "$U_3$", "#1f4e79"), (7.0, "$U_2$", "#c9563c")):
    bx.plot([x, x], [0.85, 1.15], color=c, lw=1.2); bx.text(x, 0.7, name, fontsize=9, color=c, ha="center", va="top")
bx.fill_between([4.0, 7.0], 1.15, 3.1, color="#e8e8e8", lw=0)
def trio(x, y, bound, pair=False):
    pts = [(x - 0.32, y - 0.15), (x + 0.32, y - 0.15), (x, y + 0.3)]
    for (px, py), c in zip(pts, ("#1f4e79", "#c9563c", "#3a8f5c")):
        bx.plot(px, py, "o", color=c, ms=5.5)
    if bound: bx.add_patch(Circle((x, y), 0.62, fill=False, ls="--", lw=0.8, color="black"))
    if pair: bx.add_patch(Ellipse((x, y - 0.15), 1.05, 0.55, fill=False, ls="--", lw=0.8, color="#c9563c"))
trio(2.2, 2.2, False); bx.text(2.2, 1.35, "free", fontsize=7, ha="center")
trio(5.5, 2.2, True); bx.text(5.5, 1.35, "trimer, no pair", fontsize=7, ha="center")
trio(8.4, 2.2, True, pair=True); bx.text(8.4, 1.35, "pair and trimer", fontsize=7, ha="center")
bx.text(0.2, 3.0, "(b)", fontsize=9, fontweight="bold", va="top")
fig.savefig("fig0.pdf"); fig.savefig("fig0.png", dpi=150)
