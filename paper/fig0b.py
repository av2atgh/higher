"""Schematic of the uniform-sector calculation: (a) the pair problem reduced
to a half line in the distance d; (b) the moves in the three-body orbit (a,b,c)."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"font.size": 9, "axes.linewidth": 0.6})
B, R, G = "#1f4e79", "#c9563c", "#3a8f5c"
fig = plt.figure(figsize=(3.4, 3.9))
# ---------------- (a) half line
ax = fig.add_axes([0.0, 0.62, 1.0, 0.38]); ax.axis("off"); ax.set_xlim(-0.6, 10.2); ax.set_ylim(-1.9, 2.4)
xs = [0.6, 2.6, 4.6, 6.6, 8.6]
for i, x in enumerate(xs[:-1]):
    ax.plot([x, xs[i + 1]], [0, 0], color="black", lw=0.9)
    ax.text(0.5 * (x + xs[i + 1]), 0.22, "$J_0$" if i == 0 else "$J$", fontsize=8, ha="center", va="bottom")
ax.plot([xs[-1], xs[-1] + 1.0], [0, 0], color="black", lw=0.9, ls=":")
for i, x in enumerate(xs):
    ax.plot(x, 0, "o", color="white", mec="black", ms=7, mew=0.8, zorder=3)
    ax.text(x, -0.45, f"$d={i}$", fontsize=7, ha="center", va="top")
nd = ["$1$", "$K{+}1$", "$(K{+}1)K$", "$(K{+}1)K^{2}$", "$(K{+}1)K^{3}$"]
for x, s in zip(xs, nd): ax.text(x, -1.05, s, fontsize=7, ha="center", va="top")
ax.text(-0.5, -1.05, "$n_d$", fontsize=7, ha="left", va="top")
ax.text(xs[0], 0.35, "$-U$", fontsize=8, ha="center", va="bottom", color=R)
# mini configurations above sites 0,1,2: a short piece of tree (path) with the two particles
def mini(xc, d):
    y = 1.55; step = 0.42
    pts = [xc + (k - d / 2) * step for k in range(d + 1)] if d > 0 else [xc]
    for k in range(len(pts) - 1): ax.plot([pts[k], pts[k + 1]], [y, y], color="#b0b0b0", lw=0.8)
    for p in pts: ax.plot(p, y, "o", color="white", mec="#808080", ms=3, mew=0.6, zorder=2)
    if d == 0:
        ax.plot(xc - 0.07, y, "o", color=B, ms=4.5, zorder=3); ax.plot(xc + 0.07, y, "o", color=R, ms=4.5, zorder=3)
    else:
        ax.plot(pts[0], y, "o", color=B, ms=4.5, zorder=3); ax.plot(pts[-1], y, "o", color=R, ms=4.5, zorder=3)
for x, d in zip(xs[:3], (0, 1, 2)): mini(x, d)
ax.text(-0.5, 2.3, "(a)", fontsize=9, fontweight="bold", va="top")
# ---------------- (b) three-body moves on a tree piece
bx = fig.add_axes([0.0, 0.0, 1.0, 0.62]); bx.set_aspect("equal"); bx.axis("off")
K, D = 2, 3
pos = {(): (0.0, 0.0)}; edges = []
def grow(path, a0, a1, depth):
    if depth == D: return
    n = K + 1 if depth == 0 else K
    for i in range(n):
        b0 = a0 + (a1 - a0) * i / n; b1 = a0 + (a1 - a0) * (i + 1) / n; am = 0.5 * (b0 + b1); r = 1.0 + 0.95 * depth
        c = path + (i,); pos[c] = (r * np.cos(am), r * np.sin(am)); edges.append((path, c)); grow(c, b0, b1, depth + 1)
grow((), np.pi / 2, np.pi / 2 + 2 * np.pi, 0)
for u, v in edges: bx.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="#b0b0b0", lw=0.8, zorder=1)
for v, (x, y) in pos.items(): bx.plot(x, y, "o", color="white", mec="#808080", ms=3.2, mew=0.6, zorder=2)
# particles: 1 at the median (a=0); 2 at distance 1 in branch 1 (b=1); 3 at distance 2 in branch 2 (c=2)
p1, p2, p3 = (), (1,), (2, 0)
m = pos[()]
bx.plot(m[0], m[1], "o", color="black", ms=4, zorder=4)
def dot(v, c, lab, dx=0, dy=0):
    x, y = pos[v]; bx.plot(x + dx, y + dy, "o", color=c, ms=7.5, zorder=5); bx.text(x + dx, y + dy, lab, color="white", fontsize=6, ha="center", va="center", zorder=6)
dot(p1, B, "1", dx=0.0, dy=0.28); dot(p2, R, "2"); dot(p3, G, "3")
bx.text(m[0] - 0.2, m[1] - 0.05, "median", fontsize=7, ha="right", va="top")
def arrow(u, v, c, shrinkA=6, shrinkB=3, rad=0.0):
    bx.annotate("", xy=pos[v], xytext=pos[u], arrowprops=dict(arrowstyle="-|>", color=c, lw=1.1, shrinkA=shrinkA, shrinkB=shrinkB, connectionstyle=f"arc3,rad={rad}"), zorder=6)
# (i) particle 3 outward: c -> c+1, K ways
arrow(p3, p3 + (0,), G); arrow(p3, p3 + (1,), G)
x, y = pos[p3]; bx.text(x + 0.05, y + 0.55, "$c\\to c{+}1$\n$K$ ways", fontsize=6.5, color=G, ha="left", va="bottom")
# (ii) particle 2 inward: b -> b-1, 1 way
arrow(p2, (), R, shrinkA=9, shrinkB=6, rad=0.45)
x, y = pos[p2]; bx.text(x - 0.45, y + 0.15, "$b\\to b{-}1$\n1 way", fontsize=6.5, color=R, ha="right", va="center")
# (iii) particle 1 from the median into a fresh branch: a -> 1, K-1 ways
arrow((), (0,), B, shrinkA=10, shrinkB=4)
x, y = pos[(0,)]; bx.text(x + 0.25, y + 0.25, "$a\\to1$\n$K{-}1$ ways", fontsize=6.5, color=B, ha="left", va="bottom")
# (iv) particle 1 towards particle 2: median moves, (0,b,c) -> (0,b-1,c+1)
bx.annotate("", xy=(pos[(1,)][0] + 0.12, pos[(1,)][1] + 0.05), xytext=(pos[()][0] + 0.12, pos[()][1] + 0.28), arrowprops=dict(arrowstyle="-|>", color=B, lw=1.1, shrinkA=6, shrinkB=9, connectionstyle="arc3,rad=-0.55"), zorder=6)
x, y = pos[(1,)]; bx.plot(x, y, "o", color="none", mec="black", ms=11, mew=1.0, ls="--", zorder=3)
bx.text(x + 0.35, y - 0.05, "new median:\n$(0,b,c)\\to(0,b{-}1,c{+}1)$", fontsize=6.5, color="black", ha="left", va="top")
bx.text(-3.4, 3.1, "(b)", fontsize=9, fontweight="bold", va="top")
bx.set_xlim(-3.5, 3.9); bx.set_ylim(-3.4, 3.2)
fig.savefig("fig0b.pdf"); fig.savefig("fig0b.png", dpi=150)
