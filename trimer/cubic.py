"""Cubic lattice, zero total momentum: pair and trimer thresholds by Birman-Schwinger
in a box |r|_inf <= R of relative coordinates. Pair: r = x1-x2, hopping 2t, edge -12t.
Trimer: r1 = x1-x3, r2 = x2-x3; hops of 1: r1+-e; of 2: r2+-e; of 3: (r1,r2)-+(e,e); edge -18t."""
import numpy as np, scipy.sparse as sp, sys, time
from bs import threshold_from
def pair(R):
    n = 2 * R + 1; N = n ** 3
    g = np.indices((n, n, n)).reshape(3, -1).T - R
    idx = {tuple(x): i for i, x in enumerate(g)}
    rows, cols = [], []
    for i, x in enumerate(g):
        for ax in range(3):
            for s in (-1, 1):
                y = x.copy(); y[ax] += s
                j = idx.get(tuple(y))
                if j is not None: rows.append(i); cols.append(j)
    H0 = sp.csr_matrix((-2.0 * np.ones(len(rows)), (rows, cols)), shape=(N, N))
    V = (np.abs(g).sum(1) == 0).astype(float)
    return threshold_from(H0, V, -12.0)
def trimer(R):
    n = 2 * R + 1; N = n ** 6
    g = np.indices((n,) * 6).reshape(6, -1).T - R
    strides = np.array([n ** k for k in range(5, -1, -1)])
    def index(y):
        ok = np.all(np.abs(y) <= R, axis=1)
        return ok, ((y + R) * strides).sum(1)
    rows, cols = [], []
    moves = []
    for ax in range(3):
        for s in (-1, 1):
            e = np.zeros(6, int); e[ax] = s; moves.append(e)               # particle 1
            e = np.zeros(6, int); e[3 + ax] = s; moves.append(e)           # particle 2
            e = np.zeros(6, int); e[ax] = s; e[3 + ax] = s; moves.append(e)  # particle 3 (moves both)
    src = np.arange(N)
    for e in moves:
        ok, j = index(g + e)
        rows.append(src[ok]); cols.append(j[ok])
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    H0 = sp.csr_matrix((-np.ones(len(rows)), (rows, cols)), shape=(N, N))
    r1, r2 = g[:, :3], g[:, 3:]
    V = ((np.abs(r1).sum(1) == 0).astype(float) + (np.abs(r2).sum(1) == 0) + (np.abs(r1 - r2).sum(1) == 0))
    return threshold_from(H0, V, -18.0)
if __name__ == "__main__":
    for R in (2, 3, 4, 5, 6, 8, 10, 14):
        t0 = time.time(); print(f"pair R={R}: U2={pair(R):.5f}  (Watson: 7.91355)  ({time.time()-t0:.0f}s)", flush=True)
    for R in (2, 3, 4, 5, 6):
        t0 = time.time(); print(f"trimer R={R}: U3={trimer(R):.5f}  ({time.time()-t0:.0f}s)", flush=True)
