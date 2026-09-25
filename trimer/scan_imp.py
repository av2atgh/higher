import numpy as np, json, time, networkx as nx, scipy.sparse as sp, scipy.sparse.linalg as spl
from tree4 import E, impurity_threshold
from tree3 import pair_threshold

# --- check vs ED on a finite Cayley tree with the impurity at the root (localized states)
def cayley(z, depth):
    G = nx.Graph(); G.add_node(0); frontier = [0]; nxt = 1
    for d in range(depth):
        new = []
        for v in frontier:
            for _ in range(z if v == 0 else z - 1):
                G.add_edge(v, nxt); new.append(nxt); nxt += 1
        frontier = new
    return G
def ed(A, m, U, eps0):
    N = A.shape[0]; I = sp.identity(N, format="csr")
    ops = []
    for p in range(m):
        f = [I] * m; f[p] = A
        T = f[0]
        for q in range(1, m): T = sp.kron(T, f[q], format="csr")
        ops.append(T)
    Hk = -sum(ops)
    g = np.indices([N] * m).reshape(m, -1).T
    V = eps0 * (g == 0).sum(1).astype(float)
    for a in range(m):
        for b in range(a + 1, m): V -= U * (g[:, a] == g[:, b])
    H = (Hk + sp.diags(V)).tocsr()
    return spl.eigsh(H, k=1, which="SA", tol=1e-9)[0][0]
K = 2
for m, U, eps0 in ((2, 2.0, -1.0), (3, 2.0, -1.0), (3, 3.0, -0.5)):
    e_red = E(m, K, U, eps0, 40)
    for depth in (4, 5):
        A = nx.to_scipy_sparse_array(cayley(K + 1, depth), format="csr", dtype=float)
        print(f"check K={K} m={m} U={U} eps0={eps0}: Cayley depth {depth} ED {ed(A, m, U, eps0):.6f}  reduced {e_red:.6f}", flush=True)
# single-particle impurity energy closed form: eps0 - E - z t^2 g(E) = 0, g = (-E - sqrt(E^2-4K))/(2K)
from scipy.optimize import brentq
for K in (2, 3):
    eps0 = -1.5 * impurity_threshold(K)
    f = lambda Ev: eps0 - Ev - (K + 1) * (-Ev - np.sqrt(Ev * Ev - 4 * K)) / (2 * K)
    print(f"check K={K} m=1 eps0={eps0:.4f}: closed form {brentq(f, -10, -2*np.sqrt(K)-1e-9):.6f}  reduced {E(1, K, 0.0, eps0, 400):.6f}", flush=True)

# --- thresholds vs impurity depth
def thr(m, K, eps0, L, tol=2e-3):
    edge = -2 * m * np.sqrt(K)
    lo, hi = 0.0, 2 * pair_threshold(K) + 2
    if E(m, K, lo, eps0, L) < edge - 1e-9: return 0.0
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if E(m, K, mid, eps0, L) < edge - 1e-9: hi = mid
        else: lo = mid
    return 0.5 * (lo + hi)
out = {}
for K in (2, 3):
    ec = impurity_threshold(K)
    print(f"\n=== K={K}: eps_c={ec:.4f}, invariant-sector U2={pair_threshold(K):.4f}", flush=True)
    out[K] = []
    for frac in (0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0):
        eps0 = -frac * ec
        row = [eps0]
        for L in (40, 60):
            t0 = time.time()
            u2 = thr(2, K, eps0, 2 * L); u3 = thr(3, K, eps0, L)
            row += [L, u2, u3]
            print(f"  eps0={eps0:8.4f} ({frac:.2f} eps_c)  L={L}: U2={u2:.4f}  U3={u3:.4f}  U3/U2={u3/u2 if u2>0 else float('nan'):.4f}  ({time.time()-t0:.0f}s)", flush=True)
        out[K].append(row)
    json.dump(out, open("scan_imp.json", "w"), indent=1)
