"""Two-body Q-sector ground state on girth>=6 random 3-regular graphs, N up to 1600:
N-scaling of the pair coincidence probability P2 and of E2Q - 2 edge1."""
import numpy as np, networkx as nx, scipy.sparse.linalg as spl, scipy.sparse as sp, sys, json, time
from rrg_girth import girth
def pair(N, seed, Us):
    G = nx.random_regular_graph(3, N, seed=seed)
    A = nx.to_scipy_sparse_array(G, format="csr", dtype=float)
    lam = spl.eigsh(A, k=3, which="LA", tol=1e-10)[0]; lam1p = -np.sort(lam)[-2]
    u = np.ones(N) / np.sqrt(N)
    c2 = np.eye(N)
    def Q2(v):
        v = v.reshape(N, N); v = v - np.outer(u, u @ v); v = v - np.outer(v @ u, u); return v.ravel()
    def H2(v, U):
        v = v.reshape(N, N); return (-(A @ v) - (A @ v.T).T - U * c2 * v).ravel()
    out = []
    for U in Us:
        op = spl.LinearOperator((N * N, N * N), matvec=lambda v: Q2(H2(Q2(v), U)), dtype=float)
        e, v = spl.eigsh(op, k=1, which="SA", v0=Q2(np.random.default_rng(seed).standard_normal(N * N)), tol=1e-8, maxiter=5000)
        p2 = (v[:, 0].reshape(N, N) ** 2).diagonal().sum()
        out.append([N, seed, U, lam1p, e[0], p2])
        print(f"N={N} seed={seed} U={U:.2f}: edge1={lam1p:.4f}  E2Q-2edge1={e[0]-2*lam1p:+.4f}  P2={p2:.4f}  N*P2={N*p2:.2f}", flush=True)
    return out
if __name__ == "__main__":
    res = []
    for N in (100, 200, 400, 800, 1600):
        found = []; seed = 0; t0 = time.time()
        while len(found) < 2:
            seed += 1
            if girth(nx.random_regular_graph(3, N, seed=seed)) >= 6: found.append(seed)
        print(f"N={N}: girth>=6 seeds {found} ({time.time()-t0:.0f}s)", flush=True)
        for s in found:
            res += pair(N, s, (0.0, 1.0, 1.4, 1.7, 2.0, 2.3)); json.dump(res, open("rrg_pair_scaling.json", "w"))
