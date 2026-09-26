"""Three distinguishable particles on finite random regular graphs, K=2 (degree 3),
in the sector orthogonal to the Perron (uniform) mode on every particle: Q = (1-P)^{x3}.
Reports bulk single-particle edge lam1' (lowest non-Perron eigenvalue), E2^Q, E3^Q."""
import numpy as np, networkx as nx, scipy.sparse as sp, scipy.sparse.linalg as spl, sys, json, time
K = 2; z = K + 1
def run(N, seed, Us):
    G = nx.random_regular_graph(z, N, seed=seed)
    A = nx.to_scipy_sparse_array(G, format="csr", dtype=float)
    lam = np.linalg.eigvalsh(A.toarray())[::-1]          # descending: lam[0] = z
    lam1p = -lam[1]; lam2p = -lam[2]                      # bulk edge energies (t=1): -lam
    u = np.ones(N) / np.sqrt(N)
    idx = np.arange(N)
    # two-body in Q
    I2 = np.stack(np.meshgrid(idx, idx, indexing="ij"), -1).reshape(-1, 2)
    c2 = (I2[:, 0] == I2[:, 1]).astype(float)
    def Q2(v):
        v = v.reshape(N, N); v = v - np.outer(u, u @ v); v = v - np.outer(v @ u, u); return v.ravel()
    def H2(v, U):
        v = v.reshape(N, N); return (-(A @ v) - (A @ v.T).T - U * c2.reshape(N, N) * v).ravel()
    # three-body in Q
    I3 = np.indices((N, N, N)).reshape(3, -1).T
    c3 = ((I3[:, 0] == I3[:, 1]) + (I3[:, 0] == I3[:, 2]) + (I3[:, 1] == I3[:, 2])).astype(float).reshape(N, N, N)
    def Q3(v):
        v = v.reshape(N, N, N)
        v = v - np.tensordot(u, np.tensordot(u, v, axes=(0, 0)), axes=0)
        v = v - np.einsum("j,ik->ijk", u, np.tensordot(u, v, axes=(0, 1)))
        v = v - np.einsum("k,ij->ijk", u, np.tensordot(u, v, axes=(0, 2)))
        return v.ravel()
    def H3(v, U):
        v = v.reshape(N, N, N)
        w = -np.tensordot(A.toarray(), v, axes=(1, 0)) - np.einsum("jl,ilk->ijk", A.toarray(), v) - np.tensordot(v, A.toarray(), axes=(2, 1))
        return (w - U * c3 * v).ravel()
    res = []
    for U in Us:
        t0 = time.time()
        op2 = spl.LinearOperator((N * N, N * N), matvec=lambda v: Q2(H2(Q2(v), U)), dtype=float)
        v0 = Q2(np.random.default_rng(seed).standard_normal(N * N))
        e2 = spl.eigsh(op2, k=1, which="SA", v0=v0, tol=1e-8)[0][0]
        op3 = spl.LinearOperator((N ** 3, N ** 3), matvec=lambda v: Q3(H3(Q3(v), U)), dtype=float)
        v0 = Q3(np.random.default_rng(seed + 1).standard_normal(N ** 3))
        e3 = spl.eigsh(op3, k=1, which="SA", v0=v0, tol=1e-7, maxiter=3000)[0][0]
        res.append([N, seed, U, lam1p, lam2p, e2, e3])
        print(f"N={N} seed={seed} U={U:.2f}: edge1={lam1p:.4f} edge2={lam2p:.4f}  E2Q={e2:.4f} (E2Q-2edge1={e2-2*lam1p:+.4f})  E3Q={e3:.4f} (E3Q-3edge1={e3-3*lam1p:+.4f}, E3Q-E2Q-edge1={e3-e2-lam1p:+.4f})  ({time.time()-t0:.0f}s)", flush=True)
    return res
if __name__ == "__main__":
    N = int(sys.argv[1]); out = []
    for seed in (1, 2, 3):
        out += run(N, seed, (0.0, 1.0, 1.4, 1.7, 2.0, 2.3))
        json.dump(out, open(f"rrg_q_N{N}.json", "w"))
