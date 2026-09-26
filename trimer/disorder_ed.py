"""Finite Cayley tree with box disorder on every site, three distinguishable particles.
Ground state of m=2 and m=3 at attraction U; coincidence probabilities:
P2 = prob(two particles on one site) in the two-body ground state,
P3tri = prob(all three on one site), P3pair = prob(exactly two on one site) in the three-body ground state.
Borromean signature: P3tri large while P2 small."""
import numpy as np, networkx as nx, scipy.sparse as sp, scipy.sparse.linalg as spl, sys, json, time
def cayley(z, depth):
    G = nx.Graph(); G.add_node(0); frontier = [0]; nxt = 1
    for d in range(depth):
        new = []
        for v in frontier:
            for _ in range(z if v == 0 else z - 1):
                G.add_edge(v, nxt); new.append(nxt); nxt += 1
        frontier = new
    return G
def ground(A, eps, m, U):
    N = A.shape[0]; I = sp.identity(N, format="csr")
    Hk = None
    for p in range(m):
        f = [I] * m; f[p] = -A
        T = f[0]
        for q in range(1, m): T = sp.kron(T, f[q], format="csr")
        Hk = T if Hk is None else Hk + T
    g = np.indices([N] * m).reshape(m, -1).T
    V = eps[g].sum(1)
    coinc = np.zeros(len(g))
    for a in range(m):
        for b in range(a + 1, m): coinc += (g[:, a] == g[:, b])
    H = (Hk + sp.diags(V - U * coinc)).tocsr()
    e, v = spl.eigsh(H, k=1, which="SA", tol=1e-8)
    w = v[:, 0] ** 2
    if m == 2: return e[0], [w[coinc > 0].sum()]
    return e[0], [w[coinc == 3].sum(), w[coinc == 1].sum()]
if __name__ == "__main__":
    K = int(sys.argv[1]); depth = int(sys.argv[2]); z = K + 1
    G = cayley(z, depth); A = nx.to_scipy_sparse_array(G, format="csr", dtype=float); N = A.shape[0]
    U2c, U3c = {2: (1.9744, 1.3793), 3: (3.2570, 2.1652)}[K]     # clean L2 thresholds at eps0 = 0
    out = []
    for W in (1.0, 2.0):
        for Uf in (0.4, 0.55, 0.7, 0.8, 0.9, 1.0, 1.15):
            U = Uf * U2c; acc = []
            for seed in range(8):
                eps = np.random.default_rng(1000 * seed + int(10 * W)).uniform(-W / 2, W / 2, N)
                e2, (p2,) = ground(A, eps, 2, U); e3, (ptri, ppair) = ground(A, eps, 3, U)
                acc.append([p2, ptri, ppair])
            a = np.array(acc); mean = a.mean(0); std = a.std(0) / np.sqrt(len(a))
            out.append([W, U, Uf] + mean.tolist() + std.tolist())
            print(f"K={K} N={N} W={W} U={U:.3f} ({Uf:.2f} U2_L2): P2={mean[0]:.3f}+-{std[0]:.3f}  P3tri={mean[1]:.3f}+-{std[1]:.3f}  P3pair={mean[2]:.3f}+-{std[2]:.3f}", flush=True)
            json.dump(out, open(f"disorder_ed_K{K}_d{depth}.json", "w"))
