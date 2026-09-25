import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, networkx as nx
from tree3 import build, ground, pair_threshold, pair_energy

def ed_graph(A, U, t=1.0):
    """Three distinguishable particles on graph with adjacency A (N x N sparse)."""
    N = A.shape[0]
    I = sp.identity(N, format="csr")
    Hk = -t * (sp.kron(sp.kron(A, I), I) + sp.kron(sp.kron(I, A), I) + sp.kron(sp.kron(I, I), A))
    i = np.arange(N)
    g = np.stack(np.meshgrid(i, i, i, indexing="ij"), -1).reshape(-1, 3)
    V = -U * ((g[:, 0] == g[:, 1]).astype(float) + (g[:, 0] == g[:, 2]) + (g[:, 1] == g[:, 2]))
    H = (Hk + sp.diags(V)).tocsr()
    return spl.eigsh(H, k=1, which="SA", tol=1e-9)[0][0]

# 1) K=1: infinite line vs ring ED (ground state is in the K=0 sector)
for U in (1.0, 3.0):
    H, S = build(1, U, 120)
    e_tree = ground(H)[0]
    A = nx.to_scipy_sparse_array(nx.cycle_graph(41), format="csr", dtype=float)
    print(f"K=1 U={U}: tree-reduced {e_tree:.8f}  ring41 ED {ed_graph(A, U):.8f}")

# 2) K=3: infinite tree vs finite random 4-regular graph, tightly bound trimer
for U in (4.0, 6.0):
    H, S = build(3, U, 60)
    e_tree = ground(H)[0]
    es = []
    for seed in (1, 2):
        G = nx.random_regular_graph(4, 40, seed=seed)
        A = nx.to_scipy_sparse_array(G, format="csr", dtype=float)
        es.append(ed_graph(A, U))
    print(f"K=3 U={U}: tree-reduced {e_tree:.6f}  RRG(4,40) ED {es}")

# 3) large U: E3 = -3U - 3 z t^2/(2U) + O(t^4/U^3)
for K in (2, 3):
    z = K + 1
    U = 30.0
    H, S = build(K, U, 40)
    print(f"K={K} U={U}: {ground(H)[0]:.6f}  2nd order {-3*U - 3*z/(2*U):.6f}")

# 4) pair: numeric half-line vs closed form
for K in (2, 3):
    J = 2*np.sqrt(K); J0 = 2*np.sqrt(K+1)
    for U in (pair_threshold(K)*1.5, 4.0):
        n = 400
        Hp = np.zeros((n, n)); Hp[0, 0] = -U
        Hp[0, 1] = Hp[1, 0] = -J0
        for d in range(1, n-1):
            Hp[d, d+1] = Hp[d+1, d] = -J
        print(f"K={K} U={U:.4f}: pair numeric {np.linalg.eigvalsh(Hp)[0]:.8f}  closed form {pair_energy(K, U):.8f}  U2={pair_threshold(K):.5f}")
